"""
analyze_screen_brightness_map.py

Screen brightness at each location of a ScreenBrightnessMap run (labpack), from one continuous
spectrometer recording made with clandinin_labpack's misc/record_spectrometer.py (--preset map).

    python analyze_screen_brightness_map.py --recording E:/ProjectorCalibration/..._map_center.npz

Run from the visanalysis repo root. The protocol's layout and phase times must match the run;
the defaults are the protocol's defaults, and --radii, --n-directions and the phase times
override them.

HOW LOCATIONS ARE FOUND. Each location ends with `black_time` of a black full field and then
`white_time` of a white one. That [dark][bright][dark] shape -- the last dark being the gap to the
next trial -- is matched against the whole recording, and the N best non-overlapping matches, in
time order, are the N locations. Nothing depends on the recording's clock matching the rig's, and
pauses between locations do not matter.

WHAT IS REPORTED, per visit (locations.csv) and per spot averaged over its visits (spots.csv),
plus summary.json and brightness_map.png:
* white: full-field brightness there, relative to the axis. By default every visit is divided by
  the AVERAGE of all the axis visits: the light source drifts little (~0.5% over a run), while
  each axis visit carries its own placement error, so averaging them is the better reference.
  --reference interp instead divides by the axis interpolated in time between the axis visits,
  for a run whose source really does drift -- at the cost of letting one mis-seated axis visit
  shift the spots measured around it;
* black: the black level there, as a percentage of that location's own white -- scattered light
  can make it vary across the bowl;
* aim: read over the last `--aim-window` seconds of the blink -- the moment the probe has to be
  in place, since the black and white phases follow directly. Only the spot is lit then, and it
  is on half the time, so the reading is doubled; a probe on the spot reads close to the white
  one. aim / white well below ~0.5 means the probe was not settled on the spot when the blink
  ended; treat that location with suspicion. Still positioning earlier in the blink is fine.
* spread: with --n-passes > 1 every spot is visited more than once, each a fresh placement of the
  probe, and the spread between visits is the hand-placement error. A spot whose visits differ by
  more than --flag-spread percent is flagged. A bias repeated at every visit is not caught.

OLDER RECORDINGS used other timings; pass them explicitly, e.g. for the 2026-09-28 single-pass
runs: --n-passes 1 --aim-time 15 --black-time 15 --white-time 15.

Each phase is read only after `--settle` seconds, because the spectrometer carries part of each
reading into the next ones (see clandinin_labpack docs/projector-photometry-brukerjr.md).
"""
import argparse
import json
import pathlib
import sys

import numpy as np

from visanalysis.util.octave_stim import _find_labpack

_find_labpack()
from labpack.visual_stim.clandinin.brightness_map import map_points, phase_times  # noqa: E402


def load(path, band=None):
    d = np.load(path, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
    darks = [x for x in (d["dark_before"], d["dark_after"]) if len(x)]
    dark = np.concatenate(darks).mean(axis=0)
    wl = d["wavelengths"]
    S = d["spectra"].astype(np.float64) - dark
    if band is None:
        total = S.sum(axis=1)
        bright = S[total >= np.quantile(total, 0.9)].mean(axis=0)
        i = int(np.argmax(bright))
        lo = hi = i
        while lo > 0 and bright[lo - 1] > 0.05 * bright[i]:
            lo -= 1
        while hi < len(bright) - 1 and bright[hi + 1] > 0.05 * bright[i]:
            hi += 1
        band = (float(wl[lo]), float(wl[hi]))
    m = (wl >= band[0]) & (wl <= band[1])
    t = d["t_perf"] - d["t_perf"][0]
    return dict(t=t, y=S[:, m].sum(axis=1), meta=meta, band=band,
                integ=meta["integration_us"] * 1e-6)


def find_measurements(t, y, n, black_time, white_time, gap):
    """Start times (of the black phase) of the n best [dark][bright][dark] matches."""
    dt = 0.05
    grid = np.arange(t[0], t[-1], dt)
    yi = np.interp(grid, t, y)
    L = int(round((black_time + white_time + gap) / dt))
    template = np.concatenate([np.zeros(int(round(black_time / dt))),
                               np.ones(int(round(white_time / dt))),
                               np.zeros(L - int(round(black_time / dt))
                                        - int(round(white_time / dt)))])
    tc = template - template.mean()
    scores = np.full(len(grid), -np.inf)
    for i in range(0, len(grid) - L):
        w = yi[i:i + L]
        wc = w - w.mean()
        den = np.sqrt((wc * wc).sum() * (tc * tc).sum())
        if den > 0:
            scores[i] = (wc * tc).sum() / den
    chosen = []
    min_sep = int(round((black_time + white_time) / dt))
    for i in np.argsort(scores)[::-1]:
        if not np.isfinite(scores[i]):
            break
        if all(abs(i - j) >= min_sep for j in chosen):
            chosen.append(i)
        if len(chosen) == n:
            break
    chosen = sorted(chosen)
    return [float(grid[i]) for i in chosen], [float(scores[i]) for i in chosen]


def window_median(t, y, integ, a, b):
    """Median of spectra whose whole integration window [t - integ, t] lies in [a, b]."""
    m = (t - integ >= a) & (t <= b)
    return (float(np.median(y[m])) if m.any() else np.nan), int(m.sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--recording", required=True)
    ap.add_argument("--radii", default="30, 55")
    ap.add_argument("--n-directions", type=int, default=4)
    ap.add_argument("--n-passes", type=int, default=2)
    ap.add_argument("--aim-time", type=float, default=15.0)
    ap.add_argument("--blink-time", type=float, default=5.0)
    ap.add_argument("--black-time", type=float, default=8.0)
    ap.add_argument("--white-time", type=float, default=8.0)
    ap.add_argument("--reference", choices=("mean", "interp"), default="mean",
                    help="divide by the mean of the axis visits, or interpolate them in time")
    ap.add_argument("--aim-window", type=float, default=2.0,
                    help="seconds at the end of the blink over which the aim check is read")
    ap.add_argument("--flag-spread", type=float, default=10.0,
                    help="flag a spot whose visits differ by more than this, in percent")
    ap.add_argument("--gap", type=float, default=1.0,
                    help="dark seconds between trials (the protocol's tail_time + pre_time)")
    ap.add_argument("--settle", type=float, default=5.0,
                    help="seconds each phase is left to settle before it is read")
    ap.add_argument("--band", type=float, nargs=2, metavar=("NM_LO", "NM_HI"))
    ap.add_argument("--out")
    args = ap.parse_args()

    ph = phase_times(args.aim_time, args.blink_time, args.black_time, args.white_time)
    radii = [float(x) for x in args.radii.split(",") if x.strip()]
    points = map_points(radii, args.n_directions, n_passes=args.n_passes)
    rec = load(args.recording, args.band)
    t, y, integ = rec["t"], rec["y"], rec["integ"]

    starts, scores = find_measurements(t, y, len(points), args.black_time, args.white_time,
                                       args.gap)
    if len(starts) < len(points):
        sys.exit(f"found only {len(starts)} of {len(points)} black/white measurements")

    rows = []
    edge = 0.15
    for p, s0, sc in zip(points, starts, scores):
        white, nw = window_median(t, y, integ, s0 + args.black_time + args.settle,
                                  s0 + args.black_time + args.white_time - edge)
        black, nb = window_median(t, y, integ, s0 + args.settle, s0 + args.black_time - edge)
        # The end of the blink, where the probe must be in place. The spot is on for half of
        # each blink cycle, so the MEAN (not the median, which would pick one state) is doubled.
        a0 = s0 - min(args.aim_window, args.blink_time)
        m = (t - integ >= a0) & (t <= s0 - edge)
        aim, na = (2.0 * float(np.mean(y[m])) if m.any() else np.nan), int(m.sum())
        rows.append(dict(p, time=s0, match=sc, white=white, black=black, aim=aim,
                         n_white=nw, n_black=nb, n_aim=na))

    axis = [r for r in rows if r["is_axis"]]
    ta = np.array([r["time"] for r in axis])
    wa = np.array([r["white"] for r in axis])
    ref = (np.full(len(rows), wa.mean()) if args.reference == "mean"
           else np.interp([r["time"] for r in rows], ta, wa))
    for r, f in zip(rows, ref):
        r["white_rel"] = r["white"] / f
        r["black_percent"] = 100.0 * r["black"] / r["white"]
        r["aim_ratio"] = r["aim"] / r["white"] if r["white"] > 0 else np.nan
    drift = 100.0 * (wa.max() - wa.min()) / wa.mean()

    # Per spot, over its visits. The axis is excluded from the spread: it is the reference each
    # visit is divided by, so its white_rel is 1 by construction.
    spots = []
    for label in dict.fromkeys(r["label"] for r in rows):
        v = [r for r in rows if r["label"] == label]
        w = np.array([r["white_rel"] for r in v])
        spread = (100.0 * (w.max() - w.min()) / w.mean()) if len(w) > 1 and not v[0]["is_axis"] \
            else np.nan
        spots.append(dict(label=label, rho=v[0]["rho"], psi=v[0]["psi"], az=v[0]["az"],
                          el=v[0]["el"], is_axis=v[0]["is_axis"], n_visits=len(v),
                          white_rel=float(w.mean()), spread_percent=float(spread),
                          visits=[float(x) for x in w],
                          black_percent=float(np.mean([r["black_percent"] for r in v])),
                          min_aim_ratio=float(np.nanmin([r["aim_ratio"] for r in v]))))

    rp = pathlib.Path(args.recording)
    out = pathlib.Path(args.out) if args.out else rp.with_name(rp.stem + "_analysis")
    out.mkdir(parents=True, exist_ok=True)
    keys = ["index", "label", "rho", "psi", "time", "match", "white_rel", "black_percent",
            "aim_ratio", "white", "black", "aim", "n_white", "n_black", "n_aim"]
    with open(out / "locations.csv", "w") as f:
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join(f"{r[k]:.6g}" if isinstance(r[k], float) else str(r[k])
                             for k in keys) + "\n")
    skeys = ["label", "rho", "psi", "n_visits", "white_rel", "spread_percent", "black_percent",
             "min_aim_ratio"]
    with open(out / "spots.csv", "w") as f:
        f.write(",".join(skeys) + "\n")
        for sp in spots:
            f.write(",".join(f"{sp[k]:.6g}" if isinstance(sp[k], float) else str(sp[k])
                             for k in skeys) + "\n")
    with open(out / "summary.json", "w") as f:
        json.dump(dict(recording=str(rp), meta=rec["meta"], band_nm=rec["band"],
                       axis_spread_percent=drift, reference=args.reference, radii=radii,
                       n_directions=args.n_directions,
                       n_passes=args.n_passes, spots=spots, locations=rows),
                  f, indent=2, default=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    ax[0, 0].plot(t, y, lw=0.4, color="k")
    for r in rows:
        ax[0, 0].axvspan(r["time"], r["time"] + args.black_time, color="C0", alpha=0.2)
        ax[0, 0].axvspan(r["time"] + args.black_time,
                         r["time"] + args.black_time + args.white_time, color="C1", alpha=0.2)
        ax[0, 0].text(r["time"], y.max(), str(r["index"]), fontsize=7)
    ax[0, 0].set(title="recording; blue = black phase, orange = white phase",
                 xlabel="s", ylabel="counts in band")
    for k, psi in enumerate(sorted(set(sp["psi"] for sp in spots if not sp["is_axis"]))):
        these = [sp for sp in spots if not sp["is_axis"] and sp["psi"] == psi]
        lab = these[0]["label"].split()[0]
        ax[0, 1].plot([0.0] + [sp["rho"] for sp in these],
                      [1.0] + [sp["white_rel"] for sp in these], "-", color=f"C{k}", label=lab)
        for sp in these:
            ax[0, 1].plot([sp["rho"]] * len(sp["visits"]), sp["visits"], "o", color=f"C{k}",
                          ms=4, alpha=0.7)
    ax[0, 1].set(title="white, relative to the axis (dots = visits, line = mean)",
                 xlabel="degrees from the axis", ylabel="relative brightness")
    ax[0, 1].legend()
    ax[1, 0].bar([r["index"] for r in rows], [r["black_percent"] for r in rows])
    ax[1, 0].set(title="black level, % of that location's white", xlabel="location")
    ax[1, 1].bar([r["index"] for r in rows], [r["aim_ratio"] for r in rows])
    ax[1, 1].axhline(0.5, color="C3", lw=0.8)
    ax[1, 1].set(title="aim check at the end of the blink / white (low = not in place)",
                 xlabel="location")
    fig.tight_layout()
    fig.savefig(out / "brightness_map.png", dpi=120)
    plt.close(fig)

    print(f"band {rec['band'][0]:.1f}-{rec['band'][1]:.1f} nm; axis visits spread "
          f"{drift:.2f}% (placement plus drift); reference: {args.reference} of the axis visits")
    for r in rows:
        flag = "" if r["aim_ratio"] >= 0.5 else "   <- aim check low"
        print(f"{r['index']:2d} {r['label']:9s} white {r['white_rel']:6.3f}  black "
              f"{r['black_percent']:5.2f}%  aim {r['aim_ratio']:.2f}  match r {r['match']:.2f}"
              + flag)
    if args.n_passes > 1:
        print("per spot, over visits:")
        for sp in spots:
            if sp["is_axis"]:
                continue
            flag = "   <- visits disagree" if sp["spread_percent"] > args.flag_spread else ""
            print(f"   {sp['label']:9s} white {sp['white_rel']:6.3f}  visits "
                  + " / ".join(f"{v:.3f}" for v in sp["visits"])
                  + f"  spread {sp['spread_percent']:5.1f}%" + flag)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
