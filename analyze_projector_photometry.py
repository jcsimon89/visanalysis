"""
analyze_projector_photometry.py

Projector black level and intensity-vs-code curve from a spectrometer recording of the
ProjectorPhotometrySweep protocol (labpack), recorded with clandinin_labpack's
misc/record_spectrometer.py.

    python analyze_projector_photometry.py --full E:/ProjectorCalibration/..._full_center.npz \
        --low E:/ProjectorCalibration/..._low_center.npz

Run from the visanalysis repo root (like the other analyze_* scripts): the sweep schedule is
rebuilt from labpack, found next to this checkout. --low is optional.

WHAT IT DOES
------------
1. Dark-subtracts every spectrum with the recording's own dark spectra (taken with the projector
   blocked, so they include stray room light) and sums the counts over the blue LED's band, found
   automatically from the brightest spectra unless --band is given.
2. Aligns the recording to the schedule from the light alone. Each end of the sweep carries a
   13-element Barker code; both are located by normalised cross-correlation against a template
   that accounts for each spectrum's integration window, and the pair gives the clock offset and
   the clock-rate ratio. No wiring between the spectrometer and the rig is needed.
3. Keeps only spectra whose whole integration window lies inside one step's MEASURED part --
   after the schedule's `settle` seconds, which the detector needs to shed its memory of the
   previous step, and `--edge` seconds clear of the step's end -- and takes the median per step.
4. Divides out drift using the reference steps (piecewise-linear in time).
5. Fits I = b + k * code. b / k is the line's intercept in code units, and I(0) / k is code 0's
   own light in code units; they differ exactly when the bottom of the curve is not straight, so
   read them together. Residuals are in code units too, so they read directly as "how far off the
   ideal line this code lands".
6. With --low: fits the low pass on its own (the black level in code units is scale-free, so it
   needs no merge), then scales it onto the full pass through the codes both contain and uses it
   for codes <= its maximum in the combined table.

WHAT IT REPORTS (summary.json, lut.csv, PNGs in --out)
-------------------------------------------------------
* Black level: b/k in codes, and as a fraction of code 255.
* Linearity: residual per code; worst below code 32 and overall.
* CORNER SQUARE: the pre_time idle is black with the corner square NOT toggling; the sweep's
  code-0 steps are black with it toggling. A difference means the square's light reaches the
  fiber and every reading carries it.
* SYMMETRIC PAIRS about 127 used by the ternary stimuli (0/254, 1/253, 16/238, 102/152, ...):
  whether each pair's mean equals I(127), i.e. whether a stimulus built on it keeps its mean.
* LogBandTernaryNoise's default ladder: the probe's log excursion at each mean from MEASURED
  intensities, and its mismatch -- the quantity its top_code search optimises in codes.
* lut.csv: relative intensity (I(255) = 1) for every measured code.
"""
import argparse
import json
import pathlib
import sys

import numpy as np

from visanalysis.util.octave_stim import _find_labpack

_find_labpack()
from labpack.visual_stim.clandinin.photometry_sweep import (  # noqa: E402
    BARKER_13, PhotometrySweepSchedule)

#: Symmetric pairs about 127 that stimuli in labpack show.
PAIRS = [(0, 254), (1, 253), (10, 244), (16, 238), (43, 211), (64, 190), (85, 169),
         (95, 159), (102, 152)]


# ── loading ─────────────────────────────────────────────────────────────────

def load(path):
    d = np.load(path, allow_pickle=False)
    meta = json.loads(str(d["meta"]))
    darks = [x for x in (d["dark_before"], d["dark_after"]) if len(x)]
    if not darks:
        raise ValueError(f"{path}: no dark spectra")
    dark = np.concatenate(darks).mean(axis=0)
    t = d["t_perf"] - d["t_perf"][0]
    return dict(path=str(path), meta=meta, wl=d["wavelengths"],
                spectra=d["spectra"].astype(np.float64) - dark, t=t,
                integ=meta["integration_us"] * 1e-6, dark=dark)


def find_band(rec, frac=0.05):
    """Contiguous region around the blue peak where the bright mean spectrum is above frac * peak."""
    s = rec["spectra"]
    total = s.sum(axis=1)
    bright = s[total >= np.quantile(total, 0.9)].mean(axis=0)
    i = int(np.argmax(bright))
    thr = frac * bright[i]
    lo = i
    while lo > 0 and bright[lo - 1] > thr:
        lo -= 1
    hi = i
    while hi < len(bright) - 1 and bright[hi + 1] > thr:
        hi += 1
    return float(rec["wl"][lo]), float(rec["wl"][hi]), bright


def band_signal(rec, band):
    m = (rec["wl"] >= band[0]) & (rec["wl"] <= band[1])
    return rec["spectra"][:, m].sum(axis=1)


# ── alignment ───────────────────────────────────────────────────────────────

def _box_template(sched, t, integ, a, rate, which):
    """Fraction of each spectrum's integration window spent at marker_high, for the marker
    starting at schedule time marker_starts[which]. Windows are [t - integ, t] on the PC clock;
    the offset `a` absorbs any constant readout latency."""
    m0 = sched.marker_starts[which]
    unit = sched.marker_unit
    edges = m0 + unit * np.arange(len(BARKER_13) + 1)
    bits = np.array(BARKER_13, dtype=float)
    lo = (t - integ - a) / rate
    hi = (t - a) / rate
    out = np.zeros_like(t)
    for k in range(len(bits)):
        if bits[k]:
            ov = np.clip(np.minimum(hi, edges[k + 1]) - np.maximum(lo, edges[k]), 0, None)
            out += ov
    return out / np.maximum(hi - lo, 1e-12)


def _corr(y, x):
    x = x - x.mean()
    y = y - y.mean()
    den = np.sqrt((x * x).sum() * (y * y).sum())
    return float((x * y).sum() / den) if den > 0 else -1.0


def _marker_corr(sched, t, y, integ, a, rate, which):
    """Correlation of the recording with marker `which` under the clock map (a, rate), over the
    marker plus half a unit of its OUTER side, which is dark (marker_low) as the template assumes.
    The inner side is a reference step at another level, so it is left out."""
    m0 = sched.marker_starts[which]
    pad = 0.5 * sched.marker_unit
    before, after = (pad, 0.0) if which == 0 else (0.0, pad)
    i0, i1 = np.searchsorted(t, [a + rate * m0 - before,
                                 a + rate * (m0 + sched.marker_duration) + after])
    if i1 - i0 < 2 * len(BARKER_13):
        return -1.0
    return _corr(y[i0:i1], _box_template(sched, t[i0:i1], integ, a, rate, which))


def align(sched, rec, y):
    """Clock map t_pc = a + rate * tau from the two Barker markers."""
    t, integ = rec["t"], rec["integ"]
    tau2 = sched.marker_starts[1]
    span = sched.marker_duration
    step = min(0.02, sched.marker_unit / 10)
    cands = np.arange(t[0] - 1.0, t[-1] - span, step)
    c1 = np.array([_marker_corr(sched, t, y, integ, a, 1.0, 0) for a in cands])
    # The second marker is tau2 later; search it only where a first-marker candidate places it.
    best = (-3.0, None, None)
    for i in np.argsort(c1)[::-1][:20]:
        a1 = cands[i]
        for rate in 1.0 + np.linspace(-5e-4, 5e-4, 101):
            c2 = _marker_corr(sched, t, y, integ, a1, rate, 1)
            if c1[i] + c2 > best[0]:
                best = (c1[i] + c2, a1, rate)
    _, a1, rate = best
    # Refine jointly on a fine grid.
    fine = []
    for da in np.arange(-0.05, 0.0501, 0.002):
        for dr in np.linspace(-2e-5, 2e-5, 21):
            c = _marker_corr(sched, t, y, integ, a1 + da, rate + dr, 0) + \
                _marker_corr(sched, t, y, integ, a1 + da, rate + dr, 1)
            fine.append((c, a1 + da, rate + dr))
    _, a, r = max(fine)
    return dict(a=float(a), rate=float(r),
                corr_first=_marker_corr(sched, t, y, integ, a, r, 0),
                corr_second=_marker_corr(sched, t, y, integ, a, r, 1))


# ── per-step values ─────────────────────────────────────────────────────────

def step_values(sched, rec, y, clock, edge):
    """Median band signal per schedule segment, from spectra wholly inside its measured part:
    the schedule's settle time after the start (markers excepted), `edge` before the end."""
    t, integ = rec["t"], rec["integ"]
    lo = (t - integ - clock["a"]) / clock["rate"]
    hi = (t - clock["a"]) / clock["rate"]
    rows = []
    for s, e, code, kind in sched.segments:
        start = s + (edge if kind == "marker" else max(edge, sched.settle))
        m = (lo >= start) & (hi <= e - edge)
        v = y[m]
        rows.append(dict(start=s, end=e, code=code, kind=kind, n=int(m.sum()),
                         median=float(np.median(v)) if len(v) else np.nan,
                         sd=float(np.std(v)) if len(v) > 1 else np.nan))
    # Idle before the first marker: black, corner square NOT toggling.
    m = (lo >= -10.0) & (hi <= -edge)
    idle = dict(n=int(m.sum()), median=float(np.median(y[m])) if m.any() else np.nan)
    return rows, idle


def drift_correct(rows):
    refs = [(0.5 * (r["start"] + r["end"]), r["median"]) for r in rows
            if r["kind"] == "reference" and np.isfinite(r["median"])]
    if len(refs) < 2:
        raise ValueError("fewer than two usable reference steps; cannot correct drift")
    tr, vr = map(np.array, zip(*refs))
    level = vr.mean()
    for r in rows:
        f = np.interp(0.5 * (r["start"] + r["end"]), tr, vr) / level
        r["value"] = r["median"] / f
        r["sd_value"] = r["sd"] / f
    return dict(times=tr.tolist(), values=(vr / level).tolist(),
                span_percent=float(100 * (vr.max() - vr.min()) / level))


def per_code(rows):
    by = {}
    for r in rows:
        if r["kind"] == "test" and np.isfinite(r["value"]):
            by[r["code"]] = (r["value"], r["sd_value"], r["n"])
    return by


def fit_line(by, codes=None):
    c = np.array(sorted(by if codes is None else [k for k in codes if k in by]))
    v = np.array([by[k][0] for k in c])
    A = np.stack([np.ones_like(c, dtype=float), c.astype(float)], axis=1)
    (b, k), *_ = np.linalg.lstsq(A, v, rcond=None)
    resid_codes = (v - (b + k * c)) / k
    return dict(b=float(b), k=float(k), black_codes=float(b / k),
                codes=c.tolist(), resid_codes=resid_codes.tolist())


# ── reports ─────────────────────────────────────────────────────────────────

def pair_report(I):
    out = []
    if 127 not in I:
        return out
    mid = I[127]
    for lo, hi in PAIRS:
        if lo in I and hi in I:
            out.append(dict(pair=[lo, hi],
                            mean_vs_127_percent=float(100 * ((I[lo] + I[hi]) / 2 / mid - 1)),
                            weber_up=float(I[hi] / mid - 1),
                            weber_down=float(1 - I[lo] / mid)))
    return out


def log_band_report(I):
    try:
        from labpack.visual_stim.clandinin.log_band_ternary import LogBandTernaryStimulus
    except Exception as e:
        return dict(error=str(e))
    s = LogBandTernaryStimulus()
    shown = [s.display_codes[c + s.span] for c in s.attainable_composites]
    if not all(c in I for c in shown):
        return dict(error=f"not every ladder code was measured: {shown}")
    measured = [I.get(c, np.nan) for c in s.display_codes]
    spans_meas = s.probe_log_spans(codes=measured)
    spans_nom = s.probe_log_spans()
    w0 = s.band_weights[0]
    return dict(codes_shown=shown,
                probe_log_span_nominal=spans_nom, probe_log_span_measured=spans_meas,
                mismatch_nominal_percent=float(100 * (max(spans_nom) - min(spans_nom))
                                               / np.mean(spans_nom)),
                mismatch_measured_percent=float(100 * (max(spans_meas) - min(spans_meas))
                                                / np.mean(spans_meas)),
                mean_step_nominal=s.mean_step,
                mean_step_measured=float(measured[w0 + s.span] / measured[-w0 + s.span]))


# ── figures ─────────────────────────────────────────────────────────────────

def figures(out, runs, fit_full, combined, band_info):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for name, run in runs.items():
        sched, rec, y, clock, rows = (run[k] for k in ("sched", "rec", "y", "clock", "rows"))
        fig, ax = plt.subplots(3, 1, figsize=(11, 8))
        tau = (rec["t"] - clock["a"]) / clock["rate"]
        ax[0].plot(tau, y, lw=0.5, color="k")
        for r in rows:
            if np.isfinite(r.get("value", np.nan)):
                ax[0].hlines(r["median"], r["start"], r["end"],
                             color="C1" if r["kind"] == "reference" else "C0", lw=1)
        ax[0].set(title=f"{name}: band signal on the stimulus clock (lines = step medians)",
                  xlabel="stimulus time (s)", ylabel="counts in band")
        for j, which in enumerate((0, 1)):
            m0 = sched.marker_starts[which]
            w = (tau > m0 - 2) & (tau < m0 + sched.marker_duration + 2)
            ax[1 + j].plot(tau[w], y[w], ".", ms=3, color="k")
            tmpl = _box_template(sched, rec["t"][w], rec["integ"], clock["a"], clock["rate"],
                                 which)
            lo_v, hi_v = np.percentile(y[w], [5, 95])
            ax[1 + j].plot(tau[w], lo_v + tmpl * (hi_v - lo_v), color="C3", lw=1)
            ax[1 + j].set(title=f"marker {which + 1}: r = "
                                f"{clock['corr_first' if which == 0 else 'corr_second']:.3f}",
                          xlabel="stimulus time (s)")
        fig.tight_layout()
        fig.savefig(out / f"{name}_alignment.png", dpi=120)
        plt.close(fig)

    codes = np.array(sorted(combined))
    rel = np.array([combined[c]["rel"] for c in codes])
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].plot(codes, rel, ".", ms=3)
    b, k = fit_full["b"], fit_full["k"]
    scale = b + k * 255
    ax[0, 0].plot(codes, (b + k * codes) / scale, color="C3", lw=0.8)
    ax[0, 0].set(xlabel="code", ylabel="intensity (code 255 = 1)",
                 title=(f"line intercept {fit_full['black_codes']:.2f} codes; code 0 "
                        f"{combined[0]['value'] / k:.2f} codes" if 0 in combined else
                        f"line intercept {fit_full['black_codes']:.2f} codes"))
    ax[0, 1].plot(fit_full["codes"], fit_full["resid_codes"], ".", ms=3)
    ax[0, 1].axhline(0, color="k", lw=0.5)
    ax[0, 1].set(xlabel="code", ylabel="residual (codes)", title="full pass: off the line")
    low = codes <= 48
    ax[1, 0].plot(codes[low], rel[low], ".", ms=4)
    ax[1, 0].plot(codes[low], (b + k * codes[low]) / scale, color="C3", lw=0.8)
    ax[1, 0].axhline(b / scale, color="C2", lw=0.8, ls="--", label="line intercept")
    ax[1, 0].set(xlabel="code", ylabel="intensity (code 255 = 1)", title="low codes")
    ax[1, 0].legend()
    wl, bright = band_info["wl"], band_info["bright"]
    ax[1, 1].semilogy(wl, np.clip(bright, 1e-3 * bright.max(), None))
    ax[1, 1].axvspan(*band_info["band"], color="C0", alpha=0.2)
    ax[1, 1].set(xlabel="wavelength (nm)", ylabel="counts", title="bright spectrum and band")
    fig.tight_layout()
    fig.savefig(out / "linearity.png", dpi=120)
    plt.close(fig)


# ── main ────────────────────────────────────────────────────────────────────

def parse_overrides(text):
    """'hold=2,settle=0' -> {'hold': 2.0, 'settle': 0.0}, for re-analysing a recording made
    with an older schedule."""
    out = {}
    for item in filter(None, (text or "").split(",")):
        k, v = item.split("=")
        out[k.strip()] = float(v)
    return out


def run_pass(path, preset, seed, band, edge, overrides=None):
    rec = load(path)
    got = rec["meta"].get("preset")
    if got and got != preset:
        print(f"note: {path} was recorded as preset {got!r}, analysed as {preset!r}")
    sched = PhotometrySweepSchedule(preset=preset, seed=seed, fps=rec["meta"].get("fps", 120.0),
                                    **(overrides or {}))
    lo, hi, bright = find_band(rec)
    band = tuple(band) if band else (lo, hi)
    y = band_signal(rec, band)
    clock = align(sched, rec, y)
    rows, idle = step_values(sched, rec, y, clock, edge)
    drift = drift_correct(rows)
    return dict(rec=rec, sched=sched, y=y, clock=clock, rows=rows, idle=idle, drift=drift,
                band=band, bright=bright)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[1])
    ap.add_argument("--full", required=True, help="recording of the 'full' preset (.npz)")
    ap.add_argument("--low", help="recording of the 'low' preset (.npz), optional")
    ap.add_argument("--seed", type=int, default=0, help="the protocol's seed (default 0)")
    ap.add_argument("--band", type=float, nargs=2, metavar=("NM_LO", "NM_HI"),
                    help="wavelength band to sum; default: found from the data")
    ap.add_argument("--edge", type=float, default=0.15,
                    help="seconds kept clear of each step's end (its start is set by the "
                         "schedule's settle time)")
    ap.add_argument("--full-schedule", default="",
                    help="schedule overrides for an older recording, e.g. 'hold=2,settle=0'")
    ap.add_argument("--low-schedule", default="", help="same, for the low pass")
    ap.add_argument("--out", help="output folder (default: <full>_analysis next to it)")
    args = ap.parse_args()

    full_path = pathlib.Path(args.full)
    out = pathlib.Path(args.out) if args.out else full_path.with_name(full_path.stem + "_analysis")
    out.mkdir(parents=True, exist_ok=True)

    runs = {"full": run_pass(full_path, "full", args.seed, args.band, args.edge,
                             parse_overrides(args.full_schedule))}
    if args.low:
        runs["low"] = run_pass(pathlib.Path(args.low), "low", args.seed,
                               args.band or runs["full"]["band"], args.edge,
                               parse_overrides(args.low_schedule))

    F = runs["full"]
    by_full = per_code(F["rows"])
    fit_full = fit_line(by_full)
    summary = dict(
        full=dict(path=F["rec"]["path"], meta=F["rec"]["meta"], band_nm=F["band"],
                  clock=F["clock"], drift=F["drift"], n_spectra=len(F["y"]),
                  black_level_codes=fit_full["black_codes"],
                  code0_codes=float(by_full[0][0] / fit_full["k"]) if 0 in by_full else None,
                  black_fraction_of_255=fit_full["b"] / (fit_full["b"] + 255 * fit_full["k"]),
                  max_abs_resid_codes=float(np.max(np.abs(fit_full["resid_codes"]))),
                  max_abs_resid_codes_below_32=float(np.max(np.abs(
                      [r for c, r in zip(fit_full["codes"], fit_full["resid_codes"])
                       if c < 32]))),
                  steps_without_samples=[r["code"] for r in F["rows"]
                                         if r["kind"] == "test" and r["n"] == 0]))

    zero = by_full.get(0, (np.nan,))[0]
    idle = F["idle"]["median"] / (np.interp(0, F["drift"]["times"], F["drift"]["values"]))
    summary["corner_square"] = dict(
        idle_black_counts=float(idle), sweep_code0_counts=float(zero),
        difference_codes=float((zero - idle) / fit_full["k"]),
        note="positive = the toggling corner square adds light at the fiber")

    combined = {c: dict(value=v, source="full") for c, (v, _, _) in by_full.items()}
    if "low" in runs:
        L = runs["low"]
        by_low = per_code(L["rows"])
        fit_low = fit_line(by_low)
        overlap = [c for c in by_low if c in by_full and c >= 8]
        s = (np.sum([by_low[c][0] * by_full[c][0] for c in overlap])
             / np.sum([by_low[c][0] ** 2 for c in overlap]))
        for c, (v, _, _) in by_low.items():
            combined[c] = dict(value=s * v, source="low")
        summary["low"] = dict(path=L["rec"]["path"], meta=L["rec"]["meta"], clock=L["clock"],
                              drift=L["drift"], black_level_codes=fit_low["black_codes"],
                              code0_codes=(float(by_low[0][0] / fit_low["k"]) if 0 in by_low
                                           else None),
                              max_abs_resid_codes=float(np.max(np.abs(fit_low["resid_codes"]))),
                              scale_onto_full=float(s), overlap_codes=overlap)
    top = combined[255]["value"] if 255 in combined else max(v["value"] for v in combined.values())
    for c in combined:
        combined[c]["rel"] = combined[c]["value"] / top
    I = {c: v["rel"] for c, v in combined.items()}
    summary["symmetric_pairs"] = pair_report(I)
    summary["log_band_default_ladder"] = log_band_report(I)

    with open(out / "summary.json", "w") as f:
        json.dump(summary, f, indent=2, default=float)
    with open(out / "lut.csv", "w") as f:
        f.write("code,relative_intensity,source\n")
        for c in sorted(combined):
            f.write(f"{c},{combined[c]['rel']:.6g},{combined[c]['source']}\n")
    figures(out, runs, fit_full, combined,
            dict(wl=F["rec"]["wl"], bright=F["bright"], band=F["band"]))

    s = summary["full"]
    print(f"band {F['band'][0]:.1f}-{F['band'][1]:.1f} nm; markers r = "
          f"{F['clock']['corr_first']:.3f} / {F['clock']['corr_second']:.3f}; clock rate "
          f"{1e6 * (F['clock']['rate'] - 1):+.1f} ppm; drift {F['drift']['span_percent']:.2f}%")
    print(f"black level: line intercept {s['black_level_codes']:.2f} codes "
          f"({100 * s['black_fraction_of_255']:.3f}% of code 255), code 0 itself "
          f"{s['code0_codes']:.2f} codes"
          + (f"; low pass: intercept {summary['low']['black_level_codes']:.2f}, code 0 "
             f"{summary['low']['code0_codes']:.2f} codes" if "low" in summary else ""))
    print(f"worst residual {s['max_abs_resid_codes']:.2f} codes overall, "
          f"{s['max_abs_resid_codes_below_32']:.2f} below code 32")
    print(f"corner square: {summary['corner_square']['difference_codes']:+.2f} codes")
    for p in summary["symmetric_pairs"]:
        print(f"pair {p['pair'][0]:3d}/{p['pair'][1]:3d}: mean vs 127 "
              f"{p['mean_vs_127_percent']:+.2f}%")
    lb = summary["log_band_default_ladder"]
    if "error" not in lb:
        print(f"LogBand ladder: contrast mismatch {lb['mismatch_nominal_percent']:.2f}% nominal, "
              f"{lb['mismatch_measured_percent']:.2f}% measured; mean step "
              f"{lb['mean_step_nominal']:.2f}x nominal, {lb['mean_step_measured']:.2f}x measured")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
