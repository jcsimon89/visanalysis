"""Model against data: run a model on what a recording showed, sample it where the data were
sampled, and report how well it predicts.

compare(model, segments)  -> dict(rows, per_state, sta_data, sta_model, ...)
    rows: every measured sample with its prediction, state, time since switch, plane.
    per_state: r between prediction and data (all rows, and steady rows only), SD of each, the
    mean level of each; the model's level error is the offset the LN did not get right.
summary_figure(...)       filters, nonlinearities, an example trace (data and prediction), the
                          switch-triggered averages of both, and the numbers.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from visanalysis.model import measurement, readouts
from visanalysis.model.models import SwitchedLN, Response
from visanalysis.model.plots import state_colors

__all__ = ["compare", "summary_figure"]


def _predict_segment(model, seg):
    r = model.respond(seg["stim"])
    times = [p["t"] for p in seg["planes"]]
    sampled = measurement.sample(r, times=times)
    return r, sampled


def compare(model, segments, steady_start=2.0, label="contrast_state", sta_pre=2.0, sta_post=None, sta_bin=0.25):
    rows = dict(t=[], y=[], pred=[], state=[], since=[], plane=[], seg=[], w=[])
    data_resp, model_resp, fine = [], [], []
    for k, seg in enumerate(segments):
        r, sampled = _predict_segment(model, seg)
        fine.append(r)
        for p, s in zip(seg["planes"], sampled):
            # sample() drops times outside the stimulus; align data to the kept times
            m = (p["t"] >= r.t[0]) & (p["t"] <= r.t[-1])
            y = np.asarray(p["y"], float)[m]
            rows["t"].append(s.t)
            rows["y"].append(y)
            rows["pred"].append(s.y)
            rows["state"].append(np.asarray(s.labels[label]))
            rows["since"].append(np.asarray(s.labels.get("since_switch", np.full(s.t.size, np.nan)), float))
            rows["plane"].append(np.full(s.t.size, p["plane"]))
            rows["seg"].append(np.full(s.t.size, k))
            rows["w"].append(np.full(s.t.size, float(p["n_roi"])))      # ROIs behind each plane mean
            data_resp.append(Response(s.t, y, labels=s.labels, name="data", plane=p["plane"], valid_from=s.valid_from))
            model_resp.append(s)
    rows = {k: np.concatenate(v) for k, v in rows.items()}
    ok = np.isfinite(rows["y"]) & np.isfinite(rows["pred"])
    per_state = {}
    for v in np.unique(rows["state"]):
        if v == 0:
            continue
        m = ok & (rows["state"] == v)
        with np.errstate(invalid="ignore"):
            ss = m & (np.nan_to_num(rows["since"], nan=np.inf) >= steady_start)
        d = dict(n=int(m.sum()), n_steady=int(ss.sum()))
        for nm, mm in (("all", m), ("steady", ss)):
            if mm.sum() > 2:
                d["r_" + nm] = float(np.corrcoef(rows["y"][mm], rows["pred"][mm])[0, 1])
                # slope of data on prediction: 1 when the model's gain is right whatever the noise
                # (a least-squares model predicts the conditional mean, so its SD is r x data SD)
                d["slope_" + nm] = float(np.polyfit(rows["pred"][mm], rows["y"][mm], 1)[0])
                d["sd_data_" + nm] = float(rows["y"][mm].std())
                d["sd_model_" + nm] = float(rows["pred"][mm].std())
                d["mean_data_" + nm] = float(rows["y"][mm].mean())
                d["mean_model_" + nm] = float(rows["pred"][mm].mean())
                d["rmse_" + nm] = float(np.sqrt(np.mean((rows["y"][mm] - rows["pred"][mm]) ** 2)))
        per_state[int(v)] = d
    r_all = float(np.corrcoef(rows["y"][ok], rows["pred"][ok])[0, 1]) if ok.sum() > 2 else np.nan
    frozen = _frozen_average(rows, segments, steady_start, ok)
    post = sta_post or min(10.0, float(np.median([s.get("dwell_s", 10.0) for s in segments])))
    sta_d = readouts.switch_triggered_average(data_resp, label, sta_pre, post, sta_bin)
    sta_m = readouts.switch_triggered_average(model_resp, label, sta_pre, post, sta_bin)
    return dict(rows=rows, per_state=per_state, r_all=r_all, sta_data=sta_d, sta_model=sta_m, fine=fine,
                steady_start=steady_start, label=label, out_of_range=fine[0].out_of_range if fine else {},
                frozen_average=frozen)


def _frozen_average(rows, segments, steady_start, ok, min_frac=0.5):
    """Frozen noise: every segment showed the same stimulus, so data can be averaged over segments at
    the same STIMULUS POSITION (time since the segment's onset, binned at the sample interval; the
    imaging clock is not locked to the stimulus, so sample indices do not line up across segments).
    The model is compared with that average, and odd/even segment halves give a noise ceiling: the
    correlation between the two half-averages, Spearman-Brown corrected to the full count, is the
    reliability of the full average, and its square root is the r a perfect model could reach."""
    seeds = [s.get("seed", s["stim"].meta.get("seed")) for s in segments]
    if len(segments) < 2 or not all(sd is not None and sd == seeds[0] for sd in seeds):
        return {}
    n_seg = len(segments)
    t0 = np.array([s.get("t0", s["stim"].t[0]) for s in segments])
    pos = rows["t"] - t0[rows["seg"]]
    # the bin is one plane's sample interval (the volume period), so each (plane, bin) gets about one
    # sample per segment; taking it across planes would give a bin a fraction of the volume period
    s0 = rows["seg"] == 0
    p0 = rows["plane"][s0][0]
    dt = float(np.median(np.diff(np.sort(rows["t"][s0 & (rows["plane"] == p0)]))))
    bins = np.floor(pos / dt).astype(int)
    key = rows["plane"].astype(np.int64) * 1_000_000 + bins
    u, inv = np.unique(key, return_inverse=True)
    def avg(sel):
        w = (ok & sel).astype(float)
        cnt = np.bincount(inv, w, u.size)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.bincount(inv, rows["y"] * w, u.size) / cnt, cnt
    ya, cnt = avg(np.ones(ok.size, bool))
    yo, co = avg(rows["seg"] % 2 == 0)
    ye, ce = avg(rows["seg"] % 2 == 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pa = np.bincount(inv, np.where(ok, rows["pred"], 0.0), u.size) / np.maximum(cnt, 1)
    # labels of a bin from the first segment holding it
    st_b = np.zeros(u.size, int)
    sn_b = np.full(u.size, np.nan)
    first = np.full(u.size, -1)
    for i in range(inv.size):
        if first[inv[i]] < 0 and ok[i]:
            first[inv[i]] = i
    good = first >= 0
    st_b[good] = rows["state"][first[good]]
    sn_b[good] = rows["since"][first[good]]
    enough = good & (cnt >= min_frac * n_seg) & (co > 0) & (ce > 0) & np.isfinite(ya) & np.isfinite(pa)
    frozen = {}
    for v in np.unique(st_b[enough]):
        if v == 0:
            continue
        with np.errstate(invalid="ignore"):
            mm = enough & (st_b == v) & (np.nan_to_num(sn_b, nan=np.inf) >= steady_start)
        if mm.sum() > 2 and pa[mm].std() > 0:
            r_oe = float(np.corrcoef(yo[mm], ye[mm])[0, 1])
            rel_half = max(r_oe, 0.0)
            rel_full = 2 * rel_half / (1 + rel_half) if rel_half > 0 else 0.0
            r_model = float(np.corrcoef(ya[mm], pa[mm])[0, 1])
            frozen[int(v)] = dict(r=r_model, slope=float(np.polyfit(pa[mm], ya[mm], 1)[0]),
                                  sd_data=float(ya[mm].std()), sd_model=float(pa[mm].std()), n=int(mm.sum()),
                                  n_trials=n_seg, r_odd_even=r_oe, ceiling=float(np.sqrt(rel_full)),
                                  fraction_of_ceiling=float(r_model / np.sqrt(rel_full)) if rel_full > 0 else np.nan,
                                  bin_s=dt)
    return frozen


def _lns(model):
    return dict(model.states) if isinstance(model, SwitchedLN) else {"all": model}


def _pooled_trace(rows, segments, seg_index, bin_s, plane=None):
    """Data and prediction pooled over all planes (and over segments when the noise was frozen),
    binned by stimulus position. Returns (t on segment seg_index's timeline, data, model, n_seg,
    n_roi) -- the planes sit at different phases of the volume, so pooling them samples the
    response finer than the volume rate."""
    seeds = [s.get("seed", s["stim"].meta.get("seed")) for s in segments]
    frozen = len(segments) > 1 and all(sd is not None and sd == seeds[0] for sd in seeds)
    t0 = np.array([s.get("t0", s["stim"].t[0]) for s in segments])
    use = np.ones(rows["t"].size, bool) if frozen else rows["seg"] == seg_index
    if plane is not None:
        use &= rows["plane"] == plane
    pos = rows["t"][use] - t0[rows["seg"][use]]
    y, pr, w = rows["y"][use], rows["pred"][use], rows["w"][use]
    ok = np.isfinite(y) & np.isfinite(pr)
    w = np.where(ok, w, 0.0)                     # a plane mean counts for the ROIs behind it
    b = np.floor(pos / bin_s).astype(int)
    u, inv = np.unique(b, return_inverse=True)
    cnt = np.bincount(inv, w, u.size)
    n_in = np.bincount(inv, ok.astype(float), u.size)
    with np.errstate(invalid="ignore", divide="ignore"):
        yb = np.bincount(inv, np.where(ok, y, 0.0) * w, u.size) / cnt
        pb = np.bincount(inv, np.where(ok, pr, 0.0) * w, u.size) / cnt
    # the imaging clock is not locked to the stimulus, so samples cluster unevenly over position
    # bins; a bin holding a few samples is noise, not data -- leave it out
    thin = n_in < 0.4 * np.median(n_in[n_in > 0])
    yb[thin] = np.nan
    pb[thin] = np.nan
    t = t0[seg_index] + (u + 0.5) * bin_s
    n_roi = sum(p["n_roi"] for p in segments[seg_index]["planes"] if plane is None or p["plane"] == plane)
    return t, yb, pb, (len(segments) if frozen else 1), n_roi


def summary_figure(model, segments, res, title=None, window_s=16.0, plane=None, seg_index=None, trace_bin_s=None):
    """plane=None pools every plane (and, for frozen noise, every segment) by stimulus position into
    bins of trace_bin_s (default: a quarter of the volume period, at least 20 ms); a plane number shows
    that plane alone at its own sample times."""
    lns = _lns(model)
    colors = state_colors(list(lns))
    label = res["label"]
    fig = plt.figure(figsize=(16, 8.5))
    gs = gridspec.GridSpec(2, 3, height_ratios=(1, 1.1), figure=fig)
    fig.suptitle(title or "{} vs data".format(model.name), fontsize=11)

    a = fig.add_subplot(gs[0, 0])
    for v, ln in lns.items():
        lag, d = ln.filter.tabulate(0.005)
        a.plot(lag, d, color=colors[v], lw=1.6, label="{} ({})".format(ln.name or v, ln.filter.kind))
    a.axhline(0, color="0.7", lw=0.5)
    a.set_xlim(0, 0.8)
    a.set_xlabel("lag (s)")
    a.set_ylabel("filter density ({} / contrast / s)".format(next(iter(lns.values())).filter.output_units))
    a.set_title("filters", fontsize=10)
    a.legend(fontsize=8)

    a = fig.add_subplot(gs[0, 1])
    for v, ln in lns.items():
        nl = ln.nonlinearity
        rg = nl.range or (-1, 1)
        x = np.linspace(rg[0], rg[1], 200)
        a.plot(x, nl(x), color=colors[v], lw=1.6, label="{} ({})".format(ln.name or v, nl.kind))
    a.set_xlabel("generator (contrast)")
    a.set_ylabel("response (dF/F)")
    a.set_title("nonlinearities (measured range)", fontsize=10)
    a.legend(fontsize=8)

    a = fig.add_subplot(gs[0, 2])
    for (f, to), d in sorted(res["sta_data"].items()):
        c = colors.get(to, "k")
        a.plot(d["lag"], d["mean"], color=c, lw=1.6, label="data {} -> {} (n={})".format(f, to, d["n"]))
        a.fill_between(d["lag"], d["mean"] - d["se"], d["mean"] + d["se"], color=c, alpha=0.2, lw=0)
    for (f, to), d in sorted(res["sta_model"].items()):
        a.plot(d["lag"], d["mean"], "--", color=colors.get(to, "k"), lw=1.2, label="model {} -> {}".format(f, to))
    a.axvline(0, color="0.6", lw=0.6)
    a.set_xlabel("time from switch (s)")
    a.set_ylabel("response (dF/F)")
    a.set_title("switch-triggered average", fontsize=10)
    a.legend(fontsize=7)

    # example trace: the plane with most ROIs, in the segment with most finite data, around a switch
    rows = res["rows"]
    if seg_index is None:
        segs, counts = np.unique(rows["seg"][np.isfinite(rows["y"])], return_counts=True)
        seg_index = int(segs[np.argmax(counts)]) if segs.size else 0
    seg = segments[seg_index]
    one_plane = (rows["seg"] == seg_index) & (rows["plane"] == (plane if plane is not None else rows["plane"][0]))
    vol = float(np.median(np.diff(np.sort(rows["t"][one_plane])))) if one_plane.sum() > 1 else 0.1
    if trace_bin_s is None:
        if plane is not None:
            trace_bin_s = vol
        else:
            # the finest bin that the sampling fills evenly: samples recur at similar stimulus positions
            # across repeats, so a bin finer than that spacing comes out alternately full and empty
            for cand in (max(vol / 4.0, 0.02), max(vol / 2.0, 0.02), vol):
                trace_bin_s = cand
                yy = _pooled_trace(rows, segments, seg_index, cand, plane)[1]
                if np.mean(~np.isfinite(yy)) < 0.05:
                    break
    t, y, pr, n_avg, n_roi = _pooled_trace(rows, segments, seg_index, trace_bin_s, plane)
    tr = seg["stim"].transitions(label)
    tr = [x for x in tr if 0 not in (x[1], x[2])]
    span = (t[np.isfinite(y)].min(), t[np.isfinite(y)].max()) if np.any(np.isfinite(y)) else (t[0], t[-1])
    window_s = min(window_s, span[1] - span[0])
    if len(tr) >= 2:
        tc = seg["stim"].t[tr[1][0]]
        w0 = tc - window_s * 0.3
    else:
        w0 = span[0]
    w0 = min(max(w0, span[0]), span[1] - window_s)
    w1 = w0 + window_s
    a = fig.add_subplot(gs[1, :2])
    ww = (t >= w0) & (t < w1)
    a.plot(t[ww], y[ww], color="0.25", lw=1.0, label="data: {} ROIs{}, {} ms bins{}".format(
        n_roi, " in plane {}".format(plane) if plane is not None else " in all {} planes".format(len(seg["planes"])),
        int(round(trace_bin_s * 1e3)), ", averaged over {} frozen-noise repeats".format(n_avg) if n_avg > 1 else ""))
    a.plot(t[ww], pr[ww], color="#ff7f0e", lw=1.3, label="model")
    sl = seg["stim"].labels[label]
    sw = (seg["stim"].t >= w0) & (seg["stim"].t < w1)
    lab, tt = sl[sw], seg["stim"].t[sw]
    if tt.size:
        ch = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1, tt.size]
        for i0, i1 in zip(ch[:-1], ch[1:]):
            v = int(lab[i0])
            if v in colors:
                a.axvspan(tt[i0], tt[min(i1, tt.size - 1)], color=colors[v], alpha=0.08, lw=0)
    rw = float(np.corrcoef(y[ww][np.isfinite(y[ww])], pr[ww][np.isfinite(y[ww])])[0, 1]) if np.isfinite(y[ww]).sum() > 2 else np.nan
    a.set_xlim(w0, w1)
    a.set_xlabel("series time (s)")
    a.set_ylabel("dF/F")
    a.set_title("example: {} (r in window {:.2f})".format(seg["stim"].name, rw), fontsize=10)
    a.legend(fontsize=8, loc="upper right")

    a = fig.add_subplot(gs[1, 2])
    a.axis("off")
    lines = ["prediction vs data, all planes, {} rows".format(int(np.isfinite(rows["y"]).sum())),
             "overall r {:.3f}".format(res["r_all"]), ""]
    for v, d in res["per_state"].items():
        nm = lns[v].name if v in lns else str(v)
        lines.append("{} (state {}):".format(nm, v))
        for k in ("all", "steady"):
            if "r_" + k in d:
                lines.append("  {:6s} r {:.3f}  slope {:.2f}  rmse {:.4f}".format(k, d["r_" + k], d["slope_" + k], d["rmse_" + k]))
                lines.append("         sd data {:.4f} model {:.4f}; mean data {:+.4f} model {:+.4f}".format(
                    d["sd_data_" + k], d["sd_model_" + k], d["mean_data_" + k], d["mean_model_" + k]))
    if res.get("frozen_average"):
        lines.append("")
        lines.append("data averaged over {} frozen-noise trials, steady rows:".format(
            next(iter(res["frozen_average"].values()))["n_trials"]))
        for v, d in res["frozen_average"].items():
            lines.append("  state {:+d}: r {:.3f}  slope {:.2f}  ceiling {:.3f} (odd/even r {:.3f}) -> {:.0%} of ceiling".format(
                v, d["r"], d["slope"], d["ceiling"], d["r_odd_even"], d["fraction_of_ceiling"]))
    lines.append("")
    lines.append("slope = data regressed on prediction (1 = gain right)")
    lines.append("steady = >= {:.1f} s after a switch".format(res["steady_start"]))
    for k, f in res["out_of_range"].items():
        if f > 0.02:
            lines.append("{}: {:.0%} of samples beyond the NL's measured range".format(k, f))
    a.text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=8, family="monospace", transform=a.transAxes)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig
