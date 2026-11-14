"""Contrast adaptation from ContrastSwitchTernaryNoise (labpack JCS_protocol).

WHAT THE STIMULUS IS. Full-field ternary noise (the PROBE) whose amplitude
switches between a low and a high contrast every `modulator_dwell` seconds at a
fixed mean ('fixed' schedule). One epoch is one trial of many dwells. The
question is what a switch does to the response: its gain, its kinetics, and how
fast each one moves; and, for a second channel such as cAMP, what the contrast
envelope itself does.

TIMING IS THE SAME AS THE OTHER NOISE PIPELINES. Stimulus time 0 is the
photodiode's first edge in the epoch (getStimulusTiming), the imaging sample of
an ROI is at time_vector_by_epoch - pre_time + roi_z_offset, frame k occupies
[k / fps, (k + 1) / fps) of stimulus time, and a dropped display frame holds
the previous frame (reconstruct_displayed_frame_indices). Two known offsets
come with that convention and are shared with every other noise analysis:
stimpack paints the corner square black on the stimulus's first frame, so the
first edge is the SECOND frame, and the display picks its frame as
round(stim_time * fps) from the stimulus PC's clock at paint time, which runs
~20 ppm fast against the acquisition clock. Together they put the modelled
stimulus about one frame (~8 ms) late at the start of a trial and ~14 ms late
at its end, so filter latencies read short by about that much.

RESPONSES. ROI traces come from bleach.roi_responses (multiplicative bleach
correction within each epoch), and dF/F is taken against each ROI's own median
in that epoch, as in fullfield_temporal.

EVERY ROI ON ITS OWN, THEN COMBINED. Each gated ROI gets its own filter per
contrast, nonlinearity, gains, non-adapting reference, gain dynamics and
stretch; nothing is averaged across ROIs before fitting. Population numbers are
medians across ROIs with 95% CIs from resampling ROIs, and every per-ROI value
is saved alongside. Per-ROI low-contrast filters are noisy (on Tm9 most do not
beat their null), which is exactly why each ROI's ratio is read against its own
non-adapting reference, which carries that ROI's filter noise. ROIs of one z
plane share a design, since they share their sample times; planes are sampled
at different times, which fills in the volume period, so the lag axis can be
finer than the volume (lag_dt). Each tap j is the displayed contrast averaged
over [t - (j + 1) lag_dt, t - j lag_dt], with t the middle of the plane's scan;
its lag is (j + 0.5) lag_dt.

LN CONVENTION: BACCUS & MEISTER (2002), as for contrast_blocks. Each contrast
state gets its OWN filter, scaled so the filtered stimulus has the variance of
the stimulus in that state, var(F * s) = var(s). The filter then carries only
shape and timing, and every gain change is in the nonlinearity. gain = slope of
the response on that generator. shared_range = the same slope over the
generator range BOTH states cover (their Fig. 9C); a static nonlinearity gives
equal slopes there whatever its shape, so a ratio away from 1 there is
adaptation, while a whole-range ratio away from 1 can be curvature alone.
Generators are out of sample: each row goes through filters fitted without its
cross-validation fold.

READ A RATIO AGAINST A NON-ADAPTING CELL, NOT AGAINST 1. A noisy low-contrast
filter biases the low gain down, and with frozen noise the low and high dwells
are different finite samples of the probe whose slow power need not match, so a
cell with no adaptation at all can read well away from 1 (1.22 on
JS256_x_JS261/fly_012). analyse_channel simulates that cell -- the high filter
in both states, at the channel's own noise -- runs it through the same
estimator, and reports every ratio divided by what it recovers
(ratio_vs_static). The model-free ratio carries the same stimulus imbalance and
is read against the exact no-adaptation value (static_reference.truth).

FROZEN NOISE. If every epoch has the same seed, the noise realisation, not the
trial, is the independent unit: confidence intervals bootstrap over dwells (the
same dwell in every trial is one realisation), and cross-validation folds hold
out the same stimulus position in every trial. The frozen seed also allows a
MODEL-FREE check: the covariance between trials at the same stimulus position
is the stimulus-locked variance alone, so its SD ratio between states compared
with the contrast ratio gives the gain ratio without any filter. With fresh
seeds, the unit is (trial, dwell) and the model-free check is skipped.

WHAT IS FITTED WHERE.
  * Rows: dwells >= first_dwell (default 2: the first dwell follows the idle
    grey, not the other contrast, and the switch out of it is the first real
    one), stimulus time within the trial.
  * Steady-state filters: rows whose WHOLE lag window lies inside the current
    state, i.e. >= steady_start s after the switch, which defaults to the
    filter length. A shorter start lets the filter's tail reach into the other
    contrast: a low-contrast row then carries high-contrast history, which
    inflates the low state's generator and its Baccus & Meister scale (on
    JS256_x_JS261/fly_012, starting at 2 s with a 3 s filter put the
    no-adaptation truth at 1.22 instead of ~1).
  * Gain after a switch: rows binned by time since the switch; in each bin the
    response is regressed on two generators, the current state's filter over
    the stimulus since the switch and the previous state's filter over the
    stimulus before it, so a bin early in the dwell is not penalised for
    history from the other contrast.
  * Null: the stimulus rolled within each trial by whole cycles (k * 2 dwell),
    which keeps the contrast state on every row and destroys only the pairing
    with the probe.

ROI GATE. Default: each ROI's own high-contrast filter must predict held-out
data better than the rolled null (BH q = 0.05). Selection uses lags >= 40 ms
only, so the display light leak at lag 0 cannot carry an ROI through; the leak
reaches every pixel, cell or not. Low contrast plays no part in choosing ROIs.
Because selection and measurement use the same high-contrast data, the gate is
also run on each half of the stimulus and the analysis repeated on the other
half (cross-fit) as a check on selection bias. A keep-list from elsewhere (for
example the series-1 1D STRF test) can be given instead.

THE SECOND CHANNEL. Its LN model is fitted and its held-out r reported against
the same null, so it is on record whether it has a filter at all; gains are
only meaningful for a channel that passes. The switch-triggered average is
reported for every channel.
"""
import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.optimize import curve_fit

from visanalysis.util import bleach, nullcal
from visanalysis.util.frame_timing import reconstruct_displayed_frame_indices
from visanalysis.util.fullfield_temporal import _parse_meta_attr

__all__ = [
    "is_contrast_switch",
    "resolve_meta",
    "load_responses",
    "build_context",
    "gate_rois",
    "analyse_channel",
    "switch_triggered",
    "model_free_gain",
    "simulate_responses",
]

FS = 1000.0                  # Hz, the grid the displayed stimulus is held on
FOLD_S = 60.0                # s of stimulus position per cross-validation fold
GAIN_BINS = np.array([0, .25, .5, 1, 1.5, 2, 3, 4, 6, 8, 10])
RIDGE_GRID = np.logspace(-3, 3, 13)
STA_BIN = 0.25
DARK_FRAC = 0.2              # an epoch whose population level falls below this fraction is dropped
FADE_FRAC = 0.1              # an ROI's epoch whose baseline is below this fraction of the channel's is dropped


# ── stimulus ────────────────────────────────────────────────────────────────

def is_contrast_switch(run_parameters):
    return "ContrastSwitch" in str(run_parameters.get("protocol_ID", ""))


def resolve_meta(epoch_params):
    meta = epoch_params.get("contrast_switch_stimulus")
    if meta is None:
        raise KeyError("epoch is missing contrast_switch_stimulus; the stimulus cannot be "
                       "reconstructed without it")
    return _parse_meta_attr(meta, "contrast_switch_stimulus")


def _generator(meta):
    from visanalysis.util.octave_stim import _find_labpack
    _find_labpack()
    from labpack.visual_stim.clandinin.contrast_switch_ternary import (
        ContrastSwitchTernaryStimulus)
    g = ContrastSwitchTernaryStimulus.from_metadata(meta)
    if not g.fixed_modulator:
        raise NotImplementedError("only the 'fixed' modulator schedule is analysed; a renewal "
                                  "schedule needs dwell bookkeeping from transitions()")
    return g


def epoch_stimulus(meta, stim_dur, pre_time, tail_time, drop_inds=(), cache=None):
    """Displayed stimulus of one epoch on the FS grid, stimulus time -pre_time .. stim_dur + tail.

    Returns dict(t0, c, cs, state, dwell, since, pos, ...). c is Weber contrast (0 in the
    grey before and after the stimulus), state -1 low / +1 high / 0 grey, dwell the index of
    the contrast dwell (-1 in the grey), since the time since that dwell began, pos the
    stimulus time. stim_dur is the photodiode-measured span of the epoch, so the frame
    count is the one the display actually showed.
    """
    key = (int(meta["seed"]), round(stim_dur, 3), tuple(int(d) for d in drop_inds))
    if cache is not None and key in cache:
        fr = cache[key]
    else:
        g = _generator(meta)
        n_fr = int(round(stim_dur * g.fps)) + 1
        idx = np.arange(n_fr)
        if len(drop_inds):
            actual = reconstruct_displayed_frame_indices(n_fr, np.asarray(drop_inds, int))
            idx = np.where(actual == -1, idx, actual)
        fr = dict(fps=float(g.fps), dw=int(g.modulator_dwell_frames), idx=idx,
                  c=g.contrast_signal(idx), state=g.band_values(0, idx).astype(np.int8),
                  realized=tuple(g.realized_contrasts))
        if cache is not None:
            cache[key] = fr
    t = -float(pre_time) + np.arange(int(np.ceil((pre_time + stim_dur + tail_time) * FS)) + 1) / FS
    k = np.floor(t * fr["fps"] + 1e-9).astype(np.int64)
    on = (k >= 0) & (k < fr["idx"].size)
    kk = np.clip(k, 0, fr["idx"].size - 1)
    c = np.where(on, fr["c"][kk], 0.0)
    state = np.where(on, fr["state"][kk], 0).astype(np.int8)
    dwell = np.where(on, fr["idx"][kk] // fr["dw"], -1)
    since = np.where(on, t - dwell * fr["dw"] / fr["fps"], np.nan)
    return dict(t0=-float(pre_time), c=c, cs=np.r_[0.0, np.cumsum(c)], state=state, dwell=dwell,
                since=since, pos=np.where(on, t, np.nan), dwell_s=fr["dw"] / fr["fps"], fps=fr["fps"],
                realized=fr["realized"])


def rolled(st, sec):
    """The stimulus rolled within the stimulus period by sec (the null); labels unchanged."""
    on = np.flatnonzero(st["dwell"] >= 0)
    c = st["c"].copy()
    c[on] = np.roll(st["c"][on], int(round(sec * FS)))
    return dict(st, c=c, cs=np.r_[0.0, np.cumsum(c)])


def taps(st, t, n_lag, lag_dt):
    """(len(t), n_lag): tap j = mean contrast over [t - (j+1) lag_dt, t - j lag_dt]."""
    edges = np.round((t[:, None] - np.arange(n_lag + 1)[None, :] * lag_dt - st["t0"]) * FS)
    edges = np.clip(edges.astype(np.int64), 0, st["cs"].size - 1)
    hi, lo = edges[:, :-1], edges[:, 1:]
    return (st["cs"][hi] - st["cs"][lo]) / np.maximum(hi - lo, 1)


def labels(st, t):
    i = np.clip(np.round((t - st["t0"]) * FS).astype(np.int64), 0, st["c"].size - 1)
    return st["state"][i], st["dwell"][i], st["since"][i], st["pos"][i]


def _state_of_dwells(ctx):
    st0 = ctx["stims"][0]
    i = np.round((np.arange(ctx["n_dwell"]) * ctx["dwell_s"] + 0.5 * ctx["dwell_s"] - st0["t0"]) * FS)
    return st0["state"][np.clip(i.astype(int), 0, st0["state"].size - 1)]


# ── responses ───────────────────────────────────────────────────────────────

def load_responses(ID, roi_set, correct_bleach="auto"):
    """dF/F per ROI per epoch against the ROI's own median in that epoch, with its timing."""
    rd = bleach.roi_responses(ID, roi_set, roi_prefix="aligned", dff="none",
                              correct_bleach=correct_bleach, quiet=True)
    er = np.asarray(rd["epoch_response"], dtype=float)
    n_roi = er.shape[0]
    off = np.asarray(rd.get("roi_z_offsets", np.zeros(n_roi)), dtype=float)
    try:
        frame_offsets = np.asarray(ID.getVolumeFrameOffsets(), dtype=float)
    except Exception:
        frame_offsets = np.zeros(1)
    plane = np.argmin(np.abs(off[:, None] - frame_offsets[None, :]), axis=1)
    # A channel that goes dark mid-series (fly_012 ch2 fell to ~0 for its last 1.25 trials while
    # ch1 carried on) has a median near zero and a dF/F of ~1e15. An epoch in which the
    # population's 5 s median falls below DARK_FRAC of the channel's typical level is dropped
    # whole: the bleach correction has already smeared the step across ~25 s on either side.
    pop = [np.nanmedian(er[:, e, :], axis=0) for e in range(er.shape[1])]
    typical = np.nanmedian(np.concatenate(pop))
    # The same for single ROIs. Some fade to background through a series -- on
    # JS256_x_JS261/fly_011 ROIs 155, 178 and 179 fell from ~500-2000 to 5-50 counts against a
    # channel median of ~2100 -- and dF/F against a baseline near zero (SD up to 127, against 0.7
    # typical) swamps any average it enters. An ROI's epoch is dropped when its median falls below
    # FADE_FRAC of the channel's median baseline. Relative to the channel, not to the ROI's own
    # level: whole-population brightness moves 20-30% between trials, which dF/F absorbs.
    with np.errstate(invalid="ignore"):
        f0_all = np.nanmedian(er, axis=2)                               # (n_roi, n_ep)
        faded = f0_all < FADE_FRAC * np.nanmedian(f0_all)
    resp, tvecs, dark = [], [], []
    for e, tv in enumerate(rd["time_vector_by_epoch"]):
        tv = np.asarray(tv, dtype=float)[:er.shape[2]]
        x = er[:, e, :tv.size]
        p = pop[e][:tv.size]
        bins = np.floor((tv - tv[0]) / 5.0).astype(int) if tv.size else np.zeros(0, int)
        lvl = np.array([np.nanmedian(p[bins == b]) for b in np.unique(bins)]) if tv.size else np.array([np.nan])
        if not np.all(lvl > DARK_FRAC * typical):
            dark.append(e)
            resp.append(np.full(x.shape, np.nan))
        else:
            with np.errstate(invalid="ignore", divide="ignore"):
                f0 = np.nanmedian(x, axis=1, keepdims=True)
                r_ = x / f0 - 1.0
            r_[faded[:, e]] = np.nan
            resp.append(r_)
        tvecs.append(tv)
    scan = float(np.median(np.diff(frame_offsets))) if frame_offsets.size > 1 else 0.0
    return dict(resp=resp, tvecs=tvecs, off=off, plane=plane, frame_offsets=frame_offsets,
                scan=scan, n_roi=n_roi, bleach_info=rd.get("bleach_info", {}), dark_epochs=dark,
                faded=[(int(r), int(e)) for r, e in zip(*np.nonzero(faded)) if e not in dark])


def build_context(ID, roi_sets, correct_bleach="auto"):
    """Everything an analysis needs that does not depend on which ROIs are kept.

    roi_sets maps a channel name ('calcium') to the ROI set holding it ('mask_ch1')."""
    _, _, sample_rate = ID.getVoltageData()
    timing = ID.getStimulusTiming()
    starts = np.asarray(timing["stimulus_start_times"], dtype=float)
    ends = np.asarray(timing["stimulus_end_times"], dtype=float)
    drops = np.asarray(timing["dropped_frame_times"], dtype=float) / sample_rate
    eps = ID.getEpochParameters()
    resp = {k: load_responses(ID, v, correct_bleach) for k, v in roi_sets.items()}
    n_ep = min(len(starts), len(eps), *[len(r["resp"]) for r in resp.values()])
    pre = np.asarray(ID.getEpochParameters("pre_time"), dtype=float)[:n_ep]
    tail = np.asarray(ID.getEpochParameters("tail_time"), dtype=float)[:n_ep]
    metas = [resolve_meta(eps[i]) for i in range(n_ep)]
    seeds = [int(m["seed"]) for m in metas]
    cache, stims, dur, n_drop = {}, [], [], []
    for i in range(n_ep):
        d = float(ends[i] - starts[i])
        fps = float(metas[i]["fps"])
        dd = drops[(drops >= starts[i]) & (drops <= ends[i])]
        n_drop.append(int(dd.size))
        di = np.round((dd - starts[i]) * fps).astype(int)
        stims.append(epoch_stimulus(metas[i], d, pre[i], tail[i], di, cache))
        dur.append(d)
    frozen = len(set(seeds)) == 1 and n_ep > 1
    g = _generator(metas[0])
    allc = np.concatenate([s["c"][s["dwell"] >= 1] for s in stims])
    alls = np.concatenate([s["state"][s["dwell"] >= 1] for s in stims])
    sd_state = {q: float(allc[alls == q].std()) for q in (-1, 1)}
    return dict(stims=stims, resp=resp, n_ep=n_ep, pre=pre, tail=tail, metas=metas, seeds=seeds,
                frozen=frozen, stim_dur=np.array(dur), n_drop=n_drop, sd_state=sd_state,
                dark={k: [e for e in r["dark_epochs"] if e < n_ep] for k, r in resp.items()},
                faded={k: [x for x in r["faded"] if x[1] < n_ep] for k, r in resp.items()},
                realized=tuple(g.realized_contrasts), contrast_ratio=float(g.contrast_ratio),
                dwell_s=stims[0]["dwell_s"], n_dwell=int(round(dur[0] / stims[0]["dwell_s"])))


# ── design: one row per imaging sample of a z plane; one response column per ROI ───────

def _nanmean0(x):
    """Column means ignoring NaN, 0 where a column is all NaN (no warning)."""
    k = np.sum(np.isfinite(x), axis=0)
    return np.nansum(x, axis=0) / np.maximum(k, 1)


def plane_rows(ctx, ch, rois, n_lag, lag_dt, first_dwell=2, roll_s=0.0, pos_range=(-np.inf, np.inf)):
    """Design for the ROIs of one z plane, which share their sampling times.

    A: (rows, n_lag) stimulus taps at each imaging sample. Y: (rows, n_roi) dF/F, NaN where an
    ROI's epoch was dropped. Nothing is averaged across ROIs: every ROI is analysed on its own
    against this A. Each ROI's response is demeaned per trial and state, because drift between
    trials and the contrast step are not stimulus-locked and the gains do not use them (Y_raw
    keeps the undemeaned values)."""
    R = ctx["resp"][ch]
    rois = np.asarray(rois, int)
    off = float(np.mean(R["off"][rois])) + R["scan"] / 2.0
    out = {k: [] for k in ("A", "Y", "state", "dwell", "since", "pos", "trial")}
    for e in range(ctx["n_ep"]):
        st = ctx["stims"][e] if not roll_s else rolled(ctx["stims"][e], roll_s)
        t = R["tvecs"][e] - ctx["pre"][e] + off
        s, d, sn, ps = labels(st, t)
        m = (d >= first_dwell) & (ps <= ctx["stim_dur"][e]) & (ps >= pos_range[0]) & (ps < pos_range[1])
        if not m.any():
            continue
        out["A"].append(taps(st, t[m], n_lag, lag_dt))
        out["Y"].append(R["resp"][e][rois].T[m])
        for k, v in (("state", s), ("dwell", d), ("since", sn), ("pos", ps)):
            out[k].append(v[m])
        out["trial"].append(np.full(m.sum(), e))
    d = {k: np.concatenate(v) for k, v in out.items()}
    d["Y_raw"] = d["Y"].copy()
    key = d["trial"] * 3 + d["state"] + 1
    for k in np.unique(key):
        r_ = key == k
        d["Y"][r_] = d["Y"][r_] - _nanmean0(d["Y"][r_])
    d["fold"] = np.clip((d["pos"] // FOLD_S).astype(int), 0, None)
    d["unit"] = d["dwell"] if ctx["frozen"] else d["trial"] * 1000 + d["dwell"]
    d["lags"] = (np.arange(n_lag) + 0.5) * lag_dt
    d["rois"] = rois
    return d


class _Stats:
    """Sufficient statistics of a design by (state, fold) over the rows one ROI can use, so a
    filter for either state with any fold held out is one small solve instead of a pass over
    the rows. ROIs of a plane that share the same usable rows share an instance."""

    def __init__(self, A, state, fold, usable):
        self.rows, self.n, self.sA, self.AtA = {}, {}, {}, {}
        for q in (-1, 1):
            for f in np.unique(fold):
                r_ = usable & (state == q) & (fold == f)
                if r_.any():
                    a = A[r_]
                    self.rows[(q, f)] = r_
                    self.n[(q, f)] = int(r_.sum())
                    self.sA[(q, f)] = a.sum(0)
                    self.AtA[(q, f)] = a.T @ a

    def keys(self, q, exclude=None):
        return [k for k in self.n if k[0] == q and k[1] != exclude]


def _stats_for(cache, d, usable):
    """The _Stats for this usable-row pattern of a plane's design, built once per pattern."""
    key = usable.tobytes()
    if key not in cache:
        cache[key] = _Stats(d["A"], d["state"], d["fold"], usable)
    return cache[key]


def _fit_stats(S, A, y, q, lam, P, exclude=None):
    """Ridge filter (second-difference penalty, lam relative to trace(A'A)/n) for state q on every
    fold but `exclude`, from the stored statistics plus this ROI's A'y. Returns the filter, the
    A and y means it was centred on, and the rows it used."""
    keys = S.keys(q, exclude)
    n = sum(S.n[k] for k in keys)
    rows = np.zeros(y.size, bool)
    for k in keys:
        rows |= S.rows[k]
    ax = sum(S.sA[k] for k in keys) / n
    ay = y[rows].sum() / n
    C = sum(S.AtA[k] for k in keys) - n * np.outer(ax, ax)
    c = A[rows].T @ y[rows] - n * ax * ay
    b = np.linalg.solve(C + lam * (np.trace(C) / C.shape[0]) * P, c)
    return b, ax, ay, rows


def _heldout_r(S, A, y, q, lam, P):
    """Held-out r of state q's filter for one ROI: each fold predicted by a fit on the others."""
    pred, obs = [], []
    for k in S.keys(q):
        b, ax, ay, _ = _fit_stats(S, A, y, q, lam, P, exclude=k[1])
        r_ = S.rows[k]
        pred.append((A[r_] - ax) @ b + ay)
        obs.append(y[r_])
    if len(pred) < 2:
        return np.nan
    pred, obs = np.concatenate(pred), np.concatenate(obs)
    return float(np.corrcoef(pred, obs)[0, 1]) if pred.size > 2 else np.nan


def penalty(n):
    D = np.diff(np.eye(n), n=2, axis=0)
    P = D.T @ D
    return P / (np.trace(P) / n) + 1e-3 * np.eye(n)


def fit(A, Y, lam, P):
    """Ridge with a second-difference penalty, relative to trace(A'A)/n so lam means the
    same thing at either contrast. Y may be 1-D or (rows, columns)."""
    ax, ay = A.mean(0), Y.mean(0)
    a = A - ax
    AtA = a.T @ a
    sc = np.trace(AtA) / AtA.shape[0]
    b = np.linalg.solve(AtA + lam * sc * P, a.T @ (Y - ay))
    return b, ax, ay


def colcorr(a, b):
    a = a - a.mean(0)
    b = b - b.mean(0)
    return (a * b).sum(0) / np.sqrt((a * a).sum(0) * (b * b).sum(0))


def cv_predict(A, Y, fold, lam, P):
    pred = np.empty(Y.shape)
    coefs = {}
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        b, ax, ay = fit(A[tr], Y[tr], lam, P)
        pred[te] = (A[te] - ax) @ b + ay
        coefs[int(f)] = b
    return pred, coefs


def ci(a, axis=0):
    return np.nanpercentile(a, [2.5, 97.5], axis=axis)


def _boot_w(n, rng):
    return np.bincount(rng.integers(0, n, n), minlength=n).astype(float)


# ── ROI gate ────────────────────────────────────────────────────────────────

def _roi_cv_r(ctx, ch, rois_by_plane, n_lag, lag_dt, lam, min_lag, steady_start, pos_range, roll_s=0.0):
    R = ctx["resp"][ch]
    r = np.full(R["n_roi"], np.nan)
    lags = (np.arange(n_lag) + 0.5) * lag_dt
    use = lags >= min_lag
    P = penalty(int(use.sum()))
    for p, rois in rois_by_plane.items():
        A, Y, F = [], [], []
        off = float(np.mean(R["off"][rois])) + R["scan"] / 2.0
        for e in range(ctx["n_ep"]):
            st = ctx["stims"][e] if not roll_s else rolled(ctx["stims"][e], roll_s)
            t = R["tvecs"][e] - ctx["pre"][e] + off
            s, d, sn, ps = labels(st, t)
            y = R["resp"][e][rois].T
            with np.errstate(invalid="ignore"):
                m = ((s > 0) & (d >= 1) & (sn >= steady_start) & (ps >= pos_range[0])
                     & (ps < pos_range[1]) & np.all(np.isfinite(y), axis=1))
            A.append(taps(st, t[m], n_lag, lag_dt)[:, use])
            Y.append(y[m])
            F.append(np.clip((ps[m] // FOLD_S).astype(int), 0, None))
        A, Y, F = np.vstack(A), np.vstack(Y), np.concatenate(F)
        if np.unique(F).size < 2:
            raise ValueError("pos_range leaves fewer than two folds")
        pred, _ = cv_predict(A, Y, F, lam, P)
        r[rois] = colcorr(pred, Y)
    return r


def gate_rois(ctx, ch="calcium", filter_length=3.0, lag_dt=0.02, min_lag=0.04, steady_start=None,
              n_rolls=14, q=0.05, pos_range=(-np.inf, np.inf), lam_grid=(10., 100., 316., 1000.),
              log=print):
    """Per-ROI high-contrast filter against the rolled-stimulus null. Returns dict with the
    BH-selected indices, p values, held-out r and the null."""
    R = ctx["resp"][ch]
    n_lag = int(round(filter_length / lag_dt))
    steady_start = filter_length if steady_start is None else steady_start
    by_plane = {int(p): np.flatnonzero(R["plane"] == p) for p in np.unique(R["plane"])}
    best = None
    for lam in lam_grid:
        r = _roi_cv_r(ctx, ch, by_plane, n_lag, lag_dt, lam, min_lag, steady_start, pos_range)
        if best is None or np.nanmedian(r) > best[1]:
            best = (lam, np.nanmedian(r), r)
    lam, _, r = best
    null = np.concatenate([_roi_cv_r(ctx, ch, by_plane, n_lag, lag_dt, lam, min_lag, steady_start,
                                     pos_range, roll_s=k * 2 * ctx["dwell_s"])
                           for k in range(1, n_rolls + 1)])
    ok = np.isfinite(r)
    p = np.ones(R["n_roi"])
    p[ok] = nullcal.pvalues(r[ok], null)[0]
    sel, _ = nullcal.benjamini_hochberg(p, q=q)
    log('  gate {} (stimulus {}..{} s): ridge {:g}, median held-out r {:.3f}, null 95th pct {:.3f}; '
        'BH q={} keeps {} of {} (p<0.05: {})'.format(ch, pos_range[0], pos_range[1], lam, np.nanmedian(r),
                                                    np.nanpercentile(null, 95), q, len(sel), R["n_roi"],
                                                    int((p < 0.05).sum())))
    return dict(keep=np.asarray(sel, int), p=p, r=r, null=null, lam=lam)


# ── LN analysis of one channel ──────────────────────────────────────────────

def _scale(A, h, sd):
    """Baccus & Meister: the factor making var(A @ h) equal the stimulus variance sd**2."""
    s = (A @ h).std()
    return sd / s if s > 0 else np.nan


def _pop(vals, rng, n_boot):
    """Median across ROIs and its 95% CI from resampling ROIs. (nan, [nan, nan], 0) if none."""
    v = np.asarray(vals, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, [np.nan, np.nan], 0
    idx = rng.integers(0, v.size, (n_boot, v.size))
    return float(np.median(v)), ci(np.median(v[idx], 1)).tolist(), int(v.size)


def _pop_ratio_of_medians(a, b, rng, n_boot):
    """median(a) / median(b) over the ROIs where both are finite, with a CI from resampling ROIs
    jointly. Used where a per-ROI quotient would be unstable (a noisy denominator near zero)."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    f = np.isfinite(a) & np.isfinite(b)
    a, b = a[f], b[f]
    if a.size == 0:
        return np.nan, [np.nan, np.nan]
    idx = rng.integers(0, a.size, (n_boot, a.size))
    with np.errstate(divide="ignore", invalid="ignore"):
        return float(np.median(a) / np.median(b)), ci(np.median(a[idx], 1) / np.median(b[idx], 1)).tolist()


def _pop_curve(M, rng, n_boot):
    """Median across ROIs of a (n_roi, n_bin) array per bin, with ROI-resampling CI; NaN bins ignored."""
    M = np.asarray(M, float)
    good = np.isfinite(M).any(1)
    M = M[good]
    if M.shape[0] == 0:
        return np.full(M.shape[1], np.nan), np.full((2, M.shape[1]), np.nan)
    idx = rng.integers(0, M.shape[0], (n_boot, M.shape[0]))
    return np.nanmedian(M, 0), ci(np.nanmedian(M[idx], 1))


def _analyse_roi(d, j, S, lam, P, sd, steady_start, n_boot, n_cal, rng):
    """One ROI (column j of a plane design) through the whole LN analysis: filters in each
    contrast, B&M generators, gain and shared-range slopes (own filter and high-filter readouts),
    the non-adapting reference, gain against time since the switch, and the stretch. Returns a
    dict, or None if either contrast has fewer than two folds of usable rows."""
    if any(len(S.keys(q)) < 2 for q in (-1, 1)):
        return None
    A, lags = d["A"], d["lags"]
    lag_dt = float(lags[1] - lags[0])
    y = d["Y"][:, j]
    ok = np.isfinite(y)
    yz = np.where(ok, y, 0.0)                    # unusable rows never enter a fit or a sum
    n = y.size
    ss = d["since"] >= steady_start
    st_m = {q: ok & ss & (d["state"] == (1 if q else -1)) for q in (0, 1)}
    is_q = {q: ok & (d["state"] == (1 if q else -1)) for q in (0, 1)}
    out = dict(cv_r=[_heldout_r(S, A, yz, -1, lam, P), _heldout_r(S, A, yz, 1, lam, P)],
               n_rows=int(ok.sum()))
    h = {q: _fit_stats(S, A, yz, (1 if q else -1), lam, P)[0] for q in (0, 1)}
    out["h"] = h
    out["h_scaled"] = {q: h[q] * _scale(A[st_m[q]], h[q], sd[q]) for q in (0, 1)}
    new = (lags[None, :] + lag_dt / 2) <= d["since"][:, None]       # tap wholly after the switch
    folds = sorted({k[1] for k in S.n})

    def generators(yv, split=True, own=True):
        """Each row through a filter fitted without the row's fold, B&M scaled in its state.
        own=True: each state's own filter (the convention). own=False: the high-contrast filter
        in both states -- immune to noise in the low-contrast filter, blind to a change of
        shape. Also returns the raw fold prediction (for residuals)."""
        g, g_new, g_old, pred = np.full(n, np.nan), np.zeros(n), np.zeros(n), np.full(n, np.nan)
        for f in folds:
            te = d["fold"] == f
            hf = {}
            for q in (0, 1):
                b, ax, ay, tr = _fit_stats(S, A, yv, (1 if q else -1), lam, P, exclude=f)
                hf[q] = (b, _scale(A[tr], b, sd[q]), ax, ay, tr)
            for q in (0, 1):
                cur = te & is_q[q]
                if not cur.any():
                    continue
                hc = hf[q][0] * hf[q][1] if own else hf[1][0] * _scale(A[hf[q][4]], hf[1][0], sd[q])
                g[cur] = A[cur] @ hc
                pred[cur] = (A[cur] - hf[q][2]) @ hf[q][0] + hf[q][3]
                if split:
                    g_new[cur] = (A[cur] * new[cur]) @ hc
                    g_old[cur] = (A[cur] * ~new[cur]) @ (hf[1 - q][0] * hf[1 - q][1])
        return g, g_new, g_old, pred

    g_own, g_new, g_old, pred = generators(yz)
    g_hf = generators(yz, split=False, own=False)[0]
    units = np.unique(d["unit"][ok])
    uidx = np.searchsorted(units, d["unit"])
    nu = units.size
    W = [np.ones(nu)] + [_boot_w(nu, rng) for _ in range(n_boot)]

    def slope_stats(m, g, yv):
        u, gm, ym = uidx[m], g[m], yv[m]
        return np.stack([np.bincount(u, minlength=nu), np.bincount(u, gm, nu), np.bincount(u, ym, nu),
                         np.bincount(u, gm * gm, nu), np.bincount(u, gm * ym, nu)], 1)

    def slope(Sx, w):
        k, sg, sy, sgg, sgy = (Sx * w[:, None]).sum(0)
        den = k * sgg - sg * sg
        return (k * sgy - sg * sy) / den if den > 0 else np.nan

    def ratios(g, yv):
        """Point estimates of (whole-range, shared-range) low/high slope ratios."""
        one = np.ones(nu)
        lim_ = np.percentile(g[st_m[0]], [1, 99])
        in_ = (g >= lim_[0]) & (g <= lim_[1])
        s = [[slope(slope_stats(st_m[q] & sel, g, yv), one) for q in (0, 1)]
             for sel in (np.ones(n, bool), in_)]
        return s[0][0] / s[0][1], s[1][0] / s[1][1]

    for sfx, g in (("", g_own), ("_high_filter", g_hf)):
        G = {q: np.array([slope(slope_stats(st_m[q], g, yz), w) for w in W]) for q in (0, 1)}
        lim = np.percentile(g[st_m[0]], [1, 99])      # the generator range the low state covers
        inr = (g >= lim[0]) & (g <= lim[1])
        SR = {q: np.array([slope(slope_stats(st_m[q] & inr, g, yz), w) for w in W]) for q in (0, 1)}
        for key, X in (("gain", G), ("shared_range", SR)):
            with np.errstate(divide="ignore", invalid="ignore"):
                out[key + sfx] = dict(low=float(X[0][0]), high=float(X[1][0]), ratio=float(X[0][0] / X[1][0]),
                                      ci=ci(X[0][1:] / X[1][1:]).tolist())
        out["shared_range" + sfx]["range"] = lim.tolist()
        out["shared_range" + sfx]["frac_high_rows_in_range"] = float(inr[st_m[1]].mean())
    out["generator_sd"] = dict(low=float(g_own[st_m[0]].std()), high=float(g_own[st_m[1]].std()))
    out["nl"] = {q: (g_own[st_m[q]], yz[st_m[q]]) for q in (0, 1)}

    # THE REFERENCE A RATIO HAS TO BE READ AGAINST: this ROI replaced by a non-adapting cell --
    # its high-contrast filter followed by its high-contrast nonlinearity as measured, the same
    # in both states -- plus its own residuals (sign-flipped per unit to vary their alignment),
    # through the same estimators. Three things move that cell's ratio away from 1: a noisy
    # filter shrinks the slope through it (the low-contrast filter is the noisier one); a static
    # nonlinearity that is not a straight line gives the two contrasts different whole-range
    # slopes; and with frozen noise the low and high dwells are different finite samples of the
    # probe whose slow power need not match. ratio_vs_static = measured / what it recovers.
    keys = ("gain", "shared_range", "gain_high_filter", "shared_range_high_filter")
    if n_cal and np.isfinite(out["cv_r"][1]):
        x = A @ h[1]
        hm = st_m[1]
        e_ = np.quantile(x[hm], np.linspace(0, 1, 21))
        b_ = np.clip(np.searchsorted(e_, x[hm], "right") - 1, 0, 19)
        cnt = np.bincount(b_, minlength=20).astype(float)
        with np.errstate(invalid="ignore", divide="ignore"):
            xs = np.bincount(b_, x[hm], 20) / cnt
            ys_ = np.bincount(b_, yz[hm], 20) / cnt
        okb = cnt > 0
        sig = np.interp(x, xs[okb], ys_[okb])             # the high state's static nonlinearity
        lin = x.copy()
        for k in (0, 2):
            r_ = ok & (d["state"] + 1 == k)
            if r_.any():
                sig[r_] -= sig[r_].mean()
                lin[r_] -= lin[r_].mean()
        resid = np.where(ok, yz - np.nan_to_num(pred), 0.0)
        truth = (lin[st_m[0]].std() / sd[0]) / (lin[st_m[1]].std() / sd[1])
        rec = []
        for _ in range(n_cal):
            flip = rng.choice([-1.0, 1.0], nu)[uidx]
            ys = np.where(ok, sig + flip * resid, 0.0)
            rec.append(ratios(generators(ys, split=False)[0], ys) + ratios(generators(ys, False, False)[0], ys))
        rec = np.array(rec)
        out["static"] = dict(truth_linear=float(truth), nonlinearity_x=xs, nonlinearity_y=ys_,
                             mean=rec.mean(0), sd=rec.std(0))
        for i, key_ in enumerate(keys):
            ref = rec[:, i].mean()
            with np.errstate(divide="ignore", invalid="ignore"):
                out[key_]["ratio_vs_static"] = float(out[key_]["ratio"] / ref)
                out[key_]["ci_vs_static"] = (np.array(out[key_]["ci"]) / ref).tolist()

    # gain against time since the switch: y ~ a g_new + b g_old + c, per bin
    nb = GAIN_BINS.size - 1
    bins = np.clip(np.searchsorted(GAIN_BINS, d["since"], "right") - 1, 0, nb - 1)
    X = np.c_[g_new, g_old, np.ones(n)]
    dyn = {}
    for q in (0, 1):
        m = is_q[q]
        Sm = np.zeros((nu, nb, 3, 3))
        B = np.zeros((nu, nb, 3))
        for b_ in range(nb):
            mb = m & (bins == b_)
            for a_ in range(3):
                B[:, b_, a_] = np.bincount(uidx[mb], X[mb, a_] * yz[mb], nu)
                for c_ in range(3):
                    Sm[:, b_, a_, c_] = np.bincount(uidx[mb], X[mb, a_] * X[mb, c_], nu)

        def gains(w):
            A_ = (Sm * w[:, None, None, None]).sum(0)
            b_ = (B * w[:, None, None]).sum(0)
            tr = np.trace(A_[:, :2, :2], axis1=1, axis2=2)[:, None, None]
            empty = A_[:, 2, 2] <= 0                     # a bin no drawn unit has rows in
            A_ = A_ + (1e-8 * tr + empty[:, None, None]) * np.eye(3)
            g_ = np.linalg.solve(A_, b_[..., None])[:, 0, 0]
            g_[empty] = np.nan
            return g_
        dyn[q] = np.array([gains(w) for w in W])
    ssb = GAIN_BINS[:-1] >= steady_start
    early = (GAIN_BINS[:-1] >= 0.25) & (GAIN_BINS[1:] <= 1.5)
    out["dyn"] = {}
    for q, nm in ((0, "to_low"), (1, "to_high")):
        with np.errstate(invalid="ignore", divide="ignore"):
            base = np.nanmean(dyn[q][:, ssb], 1, keepdims=True)
            rel = dyn[q] / base
            e = np.nanmean(dyn[q][:, early], 1) / base[:, 0]
        out["dyn"][nm] = dict(gain=dyn[q][0], relative=rel[0], early_rel=float(e[0]),
                              early_rel_ci=ci(e[1:]).tolist())

    # filter stretch, steady low against steady high (short filters refitted), bootstrap over units
    short = lags < 0.6
    ns = int(short.sum())
    XX = np.zeros((2, nu, ns, ns))
    XY = np.zeros((2, nu, ns))
    for q in (0, 1):
        m = st_m[q]
        A_ = A[m][:, short] - A[m][:, short].mean(0)
        y_ = yz[m] - yz[m].mean()
        um = uidx[m]
        for k in np.unique(um):
            uu = um == k
            XX[q, k] = A_[uu].T @ A_[uu]
            XY[q, k] = A_[uu].T @ y_[uu]
    Ps = penalty(ns)

    def stretch(w):
        hs = []
        for q in (0, 1):
            A_ = (XX[q] * w[:, None, None]).sum(0)
            if not np.trace(A_) > 0:                 # a draw with no rows in this state
                return np.nan
            hs.append(np.linalg.solve(A_ + lam * np.trace(A_) / ns * Ps, (XY[q] * w[:, None]).sum(0)))
        return stretch_fit(hs[0], hs[1], lags[short])
    sv = np.array([stretch(w) for w in W[:101]])
    out["stretch"] = dict(value=float(sv[0]), ci=ci(sv[1:]).tolist(),
                          at_bound=bool(sv[0] <= 0.5 or sv[0] >= 2.0))
    out["nu"] = int(nu)
    return out


def analyse_channel(ctx, ch, keep, filter_length=3.0, lag_dt=0.02, steady_start=None, first_dwell=2,
                    n_rolls=14, n_boot=1000, ridge=None, seed=0, pos_range=(-np.inf, np.inf), n_cal=8,
                    log=print):
    """LN analysis of every kept ROI on its own, then combined across ROIs.

    Each ROI gets its own filter per contrast, nonlinearity, gain and shared-range slopes (own
    filter and high-filter readouts), non-adapting reference, gain dynamics and stretch. The
    population numbers are medians across ROIs with 95% CIs from resampling ROIs; every per-ROI
    value is returned in the arrays, with its own dwell-bootstrap CI where one is computed. One
    ridge value per channel, chosen by the median held-out r across ROIs.

    pos_range restricts every row to that span of stimulus time (the cross-fit check).
    Returns (results dict, arrays dict for figures)."""
    rng = np.random.default_rng(seed)
    n_lag = int(round(filter_length / lag_dt))
    steady_start = filter_length if steady_start is None else steady_start
    P = penalty(n_lag)
    sd = {0: ctx["sd_state"][-1], 1: ctx["sd_state"][1]}
    R = ctx["resp"][ch]
    keep = np.asarray(keep, int)
    planes = [(int(p), keep[R["plane"][keep] == p]) for p in np.unique(R["plane"][keep])]

    def designs_at(roll_s=0.0):
        return [plane_rows(ctx, ch, rois, n_lag, lag_dt, first_dwell, roll_s=roll_s, pos_range=pos_range)
                for _, rois in planes]

    def usable(d, j):
        return np.isfinite(d["Y"][:, j]) & (d["since"] >= steady_start)

    def per_roi_r(dlist, lam):
        """Held-out r (low, high) of every kept ROI at this ridge, in kept order."""
        rr = []
        for d in dlist:
            cache = d.setdefault("_cache", {})
            for j in range(d["rois"].size):
                S = _stats_for(cache, d, usable(d, j))
                yz = np.where(np.isfinite(d["Y"][:, j]), d["Y"][:, j], 0.0)
                if any(len(S.keys(q)) < 2 for q in (-1, 1)):
                    rr.append([np.nan, np.nan])
                    continue
                rr.append([_heldout_r(S, d["A"], yz, -1, lam, P), _heldout_r(S, d["A"], yz, 1, lam, P)])
        return np.array(rr).reshape(-1, 2)

    designs = designs_at()
    rois = np.concatenate([d["rois"] for d in designs])
    grid = RIDGE_GRID if ridge is None else np.array([float(ridge)])
    scores = [np.nanmedian(per_roi_r(designs, l).sum(1)) for l in grid]
    lam = float(grid[int(np.nanargmax(scores))])
    r_roi = per_roi_r(designs, lam)                                     # (n_roi, 2)
    null = np.array([per_roi_r(designs_at((k + 1) * 2 * ctx["dwell_s"]), lam) for k in range(n_rolls)]
                    ).reshape(n_rolls, rois.size, 2) if n_rolls else np.full((0, rois.size, 2), np.nan)
    with np.errstate(invalid="ignore"):
        null_max = np.nanmax(null, 0) if n_rolls else np.full((rois.size, 2), np.nan)   # per ROI, per state
    passes = r_roi > null_max                                            # (n_roi, 2); False where NaN

    # the LN analysis of every ROI
    res = []
    for d in designs:
        for j in range(d["rois"].size):
            S = _stats_for(d["_cache"], d, usable(d, j))
            res.append(_analyse_roi(d, j, S, lam, P, sd, steady_start, n_boot, n_cal, rng))
    good = np.array([r_ is not None for r_ in res])

    def col(fn, default=np.nan):
        return np.array([fn(r_) if r_ is not None else default for r_ in res], float)

    out = dict(n_roi=int(rois.size), n_roi_analysed=int(good.sum()), n_rows=int(col(lambda r_: r_["n_rows"], 0).sum()),
               n_planes=len(planes), lag_dt=lag_dt, filter_length=filter_length, steady_start=steady_start,
               first_dwell=first_dwell, frozen=bool(ctx["frozen"]), pos_range=[float(x) for x in pos_range],
               ridge=lam, combine="median across ROIs; CIs from resampling ROIs",
               cv_r_low=float(np.nanmedian(r_roi[:, 0])), cv_r_high=float(np.nanmedian(r_roi[:, 1])),
               null_max_low=float(np.nanmedian(null_max[:, 0])), null_max_high=float(np.nanmedian(null_max[:, 1])),
               n_roi_filter_low=int(passes[:, 0].sum()), n_roi_filter_high=int(passes[:, 1].sum()),
               n_roi_filter_both=int(passes.all(1).sum()))
    with np.errstate(invalid="ignore"):
        out["p_low"] = float(np.nanmedian((1 + (null[:, :, 0] >= r_roi[None, :, 0]).sum(0)) / (1.0 + n_rolls))) if n_rolls else np.nan
        out["p_high"] = float(np.nanmedian((1 + (null[:, :, 1] >= r_roi[None, :, 1]).sum(0)) / (1.0 + n_rolls))) if n_rolls else np.nan
    # None when no null was run (the cross-fit repeats), not False. The population has a filter in a
    # state when the median ROI beats the median ROI's null there.
    out["has_filter_both_states"] = (bool(out["cv_r_high"] > out["null_max_high"] and out["cv_r_low"] > out["null_max_low"])
                                     if n_rolls else None)
    log('  {}: {} ROIs in {} planes ({} analysed, {} rows); ridge {:g}; median held-out r low {:.3f} high {:.3f}, '
        'median null max {:.3f} {:.3f}; ROIs beating their null: low {}, high {}, both {}'.format(
            ch, out["n_roi"], out["n_planes"], out["n_roi_analysed"], out["n_rows"], lam, out["cv_r_low"], out["cv_r_high"],
            out["null_max_low"], out["null_max_high"], out["n_roi_filter_low"], out["n_roi_filter_high"], out["n_roi_filter_both"]))

    # population filters: mean of the B&M-scaled per-ROI filters (sign preserved)
    lags = designs[0]["lags"]
    H = {q: np.array([r_["h_scaled"][q] if r_ is not None else np.full(n_lag, np.nan) for r_ in res]) for q in (0, 1)}
    h = {q: np.nanmean(H[q], 0) if good.any() else np.full(n_lag, np.nan) for q in (0, 1)}
    short = lags < 0.6
    for q, nm in ((0, "low"), (1, "high")):
        hs = h[q][short]
        if np.all(np.isfinite(hs)):
            i, j0, j1 = int(np.argmax(np.abs(hs))), int(np.argmin(hs)), int(np.argmax(hs))
            out["filter_" + nm] = dict(peak_abs_lag=float(lags[i]), trough_lag=float(lags[j0]), trough=float(hs[j0]),
                                       peak_lag=float(lags[j1]), peak=float(hs[j1]))
        else:
            out["filter_" + nm] = dict(peak_abs_lag=np.nan, trough_lag=np.nan, trough=np.nan, peak_lag=np.nan, peak=np.nan)

    # gains: per-ROI ratios combined
    keys = ("gain", "shared_range", "gain_high_filter", "shared_range_high_filter")
    roi_arrays = dict(roi=rois, roi_cv_r=r_roi, roi_null_max=null_max, roi_pass=passes.astype(int))
    for key in keys:
        v, c, nk = _pop(col(lambda r_: r_[key]["ratio"]), rng, n_boot)
        lo, lo_ci, _ = _pop(col(lambda r_: r_[key]["low"]), rng, n_boot)
        hi, hi_ci, _ = _pop(col(lambda r_: r_[key]["high"]), rng, n_boot)
        out[key] = dict(low=lo, high=hi, ratio=v, ci=c, low_ci=lo_ci, high_ci=hi_ci, n_roi=nk,
                        roi_iqr=np.nanpercentile(col(lambda r_: r_[key]["ratio"]), [25, 75]).tolist() if nk else [np.nan, np.nan])
        roi_arrays["roi_" + key + "_ratio"] = col(lambda r_: r_[key]["ratio"])
        roi_arrays["roi_" + key + "_ci"] = np.array([r_[key]["ci"] if r_ is not None else [np.nan, np.nan] for r_ in res])
        has_static = col(lambda r_: 1.0 if "ratio_vs_static" in r_[key] else 0.0, 0.0) > 0
        if has_static.any():
            # Measured and reference are both per-ROI quantities; the population comparison is
            # median(measured) / median(reference) over ROIs, not the median of per-ROI quotients:
            # a reference near zero (a noisy low-contrast filter) makes the quotient unstable.
            vs = col(lambda r_: r_[key].get("ratio_vs_static", np.nan))
            meas = col(lambda r_: r_[key]["ratio"] if "static" in r_ else np.nan)
            ref = col(lambda r_: r_["static"]["mean"][keys.index(key)] if "static" in r_ else np.nan)
            v2, c2 = _pop_ratio_of_medians(meas, ref, rng, n_boot)
            out[key]["ratio_vs_static"] = v2
            out[key]["ci_vs_static"] = c2
            out[key]["roi_median_ratio_vs_static"] = _pop(vs, rng, n_boot)[0]
            roi_arrays["roi_" + key + "_ratio_vs_static"] = vs
            roi_arrays["roi_" + key + "_static_ref"] = ref
    lim = np.nanmedian(np.array([r_["shared_range"]["range"] if r_ is not None else [np.nan, np.nan] for r_ in res]), 0)
    out["shared_range"]["range"] = lim.tolist()
    out["shared_range"]["frac_high_rows_in_range"] = float(np.nanmedian(col(lambda r_: r_["shared_range"]["frac_high_rows_in_range"])))
    out["shared_range_high_filter"]["range"] = np.nanmedian(
        np.array([r_["shared_range_high_filter"]["range"] if r_ is not None else [np.nan, np.nan] for r_ in res]), 0).tolist()
    out["generator_sd"] = dict(low=float(np.nanmedian(col(lambda r_: r_["generator_sd"]["low"]))),
                               high=float(np.nanmedian(col(lambda r_: r_["generator_sd"]["high"]))))

    # the non-adapting reference, per ROI then combined
    with_static = [r_ for r_ in res if r_ is not None and "static" in r_]
    if with_static:
        rec = np.array([r_["static"]["mean"] for r_ in with_static])
        truth = np.array([r_["static"]["truth_linear"] for r_ in with_static])
        out["static_reference"] = dict(
            truth_linear=float(np.nanmedian(truth)), n=int(len(with_static)), n_draws=int(n_cal),
            note="per ROI: a non-adapting cell = that ROI's high filter and high nonlinearity in both states, "
                 "with its own residuals. truth_linear = its exact B&M ratio if the nonlinearity were a straight "
                 "line; the others = what each estimator recovers for it (medians across ROIs; _sd = median of "
                 "each ROI's sd over its draws)",
            **{k: float(np.nanmedian(rec[:, i])) for i, k in enumerate(keys)},
            **{k + "_sd": float(np.nanmedian([r_["static"]["sd"][i] for r_ in with_static])) for i, k in enumerate(keys)})
        roi_arrays["roi_truth_linear"] = col(lambda r_: r_["static"]["truth_linear"] if "static" in r_ else np.nan)
        log('  {}: a non-adapting cell reads {:.3f} if linear (median ROI); through the estimators {}'.format(
            ch, out["static_reference"]["truth_linear"],
            ', '.join('{} {:.2f}'.format(k, out["static_reference"][k]) for k in keys)))
        log('  {}: RELATIVE TO A NON-ADAPTING CELL (median ROI, CI over ROIs): '.format(ch) + '; '.join(
            '{} {:.2f} [{:.2f}, {:.2f}]'.format(k, out[k]["ratio_vs_static"], *out[k]["ci_vs_static"])
            for k in keys if "ratio_vs_static" in out[k]))
    if out["has_filter_both_states"] is False:
        log('  WARNING {}: the median ROI\'s {} filter does not beat its null, so no gain ratio here is '
            'interpretable'.format(ch, "low-contrast" if out["cv_r_high"] > out["null_max_high"] else "high-contrast"))
    log('  {}: steady gain low {:.4f} high {:.4f}, ratio {:.2f} [{:.2f}, {:.2f}] (ROI IQR {:.2f}-{:.2f}); shared-range '
        'ratio {:.2f} [{:.2f}, {:.2f}]'.format(ch, out["gain"]["low"], out["gain"]["high"], out["gain"]["ratio"],
                                               *out["gain"]["ci"], *out["gain"]["roi_iqr"], out["shared_range"]["ratio"],
                                               *out["shared_range"]["ci"]))

    # nonlinearities for the figure: every ROI's (generator, response) pairs together, 12 quantile bins
    nl = {}
    for q in (0, 1):
        pairs = [r_["nl"][q] for r_ in res if r_ is not None]
        if not pairs:
            nl[q] = np.full((12, 3), np.nan)
            continue
        g = np.concatenate([p_[0] for p_ in pairs])
        yv = np.concatenate([p_[1] for p_ in pairs])
        e = np.quantile(g, np.linspace(0, 1, 13))
        b = np.clip(np.searchsorted(e, g, "right") - 1, 0, 11)
        nl[q] = np.array([[g[b == i].mean(), yv[b == i].mean(), yv[b == i].std() / np.sqrt(max((b == i).sum(), 1))]
                          for i in range(12)])

    # gain against time since the switch: per ROI, then the median ROI per bin
    out["gain_after_switch"] = dict(bins=GAIN_BINS.tolist())
    dyn_arr = {}
    for q, nm in ((0, "to_low"), (1, "to_high")):
        G = np.array([r_["dyn"][nm]["gain"] if r_ is not None else np.full(GAIN_BINS.size - 1, np.nan) for r_ in res])
        Rl = np.array([r_["dyn"][nm]["relative"] if r_ is not None else np.full(GAIN_BINS.size - 1, np.nan) for r_ in res])
        gm, gc = _pop_curve(G, rng, n_boot)
        rm, rc = _pop_curve(Rl, rng, n_boot)
        ev, ec, _ = _pop(col(lambda r_: r_["dyn"][nm]["early_rel"]), rng, n_boot)
        out["gain_after_switch"][nm] = dict(gain=gm.tolist(), gain_ci=gc.tolist(), relative=rm.tolist(),
                                            relative_ci=rc.tolist(), early_rel=ev, early_rel_ci=ec)
        dyn_arr[q] = Rl
        roi_arrays["roi_early_rel_" + nm] = col(lambda r_: r_["dyn"][nm]["early_rel"])

    sv, sc_, _ = _pop(col(lambda r_: r_["stretch"]["value"]), rng, n_boot)
    out["stretch_low_vs_high"] = dict(value=sv, ci=sc_,
                                      frac_at_bound=float(np.nanmean(col(lambda r_: 1.0 if r_["stretch"]["at_bound"] else 0.0))))
    roi_arrays["roi_stretch"] = col(lambda r_: r_["stretch"]["value"])

    # bm_scale: the per-state factor that puts h_low / h_high in the Baccus & Meister convention.
    # Every per-ROI filter is scaled before averaging, so the population filters already are, and
    # the factor is 1; visanalysis.model.export reads it together with sd_state.
    arrays = dict(lags=lags, h_low=h[0], h_high=h[1], H_low=H[0], H_high=H[1], nl_low=nl[0], nl_high=nl[1],
                  gain_dyn_low=dyn_arr[0], gain_dyn_high=dyn_arr[1], null=null,
                  bm_scale=np.array([1.0, 1.0]), sd_state=np.array([sd[0], sd[1]]), **roi_arrays)
    return out, arrays



def stretch_fit(h, h0, lags):
    """Time stretch s that best maps the reference h0 onto h, with a free amplitude."""
    best, bs = np.nan, -np.inf
    for s in np.arange(0.5, 2.01, 0.02):
        q = np.interp(lags / s, lags, h0, right=0.0)
        den = q @ q
        if den <= 0:
            continue
        a = (q @ h) / den
        sc = -np.sum((h - a * q) ** 2)
        if sc > bs:
            best, bs = s, sc
    return best


# ── model-free gain (frozen seed only) ──────────────────────────────────────

def _binned_by_dwell(ctx, ch, keep, bin_s, lo, hi):
    """(n_roi, n_ep, n_dwell, n_bin) dF/F of each kept ROI, binned by time from each dwell's
    start over [lo, hi); NaN where an ROI has no usable sample in a bin."""
    R = ctx["resp"][ch]
    keep = np.asarray(keep, int)
    dwell_s, n_dw = ctx["dwell_s"], ctx["n_dwell"]
    edges = np.arange(lo, hi + 1e-9, bin_s)
    nb = edges.size - 1
    S = np.zeros((keep.size, ctx["n_ep"], n_dw, nb))
    W = np.zeros_like(S)
    where = {int(r): i for i, r in enumerate(keep)}
    for p in np.unique(R["plane"][keep]):
        rois = keep[R["plane"][keep] == p]
        idx = np.array([where[int(r)] for r in rois])
        off = float(np.mean(R["off"][rois])) + R["scan"] / 2.0
        for e in range(ctx["n_ep"]):
            t = R["tvecs"][e] - ctx["pre"][e] + off
            Y = R["resp"][e][rois]                       # (n_roi_plane, time)
            fin = np.isfinite(Y)
            Yz = np.where(fin, Y, 0.0)
            for k in range(n_dw):
                rel = t - k * dwell_s
                m = (rel >= edges[0]) & (rel < edges[-1])
                if not m.any():
                    continue
                b = np.searchsorted(edges, rel[m], "right") - 1
                ii = np.repeat(idx[:, None], m.sum(), 1).ravel()
                bb = np.tile(b, idx.size)
                np.add.at(S[:, e, k], (ii, bb), Yz[:, m].ravel())
                np.add.at(W[:, e, k], (ii, bb), fin[:, m].ravel().astype(float))
    with np.errstate(invalid="ignore"):
        return S / np.where(W > 0, W, np.nan), 0.5 * (edges[1:] + edges[:-1])


def model_free_gain(ctx, ch, keep, window=(2.0, 10.0), first_dwell=2, bin_s=0.25, n_boot=5000, seed=4,
                    detrend=True):
    """Stimulus-locked SD per state for each kept ROI, from the covariance between trials at the
    same stimulus position; the population ratio is the median ROI, with a CI from resampling
    ROIs. Needs every trial to show the same noise, so it returns None otherwise."""
    if not ctx["frozen"]:
        return None
    rng = np.random.default_rng(seed)
    M, cen = _binned_by_dwell(ctx, ch, keep, bin_s, 0.0, ctx["dwell_s"])
    w = (cen >= window[0]) & (cen < window[1])
    n_dw = ctx["n_dwell"]
    state_of = _state_of_dwells(ctx)
    hi = np.array([k for k in range(first_dwell, n_dw) if state_of[k] > 0])
    lo = np.array([k for k in range(first_dwell, n_dw) if state_of[k] < 0])
    n_roi = M.shape[0]
    # With a frozen seed, "the same stimulus position" is also "the same time since the switch".
    # Each ROI's step and its relaxation repeat in every dwell of a state and in every trial, so
    # they would count as stimulus-locked variance, most of all at low contrast where the probe
    # response is small. Remove, per ROI and trial, the mean time course over that state's
    # dwells. Each dwell loses 1/n of its own probe signal in the process, equally in both states
    # when their dwell counts match, so the ratio is unaffected.
    for ks in ((hi, lo) if detrend else ()):
        with np.errstate(invalid="ignore"):
            mu = np.nansum(M[:, :, ks, :], axis=2, keepdims=True) / np.maximum(
                np.isfinite(M[:, :, ks, :]).sum(axis=2, keepdims=True), 1)
        M[:, :, ks, :] = M[:, :, ks, :] - mu
    lk, tot = np.full((n_roi, n_dw), np.nan), np.full((n_roi, n_dw), np.nan)
    for i in range(n_roi):
        for k in range(first_dwell, n_dw):
            seg = M[i][:, k][:, w]                               # (trials, bins); not M[i, :, k, w],
            # which puts the advanced axes first and hands back (bins, trials)
            seg = seg[np.isfinite(seg).mean(1) > 0.9]          # trials dropped as dark are absent
            seg = seg[:, np.all(np.isfinite(seg), 0)]
            if seg.shape[0] < 2 or seg.shape[1] < 2:
                continue
            x = seg - seg.mean(1, keepdims=True)
            n = x.shape[0]
            c = x @ x.T / x.shape[1]
            lk[i, k] = (c.sum() - np.trace(c)) / (n * (n - 1))
            tot[i, k] = np.trace(c) / n
    # Per ROI, the stimulus-locked variance in each state (a between-trial covariance, so noise
    # does not enter it in expectation). The population SD ratio averages these covariances
    # across ROIs BEFORE the square root: sqrt of a single ROI's covariance clipped at zero is
    # biased upward wherever the locked variance is small against the noise, which is the low
    # contrast state, and that bias would read as adaptation. Per-ROI ratios are returned too.
    with np.errstate(invalid="ignore", divide="ignore"):
        lk_h, lk_l = np.nanmean(lk[:, hi], 1), np.nanmean(lk[:, lo], 1)
        sh_roi, sl_roi = np.sqrt(np.maximum(lk_h, 0)), np.sqrt(np.maximum(lk_l, 0))
        ratio = np.where(sl_roi > 0, sh_roi / sl_roi, np.nan)
    good = np.isfinite(lk_h) & np.isfinite(lk_l)
    if not good.any():
        return None
    gh, gl = lk_h[good], lk_l[good]

    def sd_ratio(ix):
        s_h, s_l = np.sqrt(max(gh[ix].mean(), 0)), np.sqrt(max(gl[ix].mean(), 0))
        return s_h, s_l, (s_h / s_l if s_l > 0 else np.nan)
    sh, sl, r = sd_ratio(np.arange(gh.size))
    idx = rng.integers(0, gh.size, (n_boot, gh.size))
    bc = ci(np.array([sd_ratio(ix)[2] for ix in idx]))
    v = ratio[np.isfinite(ratio)]
    # The contrast ratio the cells actually saw in these windows. With a frozen seed the probe's
    # sample variance differs between the high and the low dwells, and that difference is not a
    # gain change.
    st0 = ctx["stims"][0]
    win = (st0["dwell"] >= first_dwell) & (st0["since"] >= window[0]) & (st0["since"] < window[1])
    cr = float(st0["c"][win & (st0["state"] > 0)].std() / st0["c"][win & (st0["state"] < 0)].std())
    with np.errstate(invalid="ignore"):
        noise_h = np.sqrt(np.maximum(np.nanmean(tot[:, hi] - lk[:, hi], 1), 0))
        noise_l = np.sqrt(np.maximum(np.nanmean(tot[:, lo] - lk[:, lo], 1), 0))
    return dict(sd_high=float(sh), sd_low=float(sl), sd_ratio=float(r), sd_ratio_ci=bc.tolist(),
                contrast_ratio=cr, contrast_ratio_nominal=ctx["contrast_ratio"],
                gain_ratio=float(cr / r), gain_ratio_ci=[float(cr / bc[1]), float(cr / bc[0])],
                noise_sd_high=float(np.nanmedian(noise_h)), noise_sd_low=float(np.nanmedian(noise_l)),
                n_roi=int(good.sum()), roi_sd_ratio=ratio.tolist(),
                roi_sd_ratio_median=float(np.median(v)) if v.size else np.nan,
                roi_sd_ratio_iqr=np.percentile(v, [25, 75]).tolist() if v.size else [np.nan, np.nan],
                combine="per-ROI between-trial covariances averaged across ROIs, then sqrt; CI from resampling ROIs",
                window=list(window), bin_s=bin_s)


# ── switch-triggered averages ───────────────────────────────────────────────

def _exp_tau(lag, m):
    s = (lag > 0.5) & (lag < 10)
    f = lambda t, a, b, tau: a + b * np.exp(-t / tau)
    try:
        p, _ = curve_fit(f, lag[s], m[s], p0=[m[s][-8:].mean(), m[s][0] - m[s][-8:].mean(), 2.0],
                         bounds=([-5, -5, 0.05], [5, 5, 60]), maxfev=4000)
        return float(p[2])
    except Exception:
        return np.nan


def switch_triggered(ctx, ch, keep, bin_s=STA_BIN, first_dwell=2, n_boot=1000, seed=2):
    """Each kept ROI's dF/F against time from each switch, by direction, one dwell either side:
    averaged over trials, then over the ROI's switches. The population curve is the mean across
    ROIs, with CIs from resampling ROIs; per-ROI steps are returned in the arrays."""
    rng = np.random.default_rng(seed)
    dwell_s = ctx["dwell_s"]
    curves, lag = _binned_by_dwell(ctx, ch, keep, bin_s, -dwell_s, dwell_s)
    with np.errstate(invalid="ignore"):
        C = np.nansum(curves, 1) / np.maximum(np.isfinite(curves).sum(1), 1)     # (n_roi, n_dwell, n_bin)
        C[np.isfinite(curves).sum(1) == 0] = np.nan
    state_of = _state_of_dwells(ctx)
    out, arr = {}, dict(lag=lag)
    pre_w = (lag > -3.0) & (lag < 0)
    late = (lag >= 0.75 * dwell_s) & (lag < dwell_s)
    n_roi = C.shape[0]
    for nm, sgn in (("to_high", 1), ("to_low", -1)):
        ks = np.array([k for k in range(first_dwell, ctx["n_dwell"]) if state_of[k] == sgn])
        with np.errstate(invalid="ignore"):
            per_roi = np.nanmean(C[:, ks], 1)                                       # (n_roi, n_bin)
            m = np.nanmean(per_roi, 0)
            ix = rng.integers(0, n_roi, (n_boot, n_roi))
            boot = np.nanmean(per_roi[ix], 1)
            step_roi = per_roi[:, late].mean(1) - per_roi[:, pre_w].mean(1)
        step_b = boot[:, late].mean(1) - boot[:, pre_w].mean(1)
        taus = np.array([_exp_tau(lag, b) for b in boot[:200]])
        out[nm] = dict(n_switch=int(ks.size), n_roi=int(n_roi), step=float(m[late].mean() - m[pre_w].mean()),
                       step_ci=ci(step_b).tolist(), step_roi_median=float(np.nanmedian(step_roi)),
                       step_roi_iqr=np.nanpercentile(step_roi, [25, 75]).tolist(),
                       tau=_exp_tau(lag, m), tau_ci=ci(taus).tolist(), frac_tau_at_bound=float(np.mean(taus > 55)),
                       level_pre=float(m[pre_w].mean()), level_late=float(m[late].mean()))
        arr[nm] = m
        arr[nm + "_se"] = np.nanstd(per_roi, 0, ddof=1) / np.sqrt(n_roi)
        arr[nm + "_roi"] = per_roi
        arr[nm + "_roi_step"] = step_roi
    # peak change after the switch (0.5 s boxcar) from the last 1 s before it, half the
    # difference between the two directions
    n = max(1, int(round(0.5 / bin_s)))
    pre1 = (lag >= -1) & (lag < 0)
    after = lag > 0
    U = uniform_filter1d(arr["to_high"], n, mode="nearest") - arr["to_high"][pre1].mean()
    D = uniform_filter1d(arr["to_low"], n, mode="nearest") - arr["to_low"][pre1].mean()
    sw = 0.5 * (U - D)
    i = int(np.argmax(np.abs(sw[after])))
    out["peak_swing"] = dict(value=float(sw[after][i]), at=float(lag[after][i]))
    out["step_half_difference"] = 0.5 * (out["to_high"]["step"] - out["to_low"]["step"])
    return out, arr



# ── synthetic ground truth ──────────────────────────────────────────────────

def simulate_responses(ctx, ch, h_true, lag_dt, gain_low=1.0, gain_high=1.0, nonlinearity=None,
                       gain_tau=0.0, roi_amp=None, noise_sd=0.05, level_step=0.0, h_low=None, seed=0):
    """A copy of ctx whose ch responses are an LN model of known gains at the real ROI timing.

    Each ROI: dF/F = amp * (N(gain(state, since) * (h * s)) + level_step * state) + noise.
    h is h_true in the high state and h_low (default h_true) in the low one; each is in dF/F
    per unit contrast per tap of lag_dt. gain_tau > 0 makes the gain relax exponentially to
    its new value after a switch. nonlinearity acts on the generator (None = linear)."""
    rng = np.random.default_rng(seed)
    R = ctx["resp"][ch]
    h_low = h_true if h_low is None else h_low
    n_lag = h_true.size
    amp = np.ones(R["n_roi"]) if roi_amp is None else np.asarray(roi_amp, float)
    gl = {-1: gain_low, 1: gain_high, 0: 0.0}
    new_resp = []
    for e in range(ctx["n_ep"]):
        st = ctx["stims"][e]
        out = np.full(R["resp"][e].shape, np.nan)
        for p in np.unique(R["plane"]):
            rois = np.flatnonzero(R["plane"] == p)
            off = float(np.mean(R["off"][rois])) + R["scan"] / 2.0
            t = R["tvecs"][e] - ctx["pre"][e] + off
            s, d, sn, ps = labels(st, t)
            T = taps(st, t, n_lag, lag_dt)
            g = np.where(s > 0, T @ h_true, T @ h_low)
            gn = np.array([gl[int(x)] for x in s])
            if gain_tau > 0:
                go = np.array([gl[-int(x)] for x in s])
                w = np.exp(-np.nan_to_num(sn, nan=1e9) / gain_tau)
                gn = gn * (1 - w) + go * w
            x = gn * g
            if nonlinearity is not None:
                x = nonlinearity(x)
            x = x + level_step * s
            out[rois] = amp[rois, None] * x[None, :] + noise_sd * rng.standard_normal((rois.size, t.size))
        new_resp.append(out)
    resp = dict(ctx["resp"])
    resp[ch] = dict(R, resp=new_resp)
    return dict(ctx, resp=resp)
