"""Temporal filter estimation for the full-field ternary noise stimulus.

WHAT THIS IS FOR. The full-field stimulus has no spatial structure, so every ROI
sees the identical signal and there is one filter per ROI rather than a filter
per cell per ROI. That makes the analysis much simpler than the octave STRF
path, and it removes the per-cell SNR division that makes a spatiotemporal STRF
expensive. What remains is a straight linear regression of response on stimulus
history, plus the checks that say whether the answer means anything.

RECONSTRUCTION, NOT REPLAY. The stimulus is regenerated from the recorded
`fullfield_stimulus` metadata -- seed, hold_taus, clip_k, p, fps -- through the
same generator that produced it. Nothing is read back from a movie, and a
dropped display frame needs no correction, because every value is a pure
function of (seed, band, frame). What the metadata cannot supply is WHICH frames
actually rendered; only the display machine knows that, and this module assumes
the nominal mapping frame = round(t * fps).

BINNING, NOT SAMPLING. The display updates at 120 fps and the imaging runs near
10-14 Hz, so each imaging sample integrates roughly ten display frames. The
stimulus is therefore AVERAGED within each imaging frame's window rather than
sampled at its center -- which is what the indicator and the scan physically do.
Point-sampling would alias the fast stimulus and bias the filter toward whatever
phase the sampling happened to land on.

THE ORDER THE CHECKS GO IN. Cross-validated predictive power comes first. A
filter that does not predict held-out epochs is not a filter, and every quantity
derived from it -- gain, shape, peak time -- is then a description of noise. Only
after that is passed does anything else get reported.
"""
import ast
import json
import re

import numpy as np

from visanalysis.util import bleach, nullcal

__all__ = [
    "is_fullfield_ternary",
    "resolve_fullfield_epoch_meta",
    "stimulus_on_imaging_grid",
    "build_design",
    "fit_filter",
    "choose_ridge",
    "cross_validated_r",
    "split_half_filter_r",
    "analyze_series",
    # per-ROI path, used by analyze_data_trf_fullfield.py
    "build_roi_design",
    "design_matrix",
    "fit_rois",
    "ridge_sweep",
    "split_half_rois",
    "filter_spectrum",
    "band_3db",
    # trial-by-trial, repeated stimuli only
    "repeat_groups",
    "fit_rois_by_epoch",
    "trial_reliability",
]

#: Matches "np.float64(0.5)" and friends, emitted by older recordings.
_NP_SCALAR = re.compile(r"\bnp\.\w+\(([^()]*)\)")


def _parse_meta_attr(v, name):
    """Read one metadata dict back out of an hdf5 attribute.

    stimpack stores these by stringifying the dict, so what comes back is a
    Python repr rather than JSON -- single-quoted keys, True rather than true --
    which rules out json.loads and leaves ast.literal_eval. Any numpy scalar
    that reached as_metadata reprs as a call expression and is a SyntaxError to
    literal_eval; unwrapping it is exact, since the argument is the literal
    value. Same tolerance as the octave path, for the same reasons.
    """
    if isinstance(v, dict):
        return v
    if isinstance(v, bytes):
        v = v.decode("utf-8")
    if not isinstance(v, str):
        raise TypeError("{} is a {}, expected a dict or a stringified one"
                        .format(name, type(v).__name__))
    for attempt in (json.loads,
                    ast.literal_eval,
                    lambda s: ast.literal_eval(_NP_SCALAR.sub(r"\1", s))):
        try:
            out = attempt(v)
        except Exception:
            continue
        if isinstance(out, dict):
            return out
    raise ValueError("could not parse {} from {!r}".format(name, v[:120]))


def is_fullfield_ternary(run_parameters):
    """True if a series ran the full-field ternary protocol.

    Keyed on protocol_ID, like octave_stim.is_octave_ternary: there is no movie
    file to name.
    """
    return "FullFieldTernary" in str(run_parameters.get("protocol_ID", ""))


def resolve_fullfield_epoch_meta(epoch_params):
    """Return the stimulus metadata dict for one epoch.

    Written by the protocol under the name it actually uses --
    `fullfield_stimulus`, matching JCS_protocol.get_trial_parameters. Nothing is
    defaulted: a default that changed later would silently rewrite the past, so
    a missing field is an error rather than something to paper over.
    """
    meta = epoch_params.get("fullfield_stimulus")
    if meta is None:
        present = sorted(str(k) for k in epoch_params
                         if "full" in str(k).lower() or "ternary" in str(k).lower())
        raise KeyError(
            "epoch is missing fullfield_stimulus. Without it the stimulus cannot "
            "be reconstructed -- check that process_data.py attached the epoch "
            "parameters for this series. Related keys that ARE present: {}"
            .format(present or "none"))
    return _parse_meta_attr(meta, "fullfield_stimulus")


def _generator(meta):
    # labpack is not installed in the analysis environment, but it is checked
    # out beside visanalysis. octave_stim already knows how to find it and the
    # rule is the same here, so borrow it rather than keeping two copies of the
    # search: there is no movie file, and reconstruction IS a call into labpack.
    from visanalysis.util.octave_stim import _find_labpack
    _find_labpack()
    from labpack.visual_stim.clandinin.fullfield_ternary import (
        FullFieldTernaryStimulus)
    return FullFieldTernaryStimulus.from_metadata(meta)


def stimulus_on_imaging_grid(meta, time_vector, pre_time, band=None):
    """Displayed contrast averaged within each imaging frame.

    `time_vector` is relative to the START OF THE EPOCH, which includes pre_time,
    so stimulus time is time_vector - pre_time. Samples before stimulus onset
    return 0 (the mid-grey the display holds during pre_time), which is the
    correct value rather than a fill: the cell really did see mid-grey then.

    band=None gives the displayed composite; an integer gives that band's raw
    ternary value, which is the regressor to use for a multi-band stimulus (the
    bands are independent, the composite is not a sum of them once clipped).

    The whole epoch is generated in ONE vectorised call and the per-sample means
    are taken from a cumulative sum. Asking frame by frame instead costs a
    backward hash scan per imaging sample -- tens of seconds per epoch, which
    the circular-shift null then multiplies -- for identical numbers.
    """
    gen = _generator(meta)
    tv = np.asarray(time_vector, dtype=float)
    dt = float(np.median(np.diff(tv)))
    t_stim = tv - float(pre_time)

    lo = np.maximum(t_stim - dt / 2.0, 0.0)
    hi = t_stim + dt / 2.0
    f0 = np.floor(lo * gen.fps).astype(np.int64)
    f1 = np.maximum(np.ceil(hi * gen.fps).astype(np.int64), f0 + 1)

    frames = np.arange(0, int(f1.max()), dtype=np.int64)
    vals = (gen.values(frames) if band is None
            else gen.band_values(int(band), frames).astype(float))
    csum = np.concatenate([[0.0], np.cumsum(np.asarray(vals, dtype=float))])

    out = (csum[f1] - csum[f0]) / (f1 - f0)
    out[hi <= 0.0] = 0.0                   # entirely within pre_time
    return out


def build_design(ID, channel, lag_s=10.0, roi_prefix="aligned", band=None,
                 average_rois=True, dff="none", correct_bleach="auto",
                 shift_s=0.0):
    """Design matrix and response for one series.

    Returns (X, y, epoch_id, dt, meta). X is (samples, n_lags) with column j
    holding the stimulus j samples in the past, so a fitted filter reads in
    causal lag order.

    Only samples whose ENTIRE lag window lies within their own epoch are used.
    That costs the first lag_s of every epoch and buys the guarantee that no row
    mixes stimulus from two different seeds, which would otherwise happen at
    every epoch boundary and is invisible once averaged.
    """
    rd = bleach.roi_responses(ID, channel, roi_prefix=roi_prefix, dff=dff,
                              correct_bleach=correct_bleach)
    er = rd["epoch_response"]
    tv = np.asarray(rd["time_vector_by_epoch"][0], dtype=float)
    dt = float(np.median(np.diff(tv)))
    n_lags = int(round(float(lag_s) / dt))

    eps = ID.getEpochParameters()
    n_ep = min(er.shape[1], len(eps))
    pre_time = float(ID.getRunParameters("pre_time"))

    # Baseline PER EPOCH, from the same array as er. bleach.correct flattens
    # each epoch but keeps its STARTING level (F / (env / env[0])), so epochs
    # recorded later, after more bleaching, sit lower. One baseline shared by
    # every epoch then scales each epoch's dF/F by its level -- the declining-
    # gain artifact the bleach correction exists to remove, reintroduced
    # between epochs. On fly_004 (JS256_x_JS261) ch1's epoch levels drifted
    # 1.04 -> 0.96, about 5 points of an apparent 20% trial-gain decline.
    # The per-epoch median is also the better baseline on its own terms: a
    # 300 s epoch of zero-mean noise response measures it far better than the
    # 5 s pre period. roi_response is not used either way: it is the
    # continuous series trace and, after a correction, on a different scale.
    with np.errstate(invalid="ignore"):
        f0 = np.nanmedian(er[:, :n_ep, :], axis=2)[:, :, None]
    with np.errstate(invalid="ignore"):
        if average_rois:
            R = np.nanmean((er[:, :n_ep, :] - f0) / f0, axis=0)
        else:
            R = (er[:, :n_ep, :] - f0) / f0

    meta0 = resolve_fullfield_epoch_meta(eps[0])
    Xs, ys, eid = [], [], []
    for i in range(n_ep):
        meta = resolve_fullfield_epoch_meta(eps[i])
        s = stimulus_on_imaging_grid(meta, tv, pre_time, band=band)
        if shift_s:
            # circular shift for the null: destroys the stimulus-response
            # pairing but keeps the stimulus autocorrelation intact
            s = np.roll(s, int(round(float(shift_s) / dt)))
        for k in range(n_lags, len(tv)):
            Xs.append(s[k - n_lags:k][::-1])
            ys.append(R[i, k] if average_rois else R[:, i, k])
            eid.append(i)

    X = np.asarray(Xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    eid = np.asarray(eid, dtype=int)

    # Scattered NaNs appear whenever every ROI is NaN at a timepoint. Left in,
    # they turn the whole normal-equation solve into NaN rather than degrading
    # it locally, which is a failure mode that looks like a code bug.
    ok = np.isfinite(y) if y.ndim == 1 else np.all(np.isfinite(y), axis=1)
    ok &= np.all(np.isfinite(X), axis=1)
    return X[ok], y[ok], eid[ok], dt, meta0


def fit_filter(X, y, ridge=1e-3):
    """Ridge-regularised least squares. ridge is relative to trace(X'X)/n_lags,
    so it means the same thing regardless of stimulus contrast or units."""
    Xc = X - X.mean(axis=0)
    A = Xc.T @ Xc
    lam = float(ridge) * np.trace(A) / X.shape[1]
    yc = y - y.mean(axis=0)
    return np.linalg.solve(A + lam * np.eye(X.shape[1]), Xc.T @ yc)


def _folds(eid, n_fold, rng):
    """Split by EPOCH, never by sample. Samples within an epoch share stimulus
    history and are not independent, so a sample-wise split leaks the training
    set into the test set and reports a reliability that does not exist."""
    ue = np.unique(eid)
    return np.array_split(rng.permutation(ue), min(n_fold, len(ue)))


def cross_validated_r(X, y, eid, ridge=1e-3, n_fold=8, seed=0):
    """Correlation between held-out response and its prediction.

    This is the gate. A filter that does not predict held-out epochs is not a
    filter, and nothing derived from it should be reported.
    """
    rng = np.random.default_rng(seed)
    pred = np.full(len(y), np.nan)
    for f in _folds(eid, n_fold, rng):
        te = np.isin(eid, f)
        if te.all() or not te.any():
            continue
        h = fit_filter(X[~te], y[~te], ridge)
        pred[te] = (X[te] - X[~te].mean(axis=0)) @ h
    ok = np.isfinite(pred)
    if ok.sum() < 10 or np.std(pred[ok]) == 0:
        return np.nan
    return float(np.corrcoef(pred[ok], y[ok])[0, 1])


def choose_ridge(X, y, eid, grid=None, n_fold=6, seed=0):
    """Pick ridge by cross-validated prediction, not by eye.

    Returns (best_ridge, {ridge: cv_r}). The grid spans six decades because the
    right value depends on how correlated the stimulus is, which depends on the
    hold time -- a value tuned at one tau is wrong at another.
    """
    if grid is None:
        grid = np.logspace(-6, 0, 13)
    scores = {}
    for g in grid:
        scores[float(g)] = cross_validated_r(X, y, eid, ridge=g,
                                             n_fold=n_fold, seed=seed)
    finite = {k: v for k, v in scores.items() if np.isfinite(v)}
    if not finite:
        return float(grid[len(grid) // 2]), scores
    return max(finite, key=finite.get), scores


def split_half_filter_r(X, y, eid, ridge=1e-3, n_rep=200, seed=0):
    """Reliability of the FILTER itself, splitting epochs.

    Reported Spearman-Brown corrected to the full dataset, so it is comparable
    to the cross-validated r above rather than describing half the data.
    """
    rng = np.random.default_rng(seed)
    ue = np.unique(eid)
    rs = []
    for _ in range(n_rep):
        p = rng.permutation(ue)
        a = fit_filter(X[np.isin(eid, p[::2])], y[np.isin(eid, p[::2])], ridge)
        b = fit_filter(X[np.isin(eid, p[1::2])], y[np.isin(eid, p[1::2])], ridge)
        if np.std(a) > 0 and np.std(b) > 0:
            rs.append(np.corrcoef(a, b)[0, 1])
    if not rs:
        return np.nan
    r = float(np.median(rs))
    return 2 * r / (1 + r) if r > -1 else np.nan


def integration_index(h):
    """DC gain normalised by filter power, in [-1, 1].

    1 is a flat, purely integrating filter; 0 is perfectly balanced biphasic,
    passing no DC. Monotonic in lobe balance, which a peak-to-trough ratio is
    not -- once the second lobe outgrows the first that ratio starts falling
    again and stops measuring what it is supposed to.
    """
    h = np.asarray(h, dtype=float)
    n = np.linalg.norm(h)
    return float(np.sum(h) / (n * np.sqrt(len(h)))) if n > 0 else np.nan


def analyze_series(ID, channel, lag_s=10.0, roi_prefix="aligned", band=None,
                   ridge=None, seed=0, correct_bleach="auto",
                   null_shifts=nullcal.DEFAULT_SHIFTS,
                   peak_window_s=nullcal.DEFAULT_PEAK_WINDOW_S):
    """Everything, in the order the checks have to happen.

    Returns a dict. `usable` is False when the cross-validated r fails to clear
    the threshold, and in that case the filter is still returned but should not
    be interpreted.
    """
    X, y, eid, dt, meta = build_design(ID, channel, lag_s=lag_s,
                                       correct_bleach=correct_bleach,
                                       roi_prefix=roi_prefix, band=band)
    if ridge is None:
        ridge, ridge_scores = choose_ridge(X, y, eid, seed=seed)
    else:
        ridge_scores = None

    cv_r = cross_validated_r(X, y, eid, ridge=ridge, seed=seed)
    h = fit_filter(X, y, ridge)
    lags = np.arange(len(h)) * dt

    peak = int(np.argmax(np.abs(h)))

    # Circular-shift null. With one averaged filter the question is not which
    # ROIs to keep but whether THIS filter is distinguishable from no coupling,
    # so the output is a p-value rather than a selection. The null runs the
    # identical fit -- same ridge, same design -- on a shifted stimulus.
    null_p, null_thr, null_z = np.nan, np.nan, np.nan
    if null_shifts and null_shifts > 0:
        # PER-EPOCH duration, not len(y)*dt. The shift is applied inside
        # build_design to one epoch's stimulus, and np.roll wraps modulo that
        # length -- so a shift derived from the pooled design (4 epochs, ~1200 s)
        # would wrap back to within a filter length of zero and leave real
        # filter in the null.
        _rp = ID.getRunParameters()
        dur = (float(_rp.get('pre_time', 0.0)) + float(_rp.get('stim_time', 0.0))
               + float(_rp.get('tail_time', 0.0)))
        try:
            # Lags searched for the peak, filter and null alike; see
            # nullcal.DEFAULT_PEAK_WINDOW_S.
            _pl = nullcal.peak_lags_for(peak_window_s, dt, len(h))
            draws = []
            for sh in nullcal.shifts_for(dur, lag_s, null_shifts):
                Xn, yn, en, _, _ = build_design(
                    ID, channel, lag_s=lag_s, roi_prefix=roi_prefix, band=band,
                    correct_bleach=correct_bleach, shift_s=sh)
                draws.append(nullcal.peak_z(fit_filter(Xn, yn, ridge)[None, :],
                                            peak_lags=_pl))
            null = np.concatenate(draws)
            null_z = float(nullcal.peak_z(h[None, :], peak_lags=_pl)[0])
            null_p = float(nullcal.pvalues(np.array([null_z]), null)[0][0])
            null_thr = float(np.percentile(null[np.isfinite(null)], 95))
        except ValueError:
            pass   # too short to place shifts clear of the filter

    return dict(
        null_p=null_p, null_threshold=null_thr, null_peak_z=null_z,
        null_shifts=int(null_shifts or 0),
        filter=h, lags=lags, dt=dt, ridge=ridge, ridge_scores=ridge_scores,
        cv_r=cv_r,
        split_half_r=split_half_filter_r(X, y, eid, ridge=ridge, seed=seed),
        n_epochs=int(len(np.unique(eid))), n_samples=int(len(y)),
        hold_taus=tuple(meta.get("hold_taus", ())),
        independent_samples=float(len(y) * dt / max(meta.get("hold_taus", [1.0]))),
        peak_lag_s=float(lags[peak]), peak_signed=float(h[peak]),
        gain=float(np.linalg.norm(h)),
        integration_index=integration_index(h),
        usable=bool(np.isfinite(cv_r) and cv_r > 0.1),
        meta=meta,
    )


# ─────────────────────────────────────────────────────────────────────────────
# per-ROI path
#
# build_design above averages the ROIs and answers "what is the filter here".
# The population question -- which ROIs are driven, and how their filters
# differ -- needs one filter per ROI, and three things build_design does not do:
#
#   * a separate stimulus grid per Z PLANE. Every ROI sees the same stimulus at
#     the same instant -- what differs is WHEN each plane is sampled. A volume
#     is acquired one plane at a time, so time_vector_by_epoch, which timestamps
#     the volume, is the time of the FIRST plane; the sample of an ROI in plane
#     z actually integrated the stimulus around that time plus
#     getVolumeFrameOffsets()[z]. Ignoring it applies plane 0's stimulus
#     alignment to every plane, which at 10-14 Hz is a sizeable fraction of a
#     lag bin and enters the answer as a per-plane latency bias. The octave path
#     already carries the same offset (roi_z_offsets added to the ROI's time
#     vector); a temporal filter is where it matters most.
#   * each epoch's OWN time vector. They carry real acquisition jitter -- tens
#     of ms, half a sample -- and using epoch 0's for all of them puts that
#     straight into the lag axis.
#   * keeping the per-epoch stimulus, so the circular-shift null re-assembles
#     without re-reading the responses.
#
# ROWS ARE STIMULUS TIME ONLY. A sample becomes a row when its whole lag window
# lies inside the stimulus period of its own epoch -- see _row_indices. The
# pre-stimulus grey is not counted toward the window and the tail is dropped
# outright, so no regressor tap is ever structurally zero and no response
# sample is a decay with nothing driving it.
#
# FOLDS ARE CONTIGUOUS BLOCKS, NOT EPOCHS. The protocol ships one 300 s trial,
# so a fold per epoch would leave nothing to hold out. Blocks of block_s give
# folds whatever the trial count, at the cost of a boundary: a test row near the
# start of its block has a lag window reaching back into the previous block,
# which may be training data. Those rows are dropped from the SCORE, not from
# the fit -- the leak is in what is being scored.
#
# When every epoch shares one seed (repeat_stimulus) the block index alone keys
# the fold, so the same stimulus never appears in both train and test.
# ─────────────────────────────────────────────────────────────────────────────


def _normal_equations(X, Y):
    """Centerd (X'X, X'Y) plus the training mean of X, which prediction needs."""
    xm = X.mean(axis=0)
    ym = Y.mean(axis=0)
    Xc = X - xm
    return Xc.T @ Xc, Xc.T @ (Y - ym), xm


def _ridge_solve(A, B, ridge):
    """Solve (A + lam I) h = B with lam relative to trace(A)/n, so `ridge`
    means the same thing regardless of stimulus contrast or units."""
    lam = float(ridge) * np.trace(A) / A.shape[0]
    return np.linalg.solve(A + lam * np.eye(A.shape[0]), B)


def build_roi_design(ID, channel, filter_length_s, roi_prefix="aligned",
                     band=None, correct_bleach="auto", block_s=60.0,
                     use_z_offsets=True):
    """Responses, per-epoch stimulus and fold structure for one series/channel.

    Returns a dict. The design matrix itself is assembled on demand by
    design_matrix(), because it differs per z plane and the null needs it again
    with the stimulus rolled.

    dF/F is taken against the per-ROI MEDIAN, not the pre-period: two seconds of
    baseline in front of a 300 s epoch is a noisy divisor, and the median over
    the epoch is the same quantity measured 150x better. Same convention as
    build_design.
    """
    rd = bleach.roi_responses(ID, channel, roi_prefix=roi_prefix, dff="none",
                              correct_bleach=correct_bleach)
    er = np.asarray(rd["epoch_response"], dtype=float)
    tvs = list(rd["time_vector_by_epoch"])
    eps = ID.getEpochParameters()
    pre_times = np.asarray(ID.getEpochParameters("pre_time"), dtype=float)
    stim_times = np.asarray(ID.getEpochParameters("stim_time"), dtype=float)

    n_roi = er.shape[0]
    n_ep = int(min(er.shape[1], len(eps), len(tvs), len(pre_times)))
    if n_ep == 0:
        raise ValueError("no epochs with both responses and parameters")

    dt = float(ID.getResponseTiming()["sample_period"])
    n_lags = max(1, int(round(float(filter_length_s) / dt)))

    # Baseline PER EPOCH, from the same array as er. bleach.correct flattens
    # each epoch but keeps its STARTING level (F / (env / env[0])), so epochs
    # recorded later, after more bleaching, sit lower. One baseline shared by
    # every epoch then scales each epoch's dF/F by its level -- the declining-
    # gain artifact the bleach correction exists to remove, reintroduced
    # between epochs. On fly_004 (JS256_x_JS261) ch1's epoch levels drifted
    # 1.04 -> 0.96, about 5 points of an apparent 20% trial-gain decline.
    # The per-epoch median is also the better baseline on its own terms: a
    # 300 s epoch of zero-mean noise response measures it far better than the
    # 5 s pre period. roi_response is not used either way: it is the
    # continuous series trace and, after a correction, on a different scale.
    with np.errstate(invalid="ignore"):
        f0 = np.nanmedian(er[:, :n_ep, :], axis=2)             # (n_roi, n_ep)
    f0 = np.where(np.isfinite(f0) & (f0 != 0), f0, np.nan)

    off = np.zeros(n_roi)
    if use_z_offsets:
        off = np.asarray(rd.get("roi_z_offsets", np.zeros(n_roi)), dtype=float)
    # Voxel-weighted offsets differ ROI to ROI, and each distinct value costs
    # its own stimulus reconstruction and fit. 1 ms is far below the 8.3 ms
    # display frame, so grouping at that resolution changes nothing real
    # and caps the groups at one per millisecond of volume.
    z_offsets, roi_group = np.unique(np.round(off, 3), return_inverse=True)

    metas = [resolve_fullfield_epoch_meta(eps[i]) for i in range(n_ep)]
    seeds = sorted({int(m["seed"]) for m in metas})
    repeated = len(seeds) == 1 and n_ep > 1

    stim, resp, tvecs = [], [], []
    for i in range(n_ep):
        tv = np.asarray(tvs[i], dtype=float)
        tv = tv[:int(min(len(tv), er.shape[2]))]
        tvecs.append(tv)
        with np.errstate(invalid="ignore", divide="ignore"):
            resp.append((er[:, i, :len(tv)] - f0[:, i, None]) / f0[:, i, None])
        # tv times the VOLUME (i.e. its first plane); an ROI in plane z was
        # sampled at tv + offset_z. Binning the stimulus on that shifted grid is
        # how the offset is carried -- the stimulus itself is not shifted.
        stim.append(np.stack([
            stimulus_on_imaging_grid(metas[i], tv + o, float(pre_times[i]),
                                     band=band)
            for o in z_offsets]))

    design = dict(stim=stim, resp=resp, tvecs=tvecs, metas=metas,
                  pre_times=pre_times[:n_ep], stim_times=stim_times[:n_ep],
                  dt=dt, n_lags=n_lags,
                  filter_length_s=n_lags * dt, n_roi=n_roi, n_ep=n_ep,
                  roi_group=roi_group, z_offsets=z_offsets,
                  block_s=float(block_s), repeated=repeated, seeds=seeds,
                  bleach_info=rd.get("bleach_info", {}), band=band,
                  n_bands=len(metas[0].get("hold_taus", (0.45,))))

    design["rows"] = _row_indices(design)
    Y, fold, row_epoch, t_in_block = _assemble_response(design)
    ok = np.all(np.isfinite(Y), axis=1)
    design.update(Y=Y[ok], fold=fold[ok], row_epoch=row_epoch[ok],
                  t_in_block=t_in_block[ok], row_ok=ok,
                  n_rows_dropped=int((~ok).sum()))
    return design


def _assemble_response(design):
    """Rows of (response, fold, epoch, time-since-block-start), epoch-major.

    Rows come from design['rows'], so this and design_matrix cannot disagree
    about which samples exist.
    """
    Y, fold, epoch, tib = [], [], [], []
    base = 0
    for i, k in enumerate(design["rows"]):
        if k.size == 0:
            continue
        Y.append(design["resp"][i][:, k].T)
        t = design["tvecs"][i][k] - design["pre_times"][i]
        b = np.floor(t / design["block_s"]).astype(int)
        fold.append(b if design["repeated"] else b + base)
        base += int(b.max()) + 1
        epoch.append(np.full(k.size, i))
        tib.append(t - b * design["block_s"])
    if not Y:
        raise ValueError(
            "no epoch has a stimulus period longer than the filter length")
    return (np.concatenate(Y, axis=0), np.concatenate(fold),
            np.concatenate(epoch), np.concatenate(tib))


def _row_indices(design):
    """Which imaging samples become rows, per epoch. STIMULUS TIME ONLY.

    A sample qualifies when its whole lag window lies inside the stimulus
    period of its own epoch:

        filter_length_s <= t - pre_time <= stim_time

    Both ends are deliberate. The lower end starts the analysis one filter
    length into the STIMULUS, not into the epoch, so no regressor tap ever
    falls in the pre-stimulus grey -- counting pre_time toward the window would
    let the first rows carry taps that are structurally zero. The upper end
    drops the tail entirely: the screen is grey there, so those rows are a
    response to nothing plus a lag window that is part stimulus and part not,
    and they contribute a decaying response with no matching regressor.

    That also makes the previous trial unreachable by construction, with the
    whole of tail_time plus the inter-trial interval to spare.
    """
    out = []
    for i, tv in enumerate(design["tvecs"]):
        t = np.asarray(tv, dtype=float) - float(design["pre_times"][i])
        k = np.flatnonzero((t >= design["filter_length_s"])
                           & (t <= float(design["stim_times"][i])))
        # index safety, independent of the timing above
        out.append(k[k >= design["n_lags"]])
    return out


def design_matrix(design, group=0, shift_s=0.0):
    """(rows, n_lags) stimulus history for one z plane.

    Column j is the stimulus j samples in the past, so column 0 is the sample
    the response frame itself integrated and a fitted filter reads in causal lag
    order with lag[j] = j * dt.

    shift_s rolls each epoch's stimulus circularly, which is the null: it
    destroys the stimulus-response pairing while preserving both the stimulus
    autocorrelation and the response's own structure.
    """
    n_lags = design["n_lags"]
    k_shift = int(round(float(shift_s) / design["dt"]))
    lags = np.arange(n_lags)
    parts = []
    for i, k in enumerate(design["rows"]):
        if k.size == 0:
            continue
        s = design["stim"][i][group]
        if k_shift:
            s = np.roll(s, k_shift)
        parts.append(s[k[:, None] - lags[None, :]])
    return np.concatenate(parts, axis=0)[design["row_ok"]]


def _group_rois(design, group):
    return np.flatnonzero(design["roi_group"] == group)


def fit_rois(design, ridge=1e-3, shift_s=0.0):
    """(n_roi, n_lags) filters, fitted per z plane because the design differs.

    Units: dF/F per unit displayed contrast, per sample. Divide by dt to read it
    as a continuous kernel; sum it for the DC gain.
    """
    h = np.full((design["n_roi"], design["n_lags"]), np.nan)
    for g in range(len(design["z_offsets"])):
        rois = _group_rois(design, group=g)
        if rois.size == 0:
            continue
        X = design_matrix(design, g, shift_s)
        A, B, _ = _normal_equations(X, design["Y"][:, rois])
        h[rois] = _ridge_solve(A, B, ridge).T
    return h


def _block_folds(fold, n_fold, seed=0):
    """Group the contiguous blocks into n_fold folds, at random."""
    ub = np.unique(fold)
    rng = np.random.default_rng(seed)
    return np.array_split(rng.permutation(ub), int(min(n_fold, len(ub))))


def ridge_sweep(design, grid=None, n_fold=6, seed=0, shift_s=0.0):
    """Cross-validated r per ROI at every ridge on the grid.

    shift_s rolls the stimulus exactly as design_matrix does, which turns
    this into the null for predictive power: same folds, same ridge, same
    estimator, pairing destroyed. The fold assignment is seeded, so a
    shifted call is scored on the same held-out rows as the unshifted one
    and the two are directly comparable.

    Returns (grid, scores) with scores (n_grid, n_roi). ONE ridge is then chosen
    for the whole series -- picking it per ROI would let each ROI choose the
    smoothing that flatters it, and the differences between ROIs are the point.

    The grid spans six decades because the right value depends on how correlated
    the stimulus is, which depends on the hold time: a value tuned at one tau is
    wrong at another.
    """
    if grid is None:
        # Up to 1e2, not 1e0: on synthetic data at realistic noise the choice
        # ran into a 1e0 ceiling, and a ridge pinned at the edge of its grid is
        # not a chosen value. The caller should say so if it lands there anyway.
        grid = np.logspace(-6, 2, 17)
    grid = np.asarray(grid, dtype=float)
    folds = _block_folds(design["fold"], n_fold, seed)
    scores = np.full((len(grid), design["n_roi"]), np.nan)

    for g in range(len(design["z_offsets"])):
        rois = _group_rois(design, group=g)
        if rois.size == 0:
            continue
        X = design_matrix(design, g, shift_s)
        Yg = design["Y"][:, rois]
        # Pearson r over the pooled held-out rows, accumulated fold by fold
        # rather than by keeping every prediction: (n_grid, n_roi, rows) would
        # be tens of MB per group for nothing.
        acc = np.zeros((len(grid), len(rois), 5))   # Sp, Sp2, Spy, Sy, Sy2
        n_sc = 0
        for f in folds:
            te = np.isin(design["fold"], f)
            if te.all() or not te.any():
                continue
            # Score only rows whose whole lag window lies inside the test
            # block; the rest reach back into training data.
            sc = te & (design["t_in_block"] >= design["filter_length_s"])
            if not sc.any():
                continue
            A, B, xm = _normal_equations(X[~te], Yg[~te])
            # One eigendecomposition serves the whole grid: A is fixed by the
            # fold and only the ridge changes.
            w, V = np.linalg.eigh(A)
            # A is PSD by construction; eigh can still return a tiny negative
            # for a collinear direction, and at the small end of the grid that
            # would make the shrinkage denominator cross zero.
            w = np.maximum(w, 0.0)
            VB = V.T @ B
            Xs = (X[sc] - xm) @ V
            y = Yg[sc]
            n_sc += int(sc.sum())
            for gi, lam in enumerate(grid):
                p = Xs @ (VB / (w + lam * w.sum() / len(w))[:, None])
                acc[gi, :, 0] += p.sum(axis=0)
                acc[gi, :, 1] += (p ** 2).sum(axis=0)
                acc[gi, :, 2] += (p * y).sum(axis=0)
                acc[gi, :, 3] += y.sum(axis=0)
                acc[gi, :, 4] += (y ** 2).sum(axis=0)
        if n_sc < 10:
            continue
        sp, sp2, spy, sy, sy2 = (acc[..., i] for i in range(5))
        den = np.sqrt((n_sc * sp2 - sp ** 2) * (n_sc * sy2 - sy ** 2))
        with np.errstate(invalid="ignore", divide="ignore"):
            r = (n_sc * spy - sp * sy) / den
        scores[:, rois] = np.where(den > 0, r, np.nan)
    return grid, scores


def split_half_rois(design, ridge=1e-3, n_rep=50, seed=0):
    """Per-ROI reliability of the FILTER, splitting blocks.

    Spearman-Brown corrected to the full dataset, so it is comparable to the
    cross-validated r rather than describing half the data.

    Every repeat reuses per-block normal equations, so the cost is one pass over
    the design however many repeats are asked for -- a split's Gram matrix is
    the sum of its blocks'.
    """
    rng = np.random.default_rng(seed)
    out = np.full((n_rep, design["n_roi"]), np.nan)
    ub = np.unique(design["fold"])

    for g in range(len(design["z_offsets"])):
        rois = _group_rois(design, group=g)
        if rois.size == 0:
            continue
        X = design_matrix(design, g)
        Yg = design["Y"][:, rois]
        # per block: raw cross-products and sums, from which any subset's
        # centerd normal equations follow exactly
        S = {}
        for b in ub:
            m = design["fold"] == b
            Xb, Yb = X[m], Yg[m]
            S[b] = (Xb.T @ Xb, Xb.T @ Yb, Xb.sum(axis=0), Yb.sum(axis=0),
                    int(m.sum()))

        def _fit(blocks):
            sxx = sum(S[b][0] for b in blocks)
            sxy = sum(S[b][1] for b in blocks)
            sx = sum(S[b][2] for b in blocks)
            sy = sum(S[b][3] for b in blocks)
            n = sum(S[b][4] for b in blocks)
            if n < X.shape[1] + 2:
                return None
            xm, ym = sx / n, sy / n
            return _ridge_solve(sxx - n * np.outer(xm, xm),
                                sxy - n * np.outer(xm, ym), ridge)

        for rep in range(n_rep):
            p = rng.permutation(ub)
            a, b = _fit(p[::2]), _fit(p[1::2])
            if a is None or b is None:
                continue
            for j, r in enumerate(rois):
                if np.std(a[:, j]) > 0 and np.std(b[:, j]) > 0:
                    out[rep, r] = np.corrcoef(a[:, j], b[:, j])[0, 1]

    with np.errstate(invalid="ignore"):
        r = np.nanmedian(out, axis=0)
        return np.where(r > -1, 2 * r / (1 + r), np.nan)


def filter_spectrum(h, dt, pad=8):
    """(freqs, |H|) for one filter or a stack of them.

    Zero-padded, which INTERPOLATES the transform rather than adding
    information: a filter that has decayed to zero inside the window has a
    compactly supported transform, so H(f) is already determined at every f.
    Check bleach.returns_to_zero before believing the low end.
    """
    h = np.asarray(h, dtype=float)
    n = int(h.shape[-1] * int(pad))
    return (np.fft.rfftfreq(n, d=float(dt)),
            np.abs(np.fft.rfft(np.nan_to_num(h), n=n, axis=-1)))


def band_3db(f, mag):
    """(f_lo, f_hi, f_peak) of the CONTIGUOUS -3 dB band around the peak.

    Contiguous, not "every frequency above the threshold": a noisy tail that
    happens to cross it does not extend the band the filter actually passes.
    """
    mag = np.asarray(mag, dtype=float)
    if not np.any(np.isfinite(mag)) or np.nanmax(mag) <= 0:
        return np.nan, np.nan, np.nan
    pk = int(np.nanargmax(mag))
    thr = mag[pk] / np.sqrt(2.0)
    lo = pk
    while lo > 0 and mag[lo - 1] >= thr:
        lo -= 1
    hi = pk
    while hi < len(mag) - 1 and mag[hi + 1] >= thr:
        hi += 1
    return float(f[lo]), float(f[hi]), float(f[pk])


def repeat_groups(design):
    """Epochs that showed frame-for-frame the same stimulus, grouped.

    Keyed on the WHOLE fullfield_stimulus metadata, not the seed alone: the
    same seed under a different hold time, p or clip is a different stimulus.
    Returns a list of epoch-index arrays, one per stimulus shown more than once,
    largest first. Empty when nothing repeats.
    """
    groups = {}
    for i, m in enumerate(design["metas"]):
        groups.setdefault(json.dumps(m, sort_keys=True, default=str), []).append(i)
    reps = [np.array(v) for v in groups.values() if len(v) > 1]
    return sorted(reps, key=len, reverse=True)


def fit_rois_by_epoch(design, ridge):
    """(n_ep, n_roi, n_lags): one filter per epoch, same lags and ridge as the
    full fit so trials are comparable to it and to each other. A single epoch
    is a fraction of the data, so these are noisier than the full filter."""
    out = np.full((design["n_ep"], design["n_roi"], design["n_lags"]), np.nan)
    for g in range(len(design["z_offsets"])):
        rois = _group_rois(design, group=g)
        if rois.size == 0:
            continue
        X = design_matrix(design, g)
        for i in range(design["n_ep"]):
            m = design["row_epoch"] == i
            if m.sum() < X.shape[1] + 2:
                continue
            A, B, _ = _normal_equations(X[m], design["Y"][m][:, rois])
            out[i][rois] = _ridge_solve(A, B, ridge).T
    return out


def _corr(a, b):
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def _filter_support(h_all, z=3.0, noise_frac=0.3, min_lags=3):
    """Lags up to the last one where the full filter clears its own noise.

    Correlating trials over the whole window would score mostly the tail,
    where every trial is noise and two good trials look no more alike than two
    bad ones. The noise scale is the longest lags' SD, as in nullcal.peak_z.
    An energy fraction was tried first and is fragile: a noisy tail accumulates
    energy over many lags, so a slightly less regularised fit of the same cell
    doubled the support and halved the trial agreement. When no lag clears the
    noise, the whole window is used -- no selection is better than a selection
    made from noise.
    """
    n = len(h_all)
    sd = np.nanstd(h_all[int(round((1.0 - noise_frac) * n)):])
    if not (np.isfinite(sd) and sd > 0):
        return n
    above = np.flatnonzero(np.abs(h_all) / sd > z)
    if not above.size:
        return n
    return min(max(int(above[-1]) + 2, min_lags), n)


def trial_reliability(design, ridge, h_all):
    """Per-trial filters and their consistency, for each repeated stimulus.

    Only meaningful when trials saw the SAME stimulus: then trial-to-trial
    differences are the cell's own, not a different stimulus sampling the
    filter differently. Returns None when nothing repeats.

    Per repeat group, per ROI:
      trial_filters      (n_trial, n_roi, n_lags)
      filter_loo_r       corr of each trial's filter with the mean of the
                         others, over support_lags: up to the last lag where
                         the full filter exceeds 3x its late-lag noise
      trial_gain         least-squares scale of each trial's filter onto the
                         full filter: 1 = typical, a trend = run-down or drift
      response_pair_r    mean correlation between single-trial RESPONSES on
                         the common stimulus-time grid -- model-free
      noise_ceiling      sqrt(response_pair_r), NaN unless response_pair_r is
                         more than 2 SE above zero: the best cv_r ANY model could
                         reach against single-trial data, since cv_r is scored
                         on single trials. cv_r / noise_ceiling is the share of
                         the reliable response the linear filter explains.
    """
    groups = repeat_groups(design)
    if not groups:
        return None
    h_ep = fit_rois_by_epoch(design, ridge)
    n_roi, dt = design["n_roi"], design["dt"]
    out = []
    for eps in groups:
        tf = h_ep[eps]
        n_tr = len(eps)
        support = np.array([_filter_support(h_all[r]) for r in range(n_roi)])
        loo = np.full((n_tr, n_roi), np.nan)
        gain = np.full((n_tr, n_roi), np.nan)
        for r in range(n_roi):
            L = support[r]
            ref = h_all[r]
            den = float(np.dot(ref, ref))
            for j in range(n_tr):
                others = np.nanmean(np.delete(tf[:, r, :], j, axis=0), axis=0)
                loo[j, r] = _corr(tf[j, r, :L], others[:L])
                if den > 0 and np.all(np.isfinite(tf[j, r])):
                    gain[j, r] = float(np.dot(tf[j, r], ref) / den)

        # Model-free response repeatability. Epochs start at arbitrary phase
        # relative to the scan, so samples are interpolated onto one
        # stimulus-time grid before trials are compared.
        t_lo = max(float(design["tvecs"][i][design["rows"][i]][0]
                         - design["pre_times"][i]) for i in eps
                   if design["rows"][i].size)
        t_hi = min(float(design["tvecs"][i][design["rows"][i]][-1]
                         - design["pre_times"][i]) for i in eps
                   if design["rows"][i].size)
        grid = np.arange(t_lo, t_hi, dt)
        R = np.full((n_tr, n_roi, len(grid)), np.nan)
        for j, i in enumerate(eps):
            k = design["rows"][i]
            if not k.size:
                continue
            t = design["tvecs"][i][k] - design["pre_times"][i]
            for r in range(n_roi):
                y = design["resp"][i][r, k]
                if np.all(np.isfinite(y)):
                    R[j, r] = np.interp(grid, t, y)
        pair = np.full(n_roi, np.nan)
        pair_se = np.full(n_roi, np.nan)
        resp_loo = np.full((n_tr, n_roi), np.nan)
        for r in range(n_roi):
            rs = np.array([_corr(R[a, r], R[b, r]) for a in range(n_tr)
                           for b in range(a + 1, n_tr)])
            rs = rs[np.isfinite(rs)]
            if rs.size:
                pair[r] = rs.mean()
                # Spread across pairs, not a formula in the sample count: the
                # samples are autocorrelated, so their nominal number overstates
                # how many independent comparisons each correlation rests on.
                pair_se[r] = (rs.std(ddof=1) / np.sqrt(rs.size)
                              if rs.size > 1 else np.nan)
            for j in range(n_tr):
                resp_loo[j, r] = _corr(
                    R[j, r], np.nanmean(np.delete(R[:, r], j, axis=0), axis=0))
        out.append(dict(
            epochs=eps, seed=int(design["metas"][eps[0]]["seed"]),
            trial_filters=tf, support_lags=support, filter_loo_r=loo,
            trial_gain=gain, response_pair_r=pair, response_pair_r_se=pair_se,
            response_loo_r=resp_loo,
            # NaN unless the repeatability is resolved (> 2 SE above zero). A
            # weak ROI's signal share per trial is ~cv_r**2 -- 0.001 at cv_r
            # 0.03 -- far below what the correlation of two noisy traces can
            # resolve, and a ceiling read from noise can come out at zero or
            # below cv_r, which is not a statement about the cell.
            noise_ceiling=np.where(pair > 2 * pair_se,
                                   np.sqrt(np.clip(pair, 0.0, None)), np.nan)))
    return out
