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
sampled at its centre -- which is what the indicator and the scan physically do.
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
    "resolve_fullfield_epoch_meta",
    "stimulus_on_imaging_grid",
    "build_design",
    "fit_filter",
    "choose_ridge",
    "cross_validated_r",
    "split_half_filter_r",
    "analyze_series",
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
    """
    gen = _generator(meta)
    tv = np.asarray(time_vector, dtype=float)
    dt = float(np.median(np.diff(tv)))
    t_stim = tv - float(pre_time)

    out = np.zeros(len(tv), dtype=float)
    for i, c in enumerate(t_stim):
        lo, hi = c - dt / 2.0, c + dt / 2.0
        if hi <= 0.0:
            continue                       # entirely within pre_time
        lo = max(lo, 0.0)
        f0 = int(np.floor(lo * gen.fps))
        f1 = int(np.ceil(hi * gen.fps))
        frames = np.arange(f0, max(f1, f0 + 1))
        vals = (gen.values(frames) if band is None
                else gen.band_values(int(band), frames).astype(float))
        out[i] = float(np.mean(vals))
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

    # f0 must come from the SAME array as er. bleach.roi_responses corrects
    # epoch_response only -- roi_response is the continuous series trace and
    # would carry a different envelope -- so taking f0 from roi_response after
    # a correction mixes two scales and puts the difference straight into dF/F.
    if rd.get("bleach_info", {}).get("applied"):
        f0 = np.nanmedian(er[:, :n_ep, :].reshape(er.shape[0], -1),
                          axis=1)[:, None, None]
    else:
        f0 = np.nanmedian(np.asarray(rd["roi_response"]), axis=1)[:, None, None]
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
                   null_shifts=nullcal.DEFAULT_SHIFTS):
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
            draws = []
            for sh in nullcal.shifts_for(dur, lag_s, null_shifts):
                Xn, yn, en, _, _ = build_design(
                    ID, channel, lag_s=lag_s, roi_prefix=roi_prefix, band=band,
                    correct_bleach=correct_bleach, shift_s=sh)
                draws.append(nullcal.peak_z(fit_filter(Xn, yn, ridge)[None, :]))
            null = np.concatenate(draws)
            null_z = float(nullcal.peak_z(h[None, :])[0])
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
