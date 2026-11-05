"""Temporal filters under switching contrast (labpack ContrastBlockUniformNoise).

WHAT THE STIMULUS IS. A full-field patch whose intensity is redrawn every
1/update_rate s from a uniform discrete distribution, in epochs of ~5 s that
alternate between a narrow range (low contrast) and a wide one (high contrast)
around the same mean. The epochs run back to back with ~0.1 s of idle grey
between them, so the fly sees one continuous stimulus in which only the
VARIANCE switches. That is the classic contrast-adaptation design, and the
question it asks is what the switch does to the linear filter: its gain, its
kinetics, and how fast each one moves.

WHY CONTINUOUS, NOT PER EPOCH. An epoch is shorter than two calcium filters, so
cutting the recording into epochs and fitting inside each one (as the full-field
ternary pipeline does) would discard the first filter length of every block --
which is exactly where adaptation happens. Here the whole series is ONE
continuous stimulus trace, rebuilt from the per-epoch intensity sequence the
protocol recorded, and every imaging sample is a row whose lag window may cross
a switch. Which rows go into which fit is then a matter of labelling rows by
the contrast they sit in and the time since the last switch.

RECONSTRUCTION. Every epoch carries its displayed sequence verbatim in its
`color` parameter (a TVPairs trajectory, blue channel), so nothing is
regenerated and there is no dependence on scipy's rv_discrete staying
bit-compatible with the rig. What the parameters cannot supply is WHEN each
epoch started; that comes from the photodiode.

EPOCH ONSETS. getStimulusTiming separates epochs by a silent gap of
0.9 * (pre_time + tail_time) -- 60 ms for this protocol, against real gaps of
~85-140 ms -- and it does occasionally miss one, merging two epochs into a
~10 s "epoch" (fly_013 series_003: 299 onsets for 300 epochs). Pairing onsets
with epoch parameters BY INDEX then mislabels every epoch after the merge,
swapping low and high contrast for the rest of the run. Onsets are therefore
matched to epochs through the recorded epoch_unix_time, and an onset the
photodiode missed is placed from its nearest detected neighbour plus the unix
interval between them.

CONTRAST UNITS. The regressor is Weber contrast about the idle (mean) level,
(I - idle) / idle, so the two conditions share units and a filter's amplitude is
dF/F per unit contrast in both. With full contrast normalisation the low-
contrast filter would be larger by the ratio of the stimulus SDs.

SUB-FRAME TIME RESOLUTION FROM JITTER. Volumes arrive every ~72 ms and the
stimulus changes every 33 ms, with no fixed phase between them, and block
onsets add their own jitter. Every imaging sample therefore sees the stimulus
at a different phase, and a filter can be estimated on a lag grid much finer
than the volume period -- the same trick as equivalent-time sampling. It also
matters HOW each sample integrated the stimulus: a z plane is scanned in one
frame period (~5.5 ms), not over the whole volume, so each sample is close to
a snapshot at its own time. build_design(lag_dt=..., sample_window_s=...)
sets the lag grid and the integration window; each regressor is the stimulus
averaged over [t - lag, t - lag + window], with t the start of that ROI's
plane. lag_dt = window = sample_period reproduces the old volume-binned
design. The calcium indicator still low-passes the answer; this removes the
extra blur the analysis itself was adding. switch_triggered_average(bin_s=...)
uses the same jitter for the switch-aligned responses.

LN CONVENTION (Baccus & Meister 2002, Neuron 36:909, Experimental
Procedures eq. 4). Each condition's filter is SCALED so
that the filtered stimulus has the same variance as the stimulus itself:
var(F * s) = var(s), with s the contrast the model sees (binned to the imaging
sample). The filter then carries only temporal SHAPE -- kinetics, biphasic
lobe, latency -- and has no gain of its own; all gain lives in the static
nonlinearity, whose slope is the gain in dF/F per unit contrast. The generator
signal of a high-contrast block spans ~4x the range of a low-contrast one, as
the stimulus does, so the two nonlinearities can be drawn on one axis and
compared directly. normalize_filters / normalize_folds do the scaling;
ln_model builds the nonlinearities.

GAIN CHANGE IS NOT ENOUGH ON ITS OWN. A compressive static nonlinearity also
makes the linear gain fall at high contrast, with no adaptation at all. Two
things separate them, and both are computed here: (1) a static nonlinearity acts
instantaneously, so the gain would be at its steady-state value in the first
bin after a switch, while adaptation moves it over time; (2) over the range of
generator signal BOTH conditions cover, a static nonlinearity has one slope,
whereas adaptation gives two.
"""
import ast

import numpy as np

from visanalysis.util import bleach

__all__ = [
    "is_contrast_blocks",
    "is_repeated",
    "epoch_onsets",
    "epoch_conditions",
    "stimulus_trace",
    "build_design",
    "design_matrix",
    "fit_filters",
    "cv_filters",
    "cv_r",
    "choose_ridge",
    "cycle_shifts",
    "gain_after_switch",
    "nonlinearity",
    "normalize_filters",
    "penalty_matrix",
    "normalize_folds",
    "ln_model",
    "window_filters",
    "switch_triggered_average",
    "predict_series",
    "fit_exp_recovery",
]


# ─────────────────────────────────────────────────────────────────────────────
# stimulus
# ─────────────────────────────────────────────────────────────────────────────

def is_contrast_blocks(run_parameters):
    """True if a series ran ContrastBlockUniformNoise."""
    return "ContrastBlock" in str(run_parameters.get("protocol_ID", ""))


def is_repeated(epoch_parameters):
    """True when every epoch used the same seed, i.e. frozen noise.

    Each condition then shows the identical sequence on every one of its
    epochs, so the series measures reliability rather than a filter.
    """
    seeds = {int(e["start_seed"]) for e in epoch_parameters}
    return len(seeds) == 1 and len(epoch_parameters) > 1


def epoch_onsets(ID, epoch_parameters, match_tol_s=0.25):
    """Photodiode onset of every epoch, matched through epoch_unix_time.

    Returns (onsets, info). onsets[i] is epoch i's first displayed frame on the
    imaging clock. info reports how many onsets were detected, which epochs had
    to be placed because the photodiode merged them into a neighbour, and the
    residual of detected onsets against the unix-time clock map, which is the
    check that the matching itself is sound.
    """
    st = ID.getStimulusTiming()
    det = np.asarray(st["stimulus_start_times"], dtype=float)
    u = np.array([float(e["epoch_unix_time"]) for e in epoch_parameters])
    n = len(u)
    if det.size == 0:
        raise ValueError("photodiode found no epoch onsets")

    # Clock map photodiode <- unix, refined once: anchor on the first onset,
    # match, fit slope + offset, rematch. The slope absorbs the ~1e-5 drift
    # between the two clocks, which over a 25 min run is ~15 ms.
    slope, off = 1.0, det[0] - u[0]
    for _ in range(2):
        pred = slope * u + off
        idx = np.array([int(np.argmin(np.abs(pred - d))) for d in det])
        ok = np.abs(det - pred[idx]) < match_tol_s
        # one detected onset per epoch: keep the closest if two claim it
        best = {}
        for j in np.flatnonzero(ok):
            e = idx[j]
            if e not in best or (abs(det[j] - pred[e])
                                 < abs(det[best[e]] - pred[e])):
                best[e] = j
        ep = np.array(sorted(best))
        dj = np.array([best[e] for e in ep])
        if len(ep) >= 2:
            slope, off = np.polyfit(u[ep], det[dj], 1)

    onsets = np.full(n, np.nan)
    onsets[ep] = det[dj]
    missing = np.flatnonzero(~np.isfinite(onsets))
    have = np.flatnonzero(np.isfinite(onsets))
    for i in missing:
        # nearest DETECTED neighbour plus the unix interval, rather than the
        # global fit: local latency is better predicted by its neighbour.
        k = have[np.argmin(np.abs(have - i))]
        onsets[i] = onsets[k] + slope * (u[i] - u[k])
    resid = det[dj] - (slope * u[ep] + off)
    info = dict(n_detected=int(det.size), n_epochs=int(n),
                placed_epochs=missing.astype(int),
                unmatched_detections=int(det.size - len(dj)),
                clock_resid_ms=np.percentile(resid * 1e3, [0, 50, 100]),
                n_dropped_frames=int(len(st.get("dropped_frame_times", []))))
    return onsets, info


def epoch_conditions(epoch_parameters, idle):
    """Condition index per epoch, ordered by RMS contrast (0 = lowest).

    A condition is one (color_min, color_max, num_colors) distribution; its RMS
    contrast is computed from the level set itself, so the ordering and the
    reported contrast do not depend on how the protocol names things.
    """
    keys = [(float(e["color_min"]), float(e["color_max"]),
             int(e["num_colors"])) for e in epoch_parameters]
    uniq = sorted(set(keys))
    rms = []
    for k in uniq:
        lv = np.linspace(k[0], k[1], k[2])
        rms.append(float(np.sqrt(np.mean(((lv - idle) / idle) ** 2))))
    order = np.argsort(rms)
    uniq = [uniq[i] for i in order]
    rms = [rms[i] for i in order]
    cond = np.array([uniq.index(k) for k in keys])
    return cond, dict(keys=uniq, rms_contrast=np.array(rms))


def _parse_color(v):
    if isinstance(v, bytes):
        v = v.decode("utf-8")
    if isinstance(v, str):
        v = ast.literal_eval(v)
    if isinstance(v, dict):
        v = v["tv_pairs"]
    return v


def stimulus_trace(epoch_parameters, onsets, t_end, idle, fs=1000.0):
    """Continuous Weber contrast on a 1/fs grid, plus per-sample labels.

    Returns dict(t, c, epoch, cond, since_switch, prev_cond, ...). Between
    epochs the display holds idle, which is contrast 0; that is what the cell
    saw, not a fill. `since_switch` is time since the start of the current run
    of same-condition epochs -- equal to time since epoch onset when conditions
    alternate, which they do in this protocol, but kept distinct so a
    randomised order does not silently reset the clock on a repeat.
    """
    cond, cinfo = epoch_conditions(epoch_parameters, idle)
    t = np.arange(0.0, float(t_end), 1.0 / fs)
    c = np.zeros_like(t)
    epoch = np.full(t.size, -1)
    for i, e in enumerate(epoch_parameters):
        tvp = _parse_color(e["color"])
        tt = np.array([p[0] for p in tvp], dtype=float)
        vv = np.array([p[1][2] for p in tvp], dtype=float)   # blue channel
        k0 = max(int(round(onsets[i] * fs)), 0)
        k1 = min(int(round((onsets[i] + float(e["stim_time"])) * fs)), t.size)
        if k1 <= k0:
            continue
        ts = t[k0:k1] - onsets[i]
        j = np.clip(np.searchsorted(tt, ts, side="right") - 1, 0, None)
        c[k0:k1] = (vv[j] - idle) / idle
        # the epoch owns everything up to the next onset, idle gap included
        kn = (int(round(onsets[i + 1] * fs)) if i + 1 < len(onsets)
              else k1 + int(round(float(e.get("tail_time", 0.0)) * fs)))
        epoch[k0:min(kn, t.size)] = i

    run_onset = np.empty(len(cond))
    prev = np.full(len(cond), -1)
    for i in range(len(cond)):
        if i == 0 or cond[i] != cond[i - 1]:
            run_onset[i] = onsets[i]
            prev[i] = cond[i - 1] if i else -1
        else:
            run_onset[i], prev[i] = run_onset[i - 1], prev[i - 1]
    lab = epoch >= 0
    since = np.full(t.size, np.nan)
    since[lab] = t[lab] - run_onset[epoch[lab]]
    cnd = np.full(t.size, -1)
    cnd[lab] = cond[epoch[lab]]
    pcnd = np.full(t.size, -1)
    pcnd[lab] = prev[epoch[lab]]
    return dict(t=t, c=c, epoch=epoch, cond=cnd, since_switch=since,
                prev_cond=pcnd, fs=float(fs), epoch_cond=cond, **cinfo)


# ─────────────────────────────────────────────────────────────────────────────
# design
# ─────────────────────────────────────────────────────────────────────────────

def build_design(ID, channel, filter_length_s, roi_prefix="aligned",
                 correct_bleach=True, block_s=60.0, fs=1000.0, lag_dt=None,
                 sample_window_s=None, penalty="auto", lag_basis="point"):
    """Continuous response, binned stimulus and row labels for one channel.

    Response: the whole-series ROI trace, bleach-corrected as ONE trace (the
    per-epoch correction the full-field pipeline uses cannot work on 5 s
    epochs), then dF/F against the per-ROI median. The correction high-passes
    below ~0.01 Hz; the block cycle is ~0.1 Hz, so it is untouched.

    Stimulus: averaged within each imaging sample's window on that ROI's z
    plane, as in fullfield_temporal -- the scan integrates, it does not sample.

    lag_dt: spacing of the filter's lag grid (default: the volume period).
    lag_basis: 'point' -- tap j is the stimulus at lag j*lag_dt, seen through
    the sample window; the filter is a comb of taps, fine when lag_dt is small.
    'box' -- tap j is the stimulus AVERAGED over the whole lag bin
    [j, j+1)*lag_dt; the filter is piecewise constant at lag_dt. Use 'box' for
    coarse bins (e.g. 40-50 ms): point taps that far apart misrepresent a
    smooth filter and fold the undetermined 30 Hz component down to low
    frequencies, whereas a box that wide has almost no response at 30 Hz, so
    there is nothing there left to be undetermined.
    sample_window_s: how long each sample integrates the stimulus, starting at
    its plane's frame time (default: the volume period, i.e. the old binned
    design; the physically right value is the frame period of one plane).

    Rows: every imaging sample with a full lag window that falls inside some
    epoch. Cross-validation folds are contiguous block_s chunks of the series,
    each holding several whole contrast cycles.
    """
    eps = ID.getEpochParameters()
    rp = ID.getRunParameters()
    idle = float(rp.get("idle_color", 0.5))
    onsets, onset_info = epoch_onsets(ID, eps)

    rd = ID.getRoiResponses(channel, roi_prefix=roi_prefix, dff="none")
    F = np.asarray(rd["roi_response"], dtype=float)
    timing = ID.getResponseTiming()
    tvol = np.asarray(timing["time_vector"], dtype=float)[:F.shape[1]]
    F = F[:, :tvol.size]
    dt = float(timing["sample_period"])
    n_roi = F.shape[0]

    Y = np.empty_like(F)
    frac = []
    for r in range(n_roi):
        f = F[r]
        if correct_bleach:
            try:
                f, _, inf = bleach.correct(f, dt)
                frac.append(inf["frac_removed"])
            except ValueError:
                pass
        with np.errstate(invalid="ignore", divide="ignore"):
            Y[r] = f / np.nanmedian(f) - 1.0

    stim = stimulus_trace(eps, onsets, tvol[-1] + 1.0, idle, fs=fs)
    cs = np.concatenate([[0.0], np.cumsum(stim["c"])])

    off = np.asarray(rd.get("roi_z_offsets", np.zeros(n_roi)), dtype=float)
    z_offsets, roi_group = np.unique(np.round(off, 3), return_inverse=True)
    lag_dt = dt if lag_dt is None else float(lag_dt)
    window = dt if sample_window_s is None else float(sample_window_s)
    n_lags = max(1, int(round(float(filter_length_s) / lag_dt)))
    span = (n_lags - 1) * lag_dt
    if penalty == "auto":
        penalty = "smooth" if lag_dt < 0.9 * dt else "ridge"
    # Row labels from plane 0; planes differ by < 20 ms, far below any window.
    k0 = np.clip(np.round(tvol * fs).astype(int), 0, stim["c"].size - 1)
    rows = np.arange(tvol.size)
    # the whole lag window must lie after the first onset
    ok = ((stim["epoch"][k0[rows]] >= 0)
          & (tvol[rows] - span - window >= onsets[0])
          & np.all(np.isfinite(Y[:, rows]), 0))
    rows = rows[ok]
    t_rows = tvol[rows]
    fold = np.floor((t_rows - t_rows[0]) / float(block_s)).astype(int)

    n_cond = len(stim["keys"])
    return dict(
        Y=Y[:, rows].T, rows=rows, t=t_rows, tvol=tvol, Y_full=Y,
        c=stim["c"], cs=cs, fs=float(fs), lag_dt=lag_dt, window=window,
        lags=np.arange(n_lags) * lag_dt, _xcache={}, lag_basis=lag_basis,
        frame=0.0 if sample_window_s is None else float(sample_window_s),
        penalty=penalty_matrix(n_lags, penalty), penalty_kind=penalty,
        roi_group=roi_group, z_offsets=z_offsets,
        cond=stim["cond"][k0[rows]], prev_cond=stim["prev_cond"][k0[rows]],
        since_switch=stim["since_switch"][k0[rows]],
        epoch=stim["epoch"][k0[rows]], fold=fold, dt=dt, n_lags=n_lags,
        filter_length_s=n_lags * lag_dt, n_roi=n_roi, onsets=onsets,
        onset_info=onset_info, epoch_cond=stim["epoch_cond"],
        cond_keys=stim["keys"], rms_contrast=stim["rms_contrast"],
        stim_sd=np.array([stim["c"][stim["cond"] == q].std()
                          for q in range(n_cond)]),
        idle=idle, repeated=is_repeated(eps),
        bleach_frac_removed=float(np.mean(frac)) if frac else np.nan)


def _full_design(design, group, shift):
    """(all rows, n_lags) for one plane, from each row's EXACT time."""
    key = (int(group), int(shift))
    cache = design["_xcache"]
    if key in cache:
        return cache[key]
    fs, w = design["fs"], design["window"]
    if shift:
        cs = np.concatenate([[0.0], np.cumsum(np.roll(
            design["c"], int(round(shift * design["dt"] * fs))))])
    else:
        cs = design["cs"]
    t0 = design["t"] + design["z_offsets"][group]         # plane frame start
    # Lags are measured from the middle of the plane's scan (t0 + frame/2).
    # With the frame window that is exactly [t0, t0 + frame] at lag 0; in
    # volume mode (frame = 0) it reproduces the original design, whose lag-0
    # bin is the volume period CENTRED on the sample.
    c0 = t0 + design.get("frame", 0.0) / 2
    if design.get("lag_basis", "point") == "box":
        # bin j holds lags [j, j+1)*lag_dt: causal by construction
        dl = design["lag_dt"]
        lo_t = c0[:, None] - design["lags"][None, :] - dl
        width = dl
    else:
        lo_t = c0[:, None] - design["lags"][None, :] - w / 2
        width = w
    lo = np.clip(np.round(lo_t * fs).astype(np.int64), 0, cs.size - 2)
    hi = np.clip(np.round((lo_t + width) * fs).astype(np.int64), 0,
                 cs.size - 1)
    hi = np.maximum(hi, lo + 1)
    X = (cs[hi] - cs[lo]) / (hi - lo)
    if not shift:                   # the null's shifted copies are not reused
        cache[key] = X
    return X


def design_matrix(design, group, mask=None, shift=0):
    """(rows, n_lags) stimulus history on one z plane; column j = lag
    design['lags'][j].

    Column j is the stimulus averaged over the sample's integration window
    moved lag_j into the past, computed from that row's own acquisition time on
    its own plane -- which is what lets jitter resolve lags finer than the
    volume period. shift rolls the stimulus by that many imaging samples'
    worth of time -- the null.
    """
    X = _full_design(design, group, shift)
    return X if mask is None else X[mask]


# ─────────────────────────────────────────────────────────────────────────────
# fitting
# ─────────────────────────────────────────────────────────────────────────────

def penalty_matrix(n_lags, kind="ridge"):
    """Regulariser for the normal equations.

    'ridge'  identity: shrinks every lag equally. Right when the lag spacing
             equals the volume period, where neighbouring lags are nearly
             independent.
    'smooth' second differences (plus 1e-3 identity so it is invertible):
             penalises CURVATURE. Needed on a sub-frame lag grid: the stimulus
             holds each value for a 33 ms frame, so it has almost no power at
             30 Hz and its harmonics, the data cannot constrain the filter
             there, and a ridge leaves a 33 ms ripple in its place. The
             curvature penalty removes exactly those components and leaves the
             low frequencies the data do constrain.
    """
    if kind == "smooth" and n_lags >= 3:
        D = np.diff(np.eye(n_lags), n=2, axis=0)
        P = D.T @ D
        P = P / (np.trace(P) / n_lags)
        return P + 1e-3 * np.eye(n_lags)
    return np.eye(n_lags)


def _solve(X, Y, ridge, P=None):
    xm, ym = X.mean(0), Y.mean(0)
    Xc = X - xm
    A = Xc.T @ Xc
    if P is None:
        P = np.eye(A.shape[0])
    lam = float(ridge) * np.trace(A) / np.trace(P)
    h = np.linalg.solve(A + lam * P, Xc.T @ (Y - ym))
    return h, ym - xm @ h                     # filters (lags, rois), intercept


def fit_filters(design, mask, ridge, shift=0):
    """Filters (n_roi, n_lags) and intercepts (n_roi,) from the rows in mask.

    Units: dF/F per unit contrast, per lag bin. The regulariser is
    design['penalty'] (see penalty_matrix), weighted relative to trace(X'X), so the same number regularises the low- and
    high-contrast fits equally even though their stimulus power differs ~17x.
    """
    h = np.full((design["n_roi"], design["n_lags"]), np.nan)
    b = np.full(design["n_roi"], np.nan)
    for g in range(len(design["z_offsets"])):
        rois = np.flatnonzero(design["roi_group"] == g)
        if rois.size == 0 or not mask.any():
            continue
        X = design_matrix(design, g, mask, shift)
        hh, bb = _solve(X, design["Y"][mask][:, rois], ridge,
                        design.get("penalty"))
        h[rois], b[rois] = hh.T, bb
    return h, b


def _fold_sets(design, n_fold, seed=0):
    ub = np.unique(design["fold"])
    rng = np.random.default_rng(seed)
    return np.array_split(rng.permutation(ub), int(min(n_fold, len(ub))))


def cv_filters(design, fit_mask, ridge, n_fold=5, shift=0, seed=0):
    """Out-of-fold filters: one fit per fold, from fit_mask rows NOT in it.

    Returns (h_folds (n_fold, n_roi, n_lags), b_folds, row_fold) where row_fold
    assigns every row -- in fit_mask or not -- to the fold whose filter did not
    see it. Anything evaluated on a row with its own fold's filter is therefore
    held out, which is what lets the gain-after-switch curve include the
    steady-state bins without scoring the filter on its own training data.
    """
    sets = _fold_sets(design, n_fold, seed)
    row_fold = np.full(design["fold"].size, -1)
    hf, bf = [], []
    for k, fs in enumerate(sets):
        te = np.isin(design["fold"], fs)
        row_fold[te] = k
        h, b = fit_filters(design, fit_mask & ~te, ridge, shift)
        hf.append(h)
        bf.append(b)
    return np.array(hf), np.array(bf), row_fold


def _generator(design, h_folds, b_folds, row_fold, mask, shift=0):
    """Out-of-fold linear prediction (rows_in_mask, n_roi)."""
    idx = np.flatnonzero(mask)
    out = np.full((idx.size, design["n_roi"]), np.nan)
    for g in range(len(design["z_offsets"])):
        rois = np.flatnonzero(design["roi_group"] == g)
        if rois.size == 0:
            continue
        X = design_matrix(design, g, mask, shift)
        for k in range(len(h_folds)):
            sel = row_fold[idx] == k
            if sel.any():
                out[np.ix_(sel, rois)] = (X[sel] @ h_folds[k][rois].T
                                          + b_folds[k][rois])
    return out


def _colcorr(a, b):
    a = a - np.nanmean(a, 0)
    b = b - np.nanmean(b, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.nansum(a * b, 0)
                / np.sqrt(np.nansum(a * a, 0) * np.nansum(b * b, 0)))


def cv_r(design, mask, ridge, n_fold=5, shift=0, seed=0):
    """Held-out prediction r per ROI on the rows in mask."""
    hf, bf, rf = cv_filters(design, mask, ridge, n_fold, shift, seed)
    pred = _generator(design, hf, bf, rf, mask, shift)
    return _colcorr(pred, design["Y"][mask])


def choose_ridge(design, masks, grid=None, n_fold=5):
    """One ridge for the channel, maximising median cv_r over ROIs AND
    conditions. Shared across conditions on purpose: a ridge chosen per
    condition would smooth the two filters differently, and the difference
    between them is the result."""
    grid = np.logspace(-4, 1, 11) if grid is None else np.asarray(grid)
    scores = np.array([[np.nanmedian(cv_r(design, m, lam, n_fold))
                        for m in masks] for lam in grid])
    med = np.nanmedian(scores, axis=1)
    best = (int(np.nanargmax(med)) if np.any(np.isfinite(med))
            else len(grid) // 2)
    return float(grid[best]), grid, scores


def cycle_shifts(design, n_shifts, min_s=30.0, seed=0):
    """Circular shifts (in samples) that are whole multiples of the contrast
    cycle, so a shifted stimulus still has the right VARIANCE on every row and
    only its pairing with the response is destroyed. A shift that moved the
    envelope too would make the null a test of the envelope as well."""
    ons = design["onsets"]
    n_cond = len(design["cond_keys"])
    cyc = float(np.median(np.diff(ons))) * n_cond
    n_cyc = int((design["tvol"][-1] - design["tvol"][0]) // cyc)
    k_min = int(np.ceil(min_s / cyc))
    ks = np.arange(k_min, max(k_min + 1, n_cyc - k_min))
    rng = np.random.default_rng(seed)
    ks = rng.choice(ks, size=min(n_shifts, ks.size), replace=False)
    return np.round(ks * cyc / design["dt"]).astype(int)


# ─────────────────────────────────────────────────────────────────────────────
# adaptation readouts
# ─────────────────────────────────────────────────────────────────────────────

def _slope(g, y):
    """Least-squares slope of y on g per column, with intercept."""
    gc = g - np.nanmean(g, 0)
    yc = y - np.nanmean(y, 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(gc * yc, 0) / np.nansum(gc * gc, 0)


def normalize_filters(design, h, mask):
    """Scale each ROI's filter so var(X h) = var(s) over the rows in mask.

    The Baccus & Meister convention: the filter keeps its shape and loses its
    gain, which then belongs entirely to the nonlinearity. s is column 0 of the
    design matrix, i.e. the stimulus at the same resolution the filter sees.
    Returns (scaled filters, scale per ROI). The sign of the filter is kept, so
    an OFF-type filter stays negative and its nonlinearity is still increasing.
    """
    out = np.full_like(h, np.nan)
    scale = np.full(h.shape[0], np.nan)
    for g in range(len(design["z_offsets"])):
        rois = np.flatnonzero(design["roi_group"] == g)
        if rois.size == 0 or not mask.any():
            continue
        X = design_matrix(design, g, mask)
        sd_s = X[:, 0].std()
        with np.errstate(invalid="ignore", divide="ignore"):
            k = sd_s / (X @ h[rois].T).std(0)
        scale[rois] = k
        out[rois] = h[rois] * k[:, None]
    return out, scale


def normalize_folds(design, folds, mask):
    """normalize_filters applied to every fold's filter; intercepts set to 0.

    Each fold is scaled on the same rows, so the out-of-fold generator has the
    stimulus's variance whichever fold a row belongs to. The intercept goes
    because it is part of the nonlinearity, not the filter.
    """
    hf, bf, rf = folds
    hn = np.array([normalize_filters(design, h, mask)[0] for h in hf])
    return hn, np.zeros_like(bf), rf


def _generator_split(design, folds_cur, folds_prev, mask):
    """Out-of-fold generator split at the switch: (g_cur, g_prev).

    g_cur applies the CURRENT condition's filter to the taps of each row's lag
    window that fall after the switch; g_prev applies the PREVIOUS condition's
    filter to the taps before it. A row a whole filter length past the switch
    has g_prev = 0. Both use each row's held-out fold.
    """
    idx = np.flatnonzero(mask)
    n_roi = design["n_roi"]
    g_cur = np.full((idx.size, n_roi), np.nan)
    g_prev = np.zeros((idx.size, n_roi))
    hc, _, rf = folds_cur
    lag_t = design["lags"]
    for g in range(len(design["z_offsets"])):
        rois = np.flatnonzero(design["roi_group"] == g)
        if rois.size == 0:
            continue
        # a tap is post-switch when its window STARTS after the switch, on
        # this plane's own clock
        since = design["since_switch"][mask] + design["z_offsets"][g]
        post = lag_t[None, :] <= since[:, None]
        X = design_matrix(design, g, mask)
        Xc, Xp = X * post, X * ~post
        for k in range(len(hc)):
            sel = rf[idx] == k
            if not sel.any():
                continue
            g_cur[np.ix_(sel, rois)] = Xc[sel] @ hc[k][rois].T
            if folds_prev is not None:
                g_prev[np.ix_(sel, rois)] = Xp[sel] @ folds_prev[0][k][rois].T
    return g_cur, g_prev


def _two_slopes(g1, g2, y):
    """Per-column OLS of y on [g1, g2, 1]; returns the g1 coefficient. Falls
    back to the single slope when g2 carries nothing (late bins)."""
    out = np.full(y.shape[1], np.nan)
    for r in range(y.shape[1]):
        a, b, yy = g1[:, r], g2[:, r], y[:, r]
        ok = np.isfinite(a) & np.isfinite(b) & np.isfinite(yy)
        if ok.sum() < 10:
            continue
        cols = [a[ok], np.ones(ok.sum())]
        if np.ptp(b[ok]) > 0:
            cols.insert(1, b[ok])
        coef = np.linalg.lstsq(np.stack(cols, 1), yy[ok], rcond=None)[0]
        out[r] = coef[0]
    return out


def gain_after_switch(design, ss_mask, ridge, bin_edges, n_fold=5,
                      normalize=True):
    """Gain on the NEW stimulus vs time since the switch, relative to each
    condition's own steady state, from OUT-OF-FOLD steady-state filters.

    Early after a switch a row's lag window straddles it: part of what the cell
    is responding to was shown under the previous contrast. Scoring the whole
    window against the new condition's filter would fold the response to that
    old stimulus into the "gain" -- after a high-to-low switch it drags the
    first bin toward zero whatever the cell is doing. So the generator is split
    at the switch (see _generator_split) and the response regressed on both
    parts; the coefficient on the post-switch part is the gain reported, and
    the pre-switch part is a nuisance regressor. Late bins have no pre-switch
    taps and reduce to a plain slope.

    1 means adapted to the current contrast. The first bins' distance from 1 is
    the gain the cell carried over from the previous block; how fast it closes
    is the adaptation time course. Both folds use the same fold assignment, so
    the two filters in a row never saw that row.

    With normalize=True (the default) the filters follow the Baccus & Meister
    convention (see normalize_filters), so the slope IS the gain of the
    nonlinearity in dF/F per unit contrast: both conditions share one axis,
    each converges to its own steady-state gain, and the first bins show where
    it started from. Dividing by the steady-state value gives the relative
    curve.

    Returns gain (n_cond, n_bins, n_roi), n_rows (n_cond, n_bins), and the
    fold filters per condition, which the other readouts reuse.
    """
    n_cond = len(design["cond_keys"])
    nb = len(bin_edges) - 1
    gain = np.full((n_cond, nb, design["n_roi"]), np.nan)
    nrow = np.zeros((n_cond, nb), dtype=int)
    folds = {}
    for q in range(n_cond):
        m = ss_mask & (design["cond"] == q)
        folds[q] = cv_filters(design, m, ridge, n_fold)
        if normalize:
            folds[q] = normalize_folds(design, folds[q], m)
    for q in range(n_cond):
        for p in [c for c in range(-1, n_cond) if c != q]:
            allq = (design["cond"] == q) & (design["prev_cond"] == p)
            if not allq.any():
                continue
            g1, g2 = _generator_split(design, folds[q],
                                      folds.get(p) if p >= 0 else None, allq)
            y = design["Y"][allq]
            ts = design["since_switch"][allq]
            for j in range(nb):
                m = (ts >= bin_edges[j]) & (ts < bin_edges[j + 1])
                if m.sum() <= 10:
                    continue
                gj = _two_slopes(g1[m], g2[m], y[m])
                # several previous conditions (>2 contrasts): row-weighted mean
                w0 = nrow[q, j]
                gain[q, j] = (gj if w0 == 0 else
                              (np.nan_to_num(gain[q, j]) * w0 + gj * m.sum())
                              / (w0 + m.sum()))
                nrow[q, j] += int(m.sum())
    return gain, nrow, folds


def window_filters(design, ridge, windows):
    """Filters fitted within time-since-switch windows, per condition.

    windows is a list of (t0, t1). A window that starts before one filter
    length has lags reaching into the previous block; its filter is still the
    cell's CURRENT kernel applied to the real stimulus history, which is the
    point -- but its estimate mixes both blocks' stimulus statistics, so read
    its shape with that in mind. Returns (n_cond, n_windows, n_roi, n_lags).
    """
    n_cond = len(design["cond_keys"])
    out = np.full((n_cond, len(windows), design["n_roi"], design["n_lags"]),
                  np.nan)
    for q in range(n_cond):
        for w, (a, b) in enumerate(windows):
            m = ((design["cond"] == q) & (design["since_switch"] >= a)
                 & (design["since_switch"] < b))
            if m.sum() > 3 * design["n_lags"]:
                out[q, w] = fit_filters(design, m, ridge)[0]
    return out


def ln_model(design, ss_mask, folds, n_bins=12):
    """Steady-state LN model per condition, Baccus & Meister convention.

    folds must be NORMALIZED (gain_after_switch returns them so). For each
    condition the generator is that condition's own out-of-fold filter applied
    to its own stimulus, so its SD equals the stimulus SD; the nonlinearity is
    the mean response in quantile bins of it. Reported per condition and ROI:

      x, y          (n_cond, n_bins, n_roi) the nonlinearity
      gain          OLS slope of response on generator over the whole
                    distribution: dF/F per unit contrast
      shared_gain   slope within the lowest condition's 2.5-97.5 percentile
                    generator range, where every condition has data. A static
                    nonlinearity would give equal shared_gain in every
                    condition; adaptation does not.
      g_sd          generator SD per condition (= stimulus SD, the check)
    """
    n_cond = len(design["cond_keys"])
    gens = []
    for q in range(n_cond):
        m = ss_mask & (design["cond"] == q)
        hf, bf, rf = folds[q]
        gens.append((_generator(design, hf, bf, rf, m), design["Y"][m]))
    lo, hi = np.nanpercentile(gens[0][0], [2.5, 97.5], axis=0)
    n_roi = design["n_roi"]
    x = np.full((n_cond, n_bins, n_roi), np.nan)
    y = np.full_like(x, np.nan)
    gain = np.full((n_cond, n_roi), np.nan)
    shared = np.full_like(gain, np.nan)
    g_sd = np.full_like(gain, np.nan)
    for q, (g, yy) in enumerate(gens):
        gain[q] = _slope(g, yy)
        g_sd[q] = np.nanstd(g, 0)
        qs = np.nanpercentile(g, np.linspace(0, 100, n_bins + 1), axis=0)
        for r in range(n_roi):
            if not np.all(np.isfinite(qs[:, r])):
                continue
            bi = np.clip(np.searchsorted(qs[:, r], g[:, r], "right") - 1,
                         0, n_bins - 1)
            for j in range(n_bins):
                sel = bi == j
                if sel.any():
                    x[q, j, r] = g[sel, r].mean()
                    y[q, j, r] = yy[sel, r].mean()
            sel = (g[:, r] >= lo[r]) & (g[:, r] <= hi[r])
            if sel.sum() > 20:
                shared[q, r] = _slope(g[sel, r:r + 1], yy[sel, r:r + 1])[0]
    return dict(x=x, y=y, gain=gain, shared_gain=shared, g_sd=g_sd,
                shared_range=np.array([lo, hi]))


def nonlinearity(design, ss_mask, folds, ref_cond, n_bins=12):
    """Static nonlinearity per condition on ONE generator axis.

    The axis is the reference condition's out-of-fold filter, scaled to unit
    norm, applied to both conditions' stimuli -- so a given x is the same
    stimulus drive in both, and a difference in slope at the same x cannot be a
    static nonlinearity. Reported:

      x, y            binned (quantiles of each condition) mean response
      slope           linear gain of each condition on the shared axis
      shared_slope    gain fitted only within the lowest condition's
                      2.5-97.5 percentile range of x, i.e. where the two
                      overlap. Adaptation: these differ. Static
                      nonlinearity: they match.
    """
    hf, _, rf = folds[ref_cond]
    nrm = np.sqrt(np.nansum(hf ** 2, axis=2, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        hu = hf / nrm
    zero_b = np.zeros(hf.shape[:2])
    n_cond = len(design["cond_keys"])
    gens, ranges = [], []
    for q in range(n_cond):
        m = ss_mask & (design["cond"] == q)
        g = _generator(design, hu, zero_b, rf, m)
        gens.append((g, design["Y"][m]))
        ranges.append(np.nanpercentile(g, [2.5, 97.5], axis=0))
    lo, hi = ranges[0]
    xs, ys, slope, shared = [], [], [], []
    for g, y in gens:
        qs = np.nanpercentile(g, np.linspace(0, 100, n_bins + 1), axis=0)
        xb = np.full((n_bins, g.shape[1]), np.nan)
        yb = np.full_like(xb, np.nan)
        sh = np.full(g.shape[1], np.nan)
        for r in range(g.shape[1]):
            if not np.all(np.isfinite(qs[:, r])):
                continue
            bi = np.clip(np.searchsorted(qs[:, r], g[:, r], "right") - 1,
                         0, n_bins - 1)
            for j in range(n_bins):
                s = bi == j
                if s.any():
                    xb[j, r] = g[s, r].mean()
                    yb[j, r] = y[s, r].mean()
            s = (g[:, r] >= lo[r]) & (g[:, r] <= hi[r])
            if s.sum() > 20:
                sh[r] = _slope(g[s, r:r + 1], y[s, r:r + 1])[0]
        xs.append(xb)
        ys.append(yb)
        slope.append(_slope(g, y))
        shared.append(sh)
    return dict(x=np.array(xs), y=np.array(ys), slope=np.array(slope),
                shared_slope=np.array(shared), shared_range=np.array([lo, hi]))


def switch_triggered_average(design, to_cond, pre_s, post_s, bin_s=None):
    """Responses aligned to every switch INTO to_cond, per ROI.

    The contrast-envelope response: what the old epoch analysis measured, and
    for a signal with no linear filter at the probe's timescale the main thing
    there is to see. Returns (lag, trials (n_switch, n_roi, n_t)).

    bin_s=None interpolates each switch onto the volume grid (smooths at the
    volume period). With bin_s, every sample is placed at its EXACT time since
    the switch -- its plane's own frame time -- and averaged in bins of bin_s:
    because switches fall at random phases of the volume clock, each bin is
    filled by different switches, and the average is resolved at bin_s rather
    than at the volume period. A switch then contributes to a bin only when one
    of its samples landed there, so trials carries NaN elsewhere; average with
    nanmean.
    """
    dt = design["dt"]
    tvol = design["tvol"]
    ec = design["epoch_cond"]
    ons = design["onsets"]
    sw = [ons[i] for i in range(1, len(ec))
          if ec[i] == to_cond and ec[i - 1] != to_cond]
    if bin_s:
        edges = np.arange(-pre_s, post_s + bin_s / 2, bin_s)
        lag = 0.5 * (edges[:-1] + edges[1:])
        Yf = design["Y_full"]
        off = design["z_offsets"][design["roi_group"]]
        trials = np.full((len(sw), Yf.shape[0], lag.size), np.nan)
        keep = np.ones(len(sw), bool)
        for i, t0 in enumerate(sw):
            if t0 - pre_s < tvol[0] or t0 + post_s > tvol[-1]:
                keep[i] = False
                continue
            for r in range(Yf.shape[0]):
                rel = tvol + off[r] - t0
                m = (rel >= edges[0]) & (rel < edges[-1])
                b = np.searchsorted(edges, rel[m], "right") - 1
                cnt = np.bincount(b, minlength=lag.size)
                tot = np.bincount(b, weights=Yf[r, m], minlength=lag.size)
                with np.errstate(invalid="ignore", divide="ignore"):
                    trials[i, r] = np.where(cnt > 0, tot / cnt, np.nan)
        return lag, trials[keep]
    lag = np.arange(-pre_s, post_s + dt / 2, dt)
    trials = []
    for t0 in sw:
        tq = t0 + lag
        if tq[0] < tvol[0] or tq[-1] > tvol[-1]:
            continue
        trials.append(np.array([np.interp(tq, tvol, y)
                                for y in design["Y_full"]]))
    return lag, np.array(trials)


def predict_series(design, filters, intercepts, gain=None, bin_edges=None):
    """Prediction of a whole series from filters fitted ELSEWHERE.

    filters[q] (n_roi, n_lags) and intercepts[q] are applied to the rows in
    condition q. With gain (n_cond, n_bins, n_roi) and bin_edges, each row's
    filter is also scaled by the gain of its time-since-switch bin -- the
    adaptation dynamics -- so comparing the two predictions asks whether the
    dynamics explain anything the steady-state filters do not.
    """
    pred = np.full(design["Y"].shape, np.nan)
    for q in range(len(filters)):
        m = design["cond"] == q
        if not m.any():
            continue
        idx = np.flatnonzero(m)
        for g in range(len(design["z_offsets"])):
            rois = np.flatnonzero(design["roi_group"] == g)
            if rois.size == 0:
                continue
            X = design_matrix(design, g, m)
            p = X @ filters[q][rois].T
            if gain is not None:
                ts = design["since_switch"][m]
                bi = np.clip(np.searchsorted(bin_edges, ts, "right") - 1,
                             0, gain.shape[1] - 1)
                gq = np.nan_to_num(gain[q][:, rois], nan=1.0)
                p = p * gq[bi]
            pred[np.ix_(idx, rois)] = p + intercepts[q][rois]
    return pred


def cycle_average(design, values, period_epochs):
    """Average a per-row quantity over repeats of a frozen cycle.

    For a repeated series every cycle of period_epochs epochs shows the same
    stimulus, so rows at the same time-in-cycle can be averaged. Returns
    (t_in_cycle grid, mean (n_t, n_roi), split-half r per ROI with
    Spearman-Brown correction -- the noise ceiling for predicting that mean).
    """
    ons = design["onsets"]
    cyc_id = design["epoch"] // period_epochs
    t_in = design["t"] - ons[cyc_id * period_epochs]
    dt = design["dt"]
    edges = np.arange(0, t_in.max() + dt, dt)
    bi = np.clip(np.searchsorted(edges, t_in, "right") - 1, 0, edges.size - 2)
    ids = np.unique(cyc_id)
    half = [ids[::2], ids[1::2]]
    means = []
    for sel in [ids] + half:
        m = np.isin(cyc_id, sel)
        acc = np.full((edges.size - 1, values.shape[1]), np.nan)
        for j in range(edges.size - 1):
            s = m & (bi == j)
            if s.sum():
                acc[j] = np.nanmean(values[s], 0)
        means.append(acc)
    r = _colcorr(means[1], means[2])
    ceiling = 2 * r / (1 + r)
    return edges[:-1], means[0], ceiling


def fit_exp_recovery(t, g):
    """g(t) = 1 + A exp(-t / tau), least squares over a tau grid.

    Coarse on purpose: with a handful of bins the time constant is resolved to
    roughly the bin width, and a nonlinear fitter would report more precision
    than the data hold. Returns (A, tau) or (nan, nan).
    """
    ok = np.isfinite(g)
    if ok.sum() < 3:
        return np.nan, np.nan
    t, y = np.asarray(t)[ok], np.asarray(g)[ok] - 1.0
    best = (np.inf, np.nan, np.nan)
    for tau in np.logspace(-1.5, 1, 200):
        e = np.exp(-t / tau)
        A = float(e @ y / (e @ e))
        sse = float(((y - A * e) ** 2).sum())
        if sse < best[0]:
            best = (sse, A, tau)
    return best[1], best[2]
