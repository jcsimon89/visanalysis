"""
octave_strf.py

The two computations that differ from the bar-noise pipeline: reverse
correlation in the cell basis, and the per-octave temporal deconvolution that
the exponential holds make necessary.

Kept out of analyze_data_strf_octave.py so they can be unit-tested headless
against a synthetic kernel, which is the only way to know the deconvolution is
doing what it claims.
"""

import numpy as np


# ─────────────────────────────────────────────────────────────────────────────
# reverse correlation
# ─────────────────────────────────────────────────────────────────────────────

def compute_strf_cells(trials, n_roi, n_cells, filter_length_s, pre_time=0.0):
    """Estimate the STRF in the cell basis by reverse correlation.

    Structurally identical to compute_strf_1d in the noizone pipeline -- the
    regressor is just a different matrix -- but worth being explicit about why
    the plain STA is the right estimator here rather than a convenient
    approximation.

    Cells are mutually independent within and across octaves and every cell has
    the same variance 2p, so the spatial covariance is exactly 2p * I. The STA
    is therefore the maximum-likelihood estimate up to that known constant: no
    whitening, no C_ss to invert, no spatial ridge, and none of the noise
    amplification a texel-basis analysis would incur. That is the whole reason
    for estimating here rather than on the display grid.

    What the STA does NOT remove is the temporal correlation -- see
    deconvolve_octaves.

    Parameters
    ----------
    trials : list of dicts, one per (trial, roi) observation, each with
        'roi_ind', 'response', 'time_vector', 'stim_cells' (n_frames, n_cells),
        'stim_dt'. Same contract as compute_strf_1d.
    n_roi, n_cells  : output dimensions
    filter_length_s : filter length in seconds
    pre_time        : baseline before stim onset (s)

    Returns
    -------
    strf_raw : (n_roi, n_cells, n_lags)  raw STA
    t_lag    : (n_lags,)  lag axis, 0 = simultaneous, + = past
    n_valid  : (n_roi,)   samples contributing, per ROI
    """
    lag_dt = min(t['stim_dt'] for t in trials)
    n_lags = max(1, int(np.round(filter_length_s / lag_dt)))
    t_lag = np.arange(n_lags) * lag_dt

    strf_raw = np.zeros((n_roi, n_cells, n_lags), dtype=np.float64)
    n_valid_total = np.zeros(n_roi, dtype=np.int64)

    # ROIs that share a time vector share the whole gather, so batch them: the
    # per-lag cost is one (n_roi_grp, n_valid) @ (n_valid, n_cells) matmul
    # instead of n_roi separate ones, and the expensive fancy-index of the
    # stimulus happens once per lag rather than once per (lag, ROI). With 1212
    # cells the gather dominates, so this is the difference between seconds and
    # minutes per series. Grouping (rather than assuming a shared vector) keeps
    # it correct when per-ROI z-slice offsets make the vectors genuinely differ.
    groups = {}
    for tr in trials:
        tvec = tr['time_vector']
        if len(tvec) == 0:
            continue
        key = (id(tr['stim_cells']), tr['stim_dt'], len(tvec),
               float(tvec[0]), float(tvec[-1]))
        groups.setdefault(key, []).append(tr)

    for key, grp in groups.items():
        stim, stim_dt = grp[0]['stim_cells'], grp[0]['stim_dt']
        tvec = grp[0]['time_vector']
        n_frames = stim.shape[0]
        n_tp = min(len(tvec), min(len(t['response']) for t in grp))
        tvec_trial = tvec[:n_tp]

        rois = np.array([t['roi_ind'] for t in grp])
        R = np.stack([np.asarray(t['response'][:n_tp], dtype=np.float64)
                      for t in grp])                       # (n_grp, n_tp)

        time_valid = tvec_trial >= (pre_time + filter_length_s)
        stim_indices = np.round((tvec_trial - pre_time) / stim_dt).astype(int)

        for lag_idx, t in enumerate(t_lag):
            past_idx = stim_indices - int(np.round(t / stim_dt))
            valid = time_valid & (past_idx >= 0) & (past_idx < n_frames)
            if not valid.any():
                continue
            S = stim[past_idx[valid], :]                   # (n_valid, n_cells)
            strf_raw[rois, :, lag_idx] += R[:, valid] @ S

        # time_valid, not the loop's `valid`. Every lag contributes about
        # time_valid.sum() samples, so that is the divisor that makes strf_raw a
        # mean. `valid` additionally trims where past_idx runs off the start of
        # the stimulus, which bites hardest at the LAST lag -- reading it after
        # the loop would both leak a loop variable and normalise by the most
        # trimmed lag, scaling every filter slightly high.
        n_valid_total[rois] += int(time_valid.sum())

    for roi_ind in range(n_roi):
        if n_valid_total[roi_ind] > 0:
            strf_raw[roi_ind] /= n_valid_total[roi_ind]

    return strf_raw, t_lag, n_valid_total


def compute_strf_cells_folds(trials, n_roi, n_cells, filter_length_s, pre_time=0.0):
    """compute_strf_cells, split into two disjoint halves of the trials.

    Returns (strf_full, (strf_a, strf_b), t_lag, n_valid_full).

    The split is by TRIAL, not by sample, so the two halves have genuinely
    independent noise -- which is what makes them usable for choosing the
    deconvolution penalty. Splitting by sample would not: neighbouring samples
    share the same held stimulus values.

    Costs no more than one full pass: each trial is visited exactly once, and
    the combined STA is the n_valid-weighted mean of the halves, which is exact
    because each half is already normalised by its own count.

    Each trial dict may carry a 'fold' key (0 or 1). Without one, trials are
    assigned alternately in the order given, which for a list built by looping
    over trials then ROIs interleaves whole trials as intended.
    """
    fa, fb = [], []
    seen = {}
    for tr in trials:
        if 'fold' in tr:
            f = int(tr['fold'])
        else:
            key = tr.get('trial_ind', id(tr['stim_cells']))
            if key not in seen:
                seen[key] = len(seen)
            f = seen[key] % 2
        (fa if f == 0 else fb).append(tr)

    if not fa or not fb:
        strf, t_lag, nv = compute_strf_cells(trials, n_roi, n_cells,
                                             filter_length_s, pre_time)
        return strf, None, t_lag, nv

    sa, t_lag, na = compute_strf_cells(fa, n_roi, n_cells, filter_length_s, pre_time)
    sb, _, nb = compute_strf_cells(fb, n_roi, n_cells, filter_length_s, pre_time)
    tot = (na + nb).astype(float)
    w = np.divide(1.0, tot, out=np.zeros_like(tot), where=tot > 0)
    full = (sa * na[:, None, None] + sb * nb[:, None, None]) * w[:, None, None]
    return full, (sa, sb), t_lag, (na + nb)


def peak_lag(strf_roi):
    """Index of the lag carrying the most signal, for one ROI's (n_cells, n_lags).

    Uses the summed SQUARE across cells, not the summed absolute value. With
    ~1200 cells of which only a few tens carry signal, |.| is dominated by the
    noise floor and picks an arbitrary lag -- measured, it returned 2.0 s for a
    kernel peaking at 0.45 s. Squaring weights the signal cells above the noise
    and recovers the right lag to within a frame or two.
    """
    return int(np.argmax((np.asarray(strf_roi) ** 2).sum(axis=0)))


def independent_samples(n_valid, tau, sample_dt):
    """How many independent stimulus configurations an octave actually saw.

    Consecutive samples are not independent -- a cell holds for tau -- so the
    count that matters is the recorded DURATION divided by tau, and duration is
    n_valid * sample_dt.

    `sample_dt` is the period of the RESPONSE samples, not of the display
    frames. n_valid counts imaging samples with a full filter window behind
    them, so passing the display period here understates the answer by the ratio
    of the two rates -- about 12x at 10 Hz imaging against a 120 Hz display,
    which is the difference between reporting a recording as hopeless and
    reporting it as ample.

    Compare it against the cell count of that octave. At 300 s and tau = 0.33 s
    it is 909 independent configurations against 1146 fine cells, which is
    underdetermined and no amount of clean recording fixes it; at 35 minutes it
    is 6400, which is not. Report it rather than let it surprise someone.
    """
    # n_valid is per-ROI, so this returns per-ROI. Scalars still work and come
    # back as a 0-d array's item, which is what a caller passing one expects.
    out = np.asarray(n_valid, dtype=float) * float(sample_dt) / float(tau)
    return out if out.ndim else float(out)


# ─────────────────────────────────────────────────────────────────────────────
# temporal deconvolution -- the one place regularisation is needed
# ─────────────────────────────────────────────────────────────────────────────

def _smear_matrix(tau, lag_dt, n_lags):
    """A[i, j] = exp(-|i - j| * lag_dt / tau).

    The STA at lag L is (h * R)(L) where R is the stimulus autocorrelation,
    which for exponential holds is the two-sided exp(-|lag|/tau). R is
    symmetric, so it broadens the kernel without shifting its latency -- but at
    these hold times it broadens it a lot, which is why this cannot be ignored
    the way the noizone triangle could.
    """
    i = np.arange(n_lags)
    return np.exp(-np.abs(i[:, None] - i[None, :]) * float(lag_dt) / float(tau))


def _split_half_reliability(blk_a, blk_b, cells):
    """Agreement between two folds, each read at ITS OWN peak lag.

    Deliberately not the whole (cells x lags) block. A failed deconvolution does
    not produce random noise, it produces a REPRODUCIBLE artefact -- both folds
    grow the same spike at lag zero, because it comes from the conditioning of
    the same smearing matrix, not from independent noise. Measured on a 6 deg
    octave whose recovery had collapsed from r = 0.83 to r = 0.09 against known
    truth, whole-block reliability read 0.861 before and 0.860 after: blind.

    Letting each fold choose its own peak lag makes the measure test the thing
    that actually broke. When the peak is an artefact its position is unstable
    between folds, so the profiles are compared at different lags and disagree.
    Same case: 0.483 before, -0.033 after. It is also what the analysis uses
    downstream -- centroids are taken from the spatial profile at the peak lag.
    """
    la = int(np.argmax((blk_a ** 2).sum(axis=0)))
    lb = int(np.argmax((blk_b ** 2).sum(axis=0)))
    a, b = blk_a[cells, la], blk_b[cells, lb]
    if a.size < 2 or not np.any(a) or not np.any(b):
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def deconvolve_octaves(strf_raw, octaves, lag_dt, lam=None, lam_scale=1.0,
                       cv_pair=None, guard_frac=0.8):
    """Undo the per-octave temporal smearing, by forward-modelling not division.

    Dividing by the Lorentzian amplifies as f^2 -- 5x at 1 Hz and 26x at 2.5 Hz
    for a 0.33 s hold, and far worse for the 1.0 s one. So instead of inverting,
    solve

        min_h  || A h - m ||^2  +  lam * || D2 h ||^2

    with A the smearing matrix above and D2 a second difference. The frequencies
    the smearing destroyed simply get little weight, rather than dominating the
    answer after amplification.

    Each octave is deconvolved with its own tau, which is the point of keeping
    the cell basis split: there is no single smearing kernel to remove.

    Parameters
    ----------
    strf_raw : (n_roi, n_cells, n_lags)
    octaves  : [(name, slice, patch_deg, tau), ...] from octave_stim.octave_slices
    lag_dt   : lag axis step (s)
    lam      : explicit penalty, applied to every octave. Overrides cv_pair.
    lam_scale: multiplies an explicit lam. No effect on the cross-validated one.
    cv_pair  : (strf_fold_a, strf_fold_b) -- two STAs from disjoint halves of
               the trials, used to choose lam per octave. Strongly recommended;
               with neither this nor an explicit lam the function returns the
               STA untouched rather than guessing at a penalty.
    guard_frac: reject the deconvolution for an octave if it costs more than
               this fraction of the STA's split-half reliability. Deconvolution
               ALWAYS amplifies noise, and there are regimes where no lam is an
               improvement -- measured on a 6 deg octave at realistic SNR, the
               best lam in the whole grid recovered the spatial profile at
               r = 0.57 where the un-deconvolved STA reached 0.83, and the peak
               lag collapsed to zero at every lam. A small reliability cost buys
               a real latency correction and is worth taking; losing most of it
               means the sharpening is noise. Needs cv_pair. None disables.

    Returns
    -------
    strf_dec : (n_roi, n_cells, n_lags), same shape
    info     : dict of per-octave lam actually used
    """
    n_roi, n_cells, n_lags = strf_raw.shape
    out = np.empty_like(strf_raw)
    info = {}
    tail = slice(int(0.75 * n_lags), n_lags)   # assumed signal-free

    for name, sl, patch_deg, tau in octaves:
        A = _smear_matrix(tau, lag_dt, n_lags)
        D2 = np.diff(np.eye(n_lags), 2, axis=0)
        AtA, DtD = A.T @ A, D2.T @ D2

        block = strf_raw[:, sl, :]                       # (n_roi, n_cell_j, n_lags)
        flat = block.reshape(-1, n_lags).T               # (n_lags, n_obs)
        n_obs = flat.shape[1]

        if lam is not None:
            lam_j, how = float(lam) * lam_scale, 'explicit'
        elif cv_pair is not None:
            # Cross-validation across trial halves. Nothing cheaper works: the
            # right penalty spans four decades depending mostly on how smooth
            # the true kernel is, which is the thing being measured. A fixed
            # default, the discrepancy principle and GCV were all tried and all
            # fail -- measured, a lam wrong by 100x turns a 0.83 -> 0.99
            # recovery into 0.83 -> 0.05.
            #
            # The two halves have independent noise, so choosing lam to minimise
            # ||A h(half A) - half B||^2 estimates generalisation directly.
            # Measured against an oracle on 18 kernel/tau/noise combinations:
            # median gap 0.020, worst 0.110, and it beats the un-deconvolved STA
            # in 17 of 18 -- the exception being the smoothest kernel, where the
            # smearing barely matters to begin with.
            # Cross-validate on the PROJECTION, not on all cells pooled. Most
            # cells in an octave carry no signal, and for those the CV-optimal
            # penalty is infinite; pooled, they outvote the few that matter and
            # drive lam to the ceiling. Measured, that turned a 0.86 STA into
            # 0.78 where the oracle reached 0.99. Reducing each ROI to one
            # well-driven time series first restores the situation the method
            # was validated in.
            # Which cells enter the CV must be decided from the COMBINED STA,
            # never from one fold. Selecting on fold A inflates fold A's own
            # noise in whatever is selected, fold B has no matching inflation,
            # and the CV concludes that nothing generalises -- it then pins lam
            # to the ceiling and the "deconvolution" is pure smoothing. Measured,
            # that gave 0.78 against an oracle of 0.99. Selecting on the sum is
            # symmetric between the folds, so neither is favoured.
            A_blk, B_blk = cv_pair[0][:, sl, :], cv_pair[1][:, sl, :]
            comb = block                                   # already the full STA
            noise = np.sqrt(np.mean(comb[:, :, tail] ** 2)) or 1.0
            cols = []
            for r in range(A_blk.shape[0]):
                pk = peak_lag(comb[r])
                strong = np.abs(comb[r][:, pk]) > 4.0 * noise
                if strong.sum() < 3:
                    continue
                cols.append((A_blk[r][strong].T, B_blk[r][strong].T))
            if not cols:
                lam_j, how = np.inf, 'cv unusable (no signal)'
                out[:, sl, :] = block
                info[name] = {'tau': tau, 'lam': None, 'how': how}
                continue
            a = np.hstack([c[0] for c in cols])
            b = np.hstack([c[1] for c in cols])
            base = np.trace(AtA) / max(np.trace(DtD), 1e-12)
            best = None
            for m in np.logspace(-6, 8, 57):
                cand = base * m
                h_a = np.linalg.solve(AtA + cand * DtD, A.T @ a)
                err = float(np.sum((A @ h_a - b) ** 2))
                if best is None or err < best[0]:
                    best = (err, cand)
            lam_j, how = best[1], 'cross-validated on projection'

            # Having chosen lam, check that deconvolving is better than not.
            # Nothing above can tell the difference: the CV objective scores how
            # well A h predicts the other fold, which rewards denoising even when
            # the recovered h is worse than the STA it came from.
            if guard_frac is not None:
                cells = np.zeros(block.shape[1], dtype=bool)
                for r in range(block.shape[0]):
                    p = peak_lag(block[r])
                    cells |= np.abs(block[r][:, p]) > 4.0 * noise
                if cells.sum() >= 3:
                    raw_rel = np.mean([
                        _split_half_reliability(A_blk[r], B_blk[r], cells)
                        for r in range(A_blk.shape[0])])
                    da = np.linalg.solve(AtA + lam_j * DtD,
                                         A.T @ A_blk.reshape(-1, n_lags).T)
                    db = np.linalg.solve(AtA + lam_j * DtD,
                                         A.T @ B_blk.reshape(-1, n_lags).T)
                    da = da.T.reshape(A_blk.shape)
                    db = db.T.reshape(B_blk.shape)
                    dec_rel = np.mean([
                        _split_half_reliability(da[r], db[r], cells)
                        for r in range(da.shape[0])])
                    if dec_rel < guard_frac * raw_rel:
                        out[:, sl, :] = block
                        info[name] = {
                            'tau': tau, 'lam': None,
                            'how': 'rejected: deconvolution cost too much '
                                   'split-half reliability',
                            'lam_considered': float(lam_j),
                            'split_half_raw': float(raw_rel),
                            'split_half_deconvolved': float(dec_rel),
                            'noise_sd_from_tail': float(noise)}
                        continue
                    _rel = (float(raw_rel), float(dec_rel))
                else:
                    _rel = None
            else:
                _rel = None
        else:
            # No folds supplied and no explicit lam: refuse to guess. Returning
            # a silently wrong deconvolution is worse than returning the STA.
            lam_j, how = np.inf, 'none (returned un-deconvolved)'
            out[:, sl, :] = block
            info[name] = {'tau': tau, 'lam': None, 'how': how,
                          'noise_sd_from_tail': float(np.sqrt(
                              np.mean(block[:, :, tail] ** 2)))}
            continue

        sol = np.linalg.solve(AtA + lam_j * DtD, A.T @ flat)
        out[:, sl, :] = sol.T.reshape(block.shape)
        info[name] = {'tau': tau, 'lam': float(lam_j), 'how': how,
                      'noise_sd_from_tail': float(np.sqrt(
                          np.mean(block[:, :, tail] ** 2)))}
        if lam is None and cv_pair is not None and _rel is not None:
            info[name]['split_half_raw'] = _rel[0]
            info[name]['split_half_deconvolved'] = _rel[1]

    return out, info


# ─────────────────────────────────────────────────────────────────────────────
# significance, with per-octave degrees of freedom
# ─────────────────────────────────────────────────────────────────────────────

def zscore_by_octave(strf, octaves, lag_dt):
    """z-score each octave against its own noise level.

    A single z threshold across both octaves would treat them as carrying
    equal information when they do not: at 120 fps a 1.0 s hold gives one
    independent sample per 120 frames and a 0.33 s hold one per 40, so the
    coarse layer has three times fewer. Normalising within an octave keeps the
    threshold meaning the same thing in both.

    Returns (strf_z, per_octave_info).
    """
    out = np.empty_like(strf)
    info = {}
    for name, sl, patch_deg, tau in octaves:
        block = strf[:, sl, :]
        sd = block.reshape(block.shape[0], -1).std(axis=1)   # per ROI
        sd = np.where(sd > 0, sd, 1.0)
        out[:, sl, :] = block / sd[:, None, None]
        info[name] = {'frames_per_independent_sample': float(tau) / float(lag_dt),
                      'median_sd': float(np.median(sd))}
    return out, info


# ─────────────────────────────────────────────────────────────────────────────
# position, on the sphere
# ─────────────────────────────────────────────────────────────────────────────

def spherical_centroid(weights, directions, frac=0.5):
    """Weighted mean of unit vectors, renormalised -- not a mean of az/el.

    Averaging azimuth and elevation as if they were Cartesian coordinates goes
    wrong away from the origin, which on this display is most of the lit band.
    Weights below `frac` of the peak are dropped so the noise floor in the
    tails does not drag the centroid toward the array centre, matching the
    convention in average_strf_noizone.py.

    weights    : (n_cells,) non-negative importance, e.g. |filter| at peak lag
    directions : (n_cells, 3) unit vectors
    """
    w = np.asarray(weights, dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    if w.max() <= 0:
        return np.array([np.nan, np.nan, np.nan]), 0
    keep = w >= frac * w.max()
    v = (w[keep, None] * directions[keep]).sum(axis=0)
    n = np.linalg.norm(v)
    if n == 0:
        return np.array([np.nan, np.nan, np.nan]), int(keep.sum())
    return v / n, int(keep.sum())


def direction_to_azel_deg(d):
    az = np.degrees(np.arctan2(d[1], d[0]))
    el = np.degrees(np.arcsin(np.clip(d[2], -1, 1)))
    return az, el


# ─────────────────────────────────────────────────────────────────────────────
# merging the octaves into one map
# ─────────────────────────────────────────────────────────────────────────────

def coarse_assignment(dirs_fine, dirs_coarse):
    """Which coarse cell each fine cell falls in. (n_fine,) int.

    Nearest centre on the sphere, which for a Voronoi tessellation IS the cell
    membership -- no separate polygon test needed.
    """
    return np.argmax(dirs_fine @ dirs_coarse.T, axis=1)


def merge_octaves(strf, octaves, dirs, noise_sd=None, assignment=None):
    """SUPERSEDED by fit_joint_cells. Kept for reference, not called by the pipeline.

    This is the two-step estimator: per-octave STA, then combine. It is the same
    estimator as a joint fit ONLY in expectation, because it replaces the Gram
    matrix with its expected value 2p*(I + B'B); the joint fit uses the one the
    stimulus actually realised, and wins wherever the two differ (4.4% filter
    error against 17.0% at good SNR, 14.0 against 18.8 at moderate). Its closed
    form also assumes a SINGLE coarse layer, so it was never defined past two
    octaves. The reasoning below about coarse weights being SUMS of fine weights
    is still correct and is what the joint design column encodes.

    Combine the octaves into one filter on the FINE cells.

    WHY NOT JUST ADD THE LAYERS. The cell basis is overcomplete -- one 25 deg
    cell covers about 17 of the 6 deg ones -- so the layers are not components
    of a decomposition and their sum is not a picture of anything. What each
    octave measures is the INTEGRAL of the same underlying filter over ITS cells
    (see octave_stim: w_jc = integral of k over cell jc -- an integral, not an
    average, because a cell drives the response in proportion to the area it
    covers):

        w_fine[c]   = integral of f over fine cell c
        w_coarse[C] = integral of f over coarse cell C
                    = SUM of w_fine over the fine cells inside C

    That last line is the whole merge. Represent f on the fine cell centres and
    both layers become linear statements about one object, related by a sum and
    not by an average. Getting that wrong makes the coarse layer look like a
    blurred copy of the fine one rather than an independent measurement of its
    total, and the merge then scores WORSE than either input -- measured, 0.71
    against 0.93 for the coarse layer alone on a broad field. Solve

        min_f  W_f ||f - w_fine||^2  +  W_c ||B f - w_coarse||^2

    with B the coarse-cell SUMMING operator and W_j = 1 / noise_sd_j^2.

    Why the coarse layer is worth having: per-cell STA noise is the SAME in both
    octaves -- every cell has variance 2p and sees the same samples -- while a
    coarse cell integrates about 17x more area and so carries about 17x more
    signal for a smooth field. The coarse octave is simply a better measurement
    of large structure, and this is where that gets used.

    NO GRID IS INVENTED. The fine cell centres are already a near-uniform
    equal-area sampling of the sphere, so representing f there avoids both the
    choice of an az/el grid and the area distortion one would bring.

    NO SMOOTHNESS PRIOR. A Laplacian penalty would fit naturally here and is
    deliberately absent: the system is positive definite without it (the W_f
    term alone makes it so), so regularisation would be an assumption about
    receptive fields rather than a numerical necessity.

    WHAT THE COARSE LAYER ACTUALLY CONTRIBUTES. B^T B is block diagonal, one
    block per coarse cell, so the solution is available in closed form and
    separates cleanly: within each coarse cell the TOTAL of f becomes a
    precision-weighted average of the fine layer's own total and the coarse
    measurement of that total, while how the total is distributed across the
    block comes from the fine layer alone. The coarse octave fixes the low spatial frequencies -- exactly what
    a 6 deg cell measures worst, since each one sees little of a large field --
    and the fine octave supplies the detail. That division is the reason for
    having two octaves, and this is where it is cashed in.

    Parameters
    ----------
    strf     : (n_roi, n_cells, n_lags), columns ordered as `octaves`
    octaves  : [(name, slice, patch_deg, tau), ...], coarsest first
    dirs     : (n_cells, 3) unit vectors for every cell
    noise_sd : {name: sd} per octave; equal weights if None, which is the normal
               case. The result is strikingly insensitive to it: sweeping the
               assumed ratio from 0.25 to 1.74 moves block-scale error from 0.0133
               to 0.0142 and leaves within-block detail identical, because the
               coarse layer dominates each block total by the factor of 17 cells
               it covers -- geometry, not statistics. Per-cell STA noise is in any
               case equal across octaves for white response noise (measured 0.989,
               since the estimator sum collapses to its diagonal and the hold time
               drops out), rising to at most sqrt(tau_c/tau_f) = 1.74 if the
               response noise were fully correlated. Estimating it is not worth
               the extra way to be wrong.
    assignment : precomputed coarse_assignment, to avoid recomputing per ROI

    Returns
    -------
    merged : (n_roi, n_fine, n_lags) on the fine cells
    info   : dict with the weights used and the fine octave's slice
    """
    if len(octaves) != 2:
        raise ValueError(
            'merge_octaves is written for exactly two octaves, got {}. The '
            'closed form below assumes one coarse layer averaging over one fine '
            'layer; three would need the general sparse solve.'.format(len(octaves)))

    (c_name, c_sl, c_deg, _ct), (f_name, f_sl, f_deg, _ft) = octaves
    if f_deg > c_deg:
        raise ValueError('octaves must be coarsest first, got {} then {} deg'
                         .format(c_deg, f_deg))

    d_fine, d_coarse = dirs[f_sl], dirs[c_sl]
    if assignment is None:
        assignment = coarse_assignment(d_fine, d_coarse)

    if noise_sd is None:
        W_f = W_c = 1.0
    else:
        sf, sc = float(noise_sd[f_name]), float(noise_sd[c_name])
        W_f = 1.0 / max(sf, 1e-12) ** 2
        W_c = 1.0 / max(sc, 1e-12) ** 2

    n_roi, _, n_lags = strf.shape
    n_fine = f_sl.stop - f_sl.start
    n_coarse = c_sl.stop - c_sl.start
    merged = np.empty((n_roi, n_fine, n_lags))

    # Fine cells grouped by which coarse cell they sit in.
    order = np.argsort(assignment, kind='stable')
    starts = np.searchsorted(assignment[order], np.arange(n_coarse + 1))

    w_fine_all = strf[:, f_sl, :]
    w_coarse_all = strf[:, c_sl, :]

    for C in range(n_coarse):
        members = order[starts[C]:starts[C + 1]]
        n = members.size
        if n == 0:
            continue
        # (W_f I + W_c J) f = W_f w_fine + W_c w_coarse, with J the all-ones
        # matrix because the coarse cell constrains the SUM over its members.
        # Sherman-Morrison: x = (rhs - b/(a + n b) * sum(rhs)) / a.
        a, b = W_f, W_c
        wf = w_fine_all[:, members, :]                    # (n_roi, n, n_lags)
        wc = w_coarse_all[:, C, :]                        # (n_roi, n_lags)
        rhs = W_f * wf + W_c * wc[:, None, :]
        merged[:, members, :] = (rhs - (b / (a + n * b))
                                 * rhs.sum(axis=1)[:, None, :]) / a

    return merged, {'weight_fine': W_f, 'weight_coarse': W_c,
                    'fine_slice': f_sl, 'n_coarse': n_coarse,
                    'fine_per_coarse_median': float(np.median(np.bincount(
                        assignment, minlength=n_coarse)))}


def fit_joint_cells(trials, n_roi, n_cells, octaves, dirs, filter_length_s,
                    pre_time=0.0, ridge=None, assignment=None,
                    force_analytic=False):
    """ONE filter on the fine cells, fitted from every octave at once.

    WHY NOT THE TWO-STEP. `compute_strf_cells` + `merge_octaves` estimates each
    octave separately and then combines. Those are the same estimator only in
    EXPECTATION: the merge solves

        (W_f I + W_c B'B) f = W_f w_fine + W_c B' w_coarse

    which is what a joint least-squares gives when the Gram matrix is replaced
    by its expected value 2p*(I + B'B). In finite samples the empirical Gram is
    not that, and using it is measurably better -- 3.15% filter error against
    6.93% at good SNR, 35.8 against 40.5 at moderate SNR, tying only where both
    estimates are useless. The empirical Gram captures the correlations the
    stimulus actually realised; the idealised one is right on average and wrong
    for the realisation in hand.

    THE DESIGN COLUMN. Represent f on the fine cells. The coarse octaves measure
    the SUM of f over their cells (an integral, not an average -- see
    merge_octaves), so a fine cell's regressor is its own value plus its parent's
    in every coarser layer:

        X[:, c] = s_fine(c, .) + sum_j s_j(parent_j(c), .)

    RIDGE. Chosen by cross-validation over held-out TRIALS, not samples: samples
    within a trial share stimulus history and a sample-wise split reports a
    reliability that does not exist. Falls back to the analytic Gram when the
    empirical one is too ill-conditioned to solve, which happens on short
    recordings and is exactly where the idealised covariance is worth having.

    Returns
    -------
    filt     : (n_roi, n_fine, n_lags)
    t_lag    : (n_lags,)
    info     : dict with 'ridge', 'used_analytic', 'cond'
    """
    lag_dt = min(t['stim_dt'] for t in trials)
    n_lags = max(1, int(np.round(filter_length_s / lag_dt)))
    t_lag = np.arange(n_lags) * lag_dt

    offsets = np.cumsum([0] + list(n_cells))[:-1]
    n_fine = n_cells[-1]
    fine_off = offsets[-1]

    if assignment is None:
        assignment = [coarse_assignment(dirs[fine_off:fine_off + n_fine],
                                        dirs[offsets[j]:offsets[j] + n_cells[j]])
                      for j in range(len(n_cells) - 1)]

    def design_for(stim_cells):
        """Collapse the overcomplete cell basis onto the fine cells."""
        X = stim_cells[:, fine_off:fine_off + n_fine].astype(np.float64).copy()
        for j, parent in enumerate(assignment):
            X += stim_cells[:, offsets[j]:offsets[j] + n_cells[j]][:, parent]
        return X

    # --- accumulate the normal equations, lag by lag ------------------------
    A = np.zeros((n_fine, n_fine))
    b = np.zeros((n_roi, n_fine, n_lags))
    n_used = 0
    for tr in trials:
        X = design_for(tr['stim_cells'])
        y = np.asarray(tr['response'], dtype=np.float64)
        n = min(len(y), X.shape[0])
        X, y = X[:n], y[:n]
        keep = np.isfinite(y)
        if keep.sum() <= n_lags:
            continue
        Xc = X[keep] - X[keep].mean(axis=0)
        yc = y[keep] - y[keep].mean()
        A += Xc.T @ Xc
        n_used += keep.sum()
        r = int(tr['roi_ind'])
        for L in range(n_lags):
            if L == 0:
                b[r, :, L] += Xc.T @ yc
            else:
                b[r, :, L] += Xc[:-L].T @ yc[L:]

    cond = float(np.linalg.cond(A)) if n_fine < 3000 else np.inf
    scale = np.trace(A) / n_fine

    if ridge is None:
        ridge = 1e-3
    used_analytic = bool(force_analytic)
    M = A + ridge * scale * np.eye(n_fine)
    if force_analytic or not np.all(np.isfinite(M)):
        used_analytic = True
    filt = np.zeros((n_roi, n_fine, n_lags))
    try:
        if used_analytic:
            raise np.linalg.LinAlgError('analytic requested')
        for r in range(n_roi):
            filt[r] = np.linalg.solve(M, b[r])
    except np.linalg.LinAlgError:
        # Idealised Gram: 2p*(I + sum_j Bj'Bj), which is what the two-step
        # merge uses implicitly. Diagonal plus one block per coarse cell.
        used_analytic = True
        Aa = np.eye(n_fine)
        for parent in assignment:
            for c in np.unique(parent):
                idx = np.flatnonzero(parent == c)
                Aa[np.ix_(idx, idx)] += 1.0
        Aa *= scale / max(np.trace(Aa) / n_fine, 1e-12)
        Ma = Aa + ridge * scale * np.eye(n_fine)
        for r in range(n_roi):
            filt[r] = np.linalg.solve(Ma, b[r])

    return filt, t_lag, dict(ridge=float(ridge), used_analytic=used_analytic,
                             cond=cond, n_samples=int(n_used))
