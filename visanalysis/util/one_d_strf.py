"""
one_d_strf.py

The 1D analogue of octave_strf.py for the single-layer bar/wedge ternary
stimulus (labpack SphericalOneDTernaryNoise).

REUSE, NOT DUPLICATION. Reverse correlation in a cell basis with a common
per-cell variance, and deconvolving a single exponential hold, have nothing
octave-specific about them -- they are the same computation with "one layer,
one slice covering every cell" instead of several nested ones. Rather than
re-derive and re-validate that numerics a second time, this module reuses it
directly from octave_strf:

    compute_strf_cells, compute_strf_cells_folds, peak_lag,
    independent_samples, deconvolve_octaves, zscore_by_octave

`one_layer(meta)` below builds the `octaves`-shaped argument
(`deconvolve_octaves`/`zscore_by_octave` expect a list of
`(name, slice, patch_deg, tau)`) those functions need for this stimulus's one
layer.

WHAT IS ACTUALLY DIFFERENT FROM OCTAVE
---------------------------------------
* One layer, so there is no merge and no fit_joint_cells. Bars/wedges are
  disjoint rings or sectors, not nested scales, so the cell basis here is NOT
  overcomplete: `compute_strf_cells` already IS the final per-cell estimate,
  with nothing left to combine across bases the way octave's coarse and fine
  layers have to be.

* Position is 1D and can wrap. A cell is a ring or a sector, not a point on
  the sphere, so there is no direction to average and
  `octave_strf.spherical_centroid` does not apply. `centroid_1d` below
  averages the BINNED COORDINATE instead -- circularly for wedges, since a
  family straddling the 0/360 seam would otherwise average toward its own
  antipode, exactly the reason `spherical_centroid` averages unit vectors
  rather than raw az/el in the octave case.
"""

import warnings

import numpy as np

from visanalysis.util.octave_strf import (  # noqa: F401 -- re-exported; see module docstring
    compute_strf_cells, compute_strf_cells_folds, peak_lag,
    independent_samples, deconvolve_octaves, zscore_by_octave,
)


def one_layer(meta, name='1d'):
    """The `octaves`-shaped list `deconvolve_octaves`/`zscore_by_octave`
    expect, for this stimulus's one layer: every cell, one hold_tau.

    `patch_deg` (the third field of each tuple) is meaningless for a bar or
    wedge and carried as None -- neither function reads it for anything but
    display, which this pipeline does differently (see `analyze_data_strf_
    one_d.py`), so nothing downstream depends on it being a real angle.
    """
    n = int(meta['n_cells'])
    return [(str(name), slice(0, n), None, float(meta['hold_tau']))]


def zscore_from_noise_cells(strf, visible_mask, noise_frac=0.3, min_noise_cells=3):
    """z-score each ROI against a noise estimate that carries no signal.

    WHY NOT octave_strf.zscore_by_octave. That takes the sd over the whole
    (n_cells x n_lags) block, which IS a clean noise estimate for the octave
    stimulus: roughly three quarters of its ~1200 cells are never displayed and
    almost none of the rest are driven, so the block is noise with a little
    signal in it. Here there are a few tens of cells and a large fraction of the
    visible ones can be driven, so the same sd absorbs the signal it is meant to
    normalise. z then reads systematically LOW, and worst exactly for the ROIs
    with the strongest filters -- the opposite of what a threshold needs.

    The never-visible cells are the natural pool: the projector could not draw
    them, so their STA is noise by construction, at the same per-cell variance
    as every other column. When there are too few of them -- a cell family that
    falls entirely inside the patch -- fall back to the late-lag tail, which is
    signal-free for any filter shorter than the lag axis. That is the same
    convention nullcal.peak_z uses, so the two agree about what noise is.

    Returns (strf_z, info); `info['how']` names the pool actually used, because
    the two regimes are not comparable and a reader cannot tell them apart from
    the sd alone.
    """
    strf = np.asarray(strf)
    n_roi, n_cells, n_lags = strf.shape
    noise_cells = ~np.asarray(visible_mask, dtype=bool)
    if int(noise_cells.sum()) >= int(min_noise_cells):
        pool, how = strf[:, noise_cells, :], 'never-visible cells'
    else:
        late = slice(int(round((1.0 - noise_frac) * n_lags)), n_lags)
        pool, how = strf[:, :, late], 'late-lag tail'
    # nanstd, not std: one NaN anywhere in an roi's pool makes the whole sd
    # NaN, which fails the `sd > 0` test below and silently substitutes 1.0 --
    # leaving that roi's 'z' as raw filter units, on the same axis and colour
    # scale as everyone else's z. Counted rather than merely handled, because a
    # substitution that nobody sees is the part that was wrong.
    # catch_warnings, not errstate: an all-NaN pool makes nanstd emit a
    # 'Degrees of freedom <= 0' RuntimeWarning through the warnings module,
    # which errstate does not cover. That case is counted and reported below,
    # so the warning is duplicate noise in the analysis log.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        sd = np.nanstd(pool.reshape(n_roi, -1), axis=1)
    usable = np.isfinite(sd) & (sd > 0)
    sd = np.where(usable, sd, 1.0)
    return strf / sd[:, None, None], {
        'how': how, 'median_sd': float(np.median(sd[usable]))
                    if np.any(usable) else float('nan'),
        'n_sd_fallback': int((~usable).sum()),
        # per-ROI, because the centroid's inclusion rule needs THIS roi's
        # noise rather than the population's
        'sd': sd,
        'n_noise_cells': int(noise_cells.sum())}


def _contiguous_lobe(w, pos, thresh, circular):
    """Cells above `thresh` in the unbroken run containing the peak.

    Walked in POSITION order, not index order: for wedges those differ whenever
    the family's seam falls inside the data. `circular` lets the run wrap, which
    costs nothing when it cannot -- the lit arc is contiguous in the screen
    frame, so a wrap only happens if the lobe really does straddle the ends.
    """
    order = np.argsort(pos)
    ws = w[order]
    n = ws.size
    pk = int(np.argmax(ws))
    keep_sorted = np.zeros(n, dtype=bool)
    if ws[pk] <= thresh:
        return np.zeros(n, dtype=bool)
    keep_sorted[pk] = True
    for step in (-1, 1):
        i = pk + step
        while 0 <= i < n and ws[i] > thresh:
            keep_sorted[i] = True
            i += step
        if circular and (i < 0 or i >= n):
            j = (i + n) % n
            while ws[j] > thresh and not keep_sorted[j]:
                keep_sorted[j] = True
                j = (j + step + n) % n
    out = np.zeros(n, dtype=bool)
    out[order] = keep_sorted
    return out


def centroid_1d(weights, positions_deg, mode, frac=0.5, cell_widths_deg=None,
                noise_sd=None, n_sigma=2.0):
    """Weighted center (and spread) of the dominant lobe along the stimulus's
    own axis.

    weights      : (n_cells,) non-negative importance, e.g. the dominant-lobe
                   filter at a ROI's own peak lag -- same "keep the polarity of
                   the strongest cell, drop the other sign" convention
                   `octave_strf.spherical_centroid` uses, and for the same
                   reason: a center-surround filter has both signs at one lag,
                   and averaging |w| would put the centroid between them,
                   which belongs to neither lobe.
    positions_deg: (n_cells,) from `one_d_stim.cell_positions_deg`.
    mode         : 'bars' or 'wedges'. Bars run pole to antipole and do not
                   wrap, so the mean is a plain weighted average. Wedges wrap
                   at 360, so the center is the angle of the weighted mean of
                   UNIT VECTORS on the circle -- averaging raw degrees would
                   be wrong exactly at the seam, the 1D version of why
                   `spherical_centroid` averages directions rather than az/el.
    cell_widths_deg : (n_cells,) each cell's extent, from
                   `one_d_stim.cell_widths_deg`. Optional, and worth passing.
                   Without it `spread_deg` is the weighted sd of the cell
                   CENTRES alone, which reads exactly 0 when a single cell
                   survives `frac` -- a receptive field cannot be narrower
                   than the cell that measured it, so 0 is not a width, it is
                   the absence of one. With it, the within-cell term is added:
                   the measured distribution is a mixture of uniforms, so its
                   variance is the variance of the centers plus the mean
                   within-cell variance w^2/12. That also corrects the
                   systematic under-report in the multi-cell case, not just
                   the single-cell one. For wedges it is added to the circular
                   variance, which is exact only in the small-angle limit and
                   close enough at any sane cell width.
    noise_sd     : this ROI's noise scale, in the same units as `weights`
                   (see zscore_from_noise_cells, which returns a per-ROI
                   'sd'). When given it REPLACES `frac`: cells are kept if
                   they clear `n_sigma * noise_sd` AND lie in the
                   contiguous run containing the peak.

                   Why not the fraction. A peak-relative cut of 0.5 leaves
                   ONE cell whenever the lobe is narrow, and a weighted mean
                   of one cell is that cell's center -- the estimate snaps to
                   the grid and sub-cell position is thrown away. Dropping the
                   cut entirely is exact when clean but worse than 0.5 once
                   the tails carry noise. A noise-relative cut adapts: it
                   keeps the whole lobe at high SNR and only the peak at low.
                   Measured against known centers on 15 deg cells, rms error
                   in degrees:

                       noise   frac0.5  frac0.2  frac0.0  noise>2sd
                        0.00      2.25     1.11     0.01       0.01
                        0.05      2.32     1.28     1.50       1.20
                        0.15      3.25     3.95     4.03       3.09
                        0.25      5.47     6.10     5.87       5.25

                   -- best or tied everywhere, which neither fixed cut is.
                   The contiguity requirement is what bounds it: a noise
                   excursion far from the peak cannot join the average even
                   when it clears the threshold.
    n_sigma      : multiples of noise_sd a cell must clear. 2.0 by default.
    frac         : used ONLY when noise_sd is None. Cells below this fraction
                   of the peak weight are dropped,
                   matching `spherical_centroid`'s convention, so the noise
                   floor in the tails does not drag the center toward the
                   array center.

    Returns (center_deg, spread_deg, n_kept). `center_deg` is in the SAME
    branch as the `positions_deg` passed in -- for wedges those span
    `phase_deg .. phase_deg + 360` rather than [0, 360), so the center is
    wrapped to match and stays directly comparable to them.
    For wedges `spread_deg` is the
    circular standard deviation (`sqrt(-2 ln R)` from the resultant length R,
    the standard wrapped-normal approximation); for bars it is the ordinary
    weighted standard deviation. Both are a size, not a direction -- the 1D
    equivalent of a receptive field width, which each octave cell does not
    need to report because it already IS a point.

    Returns (nan, nan, 0) if every weight is zero or non-finite.
    """
    w = np.asarray(weights, dtype=float)
    w = np.where(np.isfinite(w), w, 0.0)
    if w.max() <= 0:
        return np.nan, np.nan, 0
    if noise_sd is not None and np.isfinite(noise_sd) and noise_sd > 0:
        keep = _contiguous_lobe(w, np.asarray(positions_deg, dtype=float),
                                float(n_sigma) * float(noise_sd),
                                circular=(str(mode) == 'wedges'))
        if not keep.any():
            keep = w >= w.max()          # nothing cleared noise: the peak alone
    else:
        keep = w >= frac * w.max()
    wk = w[keep]
    pk = np.asarray(positions_deg, dtype=float)[keep]
    if wk.sum() <= 0:
        return np.nan, np.nan, int(keep.sum())

    # Mean within-cell variance of a uniform over each kept cell, w^2/12.
    within_var = 0.0
    if cell_widths_deg is not None:
        wid = np.asarray(cell_widths_deg, dtype=float)[keep]
        within_var = float(np.sum(wk * wid ** 2) / wk.sum() / 12.0)

    if str(mode) == 'wedges':
        theta = np.radians(pk)
        c = np.array([np.sum(wk * np.cos(theta)), np.sum(wk * np.sin(theta))])
        r = np.linalg.norm(c)
        if r == 0:
            return np.nan, np.nan, int(keep.sum())
        # Wrap into the same branch as `positions_deg`, NOT into [0, 360).
        # one_d_tessellation.cell_bounds returns wedge bounds spanning
        # phase_deg .. phase_deg + 360 un-wrapped, so a bare % 360 puts the center
        # in a different frame from the positions it is reported beside: at
        # phase_deg = 82.5 the positions run 90..435 and a lobe on the last cells
        # comes back as 68.8 -- off the plotted axis, and not comparable to the
        # positions_deg saved next to it.
        base = float(np.min(np.asarray(positions_deg, dtype=float)))
        center = base + (float(np.degrees(np.arctan2(c[1], c[0]))) - base) % 360.0
        R = np.clip(r / wk.sum(), 1e-12, 1.0)
        circ_var = float(np.degrees(np.sqrt(-2.0 * np.log(R)))) ** 2
        spread = float(np.sqrt(circ_var + within_var))
    else:
        center = float(np.sum(wk * pk) / wk.sum())
        between_var = float(np.sum(wk * (pk - center) ** 2) / wk.sum())
        spread = float(np.sqrt(between_var + within_var))
    return center, spread, int(keep.sum())


# ─────────────────────────────────────────────────────────────────────────────
# repeated trials: is the response the same when the stimulus is?
# ─────────────────────────────────────────────────────────────────────────────

def _corr(a, b):
    """Pearson r over the finite entries, NaN when it would be meaningless."""
    a, b = np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 3 or np.std(a[ok]) == 0 or np.std(b[ok]) == 0:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def repeat_groups(trials):
    """Trials that showed frame-for-frame the same stimulus, grouped by seed.

    Returns [(seed, [trial_ind, ...]), ...] for every seed shown more than once,
    largest group first; empty when nothing repeats.

    Keyed on the SEED ALONE, unlike fullfield_temporal.repeat_groups which keys
    on the whole metadata dict. That is safe here and would not be there: these
    trials have already been grouped by the caller's geometry key, so every
    field of the stimulus except the seed is identical by construction, and one
    seed therefore IS one stimulus frame for frame. Repeats arise when the
    protocol runs with repeat_stimulus=True, which makes trial_seed() hand back
    the run seed on every trial instead of deriving a fresh one.
    """
    seen = {}
    for t in trials:
        ti = int(t['trial_ind'])
        if ti not in seen:
            seen[ti] = int(t['seed'])
    by_seed = {}
    for ti, seed in sorted(seen.items()):
        by_seed.setdefault(seed, []).append(ti)
    reps = [(s, v) for s, v in by_seed.items() if len(v) > 1]
    return sorted(reps, key=lambda sv: len(sv[1]), reverse=True)


def _filter_support(block, z=3.0, noise_frac=0.3, min_lags=3):
    """Lags up to the last one where a (n_cells, n_lags) filter clears its noise.

    Correlating trials over the whole lag window would score mostly the tail,
    where every trial is noise and two good trials look no more alike than two
    bad ones. Noise scale is the longest lags' sd, as in nullcal.peak_z and
    fullfield_temporal._filter_support -- the 1D difference is only that the
    statistic is summed across cells first, so one strong cell does not have to
    carry the whole support on its own.
    """
    block = np.asarray(block, dtype=float)
    n_lags = block.shape[-1]
    late = slice(int(round((1.0 - noise_frac) * n_lags)), n_lags)
    sd = np.nanstd(block[..., late])
    if not np.isfinite(sd) or sd <= 0:
        return n_lags
    power = np.sqrt(np.nansum(block ** 2, axis=0))          # across cells
    above = np.flatnonzero(power > z * sd * np.sqrt(block.shape[0]))
    return int(max(min_lags, above[-1] + 1)) if above.size else int(min_lags)


def trial_reliability(trials, n_roi, n_cells, strf_full, filter_length_s,
                      pre_time=0.0, visible_mask=None, lag_dt=None):
    """Per-ROI consistency across trials that saw the SAME stimulus.

    Only meaningful when the stimulus repeats: then trial-to-trial differences
    are the cell's own -- drift, run-down, state -- rather than a different
    noise draw sampling the filter differently. Returns None when nothing
    repeats, so a caller can simply skip the section.

    Measures and names mirror fullfield_temporal.trial_reliability so the two
    pipelines report the same quantities:

      filter_loo_r      corr of each trial's STRF with the mean of the others,
                        over the visible cells and the support lags
      trial_gain        least-squares scale of each trial's STRF onto the full
                        one: 1 is typical, a trend across trials is run-down
      response_pair_r   mean correlation between single-trial RESPONSES on a
                        common stimulus-time grid -- model-free, so it stands
                        even if the linear filter is a poor description
      noise_ceiling     sqrt(response_pair_r), NaN unless response_pair_r is
                        more than 2 SE above zero: the most any model could
                        explain of single-trial data
      trial_time_courses / full_time_course
                        each trial's STRF, and the all-trial STRF, projected
                        onto that ROI's own spatial profile -- one time course
                        per trial, which is what the repeatability panels draw

    Per-trial STRFs are NOT returned. The fullfield version returns its
    per-trial filters because they are one vector per ROI; here each carries a
    cell axis as well, so a 10-trial group at 1200 lags is tens of MB per
    series and the summary statistics are what gets read anyway.
    """
    groups = repeat_groups(trials)
    if not groups:
        return None

    if visible_mask is None:
        visible_mask = np.ones(n_cells, dtype=bool)
    visible_mask = np.asarray(visible_mask, dtype=bool)

    by_trial = {}
    for t in trials:
        by_trial.setdefault(int(t['trial_ind']), []).append(t)

    out = []
    for seed, tis in groups:
        n_tr = len(tis)

        # One STA per trial, same lags and same estimator as the full fit, so
        # the trials are comparable to it and to each other.
        tf = np.full((n_tr, n_roi, int(visible_mask.sum()),
                      strf_full.shape[-1]), np.nan)
        for j, ti in enumerate(tis):
            # same lag axis as the full fit, or the per-trial filters would
            # not be comparable to it or to each other
            s = compute_strf_cells(by_trial[ti], n_roi, n_cells,
                                   filter_length_s, pre_time, lag_dt)[0]
            tf[j] = s[:, visible_mask, :]

        # One time course per trial per ROI: the trial's STRF projected onto the
        # full STRF's spatial profile at its peak lag. Stored INSTEAD of the
        # per-trial (cell x lag) blocks -- it is the object the repeatability
        # panels actually draw, it matches the shape the full-field pipeline
        # plots, and it is the same reduction the deconvolution CV performs
        # (one well-driven time series per ROI), so the two agree about what
        # 'this ROI's filter' means.
        n_lags = strf_full.shape[-1]
        tc = np.full((n_tr, n_roi, n_lags), np.nan)
        full_tc = np.full((n_roi, n_lags), np.nan)
        for r in range(n_roi):
            ref = strf_full[r][visible_mask]
            prof = ref[:, int(np.argmax((ref ** 2).sum(axis=0)))]
            nrm = float(np.linalg.norm(prof))
            if nrm == 0 or not np.isfinite(nrm):
                continue
            # Oriented so the DOMINANT CELL is positive, which leaves the
            # projection carrying the roi's own polarity: an OFF cell reads as a
            # negative-going time course instead of being silently flipped to
            # look like an ON one. Projecting onto an unoriented profile makes
            # the value at the peak lag positive by construction, for every roi.
            # This is per-roi orientation for DISPLAY, not sign-alignment before
            # averaging -- see util.nullcal on why the latter manufactures peaks.
            prof = prof / nrm
            prof = prof * np.sign(prof[int(np.argmax(np.abs(prof)))])
            full_tc[r] = prof @ ref
            for j in range(n_tr):
                tc[j, r] = prof @ tf[j, r]

        support = np.array([_filter_support(strf_full[r][visible_mask])
                            for r in range(n_roi)])
        loo = np.full((n_tr, n_roi), np.nan)
        gain = np.full((n_tr, n_roi), np.nan)
        for r in range(n_roi):
            L = int(support[r])
            ref = strf_full[r][visible_mask][:, :L]
            den = float(np.nansum(ref * ref))
            for j in range(n_tr):
                others = np.nanmean(np.delete(tf[:, r, :, :L], j, axis=0), axis=0)
                loo[j, r] = _corr(tf[j, r, :, :L], others)
                if den > 0 and np.all(np.isfinite(tf[j, r, :, :L])):
                    gain[j, r] = float(np.nansum(tf[j, r, :, :L] * ref) / den)

        # Model-free response repeatability. Trials start at an arbitrary phase
        # relative to the imaging clock, so samples go onto one stimulus-time
        # grid before trials are compared -- and the grid is built per ROI,
        # since a per-ROI z-offset shifts its time vector.
        resp = {}
        for j, ti in enumerate(tis):
            for t in by_trial[ti]:
                resp[j, int(t['roi_ind'])] = (
                    np.asarray(t['time_vector'], dtype=float) - pre_time,
                    np.asarray(t['response'], dtype=float))

        pair = np.full(n_roi, np.nan)
        pair_se = np.full(n_roi, np.nan)
        resp_loo = np.full((n_tr, n_roi), np.nan)
        for r in range(n_roi):
            have = [j for j in range(n_tr) if (j, r) in resp]
            if len(have) < 2:
                continue
            t_lo = max(resp[j, r][0][0] for j in have)
            t_hi = min(resp[j, r][0][-1] for j in have)
            if not np.isfinite(t_lo) or not np.isfinite(t_hi) or t_hi <= t_lo:
                continue
            dts = [np.median(np.diff(resp[j, r][0])) for j in have]
            grid = np.arange(t_lo, t_hi, float(np.median(dts)))
            if grid.size < 3:
                continue
            R = np.full((n_tr, grid.size), np.nan)
            for j in have:
                tv, y = resp[j, r]
                n = min(len(tv), len(y))
                ok = np.isfinite(y[:n])
                if ok.sum() >= 3:
                    R[j] = np.interp(grid, tv[:n][ok], y[:n][ok])
            rs = np.array([_corr(R[a], R[b]) for a in have
                           for b in have if b > a])
            rs = rs[np.isfinite(rs)]
            if rs.size:
                pair[r] = rs.mean()
                # Spread across pairs, not a formula in the sample count: the
                # samples are autocorrelated, so their nominal number overstates
                # how many independent comparisons each correlation rests on.
                pair_se[r] = (rs.std(ddof=1) / np.sqrt(rs.size)
                              if rs.size > 1 else np.nan)
            for j in have:
                resp_loo[j, r] = _corr(
                    R[j], np.nanmean(np.delete(R, j, axis=0), axis=0))

        out.append(dict(
            seed=int(seed), trials=np.asarray(tis, dtype=int),
            trial_time_courses=tc, full_time_course=full_tc,
            support_lags=support, filter_loo_r=loo, trial_gain=gain,
            response_pair_r=pair, response_pair_r_se=pair_se,
            response_loo_r=resp_loo,
            # NaN unless the repeatability is resolved (> 2 SE above zero). A
            # ceiling read off noise can come out at or below the model's own
            # score, which is not a statement about the cell.
            noise_ceiling=np.where(pair > 2 * pair_se,
                                   np.sqrt(np.clip(pair, 0.0, None)), np.nan)))
    return out
