# -*- coding: utf-8 -*-
"""Null-calibrated inclusion threshold for noise-stimulus filter analyses.

Applies to every noise analysis, not just the 1D bar STRF: the octave cell
basis and the full-field temporal filter have the same problem and the same
fix. Nothing here knows what the stimulus is -- the caller supplies a function
that recomputes its own filter with the stimulus shifted.

WHY NOT A FIXED z. A 1D STRF at 22 bar positions x 1200 lags is 26400 values
per ROI, and the reported statistic is the MAXIMUM over all of them. The max of
26400 standard normals sits near 4.0 and routinely exceeds 4.5, so a fixed z=3
admits essentially every ROI regardless of whether it carries signal. On Mi1
(fly_011 series_003) that turned a uniformly OFF-like cAMP population into an
apparently mixed one: every ROI whose peak came out positive had |z| below the
null level, while half the negative ones cleared it.

The right threshold is also not a constant across analyses -- the number of
comparisons behind the maximum depends on lag count, cell or position count,
and smoothing -- so it has to be measured per analysis.

HOW: CIRCULAR SHIFT. Roll the stimulus relative to the response by more than
the filter length. That destroys the stimulus-response pairing while preserving
the stimulus autocorrelation and the response's own temporal structure, neither
of which a plain shuffle keeps.

WHY NOT TAKE THE NULL FROM FAR LAGS INSTEAD. It is tempting, because far lags
are already computed and match the shift null closely (within 1-3% on Mi1).
But compute_strf_1d masks samples GLOBALLY -- `tvec >= pre_time +
filter_length_s`, so every lag shares one sample set and z stays comparable
across the lag axis -- which means extending the lag range to reach clean far
lags discards that much from the front of every epoch. A 30 s gap costs ~8% of
SNR on 300 s epochs. Shifts cost compute instead: (n_shifts + 1)x the filter
computation, offline and one-time. Compute is cheap and data is not, so shifts
are the better trade even though they look more expensive.

The threshold is an OUTPUT. Record it with the analysis, but never carry a
value into another one.

A NOTE ON SIGN. Do not sign-align ROIs to their own peak before averaging.
Aligning each ROI to its own maximum makes noise add constructively there, so
the average acquires a peak even when no ROI carries signal, and the bias grows
as the threshold gets looser. If a population is genuinely mixed in polarity,
split it and report the groups rather than folding them together.
"""
import numpy as np

#: Default shift count. min_detectable_k ~= 1/(q*n_shifts), so 4 shifts at
#: q=0.05 resolves any responsive set of >=5 ROIs. Costs 5x the filter compute.
DEFAULT_SHIFTS = 4


def peak_z(filt, noise_frac=0.3):
    """Per-ROI peak |z|, noise scale from the longest lags.

    filt : (n_roi, ..., n_lags). Works for the 1D bar STRF (n_roi, n_pos,
    n_lags), the octave cell basis (n_roi, n_cells, n_lags) and a plain
    temporal filter (n_roi, n_lags) alike -- only the first and last axes are
    interpreted.
    """
    filt = np.asarray(filt)
    n_lags = filt.shape[-1]
    late = slice(int(round((1.0 - noise_frac) * n_lags)), n_lags)
    out = np.full(filt.shape[0], np.nan)
    for r in range(filt.shape[0]):
        sd = np.nanstd(filt[r][..., late])
        if sd > 0 and np.isfinite(sd):
            out[r] = np.nanmax(np.abs(filt[r] / sd))
    return out


def shifts_for(duration_s, filter_length_s, n_shifts=DEFAULT_SHIFTS):
    """Evenly spaced shifts, all clear of the filter and of the wrap point.

    Every shift must exceed filter_length_s or real filter survives into the
    "null" and the threshold comes out too high; each must also stay that far
    from the full duration, since a circular roll makes the far end adjacent to
    zero again.
    """
    lo, hi = float(filter_length_s), float(duration_s) - float(filter_length_s)
    if hi <= lo:
        raise ValueError(
            'cannot place shifts: duration {:.1f}s leaves no room clear of a '
            '{:.1f}s filter at both ends'.format(duration_s, filter_length_s))
    return list(np.linspace(lo, hi, int(n_shifts) + 2)[1:-1])


def calibrate(make_filter, shifts, percentile=95.0, noise_frac=0.3):
    """Null from circular shifts of the stimulus.

    make_filter : callable(shift_s) -> (n_roi, ..., n_lags)
        Recomputes the caller's own filter with the stimulus rolled by
        shift_s seconds. Each pipeline supplies its own; this module does not
        need to know the stimulus representation.
    shifts : sequence of float
        Use shifts_for() unless there is a reason not to.
    """
    shifts = list(shifts)
    if not shifts:
        raise ValueError('need at least one shift; several is better')
    draws = [peak_z(make_filter(s), noise_frac) for s in shifts]
    per_shift = [dict(shift_s=float(s), n=int(np.isfinite(d).sum()),
                      p95=float(np.nanpercentile(d, 95)))
                 for s, d in zip(shifts, draws)]
    null = np.concatenate(draws)
    null = null[np.isfinite(null)]
    if null.size == 0:
        raise ValueError('null contained no usable ROIs')
    return dict(threshold=float(np.percentile(null, percentile)), null=null,
                shifts=shifts, n_shifts=len(shifts), per_shift=per_shift,
                percentile=float(percentile), n_null=int(null.size),
                noise_frac=float(noise_frac))


def pvalues(observed, null):
    """Upper-tail p per ROI against the pooled null, empirically.

    p cannot go below the floor 1/(n_null + 1). That only bites in the SPARSE
    regime: BH rejects the largest k with p_(k) <= (k/m) q and that threshold
    grows with k, so 40 of 59 ROIs at a floor of 0.0042 still pass at k = 40.
    Detecting a handful out of many is what the floor blocks -- see
    min_detectable_k, and add shifts if it matters.
    """
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if null.size == 0:
        raise ValueError('empty null')
    p = np.array([(np.sum(null >= o) + 1.0) / (null.size + 1.0)
                  for o in observed])
    floor = 1.0 / (null.size + 1.0)
    return p, dict(floor=float(floor), n_null=int(null.size),
                   min_detectable_k=int(np.ceil(observed.size * floor / 0.05)))


def benjamini_hochberg(p, q=0.05):
    """Indices surviving BH at false-discovery rate q, and the cutoff used."""
    p = np.asarray(p, dtype=float)
    m = p.size
    order = np.argsort(p)
    passing = p[order] <= (np.arange(1, m + 1) / float(m)) * q
    if not passing.any():
        return np.array([], dtype=int), 0.0
    k = int(np.max(np.flatnonzero(passing)))
    return np.sort(order[:k + 1]), float(p[order][k])


def select_rois(filt, cal, method='auto', q=0.05, noise_frac=0.3):
    """Apply the null: FDR by default, percentile only when FDR cannot work.

    There is deliberately no minimum-ROI rule. The resolution limit is
    k >= m/((n_null + 1) q), and with the null pooled over ROIs the m cancels,
    leaving min_detectable_k ~= 1/(q * n_shifts) -- independent of ROI count,
    so FDR is as usable at 10 ROIs as at 59. 'auto' falls back to the
    percentile only when that exceeds the ROI count, i.e. when BH could not
    reject anything however strong the effect.

    Returns (keep, signs, info). `signs` is for CHECKING polarity uniformity,
    never for flipping traces before averaging.
    """
    z = peak_z(filt, noise_frac)
    ok = np.isfinite(z)
    m = int(ok.sum())
    n_shifts = int(cal.get('n_shifts', len(cal.get('shifts', [1]))))
    resolvable = int(np.ceil(1.0 / (q * max(n_shifts, 1))))
    use_fdr = method == 'fdr' or (method == 'auto' and resolvable <= m)
    info = dict(n_roi=m, n_shifts=n_shifts, min_detectable_k=resolvable,
                q=float(q), method_requested=method)

    if use_fdr:
        p, meta = pvalues(z[ok], cal['null'])
        keep_ok, cut = benjamini_hochberg(p, q=q)
        keep = np.flatnonzero(ok)[keep_ok] if len(keep_ok) else np.array([], int)
        info.update(meta, method_used='fdr', p_cutoff=float(cut),
                    threshold=float(np.min(z[keep])) if len(keep) else np.nan)
        if len(keep) <= resolvable:
            info['warning'] = (
                'surviving set ({}) is at or below what {} shift(s) can resolve '
                '({}); add shifts before trusting this'
                .format(len(keep), n_shifts, resolvable))
    else:
        keep = np.flatnonzero(ok & (z >= cal['threshold']))
        info.update(method_used='percentile', threshold=float(cal['threshold']),
                    reason='min_detectable_k={} > n_roi={}'.format(resolvable, m))

    filt = np.asarray(filt)
    n_lags = filt.shape[-1]
    late = slice(int(round((1.0 - noise_frac) * n_lags)), n_lags)
    signs = np.zeros(len(keep))
    for i, r in enumerate(keep):
        zz = filt[r] / np.nanstd(filt[r][..., late])
        signs[i] = np.sign(zz.ravel()[int(np.nanargmax(np.abs(zz)))])
    return keep, signs, info
