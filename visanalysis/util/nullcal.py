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

WHICH SIGN. `peak_signs` reads the sign at the |z| extremum, which is the
detection statistic and not the polarity: see the POLARITY section at the end
for why the FIRST lobe is the polarity and `first_lobe_polarity` for how it
is read.
"""
import numpy as np

#: Default shift count, shared by every noise analysis (jcsimon, 2026-11-16).
#: The pooled null has n_roi * n_shifts draws and p cannot go below
#: 1/(n_null + 1), so with few shifts BH cannot start on a family with only a
#: handful of responsive ROIs: at 10 shifts and 213 ROIs the floor is 0.00047
#: and the first ROI would need p <= 0.00023. 20 halves the floor. Costs 21x
#: the filter compute; past ~20 the shifts of a 10 s filter in a 10 min
#: recording start to sample overlapping stretches and stop being independent.
DEFAULT_SHIFTS = 20


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
                    threshold=float(np.min(z[keep])) if len(keep) else np.nan,
                    threshold_kind='min z of kept ROIs; descriptive only, '
                                   'the cutoff is on p, not on z')
        if len(keep) <= resolvable:
            info['warning'] = (
                'surviving set ({}) is at or below what {} shift(s) can resolve '
                '({}); add shifts before trusting this'
                .format(len(keep), n_shifts, resolvable))
    else:
        keep = np.flatnonzero(ok & (z >= cal['threshold']))
        info.update(method_used='percentile', threshold=float(cal['threshold']),
                    threshold_kind='{:g}th percentile of the null; a real '
                                   'cutoff on z'.format(
                                       cal.get('percentile', 95.0)),
                    reason='min_detectable_k={} > n_roi={}'.format(resolvable, m))

    return keep, peak_signs(filt, keep, noise_frac), info


def peak_signs(filt, keep, noise_frac=0.3):
    """Sign of each kept ROI's filter at its own |z| extremum.

    Split out of select_rois so a kept set decided on ANOTHER channel can be
    given this channel's own polarities: the set says which cells responded,
    the sign has to be read from the indicator actually being reported. Same
    noise estimate as peak_z (the late lags), same rule as before.
    """
    filt = np.asarray(filt)
    n_lags = filt.shape[-1]
    late = slice(int(round((1.0 - noise_frac) * n_lags)), n_lags)
    keep = np.asarray(keep, dtype=int)
    signs = np.zeros(len(keep))
    for i, r in enumerate(keep):
        zz = filt[r] / np.nanstd(filt[r][..., late])
        signs[i] = np.sign(zz.ravel()[int(np.nanargmax(np.abs(zz)))])
    return signs


# ─────────────────────────────────────────────────────────────────────────────
# POLARITY: which lobe comes FIRST, not which is biggest
# ─────────────────────────────────────────────────────────────────────────────
# peak_signs above reads the sign at the global |z| extremum. That is the right
# statistic for DETECTION (it is what the shift null calibrates) and the wrong
# one for POLARITY: an ON cell's calcium filter is a fast positive lobe followed
# by a slower negative rebound, and whichever of the two is larger after the
# stimulus hold has smeared them wins the extremum. At hold_tau 0.4 s the smear
# favours the slow lobe, and on JS140_x_JS261 fly_022 the extremum rule called
# 21 and 28 of 75 and 72 kept Mi1 rois negative; the same rois in the OTHER pole
# family (independent data, same cells) disagreed with themselves 8 times in 61.
# Polarity is a property of the ONSET, so it is read from the first lobe, and
# the registration lag -- where the spatial profile is centroided -- is the
# first lobe's peak, because that is the direct response; a delayed surround or
# an adaptation rebound can have a different spatial footprint.
#
# HOW, and why not simply "first run above threshold". Every per-roi onset
# detector tried on fly_022 sat on one trade-off: strict enough to ignore a
# noise blip at lag 0 means missing a weak first lobe and flipping the sign
# systematically (that is the extremum rule again); permissive enough to catch
# it means a noise blip can take the call (a 2.5 sd blip at lag 0 flipped 17% of
# rois). What has the SNR to decide is the whole waveform: the family's kept
# rois share one kinetics, so the sign-aligned mean trace of the family is a
# template, the projection of each roi's trace onto it is a whole-filter
# likelihood ratio between the two orientations, and the template's own first
# lobe -- at family SNR, trivially resolved -- orients the result. Measured on
# fly_022: 61 of 61 cross-family sign agreements (extremum rule 53), 1-2
# minority rois per Mi1 family (21-28), 3 of 119 blip flips, and the same
# answer from either initialisation. The per-roi lobe detector survives as the
# INITIALISATION, as the LAG finder inside the template's first-lobe window,
# and as the fallback for families too small to carry a template.
#
# On sign-aligning: the module docstring says never to sign-align rois before
# averaging, because noise then adds constructively at every peak. The template
# IS a sign-aligned mean -- but it is a classifier built from rois the null
# already kept, used only to decide each roi's orientation, and never drawn or
# reported as a population filter. A family whose kept set is noise gets a
# noise template and random signs, which is what noise deserves; the registered
# means downstream are still built un-flipped, split by polarity group.


def _smooth(x, width):
    """Centred boxcar, same length; width <= 1 is a no-op."""
    x = np.asarray(x, dtype=float)
    width = int(width)
    if width <= 1:
        return x
    return np.convolve(x, np.ones(width) / width, mode='same')


def _runs_of_sign(x):
    """[start, stop) of every maximal run of constant non-zero sign."""
    s = np.sign(x)
    out, start = [], 0
    for i in range(1, x.size + 1):
        if i == x.size or s[i] != s[start]:
            if s[start] != 0:
                out.append((start, i))
            start = i
    return out


def _area_sd(late):
    """sd of the sum of m consecutive samples of noise with the autocovariance
    of `late`, as a function of m.

    Not sd * sqrt(m): the filter noise is correlated across lags -- through the
    stimulus hold, which puts the stimulus autocorrelation on every STA column,
    and through any smoothing applied before the runs are found. At hold_tau
    0.4 s and a 5-sample boxcar the sd of a 15-sample area is 3.6x the
    white-noise figure, so a white threshold at 3 sd is a real threshold at
    under 1 sd. The autocovariance is taken from the late lags, which are
    signal-free for any filter shorter than the lag axis, and truncated at a
    third of their length.
    """
    x = np.asarray(late, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 4:
        return lambda m: np.nan
    x = x - x.mean()
    max_lag = max(1, n // 3)
    acov = np.array([np.dot(x[:n - k], x[k:]) / n for k in range(max_lag)])
    # Bartlett taper. The raw truncated autocovariance is not a valid one --
    # its implied variance of a long sum can come out NEGATIVE (it did, on a
    # 230-sample run: sd 1e-150, z 1e147, and that run became the first lobe).
    # Tapered by (1 - k/max_lag) it is positive semidefinite, so the variance
    # of any sum is >= 0, and it grows linearly past max_lag as it should.
    acov = acov * (1.0 - np.arange(max_lag) / float(max_lag))

    def f(m):
        m = int(m)
        ks = np.arange(1, min(m, acov.size))
        v = m * acov[0] + 2.0 * np.sum((m - ks) * acov[ks])
        # never below one sample's own variance, whatever the taper left
        return float(np.sqrt(max(v, acov[0], 1e-300)))
    return f


def energy_peak_lag(block):
    """Lag of greatest summed squared value across cells, for one
    (n_cells, n_lags) filter. Squared, not absolute: with many pure-noise cells
    the |.| sum is dominated by the floor and picks an arbitrary lag."""
    block = np.asarray(block, dtype=float)
    e = np.nansum(block ** 2, axis=tuple(range(block.ndim - 1)))
    return int(np.nanargmax(e)) if np.any(np.isfinite(e)) else 0


def lobe_trace(block, lag_ind=None):
    """One temporal trace for a (n_cells, n_lags) filter: the cell time courses
    summed with NON-NEGATIVE weights, the weights being the dominant-sign lobe
    of the spatial profile at `lag_ind` (default: the energy peak), unit norm.

    Why this and not the single strongest cell: the receptive field spans a
    few cells and pooling them raises the trace SNR by about the square root
    of that count; on fly_022 it alone took the cross-family agreement of every
    rule tried from 51-53 to 56-60 of 61. Why not a rank-1 SVD: with 15 cells
    and 600 lags of mostly noise the top singular vector locks onto noise and
    its 'first lobe' lands at random lags. Why non-negative weights: the trace
    then has the SIGN OF THE FILTER at the receptive-field cells at every lag --
    nothing is flipped -- so its first lobe's sign is the cell's polarity.
    Returns (trace, weights); the trace is all-zero if the slice is.
    """
    block = np.asarray(block, dtype=float)
    if lag_ind is None:
        lag_ind = energy_peak_lag(block)
    w = np.nan_to_num(block[:, int(lag_ind)])
    if not np.any(w):
        return np.zeros(block.shape[-1]), np.zeros(block.shape[0])
    sgn = np.sign(w[int(np.argmax(np.abs(w)))])
    lobe = np.clip(w * sgn, 0.0, None)
    lobe = lobe / np.sqrt(np.sum(lobe ** 2))
    return lobe @ np.nan_to_num(block), lobe


def lobes(trace, noise_frac=0.3, smooth=5):
    """Sign-constant runs of the smoothed trace, with an area z for each.

    Returns (runs, smoothed): runs is a list of dicts with start, stop
    (half-open), peak (index of the run's extremum), sign, and z = area / sd of
    the sum of that many noise samples (see _area_sd), so z is comparable
    between a short sharp lobe and a long shallow one.

    ZERO IS THE LATE-LAG MEAN. A reverse-correlation trace carries a constant
    offset -- the response mean times the stimulus's finite-sample mean, the
    same at every lag -- and on fly_022 it exceeded two standard errors in a
    quarter to a half of the kept rois. Left in, it makes the whole tail one
    run of its sign with an enormous area, which then outranks every real lobe
    (the cAMP template's 'first lobe' ran from 0 to 9 s). The late lags are
    signal-free for any filter shorter than the lag axis, the same assumption
    peak_z makes for the noise sd, so their mean is the trace's zero.
    """
    ts = _smooth(trace, smooth)
    n = ts.size
    late_sl = slice(int(round((1.0 - float(noise_frac)) * n)), n)
    ts = ts - np.nanmean(ts[late_sl])
    late = ts[late_sl]
    sd_area = _area_sd(late)
    out = []
    for a, b in _runs_of_sign(ts):
        seg = ts[a:b]
        sd = sd_area(b - a)
        z = float(abs(seg.sum()) / sd) if sd and np.isfinite(sd) and sd > 0 else np.nan
        out.append(dict(start=int(a), stop=int(b), peak=int(a + np.argmax(np.abs(seg))),
                        sign=int(np.sign(seg[0])), z=z))
    return out, ts


def first_lobe(trace, noise_frac=0.3, smooth=5, z_min=2.0):
    """Per-trace rule: the EARLIER of the two most significant lobes.

    Ranking by significance and then taking the earlier one is what makes this
    indifferent to which lobe is bigger -- a rebound twice the first lobe's
    height still comes second -- and robust to a brief blip at lag 0, which is
    never one of the two most significant runs. It is still a single-trace
    onset test and shares the SNR limit of every such test (see the section
    comment); first_lobe_polarity uses it to initialise and as a fallback.
    Returns a dict with resolved, sign, lag (peak index), start, stop, z, and
    amp_ratio = |smoothed trace at lag| / max |smoothed trace| (1.0 when the
    first lobe is also the biggest, < 1 when a later one is).
    """
    runs, ts = lobes(trace, noise_frac, smooth)
    cands = [r for r in runs if np.isfinite(r['z']) and r['z'] >= float(z_min)]
    if not cands:
        return dict(resolved=False, sign=0, lag=None, start=None, stop=None,
                    z=np.nan, amp_ratio=np.nan)
    top = sorted(cands, key=lambda r: -r['z'])[:2]
    first = min(top, key=lambda r: r['start'])
    peak = np.nanmax(np.abs(ts))
    return dict(resolved=True, sign=first['sign'], lag=first['peak'],
                start=first['start'], stop=first['stop'], z=first['z'],
                amp_ratio=float(abs(ts[first['peak']]) / peak) if peak > 0 else np.nan)


def first_lobe_polarity(filt, keep, noise_frac=0.3, smooth=5, z_min=2.0,
                        min_roi=3, n_iter=20):
    """Polarity and first-lobe lag of each kept ROI, by template projection.

    filt : (n_roi, n_cells, n_lags), the VISIBLE cells only (as for peak_signs).
    keep : indices of the kept ROIs -- the family the template is built from.

    Steps. (1) One lobe_trace per kept roi. (2) Initial signs from first_lobe
    on each trace (extremum sign where that is unresolved). (3) Iterate: the
    template is the mean of the traces each multiplied by its sign; each sign
    becomes the sign of its trace's projection onto the template; stop at a
    fixed point. (4) Orient: the template's own first lobe (first_lobe on the
    template; its extremum if even that is unresolved) says which template
    sign is 'positive-first', and every roi's sign is multiplied by it.
    (5) Lag: each roi's own run of its sign that overlaps the template's
    first-lobe run, the largest if several, its extremum; the template's peak
    when the roi has none there. (6) Confidence: |cos| between trace and
    template -- low values mark rois whose waveform does not resemble the
    family's, which is where to look if a sign seems wrong.

    Families with fewer than `min_roi` kept rois get the per-roi rule instead
    (rule_used = 'own_lobe'), extremum sign and energy-peak lag where that is
    unresolved; the result says so per roi in `lag_source`.

    Returns a dict: signs, lag_ind, confidence, amp_ratio (per kept roi, in
    keep order), runs ((n_kept, 2) [start, stop) of the lobe the lag came
    from), traces ((n_kept, n_lags) the smoothed, late-lag-zeroed pooled
    traces, for callers that weight lags by them), extremum_signs
    (peak_signs, for comparison), own_lobe_signs (step 2), lag_source
    ('own_lobe' | 'template' | 'energy_peak'), template (n_lags),
    template_lobe (start, stop, peak) or None, rule_used, n_iter_used,
    n_flipped_vs_extremum, smooth, z_min.
    """
    filt = np.asarray(filt, dtype=float)
    keep = np.asarray(keep, dtype=int)
    n_k, n_lags = len(keep), filt.shape[-1]
    ext = np.asarray(peak_signs(filt, keep, noise_frac), dtype=int) if n_k else np.zeros(0, int)
    out = dict(signs=np.zeros(n_k, int), lag_ind=np.zeros(n_k, int),
               confidence=np.full(n_k, np.nan), amp_ratio=np.full(n_k, np.nan),
               runs=np.zeros((n_k, 2), int), traces=np.zeros((n_k, n_lags)),
               extremum_signs=ext, own_lobe_signs=np.zeros(n_k, int),
               lag_source=np.array(['energy_peak'] * n_k, dtype=object),
               template=np.full(n_lags, np.nan), template_lobe=None,
               rule_used='none', n_iter_used=0, n_flipped_vs_extremum=0,
               smooth=int(smooth), z_min=float(z_min))
    if n_k == 0:
        return out

    traces = np.array([lobe_trace(filt[r])[0] for r in keep])
    own = [first_lobe(t, noise_frac, smooth, z_min) for t in traces]
    own_signs = np.array([o['sign'] for o in own], dtype=int)
    own_lags = np.array([o['lag'] if o['resolved'] else energy_peak_lag(filt[r])
                         for o, r in zip(own, keep)], dtype=int)
    own_ratio = np.array([o['amp_ratio'] for o in own])
    own_runs = np.array([(o['start'], o['stop']) if o['resolved'] else (l, l + 1)
                         for o, l in zip(own, own_lags)], dtype=int).reshape(n_k, 2)
    out['own_lobe_signs'] = own_signs
    out['traces'] = np.array([lobes(t, noise_frac, smooth)[1] for t in traces])
    init = np.where(own_signs != 0, own_signs, ext)
    init = np.where(init != 0, init, 1)

    if n_k < int(min_roi):
        out.update(signs=init.astype(int), lag_ind=own_lags, amp_ratio=own_ratio,
                   runs=own_runs, rule_used='own_lobe',
                   lag_source=np.array(['own_lobe' if o['resolved'] else 'energy_peak'
                                        for o in own], dtype=object))
        out['n_flipped_vs_extremum'] = int(np.sum(out['signs'] != ext))
        return out

    # Smoothed and zeroed on the late lags exactly as lobes() does, so the
    # template, the projections and the per-roi runs below share one baseline.
    T = out['traces']
    s = init.astype(float)
    it = 0
    for it in range(1, int(n_iter) + 1):
        templ = (T * s[:, None]).mean(axis=0)
        s_new = np.sign(T @ templ)
        s_new = np.where(s_new == 0, s, s_new)
        if np.array_equal(s_new, s):
            break
        s = s_new
    templ = (T * s[:, None]).mean(axis=0)
    tl = first_lobe(templ, noise_frac, smooth, z_min)
    if tl['resolved']:
        orient, a0, b0, p0 = tl['sign'], tl['start'], tl['stop'], tl['lag']
        rule = 'template'
    else:
        p0 = int(np.nanargmax(np.abs(templ)))
        orient, a0, b0 = int(np.sign(templ[p0])), p0, p0 + 1
        rule = 'template (extremum-oriented: its first lobe was unresolved)'
    signs = (s * orient).astype(int)
    norms = np.linalg.norm(T, axis=1) * np.linalg.norm(templ)
    with np.errstate(invalid='ignore', divide='ignore'):
        conf = np.abs(T @ templ) / np.where(norms > 0, norms, np.nan)

    lags = np.zeros(n_k, int)
    runs = np.zeros((n_k, 2), int)
    src = np.array(['template'] * n_k, dtype=object)
    ratio = np.full(n_k, np.nan)
    for i in range(n_k):
        tc = T[i]
        cands = [(a, b) for a, b in _runs_of_sign(tc)
                 if np.sign(tc[a]) == signs[i] and a < b0 and b > a0]
        if cands:
            a, b = max(cands, key=lambda ab: np.max(np.abs(tc[ab[0]:ab[1]])))
            lags[i] = a + int(np.argmax(np.abs(tc[a:b])))
            runs[i] = (a, b)
            src[i] = 'own_lobe'
        else:
            lags[i] = p0
            runs[i] = (a0, b0)
        peak = np.nanmax(np.abs(tc))
        ratio[i] = float(abs(tc[lags[i]]) / peak) if peak > 0 else np.nan

    out.update(signs=signs, lag_ind=lags, confidence=conf, amp_ratio=ratio,
               runs=runs, lag_source=src, template=templ,
               template_lobe=(int(a0), int(b0), int(p0)),
               rule_used=rule, n_iter_used=int(it),
               n_flipped_vs_extremum=int(np.sum(signs != ext)))
    return out
