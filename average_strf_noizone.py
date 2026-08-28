"""
average_strf_noizone.py

Register per-ROI STRFs onto a common centred axis and average them.

A SEPARATE, FINAL STEP. It reads the STRFs that analyze_data_strf_noizone.py
already wrote into the hdf5 and touches nothing else, so it can be re-run as
often as you like while tuning thresholds or grids without recomputing filters
or regenerating a single per-ROI figure.

    python average_strf_noizone.py \\
        --experiment_file_directory "path/to/fly_folder" \\
        --tag final --align_channel ch2 --z_threshold "ch1:3.0,ch2:4.5"

WHY REGISTER AT ALL
-------------------
Each ROI sees a different part of the world, so a naive average of many sharp
filters is one broad smear. Shifting each to a common centre first makes the
average an estimate of the SHAPE of a typical filter. The cost is that
retinotopic POSITION is deliberately discarded -- per-ROI centroids are saved
separately so that can be analysed on its own.

DESIGN
------
* Centroid from ONE channel. --align_channel supplies the shift, and that same
  shift is applied to every channel, so cross-channel comparisons survive.
  Aligning each channel independently would destroy them.

* Inclusion needs EVERY channel to pass. --z_threshold takes either one number
  for all channels or a per-channel spec ("ch1:3.0,ch2:4.5"), and an ROI is
  included only if its peak |z| clears the bar in ALL of them.

* Centroid, not peak. The peak is one noisy sample. This is a power-weighted
  (f^2) centre of mass over cells above --centroid_frac of that ROI's own
  maximum, at its own peak lag; restricting to the significant region stops the
  noise floor in the tails dragging the centroid toward the array centre.

* NO polarity flipping. Filter sign is a property of the cell type and the
  channel, so within one experiment it is constant across ROIs and flipping
  per-ROI would be wrong. Each channel is averaged as-is, preserving its own
  sign, which may legitimately differ between channels. The per-channel sign IS
  measured and reported, together with how consistent it is across ROIs -- a
  channel whose ROIs disagree on sign is evidence that something is off
  (mixed cell types, or filters too weak to have a real sign), and that is
  worth seeing rather than averaging away.
* Values are BINNED, never interpolated, and the grid is finer than the bars.
  Each ROI's centroid sits at an arbitrary sub-bar phase, so across ROIs the
  bars sample the filter at many different offsets. Assigning each measured bar
  value to the bin its shifted position falls in preserves that; interpolating
  each ROI onto a shared coarse grid would average it away. Measured on a known
  filter with 5 / 20 / 60 ROIs, this recovers the shape 1.7x / 2.0x / 2.6x more
  accurately -- and unlike the coarse estimate, it keeps improving with N.

* Thin edges are masked. After shifting, ROIs cover different parts of the
  common grid, so points near the edges rest on fewer ROIs.
  --min_roi_per_area blanks those rather than reporting them as if they
  were as well estimated as the centre.
"""

import argparse
import os
import sys
import warnings

import h5py
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

from visanalysis.util.strf_plot import (
    _pole_label, _panel_size, _grid_shape, _tighten,
)

plt.ioff()

# point matplotlib at the ffmpeg bundled with this conda env, in case it is not
# on PATH when this is invoked from the shell wrapper. Duplicated from
# analyze_data_strf_noizone.py deliberately: this script no longer imports that
# module (so it can run without the heavy visanalysis stack), and the movie
# writer is unavailable without it.
_ffmpeg_path = os.path.join(os.path.dirname(sys.executable), 'Library', 'bin', 'ffmpeg.exe')
if os.path.exists(_ffmpeg_path):
    plt.rcParams['animation.ffmpeg_path'] = _ffmpeg_path

INK, MUTED = '#1a1a1a', '#6b6b6b'

# Grid points masked for thin coverage must not be confusable with a zero
# response: RdBu_r's midpoint at a symmetric vmin/vmax is white, which is
# exactly what 'no data' would otherwise look like.
# Fixed hue order for channels: assigned in sorted channel order and never
# cycled, so ch1 is the same colour in every figure. Blue/orange first --
# the standard CVD-safe pair, separated in lightness as well as hue.
CH_COLORS = ['#2f6ba8', '#d2691e', '#4c9a52', '#8a4b9c']
NO_DATA = '#d8d8d8'
CMAP = matplotlib.colormaps['RdBu_r'].with_extremes(bad=NO_DATA)
MOVIE_FPS = 30


# ────────────────────────────────────────────────────────────────────────────
# reading
# ────────────────────────────────────────────────────────────────────────────

def load_strf_groups(path, tag):
    """Every 1D and 2D STRF in one fly's hdf5, keyed for pairing across channels."""
    oned, twod = {}, {}
    with h5py.File(path, 'r') as f:
        root = f.get(f'STRF_{tag}')
        if root is None:
            raise KeyError(f'{path} has no /STRF_{tag} group -- run '
                           f'analyze_data_strf_noizone.py first.')
        for sn in root:
            if sn in ('combined_2d', 'averaged'):
                continue
            for ch in root[sn]:
                for rtz_key in root[sn][ch]:
                    for nt_key in root[sn][ch][rtz_key]:
                        g = root[sn][ch][rtz_key][nt_key]
                        oned[(sn, ch, float(g.attrs['rtz']),
                              float(g.attrs['noisetau']))] = {
                            'strf': g['strf'][()],
                            't_lag': g['t_lag'][()],
                            'bar_pos_deg': g['bar_pos_deg'][()],
                        }
        c2 = root.get('combined_2d')
        if c2 is not None:
            for pair in c2:
                for ch in c2[pair]:
                    for nt_key in c2[pair][ch]:
                        for rtz_key in c2[pair][ch][nt_key]:
                            g = c2[pair][ch][nt_key][rtz_key]
                            twod[(pair, ch, nt_key, rtz_key)] = {
                                'strf_2d': g['strf_2d'][()],
                                't_lag': g['t_lag'][()],
                                'bar_pos_h_deg': g['bar_pos_h_deg'][()],
                                'bar_pos_v_deg': g['bar_pos_v_deg'][()],
                                'rtz_h': float(g.attrs['rtz_h']),
                                'rtz_v': float(g.attrs['rtz_v']),
                            }
    return oned, twod


def parse_thresholds(spec, channels, default=4.0):
    """'4.0' -> every channel; 'ch1:3.0,ch2:4.5' -> per channel."""
    if spec is None:
        return {c: default for c in channels}
    spec = str(spec).strip()
    if ':' not in spec and '=' not in spec:
        return {c: float(spec) for c in channels}
    out = {}
    for part in spec.replace(';', ',').split(','):
        if not part.strip():
            continue
        k, v = part.replace('=', ':').split(':')
        out[k.strip()] = float(v)
    missing = [c for c in channels if c not in out]
    if missing:
        raise ValueError(f'--z_threshold gives no value for {missing}; '
                         f'channels present are {channels}')
    return out


# ────────────────────────────────────────────────────────────────────────────
# per-ROI summary
# ────────────────────────────────────────────────────────────────────────────

def _peak_lag(f):
    """Lag index of the largest |filter| for one ROI. f is (..., n_lags)."""
    return int(np.argmax(np.abs(f).reshape(-1, f.shape[-1]).max(axis=0)))


def _weighted_centroid(vals, coords, frac):
    """Power-weighted centre of mass over cells above frac * max|vals|."""
    a = np.abs(vals)
    peak = a.max()
    if not np.isfinite(peak) or peak <= 0:
        return tuple(np.nan for _ in coords), 0.0
    w = np.where(a >= frac * peak, vals ** 2, 0.0)
    tot = w.sum()
    if tot <= 0:
        return tuple(np.nan for _ in coords), 0.0
    return tuple(float((w * c).sum() / tot) for c in coords), float(peak)


def peak_z_and_sign(strf):
    """Per-ROI peak |z| and the sign at that peak. strf is (n_roi, ..., n_lags).

    Thresholding this (see --z_threshold, applied in main) is a de facto
    significance test, so it is worth recording what its null actually is.

    The lag axis is oversampled. compute_strf_1d builds it at the raw display
    frame period, but a noizone image is held for noisetau ms, so the number of
    INDEPENDENT lags is filter_length_s / noisetau, not n_lags -- typically a
    factor of ~5-10 fewer. That does not bias z itself: the marginal variance at
    each lag is unaffected by whether neighbouring lags are correlated. What it
    changes is the multiple-comparisons burden behind `peak`, since the expected
    maximum over M correlated entries is smaller than over M independent ones.

    So the oversampling makes a fixed threshold STRICTER than the array size
    suggests, which is the safe direction -- decimating the lag axis to noisetau
    resolution would make the same nominal threshold more permissive, not less.
    Within one noisetau nothing here needs fixing.

    The one case that IS confounded is comparing ACROSS noisetau. n_lags is the
    same for every group (filter_length_s / raw display period, independent of
    noisetau) but the effective entry count is not, so a longer hold has fewer
    independent entries at identical array size and therefore a lower null
    maximum. The same nominal threshold is then a different test in each group,
    and the long-hold condition will look more selective purely from smoothing.
    If noisetau ever becomes something we compare across rather than just a
    grouping key, calibrate per-group against a shuffled null or normalise by
    effective entry count.

    Separately, note the denominator: compute_strf_1d divides by the std of the
    whole (n_positions, n_lags) array, which for a strongly responding ROI
    includes its own signal. That deflates z for the best cells -- conservative
    again, but it means z is not comparable to a noise-only z-score.
    """
    n_roi = strf.shape[0]
    pk = np.zeros(n_roi)
    sgn = np.zeros(n_roi)
    for r in range(n_roi):
        f = strf[r]
        s = f[..., _peak_lag(f)].ravel()
        i = int(np.argmax(np.abs(s)))
        pk[r] = abs(s[i])
        sgn[r] = np.sign(s[i]) or 1.0
    return pk, sgn


def centroids_1d(strf, bar_pos_deg, centroid_frac):
    cen = np.full(strf.shape[0], np.nan)
    for r in range(strf.shape[0]):
        f = strf[r, :, _peak_lag(strf[r])]
        (c,), _ = _weighted_centroid(f, (bar_pos_deg,), centroid_frac)
        cen[r] = c
    return cen


def centroids_2d(strf_2d, pos_h, pos_v, centroid_frac):
    cen = np.full((strf_2d.shape[0], 2), np.nan)
    H = pos_h[:, None] * np.ones((1, len(pos_v)))
    V = np.ones((len(pos_h), 1)) * pos_v[None, :]
    for r in range(strf_2d.shape[0]):
        f = strf_2d[r, :, :, _peak_lag(strf_2d[r])]
        (ch_, cv), _ = _weighted_centroid(f, (H, V), centroid_frac)
        cen[r] = (ch_, cv)
    return cen


def summarise_polarity(signs, include):
    """(modal sign, fraction agreeing) among included ROIs.

    Sign is set by cell type and channel, so within one experiment it should be
    the same for every ROI in a channel. A low agreeing-fraction means it isn't,
    which is worth knowing.
    """
    s = signs[include]
    if s.size == 0:
        return 0.0, 0.0
    modal = 1.0 if (s > 0).sum() >= (s < 0).sum() else -1.0
    return modal, float((s == modal).mean())


# ────────────────────────────────────────────────────────────────────────────
# registration + averaging
# ────────────────────────────────────────────────────────────────────────────

def auto_grid_step(pos, n_roi, min_roi_per_area, cap=4):
    """Grid step fine enough to exploit centroid jitter, coarse enough to fill.

    Two limits, whichever binds first.

    HOW MANY ROIs YOU HAVE. With a step of bar/k each ROI lands in one bin out
    of every k, so a bin collects roughly n_roi/k ROIs. Asking for more than
    n_roi/min_roi_per_area subdivisions leaves most bins under-populated and
    masked, so k stops there -- with a factor of 2 margin, because the sub-bar
    phases are random rather than evenly spread.

    WHAT THE STIMULUS CAN RESOLVE. Every sample is the filter integrated over
    one bar, so the measurement is already low-passed by a box of that width,
    whose first spectral null sits at 1/bar. A step of bar/2 therefore already
    matches Nyquist for the box-limited signal and bar/4 gives 2x margin.
    Measured on a known filter with 60 ROIs, the shape error falls
    0.0108 -> 0.0070 -> 0.0063 for bar/1 -> bar/2 -> bar/4, then flattens and
    slowly REVERSES (0.0065 at bar/12, 0.0066 at bar/20) as bins get too thin
    to average well, while coverage drops 91% -> 38%. Hence cap=4: beyond it
    costs coverage and buys nothing.

    Returns (step_deg, bar_spacing_deg, k).
    """
    bar = float(np.median(np.diff(np.sort(pos))))
    k = int(np.clip(n_roi // max(2 * min_roi_per_area, 1), 1, cap))
    return bar / k, bar, k


def _accumulate(shape_tail, idx_list, values, n_bins):
    """Per-ROI bin means, then across-ROI sum / sum-of-squares / ROI count.

    Averaging within an ROI first means a bin that happens to catch two of the
    same ROI's bars counts as ONE roi, not two -- so min_roi_per_area really
    counts ROIs.
    """
    total = np.zeros((n_bins,) + shape_tail)
    total_sq = np.zeros((n_bins,) + shape_tail)
    n_contrib = np.zeros(n_bins, dtype=int)
    for idx, vals in zip(idx_list, values):
        rs = np.zeros((n_bins,) + shape_tail)
        rn = np.zeros(n_bins)
        np.add.at(rs, idx, vals)
        np.add.at(rn, idx, 1.0)
        hit = rn > 0
        rm = rs[hit] / rn[hit].reshape((-1,) + (1,) * len(shape_tail))
        total[hit] += rm
        total_sq[hit] += rm ** 2
        n_contrib[hit] += 1
    return total, total_sq, n_contrib


def _finish(total, total_sq, n_contrib, min_roi_per_area, extra):
    n = n_contrib.reshape((-1,) + (1,) * (total.ndim - 1))
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(n > 0, total / np.maximum(n, 1), np.nan)
        var = np.where(n > 1, total_sq / np.maximum(n, 1) - mean ** 2, np.nan)
        sem = np.sqrt(np.maximum(var, 0.0)) / np.sqrt(np.maximum(n, 1))
    thin = (n_contrib < min_roi_per_area)
    mean[thin] = np.nan
    sem[thin] = np.nan
    out = {'mean': mean, 'sem': sem, 'count': n_contrib,
           'min_roi_per_area': min_roi_per_area}
    out.update(extra)
    return out


def register_average_1d(strf, bar_pos_deg, centroids, include,
                        min_roi_per_area=3, grid_step_deg=None):
    """Bin each included ROI's samples by shifted position and average.

    Values are never interpolated -- each measured bar value is assigned to the
    grid bin its shifted position falls in. Because ROI centroids sit at
    arbitrary sub-bar phases, different ROIs fall in different bins, and the
    population resolves the filter more finely than the bar spacing.
    """
    sel = np.flatnonzero(include)
    step, bar, k = auto_grid_step(bar_pos_deg, len(sel), min_roi_per_area)
    if grid_step_deg:
        step = float(grid_step_deg)
    lo = float(np.nanmin(bar_pos_deg[0] - centroids[sel]))
    hi = float(np.nanmax(bar_pos_deg[-1] - centroids[sel]))
    n_bins = max(3, int(np.floor((hi - lo) / step)) + 1)
    grid = lo + step * np.arange(n_bins)

    idx_list, values = [], []
    for r in sel:
        i = np.rint(((bar_pos_deg - centroids[r]) - lo) / step).astype(int)
        ok = (i >= 0) & (i < n_bins)
        idx_list.append(i[ok])
        values.append(strf[r, ok, :])
    tot, tsq, cnt = _accumulate((strf.shape[2],), idx_list, values, n_bins)
    return _finish(tot, tsq, cnt, min_roi_per_area,
                   {'grid_deg': grid, 'n_roi': len(sel),
                    'grid_step_deg': step, 'bar_spacing_deg': bar,
                    'subdivisions': int(round(bar / step))})


def register_average_2d(strf_2d, pos_h, pos_v, centroids, include,
                        min_roi_per_area=3, grid_step_deg=None):
    """2D analogue -- bin in both axes, still no interpolation of values."""
    sel = np.flatnonzero(include)
    step_h, bar_h, _ = auto_grid_step(pos_h, len(sel), min_roi_per_area)
    step_v, bar_v, _ = auto_grid_step(pos_v, len(sel), min_roi_per_area)
    if grid_step_deg:
        step_h = step_v = float(grid_step_deg)
    lo_h = float(np.nanmin(pos_h[0] - centroids[sel, 0]))
    hi_h = float(np.nanmax(pos_h[-1] - centroids[sel, 0]))
    lo_v = float(np.nanmin(pos_v[0] - centroids[sel, 1]))
    hi_v = float(np.nanmax(pos_v[-1] - centroids[sel, 1]))
    nh = max(3, int(np.floor((hi_h - lo_h) / step_h)) + 1)
    nv = max(3, int(np.floor((hi_v - lo_v) / step_v)) + 1)
    gh = lo_h + step_h * np.arange(nh)
    gv = lo_v + step_v * np.arange(nv)

    idx_list, values = [], []
    for r in sel:
        ih = np.rint(((pos_h - centroids[r, 0]) - lo_h) / step_h).astype(int)
        iv = np.rint(((pos_v - centroids[r, 1]) - lo_v) / step_v).astype(int)
        okh, okv = (ih >= 0) & (ih < nh), (iv >= 0) & (iv < nv)
        IH, IV = np.meshgrid(ih[okh], iv[okv], indexing='ij')
        idx_list.append((IH * nv + IV).ravel())          # flattened bin index
        values.append(strf_2d[r][np.ix_(okh, okv)].reshape(-1, strf_2d.shape[3]))
    tot, tsq, cnt = _accumulate((strf_2d.shape[3],), idx_list, values, nh * nv)
    res = _finish(tot, tsq, cnt, min_roi_per_area,
                  {'grid_h_deg': gh, 'grid_v_deg': gv, 'n_roi': len(sel),
                   'grid_step_deg': (step_h, step_v),
                   'bar_spacing_deg': (bar_h, bar_v),
                   'subdivisions': (int(round(bar_h / step_h)),
                                    int(round(bar_v / step_v)))})
    res['mean'] = res['mean'].reshape(nh, nv, -1)
    res['sem'] = res['sem'].reshape(nh, nv, -1)
    res['count'] = res['count'].reshape(nh, nv)
    return res


def _smooth_box(a, width):
    """Box-average a 1D profile, edges handled by edge-padding."""
    w = max(1, int(width))
    if w == 1:
        return np.asarray(a, dtype=float)
    pad = w // 2
    return np.convolve(np.pad(np.asarray(a, dtype=float), pad, mode='edge'),
                       np.ones(w) / w, mode='same')[pad:pad + len(a)]


def _contiguous_run(profile, thresh):
    """Slice of the contiguous run through the best point of `profile`.

    A plain bounding box is not enough: with random sub-bar phases an isolated
    bin far out can clear the threshold alone and drag the axis limit to it,
    leaving the filter in a quarter of the figure. Growing outward from the
    peak stops at the first genuine gap instead.

    `profile` must already be smoothed -- see coverage_bounds.
    """
    ok = profile >= thresh
    if not ok.any():
        return slice(int(np.argmax(profile)), int(np.argmax(profile)) + 1)
    c = int(np.argmax(profile))
    lo = c
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    hi = c
    while hi < len(profile) - 1 and ok[hi + 1]:
        hi += 1
    return slice(lo, hi + 1)


def coverage_bounds(count, min_roi_per_area, frac_of_max=0.5, smooth=1):
    """Index bounds of the well-supported region, for DISPLAY only.

    Two different jobs, deliberately separated:

      * min_roi_per_area decides what is REPORTED. Points below it are NaN in
        the saved arrays -- a validity question.
      * this decides what is SHOWN. The common grid is the union of every
        shifted ROI's reach, so its edges rest on whichever one or two ROIs
        happen to extend furthest.

    The coverage profile is box-smoothed over one bar first, because raw
    bin-level coverage is spiky: with the grid subdivided k-fold each bin
    catches a different random subset of ROIs, and growing a run through that
    stops at the first unlucky bin. The threshold is then taken as
    frac_of_max of the SMOOTHED profile's own maximum -- taking it from the raw
    maximum instead made high fractions unreachable, so the run search fell
    through to "show everything" and the crop widened as the fraction rose.

    frac_of_max = 0 crops only at min_roi_per_area (widest); 0.5 is the
    default; approaching 1 keeps just the best-covered core.

    Returns a tuple of slices, one per spatial axis of `count`.
    """
    c = np.asarray(count)
    if c.max() <= 0:
        return tuple(slice(None) for _ in range(c.ndim))
    sm = smooth if isinstance(smooth, (tuple, list)) else (smooth,) * c.ndim
    out = []
    for ax in range(c.ndim):
        other = tuple(a for a in range(c.ndim) if a != ax)
        prof = _smooth_box(c.max(axis=other) if other else c, sm[ax])
        thresh = max(min_roi_per_area, frac_of_max * prof.max())
        out.append(_contiguous_run(prof, thresh))
    return tuple(out)


# ────────────────────────────────────────────────────────────────────────────
# figures
# ────────────────────────────────────────────────────────────────────────────

def _style(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=2, colors=MUTED)


def figure_1d(res, key, out_dir, tag, fmt='.pdf', crop_frac=0.5):
    sn, ch, rtz, nt = key
    (ks,) = coverage_bounds(res['count'], res.get('min_roi_per_area', 1),
                                frac_of_max=crop_frac,
                                smooth=res.get('subdivisions', 1))
    grid = res['grid_deg'][ks]
    mean, sem = res['mean'][ks], res['sem'][ks]
    fh, axes = plt.subplots(1, 2, figsize=(9, 3.6), constrained_layout=True,
                            gridspec_kw={'width_ratios': [2, 1]})
    vmax = np.nanmax(np.abs(mean)) or 1.0
    im = axes[0].imshow(mean, aspect='auto', origin='lower', cmap=CMAP,
                        vmin=-vmax, vmax=vmax,
                        extent=[res['t_lag'][0] * 1000, res['t_lag'][-1] * 1000,
                                grid[0], grid[-1]])
    axes[0].set_xlabel('Lag (ms)')
    axes[0].set_ylabel(f'{_pole_label(rtz)}\nrelative to RF centre (deg)')
    axes[0].axhline(0, color='k', lw=0.6, ls=':')
    plt.colorbar(im, ax=axes[0], label='mean z-score')
    _style(axes[0])

    li = int(np.nanargmax(np.nanmax(np.abs(mean), axis=0)))
    m, e = mean[:, li], sem[:, li]
    axes[1].fill_betweenx(grid, m - e, m + e, color='#2f6ba8', alpha=0.25, lw=0)
    axes[1].plot(m, grid, color='#2f6ba8', lw=2)
    axes[1].axvline(0, color='#c0c0c0', lw=1)
    axes[1].axhline(0, color='k', lw=0.6, ls=':')
    axes[1].set_xlabel('mean z-score')
    axes[1].set_title(f'peak lag {res["t_lag"][li]*1000:.0f} ms', fontsize=9)
    _style(axes[1])

    pol = res.get('polarity_sign', 0)
    fh.suptitle(f'Registered mean 1D STRF   {sn} {ch}  rtz={rtz:.0f}  '
                f'tau={nt:.0f}ms   n={res["n_roi"]} ROIs   '
                f'channel sign {pol:+.0f} ({100*res.get("polarity_frac", 0):.0f}% agree)   '
                f'{tag}', fontsize=10)
    name = f'strf_1d_MEAN_rtz{int(rtz)}_noisetau{int(nt)}_{ch}_{sn}_{tag}_'
    fh.savefig(os.path.join(out_dir, name + fmt), dpi=200, transparent=True,
               bbox_inches='tight')
    plt.close(fh)


def figure_1d_compare(res_by_ch, key, out_dir, tag, fmt='.pdf', crop_frac=0.5):
    """All channels' mean 1D filters in ONE figure, with marginals on both axes.

    Left column: the (position x lag) map for each channel, on a shared pair of
    axes so features line up row to row.

    Right column: the two marginals, every channel overlaid.
        top    -- SPATIAL, the profile across position at that channel's own
                  peak lag: the receptive field's shape.
        bottom -- TEMPORAL, the profile across lag at that channel's own peak
                  position: the filter's kinetics.

    Each channel is taken at ITS OWN peak, not a shared one, because a genuine
    difference in kinetics between indicators would otherwise be read as a
    difference in shape.

    Marginals are divided by their own peak MAGNITUDE, so shapes overlay
    directly while the sign is preserved -- an OFF channel stays a
    negative-going curve rather than being flipped to match an ON one. The
    legend carries each channel's true peak z, peak lag and peak position, so
    normalising for shape does not hide the scale.
    """
    sn, rtz, nt = key
    chans = sorted(res_by_ch)
    n_ch = len(chans)

    fh = plt.figure(figsize=(11, 2.4 * max(n_ch, 2) + 0.6), constrained_layout=True)
    gs = fh.add_gridspec(max(n_ch, 2), 2, width_ratios=[1.7, 1])
    ax_img = [fh.add_subplot(gs[i, 0]) for i in range(n_ch)]
    ax_sp = fh.add_subplot(gs[0, 1])
    ax_tm = fh.add_subplot(gs[1, 1])

    legend_bits = []
    for i, ch in enumerate(chans):
        res = res_by_ch[ch]
        colour = CH_COLORS[i % len(CH_COLORS)]
        (ks,) = coverage_bounds(res['count'], res.get('min_roi_per_area', 1),
                                frac_of_max=crop_frac,
                                smooth=res.get('subdivisions', 1))
        grid = res['grid_deg'][ks]
        mean, sem = res['mean'][ks], res['sem'][ks]
        t_ms = res['t_lag'] * 1000

        vmax = np.nanmax(np.abs(mean)) or 1.0
        ax = ax_img[i]
        im = ax.imshow(mean, aspect='auto', origin='lower', cmap=CMAP,
                       vmin=-vmax, vmax=vmax,
                       extent=[t_ms[0], t_ms[-1], grid[0], grid[-1]])
        ax.axhline(0, color='k', lw=0.6, ls=':')
        ax.set_ylabel(f'{ch}\n{_pole_label(rtz)}\nrel. RF centre (deg)', fontsize=8)
        if i == n_ch - 1:
            ax.set_xlabel('Lag (ms)')
        else:
            ax.set_xticklabels([])
        plt.colorbar(im, ax=ax, label='mean z', pad=0.01)
        _style(ax)

        # each channel at its OWN peak, in both directions
        li = int(np.nanargmax(np.nanmax(np.abs(mean), axis=0)))
        pi = int(np.nanargmax(np.abs(mean[:, li])))
        scale = np.abs(mean[pi, li]) or 1.0

        sp, sp_e = mean[:, li] / scale, sem[:, li] / scale
        ax_sp.fill_betweenx(grid, sp - sp_e, sp + sp_e, color=colour, alpha=0.2, lw=0)
        ax_sp.plot(sp, grid, color=colour, lw=2, label=ch)

        tm, tm_e = mean[pi, :] / scale, sem[pi, :] / scale
        ax_tm.fill_between(t_ms, tm - tm_e, tm + tm_e, color=colour, alpha=0.2, lw=0)
        ax_tm.plot(t_ms, tm, color=colour, lw=2, label=ch)

        legend_bits.append(f'{ch}: peak z {mean[pi, li]:+.2f} at '
                           f'{t_ms[li]:.0f} ms, {grid[pi]:+.1f} deg')

    ax_sp.axvline(0, color='#c0c0c0', lw=1)
    ax_sp.axhline(0, color='k', lw=0.6, ls=':')
    ax_sp.set_xlabel('normalised amplitude')
    ax_sp.set_ylabel('rel. RF centre (deg)', fontsize=8)
    ax_sp.set_title('spatial marginal (at each channel\'s peak lag)', fontsize=8)
    ax_sp.legend(fontsize=8, frameon=False)

    ax_tm.axhline(0, color='#c0c0c0', lw=1)
    ax_tm.set_xlabel('Lag (ms)')
    ax_tm.set_ylabel('normalised amplitude', fontsize=8)
    ax_tm.set_title('temporal marginal (at each channel\'s peak position)', fontsize=8)
    ax_tm.legend(fontsize=8, frameon=False)
    for a in (ax_sp, ax_tm):
        _style(a)
        a.grid(True, color='#f0f0f0', lw=0.6)
        a.set_axisbelow(True)

    n_roi = res_by_ch[chans[0]]['n_roi']
    fh.suptitle(f'Registered mean 1D STRF, channels compared   {sn}  rtz={rtz:.0f}  '
                f'tau={nt:.0f}ms   n={n_roi} ROIs   {tag}\n'
                + '     '.join(legend_bits), fontsize=9)
    name = f'strf_1d_MEAN_compare_rtz{int(rtz)}_noisetau{int(nt)}_{sn}_{tag}_'
    fh.savefig(os.path.join(out_dir, name + fmt), dpi=200, transparent=True,
               bbox_inches='tight')
    plt.close(fh)


def figure_2d(res, key, out_dir, tag, n_display_lags=24, fmt='.pdf', crop_frac=0.5):
    pair, ch, nt_key, rtz_key = key
    kh, kv = coverage_bounds(res['count'], res.get('min_roi_per_area', 1),
                                frac_of_max=crop_frac,
                                smooth=res.get('subdivisions', 1))
    gh, gv = res['grid_h_deg'][kh], res['grid_v_deg'][kv]
    extent = [gv[0], gv[-1], gh[0], gh[-1]]
    mean, t_lag = res['mean'][kh][:, kv], res['t_lag']
    vmax = np.nanmax(np.abs(mean)) or 1.0

    fh, ax = plt.subplots(figsize=_panel_size(extent, 4.2), constrained_layout=True)
    li = int(np.nanargmax(np.nanmax(np.abs(mean), axis=(0, 1))))
    ax.imshow(mean[:, :, li], aspect='equal', origin='lower', extent=extent,
              cmap=CMAP, vmin=-vmax, vmax=vmax)
    ax.axhline(0, color='k', lw=0.6, ls=':')
    ax.axvline(0, color='k', lw=0.6, ls=':')
    ax.set_xlabel(f'{_pole_label(res["rtz_v"])}\nrel. RF centre (deg)')
    ax.set_ylabel(f'{_pole_label(res["rtz_h"])}\nrel. RF centre (deg)')
    _style(ax)
    fh.suptitle(f'Registered mean 2D STRF, peak lag {t_lag[li]*1000:.0f} ms\n'
                f'{pair} {ch} {nt_key} {rtz_key}  n={res["n_roi"]}  {tag}', fontsize=9)
    name = f'strf_2d_MEAN_peak_lag_{rtz_key}_{nt_key}_{pair}_{ch}_{tag}_'
    fh.savefig(os.path.join(out_dir, name + fmt), dpi=200, transparent=True,
               bbox_inches='tight')
    plt.close(fh)

    idxs = np.unique(np.linspace(0, mean.shape[2] - 1,
                                 min(n_display_lags, mean.shape[2])).round().astype(int))
    n_rows, n_cols = _grid_shape(len(idxs), extent)
    pw, ph = _panel_size(extent, 1.7)
    fh, axes = plt.subplots(n_rows, n_cols, figsize=(pw * n_cols, ph * n_rows + 0.35),
                            constrained_layout=True)
    _tighten(fh)
    A = np.array(axes).flatten()
    for i, li in enumerate(idxs):
        ax = A[i]
        ax.imshow(mean[:, :, li], aspect='equal', origin='lower', extent=extent,
                  cmap=CMAP, vmin=-vmax, vmax=vmax)
        ax.text(0.03, 0.97, f'{t_lag[li]*1000:.0f} ms', transform=ax.transAxes,
                ha='left', va='top', fontsize=7,
                bbox=dict(fc='white', ec='none', alpha=0.75, pad=1))
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
    for ax in A[len(idxs):]:
        ax.set_visible(False)
    fh.suptitle(f'Registered mean 2D STRF  {pair} {ch} {rtz_key}  '
                f'n={res["n_roi"]}  {tag}', fontsize=10)
    name = f'strf_2d_MEAN_lags_{rtz_key}_{nt_key}_{pair}_{ch}_{tag}_'
    fh.savefig(os.path.join(out_dir, name + fmt), dpi=150, transparent=True,
               bbox_inches='tight')
    plt.close(fh)


def movie_2d(res, key, out_dir, tag, fps=MOVIE_FPS, crop_frac=0.5):
    """Averaged 2D STRF across every lag, one video per channel.

    Mirrors the per-ROI movie in analyze_data_strf_noizone.py: one colour scale
    for the whole sequence, so a weak lag is not rendered as if it were as
    strong as the peak.
    """
    pair, ch, nt_key, rtz_key = key
    kh, kv = coverage_bounds(res['count'], res.get('min_roi_per_area', 1),
                                frac_of_max=crop_frac,
                                smooth=res.get('subdivisions', 1))
    gh, gv = res['grid_h_deg'][kh], res['grid_v_deg'][kv]
    extent = [gv[0], gv[-1], gh[0], gh[-1]]
    mean, t_lag = res['mean'][kh][:, kv], res['t_lag']
    n_lags = mean.shape[2]
    vmax = np.nanmax(np.abs(mean)) or 1.0

    fh, ax = plt.subplots(figsize=_panel_size(extent, 4.2), constrained_layout=True)
    im = ax.imshow(mean[:, :, 0], aspect='equal', origin='lower', extent=extent,
                   cmap=CMAP, vmin=-vmax, vmax=vmax)
    ax.axhline(0, color='k', lw=0.6, ls=':')
    ax.axvline(0, color='k', lw=0.6, ls=':')
    ax.set_xlabel(f'{_pole_label(res["rtz_v"])}\nrel. RF centre (deg)')
    ax.set_ylabel(f'{_pole_label(res["rtz_h"])}\nrel. RF centre (deg)')
    plt.colorbar(im, ax=ax, label='mean z-score')
    _style(ax)
    title = ax.set_title('')

    def _update(li, _im=im, _title=title):
        _im.set_data(mean[:, :, li])
        _title.set_text(f'{pair} {ch} {rtz_key}  n={res["n_roi"]}  '
                        f'lag = {t_lag[li]*1000:.1f} ms')
        return _im, _title

    ani = animation.FuncAnimation(fh, _update, frames=n_lags, blit=False)
    name = f'strf_2d_MEAN_movie_{rtz_key}_{nt_key}_{pair}_{ch}_{tag}_.mp4'
    ani.save(os.path.join(out_dir, name), writer='ffmpeg', fps=fps, dpi=150)
    plt.close(fh)


def figure_inclusion(summary, thresholds, out_dir, tag, fmt='.pdf'):
    """Which ROIs made the cut, and where their receptive fields sit.

    A receptive field centre is a 2D point, but each 1D STRF constrains only
    ONE axis: bars at rtz=0 say nothing about dorsal-ventral position, and
    bars at rtz=90 say nothing about anterior-posterior. So the centroid this
    script computes per orientation is a scalar -- the projection onto that
    orientation's psi axis. The 2D position appears only when the two
    orientations are paired, which is what the middle panel does.

    Panels:
      left   -- peak |z| per channel, sorted, with each channel's threshold.
      middle -- RF centres in 2D where both orientations exist, from pairing
                the two 1D centroids by ROI. Falls back to a 1D strip if only
                one orientation was run.
      right  -- centroid distribution per orientation against the extent of
                the field the bars actually covered. A distribution pressed up
                against an edge means those RFs run off the screen, so their
                filters are truncated and their centroids biased inward.
    """
    fh, axes = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)

    # ---- left: the inclusion test ----------------------------------------
    for i, ch in enumerate(sorted(thresholds)):
        c = CH_COLORS[i % len(CH_COLORS)]
        first = True
        for s_ in summary.values():
            if ch not in s_['peak_z']:
                continue
            axes[0].plot(np.sort(s_['peak_z'][ch])[::-1], lw=1.5, color=c,
                         label=ch if first else None)
            first = False
        axes[0].axhline(thresholds[ch], color=c, lw=1.1, ls='--')
    axes[0].set_xlabel('ROI (sorted)')
    axes[0].set_ylabel('peak |z|')
    axes[0].set_title('inclusion (dashed = per-channel threshold)', fontsize=9)
    axes[0].legend(fontsize=8, frameon=False)

    # ---- middle: RF centres, 2D when both orientations are available ------
    by_rtz = {}
    for (sn, ch, rtz, nt), s_ in summary.items():
        by_rtz.setdefault(round(rtz % 180, 3), []).append(s_)
    angles = sorted(by_rtz)
    pair = None
    for a in angles:
        b = round((a + 90) % 180, 3)
        if b in by_rtz and b != a:
            pair = (a, b)
            break

    if pair is not None:
        sa, sb = by_rtz[pair[0]][0], by_rtz[pair[1]][0]
        n = min(len(sa['centroid']), len(sb['centroid']))
        xa, xb = sa['centroid'][:n], sb['centroid'][:n]
        both = sa['include'][:n] & sb['include'][:n]
        axes[1].scatter(xa[~both], xb[~both], s=18, facecolor='none',
                        edgecolor='#b0b0b0', lw=0.8, label='excluded')
        axes[1].scatter(xa[both], xb[both], s=26, color=CH_COLORS[0],
                        alpha=0.8, label=f'included (n={int(both.sum())})')
        axes[1].set_xlabel(f'{_pole_label(pair[0])} (deg)', fontsize=8)
        axes[1].set_ylabel(f'{_pole_label(pair[1])} (deg)', fontsize=8)
        axes[1].set_title('RF centres, both orientations paired', fontsize=9)
        axes[1].axhline(0, color='#e0e0e0', lw=1, zorder=0)
        axes[1].axvline(0, color='#e0e0e0', lw=1, zorder=0)
        axes[1].legend(fontsize=7, frameon=False)
    else:
        for (sn, ch, rtz, nt), s_ in summary.items():
            axes[1].plot(s_['centroid'][s_['include']],
                         s_['peak_z'][s_['align_channel']][s_['include']],
                         'o', ms=4, alpha=0.7, label=f'rtz={rtz:.0f}')
        axes[1].set_xlabel('RF centroid (deg)')
        axes[1].set_ylabel('peak |z| (align channel)')
        axes[1].set_title('only one orientation -- no 2D centre available', fontsize=9)
        axes[1].legend(fontsize=7, frameon=False)

    # ---- right: centroids vs the field that was actually covered ---------
    for i, ((sn, ch, rtz, nt), s_) in enumerate(sorted(summary.items())):
        c = CH_COLORS[i % len(CH_COLORS)]
        vals = s_['centroid'][s_['include']]
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue
        axes[2].hist(vals, bins=18, histtype='step', lw=2, color=c,
                     label=f'rtz={rtz:.0f}  (n={vals.size})')
        fr = s_.get('field_range')
        if fr is not None:
            for edge in fr:
                axes[2].axvline(edge, color=c, lw=1.1, ls='--', alpha=0.7)
    axes[2].set_xlabel('RF centroid (deg)   dashed = edge of the covered field')
    axes[2].set_ylabel('ROIs')
    axes[2].set_title('centroids vs field extent'
                      '  (piled at an edge -> truncated RFs)', fontsize=9)
    axes[2].legend(fontsize=7, frameon=False)

    for ax in axes:
        _style(ax)
        ax.grid(True, color='#f0f0f0', lw=0.6)
        ax.set_axisbelow(True)
    fh.savefig(os.path.join(out_dir, f'strf_MEAN_inclusion_{tag}_' + fmt),
               dpi=200, transparent=True, bbox_inches='tight')
    plt.close(fh)


# ────────────────────────────────────────────────────────────────────────────
# main
# ────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--experiment_file_directory', nargs='?', required=True)
    ap.add_argument('--tag', nargs='?', default='final', choices=['raw', 'final'])
    ap.add_argument('--align_channel', nargs='?', default=None,
                    help='channel supplying the centroid; its shift is applied '
                         'to every channel. Default: first channel present.')
    ap.add_argument('--z_threshold', nargs='?', default=None,
                    help='one number for all channels, or per channel as '
                         '"ch1:3.0,ch2:4.5". EVERY channel must pass.')
    ap.add_argument('--centroid_frac', nargs='?', type=float, default=0.5)
    ap.add_argument('--min_roi_per_area', nargs='?', type=int, default=3,
                    help='grid points backed by fewer ROIs than this are masked')
    ap.add_argument('--grid_step_deg', nargs='?', type=float, default=None,
                    help='output grid spacing. Default: chosen from the ROI '
                         'count so bins stay populated -- see auto_grid_step. '
                         'Values are binned, never interpolated.')
    ap.add_argument('--crop_frac_of_max', nargs='?', type=float, default=0.5,
                    help='DISPLAY ONLY. Figures are cropped to the contiguous '
                         'region whose ROI coverage reaches this fraction of '
                         'the best coverage achieved. 0 = crop only at '
                         'min_roi_per_area (widest), 0.5 = default, 0.9 = the '
                         'very core. Does not change the saved arrays.')
    ap.add_argument('--save_figs', nargs='?', default='True')
    args = ap.parse_args()

    save_figs = args.save_figs == 'True'
    fname = 'fly.hdf5' if args.tag == 'raw' else 'fly_final.hdf5'
    path = os.path.join(args.experiment_file_directory, fname)
    figs_dir = os.path.join(args.experiment_file_directory, f'{args.tag}_roi_figs')
    os.makedirs(figs_dir, exist_ok=True)

    oned, twod = load_strf_groups(path, args.tag)
    if not oned:
        print('no 1D STRFs found -- nothing to average')
        sys.exit(0)

    channels = sorted({k[1] for k in oned})
    align_ch = args.align_channel or channels[0]
    if align_ch not in channels:
        raise ValueError(f'--align_channel {align_ch} not in {channels}')
    thr = parse_thresholds(args.z_threshold, channels)

    print(f'channels: {channels}    centroid from: {align_ch}')
    print('thresholds (ALL must pass): '
          + ', '.join(f'{c} >= {thr[c]:g}' for c in channels) + '\n')

    out_1d, out_2d, summary = {}, {}, {}

    # ---- 1D ---------------------------------------------------------------
    for (sn, ch, rtz, nt) in sorted(k for k in oned if k[1] == align_ch):
        ref = oned[(sn, align_ch, rtz, nt)]
        cen = centroids_1d(ref['strf'], ref['bar_pos_deg'], args.centroid_frac)

        pk, sg = {}, {}
        include = np.isfinite(cen)
        for c in channels:
            src = oned.get((sn, c, rtz, nt))
            if src is None:
                continue
            pk[c], sg[c] = peak_z_and_sign(src['strf'])
            # this is the de facto significance test -- see peak_z_and_sign for
            # why the threshold means something different at each noisetau, and
            # why that only matters if we start comparing across nt rather than
            # treating it as a grouping key
            include &= pk[c] >= thr[c]

        print(f'  {sn} rtz={rtz:>3.0f} tau={nt:>4.0f}ms: '
              f'{int(include.sum())}/{len(cen)} ROIs pass')
        for c in sorted(pk):
            n_c = int((pk[c] >= thr[c]).sum())
            print(f'      {c}: {n_c}/{len(pk[c])} over {thr[c]:g} '
                  f'(median |z| {np.median(pk[c]):.1f}, max {pk[c].max():.1f})')
        if not include.any():
            print('      -> nothing to average')
            continue
        print(f'      centroids {np.nanmin(cen[include]):+.0f}..'
              f'{np.nanmax(cen[include]):+.0f} deg')

        summary[(sn, align_ch, rtz, nt)] = {
            'centroid': cen, 'peak_z': pk, 'sign': sg, 'include': include,
            'align_channel': align_ch,
            'field_range': (float(ref['bar_pos_deg'].min()),
                            float(ref['bar_pos_deg'].max()))}

        for c in channels:
            src = oned.get((sn, c, rtz, nt))
            if src is None:
                continue
            res = register_average_1d(src['strf'], src['bar_pos_deg'], cen,
                                      include, args.min_roi_per_area,
                                      args.grid_step_deg)
            res['t_lag'] = src['t_lag']
            pol, frac = summarise_polarity(sg[c], include)
            res['polarity_sign'], res['polarity_frac'] = pol, frac
            print(f'      {c} channel sign {pol:+.0f} '
                  f'({100*frac:.0f}% of included ROIs agree)'
                  + ('   <- MIXED, check for more than one cell type'
                     if frac < 0.9 else ''))
            out_1d[(sn, c, rtz, nt)] = res

    # ---- 2D ---------------------------------------------------------------
    for (pair, ch, nt_key, rtz_key) in sorted(k for k in twod if k[1] == align_ch):
        ref = twod[(pair, align_ch, nt_key, rtz_key)]
        cen = centroids_2d(ref['strf_2d'], ref['bar_pos_h_deg'],
                           ref['bar_pos_v_deg'], args.centroid_frac)
        include = np.isfinite(cen).all(axis=1)
        for c in channels:
            src = twod.get((pair, c, nt_key, rtz_key))
            if src is None:
                continue
            pk_c, _ = peak_z_and_sign(src['strf_2d'])
            include &= pk_c >= thr[c]
        print(f'  {pair} {rtz_key} {nt_key}: {int(include.sum())}/{len(cen)} ROIs pass')
        if not include.any():
            continue
        for c in channels:
            src = twod.get((pair, c, nt_key, rtz_key))
            if src is None:
                continue
            res = register_average_2d(src['strf_2d'], src['bar_pos_h_deg'],
                                      src['bar_pos_v_deg'], cen, include,
                                      args.min_roi_per_area, args.grid_step_deg)
            res.update(t_lag=src['t_lag'], rtz_h=src['rtz_h'], rtz_v=src['rtz_v'])
            out_2d[(pair, c, nt_key, rtz_key)] = res

    # ---- save -------------------------------------------------------------
    with h5py.File(path, 'a') as f:
        root = f[f'STRF_{args.tag}']
        if 'averaged' in root:
            del root['averaged']
        g = root.create_group('averaged')
        g.attrs.update(align_channel=align_ch,
                       z_threshold=str(thr),
                       centroid_frac=args.centroid_frac,
                       min_roi_per_area=args.min_roi_per_area,
                       polarity_flipping='none - filter sign preserved per channel',
                       grid_step_deg=str(args.grid_step_deg),
                       crop_frac_of_max=args.crop_frac_of_max)
        for (sn, ch, rtz, nt), res in out_1d.items():
            sg = g.require_group(f'oned/{sn}/{ch}/rtz_{rtz}/noisetau_{nt}')
            for k in ('grid_deg', 'mean', 'sem', 'count', 't_lag'):
                sg.create_dataset(k, data=res[k], compression='gzip')
            sg.attrs.update(n_roi=res['n_roi'],
                            polarity_sign=res['polarity_sign'],
                            polarity_frac=res['polarity_frac'])
        for (pair, ch, nt_key, rtz_key), res in out_2d.items():
            sg = g.require_group(f'twod/{pair}/{ch}/{nt_key}/{rtz_key}')
            for k in ('grid_h_deg', 'grid_v_deg', 'mean', 'sem', 'count', 't_lag'):
                sg.create_dataset(k, data=res[k], compression='gzip')
            sg.attrs['n_roi'] = res['n_roi']
        for key, s in summary.items():
            sg = g.require_group(
                f'per_roi/{key[0]}/{key[1]}/rtz_{key[2]}/noisetau_{key[3]}')
            sg.create_dataset('centroid', data=s['centroid'])
            sg.create_dataset('include', data=s['include'])
            for c, v in s['peak_z'].items():
                sg.create_dataset(f'peak_z_{c}', data=v)
            for c, v in s['sign'].items():
                sg.create_dataset(f'sign_{c}', data=v)
    print(f'\nsaved to {path} under /STRF_{args.tag}/averaged')

    if save_figs:
        for key, res in out_1d.items():
            figure_1d(res, key, figs_dir, args.tag, crop_frac=args.crop_frac_of_max)
        # one extra figure per group with every channel side by side, so the
        # shapes can be compared without flipping between files
        groups = {}
        for (sn_, ch_, rtz_, nt_), res in out_1d.items():
            groups.setdefault((sn_, rtz_, nt_), {})[ch_] = res
        for gkey, by_ch in groups.items():
            figure_1d_compare(by_ch, gkey, figs_dir, args.tag,
                              crop_frac=args.crop_frac_of_max)
        for key, res in out_2d.items():
            figure_2d(res, key, figs_dir, args.tag,
                      crop_frac=args.crop_frac_of_max)   # peak lag + lag grid
            movie_2d(res, key, figs_dir, args.tag,
                     crop_frac=args.crop_frac_of_max)    # one video per channel
        if summary:
            figure_inclusion(summary, thr, figs_dir, args.tag)
        print(f'figures -> {figs_dir}')
