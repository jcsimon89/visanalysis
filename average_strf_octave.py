"""
average_strf_octave.py

Registers each ROI's merged octave-ternary filter to its own receptive-field
centroid and averages across ROIs, giving a mean RF for a cell type.

Standalone final step, like average_strf_noizone.py: it reads the filters
already written into the hdf5 by analyze_data_strf_octave.py, so it can be
re-run on its own while tuning thresholds without recomputing anything.

    python average_strf_octave.py \
        --experiment_file_directory "path/to/fly_folder" \
        --tag final \
        --align_channel ch1 \
        --z_threshold "ch1:6,ch2:3.5"

Results go to /STRF_OCTAVE_AVG_{tag}, figures to {tag}_roi_figs.

WHY THIS IS NOT average_strf_noizone WITH DIFFERENT COLUMNS
-----------------------------------------------------------
* Cells live on a SPHERE, not on a bar axis. Registration is a rotation that
  carries each ROI's centroid onto a common reference direction, not a
  subtraction of positions. Subtracting azimuths would be wrong by cos(elevation)
  and the lit band reaches 51 degrees from the horizon.

* Cell centres are IRREGULAR and differ between flies, because the tessellation
  is rotated per run. There is no shared sampling grid to average on, so values
  are binned into a common RF-centred grid -- never interpolated, exactly as the
  noizone path bins bar samples.

* The filter being averaged is the MERGED one, which already folds both octaves
  into a single estimate on the fine cells, so there is one map rather than two
  to reconcile. The per-octave filters stay in the hdf5 alongside it if the two
  layers ever need to be compared directly.

WHAT THE AVERAGE IS AND IS NOT SAFE FOR
---------------------------------------
Position is the robust part. Centroids survive even the pathological
non-separable case to within 0.3 degrees, because centre and surround are
concentric, so alignment is trustworthy. Polarity and apparent size read from a
single lag are NOT: where a filter is non-separable within an octave the smear
can change which lobe owns the peak, flipping sign and inflating size threefold.
This script therefore reports the lag profile alongside the map, and reports a
separability statistic per ROI -- how much of the filter one spatial profile
times one time course accounts for. Read it comparatively. It is SNR-limited, so
an exactly separable filter at realistic noise scores about 0.5, and the useful
signal is which ROIs sit low relative to the rest rather than the absolute value.
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

plt.ioff()


# ─────────────────────────────────────────────────────────────────────────────
# loading
# ─────────────────────────────────────────────────────────────────────────────

def load_octave_strf(path, tag):
    """Everything analyze_data_strf_octave wrote, keyed (sn, ch, fib)."""
    group = 'STRF_OCTAVE_' + tag
    out = {}
    with h5py.File(path, 'r') as f:
        if group not in f:
            raise KeyError(
                'no /{} in {}. Run analyze_data_strf_octave.py first.'.format(group, path))
        g = f[group]
        attrs = dict(g.attrs)
        for sn in g:
            if sn == 'combined':
                continue
            for ch in g[sn]:
                for fib in g[sn][ch]:
                    sg = g[sn][ch][fib]
                    if 'strf_merged' not in sg:
                        continue
                    rec = {'strf_merged': sg['strf_merged'][()],
                           'azel': sg['merged_azel_deg'][()],
                           't_lag': sg['t_lag'][()],
                           'peak_lag_ind': sg['peak_lag_ind'][()],
                           'octaves': {}}
                    for key in sg:
                        if key.startswith('octave_'):
                            og = sg[key]
                            rec['octaves'][key[len('octave_'):]] = {
                                'centroid': og['centroid_azel_deg'][()],
                                'n_cells_in_centroid': (
                                    og['n_cells_in_centroid'][()]
                                    if 'n_cells_in_centroid' in og else None),
                                'patch_deg': float(og.attrs['patch_deg']),
                                'cell_start': int(og.attrs['cell_start']),
                                'cell_stop': int(og.attrs['cell_stop'])}
                    out[sn, ch, fib] = rec
    return out, attrs


def parse_thresholds(spec, channels, default=4.0):
    """"ch1:6,ch2:3.5" or a single number for all channels."""
    if spec is None or spec == '':
        return {c: default for c in channels}
    if ':' not in spec:
        return {c: float(spec) for c in channels}
    out = {c: default for c in channels}
    for part in spec.split(','):
        k, v = part.split(':')
        out[k.strip()] = float(v)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# geometry
# ─────────────────────────────────────────────────────────────────────────────

def azel_to_vec(az_deg, el_deg):
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], axis=-1)


def rotation_carrying(a, b):
    """Rotation matrix taking unit vector a onto unit vector b.

    Rodrigues about their common perpendicular -- the minimal rotation, so it
    introduces no spin about the receptive-field axis that would smear the
    average azimuthally.
    """
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-12:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def register_to_centroid(azel, centroid_azel):
    """Cell directions rotated so the ROI's centroid sits at (0, 0).

    Returns RF-centred (az, el) in degrees. Near the origin this is a faithful
    local coordinate; the distortion that would matter far from it is the reason
    the rotation is done on the sphere rather than by subtracting angles.
    """
    dirs = azel_to_vec(azel[:, 0], azel[:, 1])
    R = rotation_carrying(azel_to_vec(*centroid_azel), np.array([1.0, 0.0, 0.0]))
    d = dirs @ R.T
    az = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    el = np.degrees(np.arcsin(np.clip(d[:, 2], -1, 1)))
    return np.stack([az, el], axis=1)


def separability(block, min_cells=8):
    """Fraction of the (cells x lags) filter captured by one spatial profile.

    Purely descriptive: no receptive-field shape is assumed, it only asks whether
    the spatial pattern keeps its shape as lag changes. Where it does not, the
    single-lag polarity and size are unreliable, because the smear can change
    which lobe owns the peak.

    Computed only on cells that carry signal. A receptive field covers on the
    order of 1% of the fine cells, so over the whole block this statistic
    measures signal-to-noise and not separability at all -- an exactly separable
    planted filter at realistic noise scored 0.047 that way, which would flag
    every ROI ever recorded. Cells are selected by MAD, which the small
    signal-carrying fraction cannot inflate, and NaN is returned when too few
    qualify to make the question meaningful.
    """
    if not np.any(block):
        return np.nan
    mad = 1.4826 * np.median(np.abs(block - np.median(block)))
    if mad <= 0:
        return np.nan
    strong = np.abs(block).max(axis=1) > 4.0 * mad
    if strong.sum() < min_cells:
        return np.nan
    sv = np.linalg.svd(block[strong], compute_uv=False)
    return float(sv[0] ** 2 / (sv ** 2).sum())


# ─────────────────────────────────────────────────────────────────────────────
# averaging
# ─────────────────────────────────────────────────────────────────────────────

def register_average(strf, azel, centroids, include, grid_step_deg,
                     half_extent_deg, min_roi_per_bin, peak_ind):
    """Bin each included ROI's cells by RF-centred position, then average.

    Values are never interpolated. Each measured cell value is assigned to the
    bin its registered position falls in, and averaged WITHIN an ROI first, so a
    bin catching two of the same ROI's cells counts as one ROI -- which is what
    makes min_roi_per_bin count ROIs rather than samples.
    """
    edges = np.arange(-half_extent_deg, half_extent_deg + grid_step_deg, grid_step_deg)
    centres = 0.5 * (edges[:-1] + edges[1:])
    n = len(centres)
    n_lags = strf.shape[2]

    total = np.zeros((n, n, n_lags))
    total_sq = np.zeros((n, n, n_lags))
    n_roi_hit = np.zeros((n, n), dtype=int)

    for r in np.flatnonzero(include):
        c = centroids[r]
        if not np.all(np.isfinite(c)):
            continue
        rc = register_to_centroid(azel, c)
        ix = np.digitize(rc[:, 0], edges) - 1
        iy = np.digitize(rc[:, 1], edges) - 1
        ok = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
        if not ok.any():
            continue
        # per-ROI bin means first
        acc = np.zeros((n, n, n_lags))
        cnt = np.zeros((n, n))
        np.add.at(acc, (iy[ok], ix[ok]), strf[r][ok, :])
        np.add.at(cnt, (iy[ok], ix[ok]), 1.0)
        hit = cnt > 0
        m = np.zeros_like(acc)
        m[hit] = acc[hit] / cnt[hit][:, None]
        total[hit] += m[hit]
        total_sq[hit] += m[hit] ** 2
        n_roi_hit[hit] += 1

    with np.errstate(invalid='ignore', divide='ignore'):
        nn = np.maximum(n_roi_hit, 1)[:, :, None]
        mean = np.where(n_roi_hit[:, :, None] > 0, total / nn, np.nan)
        var = np.where(n_roi_hit[:, :, None] > 1, total_sq / nn - mean ** 2, np.nan)
        sem = np.sqrt(np.maximum(var, 0.0)) / np.sqrt(nn)
    thin = n_roi_hit < min_roi_per_bin
    mean[thin] = np.nan
    sem[thin] = np.nan
    return {'mean': mean, 'sem': sem, 'count': n_roi_hit, 'grid_deg': centres,
            'n_roi_included': int(include.sum())}


# ─────────────────────────────────────────────────────────────────────────────
# figures
# ─────────────────────────────────────────────────────────────────────────────

def figure_average(res, key, out_dir, tag, t_lag, fmt='.pdf'):
    sn, ch, fib = key
    mean, cnt, grid = res['mean'], res['count'], res['grid_deg']
    finite = np.isfinite(mean)
    if not finite.any():
        return
    power = np.nansum(mean ** 2, axis=(0, 1))
    li = int(np.argmax(power))
    ext = [grid[0], grid[-1], grid[0], grid[-1]]

    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.4))
    m = mean[:, :, li]
    lim = np.nanmax(np.abs(m)) or 1.0
    im = ax[0].imshow(m, origin='lower', extent=ext, cmap='RdBu_r',
                      vmin=-lim, vmax=lim, interpolation='nearest')
    ax[0].set_title('mean RF @ lag {:.2f}s   n = {} roi'.format(
        t_lag[li], res['n_roi_included']))
    ax[0].set_xlabel('deg from RF centre (az)')
    ax[0].set_ylabel('deg from RF centre (el)')
    fig.colorbar(im, ax=ax[0], fraction=0.046)

    im = ax[1].imshow(cnt, origin='lower', extent=ext, cmap='viridis',
                      interpolation='nearest')
    ax[1].set_title('ROIs per bin')
    ax[1].set_xlabel('deg from RF centre (az)')
    fig.colorbar(im, ax=ax[1], fraction=0.046)

    # radial profile: the summary that does not depend on choosing a lag well
    yy, xx = np.meshgrid(grid, grid, indexing='ij')
    rad = np.sqrt(xx ** 2 + yy ** 2)
    bins = np.arange(0, grid[-1], max(grid[1] - grid[0], 1.0))
    prof, err = [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        sel = (rad >= lo) & (rad < hi) & np.isfinite(mean[:, :, li])
        prof.append(np.nanmean(mean[:, :, li][sel]) if sel.any() else np.nan)
        err.append(np.nanstd(mean[:, :, li][sel]) / max(np.sqrt(sel.sum()), 1)
                   if sel.any() else np.nan)
    mid = 0.5 * (bins[:-1] + bins[1:])
    prof, err = np.array(prof), np.array(err)
    ax[2].plot(mid, prof, 'k', lw=1.8)
    ax[2].fill_between(mid, prof - err, prof + err, color='0.7', alpha=0.5)
    ax[2].axhline(0, color='0.6', lw=0.8)
    ax[2].set_xlabel('deg from RF centre')
    ax[2].set_ylabel('mean weight')
    ax[2].set_title('radial profile')

    fig.suptitle('{} {} {}  merged octaves, registered mean'.format(sn, ch, fib))
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, 'avg_strf_octave_{}_{}_{}{}'.format(
        sn, ch, fib, fmt)))
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# main
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--experiment_file_directory', nargs='?', required=True)
    ap.add_argument('--tag', nargs='?', default='final', choices=['raw', 'final'])
    ap.add_argument('--align_channel', nargs='?', default=None,
                    help='channel supplying the RF centroid; its registration is '
                         'applied to every channel, so the channels stay comparable')
    ap.add_argument('--align_octave', nargs='?', default=None,
                    help='which octave supplies the centroid. Default: the finest, '
                         'which localises best')
    ap.add_argument('--z_threshold', nargs='?', default=None,
                    help='one number for all channels or "ch1:6,ch2:3.5". ALL must pass')
    ap.add_argument('--min_roi_per_bin', nargs='?', type=int, default=3)
    ap.add_argument('--grid_step_deg', nargs='?', type=float, default=3.0,
                    help='output bin size. Values are binned, never interpolated; a '
                         'grid finer than the cells is resolvable because ROI '
                         'centroids sit at different sub-cell offsets')
    ap.add_argument('--half_extent_deg', nargs='?', type=float, default=60.0)
    ap.add_argument('--separability_min', nargs='?', type=float, default=0.0,
                    help='drop ROIs whose filter is less separable than this. 0, the '
                         'default, keeps everything and only reports. Use it to compare '
                         'ROIs within a dataset, not against an absolute standard: the '
                         'statistic is SNR-limited, and an exactly separable filter at '
                         'realistic noise scores around 0.5 rather than 1')
    ap.add_argument('--save_figs', nargs='?', default='True')
    args = ap.parse_args()

    path = os.path.join(args.experiment_file_directory,
                        'fly.hdf5' if args.tag == 'raw' else 'fly_' + args.tag + '.hdf5')
    if not os.path.exists(path):
        print('No such file: ' + path)
        sys.exit(1)

    data, attrs = load_octave_strf(path, args.tag)
    if not data:
        print('No merged octave STRFs found. Did analyze_data_strf_octave.py run?')
        sys.exit(0)

    channels = sorted({k[1] for k in data})
    thresholds = parse_thresholds(args.z_threshold, channels)
    print('found {} (series, channel, geometry) groups; channels {}'.format(
        len(data), channels))
    print('z thresholds: {}'.format(thresholds))
    if not attrs.get('deconvolved', False):
        print('NOTE: filters were not deconvolved, so peak lags read late by '
              'roughly 0.1-0.3 s. Positions are unaffected.')

    # ── per-ROI inclusion, decided once per (series, geometry) ───────────────
    results = {}
    for (sn, ch, fib), rec in sorted(data.items()):
        strf = rec['strf_merged']
        n_roi = strf.shape[0]
        oct_names = sorted(rec['octaves'], key=lambda k: rec['octaves'][k]['patch_deg'])
        align_oct = args.align_octave or oct_names[0]      # finest by default
        if align_oct not in rec['octaves']:
            print('  {} {} {}: no octave {!r}, skipping'.format(sn, ch, fib, align_oct))
            continue

        # z of the merged filter, per ROI, against its own scatter
        sd = strf.reshape(n_roi, -1).std(axis=1)
        sd = np.where(sd > 0, sd, 1.0)
        peak_z = np.array([np.abs(strf[r]).max() / sd[r] for r in range(n_roi)])
        sep = np.array([separability(strf[r]) for r in range(n_roi)])
        centroids = rec['octaves'][align_oct]['centroid']

        include = (peak_z >= thresholds.get(ch, 4.0))
        include &= np.all(np.isfinite(centroids), axis=1)
        if args.separability_min > 0:
            include &= (sep >= args.separability_min)

        print('  {} {} {}: {} of {} roi pass z >= {:.1f}'.format(
            sn, ch, fib, int(include.sum()), n_roi, thresholds.get(ch, 4.0)))
        _s = sep[np.isfinite(sep)]
        if _s.size:
            print('      centroid from the {} octave; separability median {:.3f}, '
                  'range {:.3f}-{:.3f}'.format(align_oct, np.median(_s), _s.min(), _s.max()))
            print('      (SNR-limited: an exactly separable planted filter at realistic '
                  'noise scores about 0.5, so a low value alone is not evidence of '
                  'non-separability -- compare ROIs against each other, not against 1)')
        else:
            print('      centroid from the {} octave; too few driven cells to judge '
                  'separability'.format(align_oct))

        if not include.any():
            continue
        res = register_average(strf, rec['azel'], centroids, include,
                               args.grid_step_deg, args.half_extent_deg,
                               args.min_roi_per_bin, rec['peak_lag_ind'])
        res['separability'] = sep
        res['peak_z'] = peak_z
        res['include'] = include
        res['t_lag'] = rec['t_lag']
        results[sn, ch, fib] = res

    if not results:
        print('\nNo ROI passed threshold in any group; nothing to average.')
        sys.exit(0)

    # ── save ─────────────────────────────────────────────────────────────────
    out_group = 'STRF_OCTAVE_AVG_' + args.tag
    with h5py.File(path, 'a') as f:
        if out_group in f:
            del f[out_group]
        g = f.create_group(out_group)
        g.attrs['grid_step_deg'] = args.grid_step_deg
        g.attrs['half_extent_deg'] = args.half_extent_deg
        g.attrs['min_roi_per_bin'] = args.min_roi_per_bin
        g.attrs['separability_min'] = args.separability_min
        for (sn, ch, fib), res in results.items():
            sg = g.require_group('{}/{}/{}'.format(sn, ch, fib))
            for k in ('mean', 'sem'):
                sg.create_dataset(k, data=res[k], compression='gzip')
            sg.create_dataset('count', data=res['count'])
            sg.create_dataset('grid_deg', data=res['grid_deg'])
            sg.create_dataset('t_lag', data=res['t_lag'])
            sg.create_dataset('separability', data=res['separability'])
            sg.create_dataset('peak_z', data=res['peak_z'])
            sg.create_dataset('included', data=res['include'])
            sg.attrs['n_roi_included'] = res['n_roi_included']
    print('\nSaved to {} under /{}'.format(path, out_group))

    if args.save_figs == 'True':
        figs_dir = os.path.join(args.experiment_file_directory, args.tag + '_roi_figs')
        os.makedirs(figs_dir, exist_ok=True)
        for key, res in results.items():
            figure_average(res, key, figs_dir, args.tag, res['t_lag'])
        print('Figures written to ' + figs_dir)
