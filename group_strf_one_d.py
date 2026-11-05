"""
group_strf_one_d.py

Pool 1D ternary STRFs across flies, and compare named groups of flies -- a
genotype against its control, say.

WHAT IS POOLED, AND WHY IT IS ROIS AND NOT FLY MEANS. Every roi is registered to
its OWN receptive-field center before averaging, so the registered coordinate
already means the same thing in every fly: degrees from that roi's center. Flies
therefore concatenate without any cross-fly alignment step -- unlike
average_strf_noizone, which has to rotate on the sphere. Pooling the rois and
registering once is not the same as averaging each fly's mean: the second
weights a fly with nine rois the same as one with two hundred, and loses the
within-fly spread that the standard error is built from.

WHICH ROIS. Whatever selection the per-fly step already made. analyze_data_strf_
one_d.py writes the chosen roi indices for every variant it builds --
`all`, `all_m10`, `all_m20_maj` and so on -- to
/STRF_ONE_D_AVG_{tag}/{sn}/{ch}/geom_{glabel}/{variant}/rois, and this reads them
back rather than re-deriving the criterion. One definition of 'central', applied
once, in the script that can see the stimulus geometry.

Note a margin means something slightly different in each fly, because it is
measured from that fly's own observed extent. That is deliberate -- the extent is
the only evidence of what each fly could see -- but it means `m10` is 'each fly's
central population', not 'one fixed region of the visual field'.

THE REPLICATE, AND WHY THIS REPORTS TWO OF THEM. Rois within a fly are not
independent: they share a preparation, a focal plane, an indicator expression
level and a registration. Pooling them is the right way to DESCRIBE the mean
filter -- more rois genuinely means a better-estimated shape. It is the wrong
denominator for TESTING a difference between genotypes, where treating 300 rois
from three flies as 300 independent samples would overstate the evidence by
roughly the number of rois per fly.

So both are reported, and they answer different questions:

  * the pooled mean, with a standard error over rois -- what the filter looks
    like
  * one kinetics number per FLY, with the flies shown individually -- whether
    the genotypes differ

With a handful of flies per group the second is a small-n comparison and is
presented as the individual values rather than as a test statistic. Three flies
against three does not support a p-value worth quoting, and printing one would
lend it an authority the design cannot carry.

Usage:
    python group_strf_one_d.py \\
        --group "control=D:/data/fly_016" \\
        --group "pdeA=D:/data/fly_008,D:/data/fly_010" \\
        --tag raw --variant all --output_directory D:/data/group
"""

import argparse
import os
import sys

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from visanalysis.util.one_d_group import (auto_grid_step, register_average,
                                          coverage_bounds)
from visanalysis.util.strf_plot import apply_figure_style

apply_figure_style()

#: Distinct, and readable next to each other on a white ground.
GROUP_COLORS = ['#1f77b4', '#d62728', '#2ca02c', '#9467bd', '#ff7f0e']


def parse_group(spec):
    """'name=dir1,dir2' -> (name, [dir1, dir2])."""
    if '=' not in spec:
        raise ValueError("--group wants 'name=dir1,dir2', got {!r}".format(spec))
    name, dirs = spec.split('=', 1)
    paths = [d.strip() for d in dirs.split(',') if d.strip()]
    if not paths:
        raise ValueError('--group {} lists no directories'.format(name))
    return name.strip(), paths


def load_fly(fly_dir, tag, variant):
    """Registered inputs for one fly, keyed (channel, geometry label).

    Returns {} when the fly has no results for this tag, so a group can name a
    fly that has not been analysed yet without the whole run failing.
    """
    path = os.path.join(fly_dir, 'fly.hdf5')
    if not os.path.exists(path):
        print('    no fly.hdf5 in {}'.format(fly_dir))
        return {}
    raw_g, avg_g = 'STRF_ONE_D_' + tag, 'STRF_ONE_D_AVG_' + tag
    out = {}
    with h5py.File(path, 'r') as f:
        if raw_g not in f or avg_g not in f:
            print('    {} has no /{} or /{}'.format(
                os.path.basename(fly_dir), raw_g, avg_g))
            return {}
        # the avg group holds both {sn}/{ch}/geom_... and a top-level {ch}/radial
        for sn in [k for k in f[avg_g] if k.startswith('sn')]:
            for ch in f[avg_g][sn]:
                for geom in f[avg_g][sn][ch]:
                    node = f[avg_g][sn][ch][geom]
                    if variant not in node or 'rois' not in node[variant]:
                        continue
                    rois = np.asarray(node[variant]['rois'][()], dtype=int)
                    if rois.size < 3:
                        continue
                    glabel = geom[len('geom_'):] if geom.startswith('geom_') else geom
                    rg = f[raw_g][sn][ch][geom]
                    vis = np.asarray(rg['visible_cells'][()], dtype=int)
                    out[ch, glabel] = dict(
                        filters=rg['strf_z'][()][np.ix_(rois, vis)],
                        positions=np.asarray(rg['screen_offset_deg'][()])[vis],
                        centers=np.asarray(rg['center_screen_deg'][()])[rois],
                        t_lag=rg['t_lag'][()],
                        cell_width=float(rg.attrs['cell_width_deg']),
                        mode=str(rg.attrs['mode']),
                        n_roi=int(rois.size), fly=os.path.basename(fly_dir))
    return out


def kinetics(grid, mean, t_lag):
    """Peak lag and the following zero crossing, read at the roi's own center.

    Both are read from the row nearest 0 degrees rather than from whichever row
    happens to be largest: the centre row is the same quantity in every group,
    which is what makes them comparable.
    """
    if grid.size == 0:
        return np.nan, np.nan
    i = int(np.nanargmin(np.abs(grid)))
    tc = mean[i]
    if not np.any(np.isfinite(tc)):
        return np.nan, np.nan
    pk = int(np.nanargmax(np.abs(tc)))
    zc = np.flatnonzero(tc[pk:] * np.sign(tc[pk]) < 0)
    return (float(t_lag[pk]),
            float(t_lag[pk + zc[0]]) if zc.size else np.nan)


def pooled(entries, step, min_roi_per_bin):
    """One registered mean over every roi of every fly in `entries`."""
    base = entries[0]
    filters = np.concatenate([e['filters'] for e in entries], axis=0)
    centers = np.concatenate([e['centers'] for e in entries], axis=0)
    out = register_average(filters, base['positions'], centers,
                           base['cell_width'], mode=base['mode'],
                           grid_step=step, min_roi_per_bin=min_roi_per_bin)
    out['t_lag'] = base['t_lag']
    out['n_fly'] = len(entries)
    out['flies'] = [e['fly'] for e in entries]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--group', action='append', required=True,
                    help="'name=dir1,dir2'. Repeat for each group. With more "
                         'than one, a comparison figure is drawn as well.')
    ap.add_argument('--tag', nargs='?', default='raw')
    ap.add_argument('--variant', nargs='?', default='all',
                    help="which per-fly selection to pool: 'all', 'all_m10', "
                         "'all_m20_maj' and so on, exactly as the per-fly step "
                         'named them. A fly missing that variant is skipped and '
                         'said so, rather than silently contributing nothing.')
    ap.add_argument('--output_directory', nargs='?', required=True)
    ap.add_argument('--grid_step_deg', nargs='?', type=float, default=None,
                    help='bin width for every group. Omit to derive it from the '
                         'largest pooled population, which is then applied to '
                         'all of them -- groups binned differently cannot be '
                         'compared.')
    ap.add_argument('--min_roi_per_bin', nargs='?', type=int, default=3)
    ap.add_argument('--mean_span_deg', nargs='?', type=float, default=75.0)
    ap.add_argument('--save_figs', nargs='?', default='True')
    args = ap.parse_args()

    save_figs = str(args.save_figs) == 'True'
    os.makedirs(args.output_directory, exist_ok=True)
    groups = [parse_group(g) for g in args.group]

    # ── load ────────────────────────────────────────────────────────────────
    loaded = {}
    for name, dirs in groups:
        print('group {!r}: {} flies'.format(name, len(dirs)))
        per_fly = []
        for d in dirs:
            got = load_fly(d, args.tag, args.variant)
            if not got:
                continue
            print('    {:<12} {} geometry/channel combination(s), {} rois'.format(
                os.path.basename(d), len(got),
                sum(v['n_roi'] for v in got.values())))
            per_fly.append((os.path.basename(d), got))
        if not per_fly:
            print('    nothing loaded for {!r}'.format(name))
        loaded[name] = per_fly

    # every (channel, geometry) any group has, and who has it
    keys = sorted({k for pf in loaded.values() for _, g in pf for k in g})
    if not keys:
        print('\nNo pooled data. Does --variant {!r} exist in these flies?'
              .format(args.variant))
        sys.exit(1)

    # ── one bin width for everything ────────────────────────────────────────
    if args.grid_step_deg is not None:
        step = float(args.grid_step_deg)
        why = 'given'
    else:
        biggest, step = None, None
        for k in keys:
            for name, pf in loaded.items():
                es = [g[k] for _, g in pf if k in g]
                if not es:
                    continue
                n = sum(e['n_roi'] for e in es)
                if biggest is None or n > biggest:
                    rel = np.concatenate(
                        [e['positions'][None, :] - e['centers'][:, None]
                         for e in es], axis=0)
                    biggest = n
                    step = auto_grid_step(es[0]['cell_width'], rel)
        why = 'from the largest pooled population ({} rois)'.format(biggest)
    print('\nbin width {:g} deg ({})'.format(step, why))

    # ── pool, and measure each fly on its own ───────────────────────────────
    results, per_fly_kin = {}, []
    for k in keys:
        for name, pf in loaded.items():
            es = [g[k] for _, g in pf if k in g]
            if not es:
                continue
            nl = min(e['t_lag'].size for e in es)
            for e in es:
                e['filters'] = e['filters'][:, :, :nl]
                e['t_lag'] = e['t_lag'][:nl]
            results[k, name] = pooled(es, step, args.min_roi_per_bin)
            # and the same measurement fly by fly, which is the unit a
            # between-group comparison is entitled to use
            for fly, g in pf:
                if k not in g:
                    continue
                one = pooled([g[k]], step, args.min_roi_per_bin)
                ks = coverage_bounds(one['n_roi_per_bin'])
                pl, zc = kinetics(one['grid_deg'][ks], one['mean'][ks],
                                  one['t_lag'])
                per_fly_kin.append(dict(ch=k[0], geom=k[1], group=name, fly=fly,
                                        n_roi=g[k]['n_roi'], peak_lag=pl,
                                        zero_cross=zc))

    # ── report ──────────────────────────────────────────────────────────────
    print('\nPOOLED (roi as replicate -- describes the filter)')
    print('{:<5} {:<26} {:<10} {:>6} {:>7} {:>10} {:>12}'.format(
        'ch', 'geometry', 'group', 'flies', 'rois', 'peak lag', 'zero cross'))
    for (k, name), out in sorted(results.items()):
        ks = coverage_bounds(out['n_roi_per_bin'])
        pl, zc = kinetics(out['grid_deg'][ks], out['mean'][ks], out['t_lag'])
        print('{:<5} {:<26} {:<10} {:>6} {:>7} {:>8.0f}ms {:>12}'.format(
            k[0], k[1][:26], name, out['n_fly'], out['n_roi'], pl * 1000,
            '{:.0f}ms'.format(zc * 1000) if np.isfinite(zc) else 'none'))

    print('\nPER FLY (fly as replicate -- the unit a group comparison may use)')
    print('{:<5} {:<26} {:<10} {:<12} {:>6} {:>10} {:>12}'.format(
        'ch', 'geometry', 'group', 'fly', 'rois', 'peak lag', 'zero cross'))
    for r in sorted(per_fly_kin, key=lambda x: (x['ch'], x['geom'], x['group'],
                                                x['fly'])):
        print('{:<5} {:<26} {:<10} {:<12} {:>6} {:>8.0f}ms {:>12}'.format(
            r['ch'], r['geom'][:26], r['group'], r['fly'], r['n_roi'],
            r['peak_lag'] * 1000,
            '{:.0f}ms'.format(r['zero_cross'] * 1000)
            if np.isfinite(r['zero_cross']) else 'none'))

    # ── figures ─────────────────────────────────────────────────────────────
    if save_figs:
        _figures(results, per_fly_kin, keys, loaded, args, step)

    # ── save ────────────────────────────────────────────────────────────────
    out_h5 = os.path.join(args.output_directory,
                          'group_strf_one_d_{}_{}.hdf5'.format(args.tag,
                                                               args.variant))
    with h5py.File(out_h5, 'w') as f:
        f.attrs['variant'] = args.variant
        f.attrs['tag'] = args.tag
        f.attrs['grid_step_deg'] = step
        f.attrs['replicate'] = ('rois pooled across flies; per-fly kinetics '
                                'stored separately for group comparison')
        for (k, name), out in results.items():
            g = f.require_group('{}/{}/{}'.format(name, k[0], k[1]))
            for key in ('grid_deg', 'mean', 'sem', 'n_roi_per_bin', 't_lag'):
                g.create_dataset(key, data=out[key])
            g.attrs['n_roi'] = out['n_roi']
            g.attrs['n_fly'] = out['n_fly']
            g.attrs['flies'] = ', '.join(out['flies'])
        kg = f.require_group('per_fly')
        for col in ('ch', 'geom', 'group', 'fly'):
            kg.create_dataset(col, data=np.array(
                [str(r[col]) for r in per_fly_kin], dtype='S64'))
        for col in ('n_roi', 'peak_lag', 'zero_cross'):
            kg.create_dataset(col, data=np.array(
                [r[col] for r in per_fly_kin], dtype=float))
    print('\nSaved {}'.format(out_h5))


def _figures(results, per_fly_kin, keys, loaded, args, step):
    span = float(args.mean_span_deg)
    names = [n for n, _ in (list(loaded.items()))]
    color = {n: GROUP_COLORS[i % len(GROUP_COLORS)] for i, n in enumerate(names)}

    # shared limits per channel, so every group and geometry is comparable
    lim = {}
    for (k, name), out in results.items():
        ks = coverage_bounds(out['n_roi_per_bin'])
        if ks.size < 3:
            continue
        g, m, e = out['grid_deg'][ks], out['mean'][ks], out['sem'][ks]
        i = np.abs(g) <= span
        g, m, e = g[i], m[i], e[i]
        cur = lim.get(k[0])
        new = (float(g[0]), float(g[-1]), float(np.nanmin(m - e)),
               float(np.nanmax(m + e)), float(np.nanmax(np.abs(m))))
        lim[k[0]] = new if cur is None else (
            min(cur[0], new[0]), max(cur[1], new[1]), min(cur[2], new[2]),
            max(cur[3], new[3]), max(cur[4], new[4]))

    n_fig = 0
    for k in keys:
        present = [n for n in names if (k, n) in results]
        if not present:
            continue
        glo, ghi, ylo, yhi, vmax = lim[k[0]]
        pad = 0.05 * (yhi - ylo) if yhi > ylo else 1.0

        # one pooled figure per group
        for name in present:
            out = results[k, name]
            ks = coverage_bounds(out['n_roi_per_bin'])
            g, m, e, t = (out['grid_deg'][ks], out['mean'][ks],
                          out['sem'][ks], out['t_lag'])
            vis = (g >= glo) & (g <= ghi)
            mv = np.where(vis[:, None], np.abs(m), np.nan)
            if not np.any(np.isfinite(mv)):
                continue
            li = int(np.nanargmax(np.nanmax(mv, axis=0)))
            pi = int(np.nanargmax(mv[:, li]))
            fig, ax = plt.subplots(1, 3, figsize=(17.5, 4.4),
                                   constrained_layout=True)
            im = ax[0].imshow(m, aspect='auto', origin='lower', cmap='RdBu_r',
                              vmin=-vmax, vmax=vmax,
                              extent=[t[0] * 1000, t[-1] * 1000,
                                      g[0] - step / 2, g[-1] + step / 2])
            ax[0].set_ylim(glo, ghi)
            ax[0].axhline(0, color='k', lw=0.8, ls=':')
            ax[0].set_xlabel('Lag (ms)')
            ax[0].set_ylabel("deg re. each roi's own RF center")
            fig.colorbar(im, ax=ax[0], label='z')
            ax[1].plot(g, m[:, li], color=color[name], lw=2)
            ax[1].fill_between(g, m[:, li] - e[:, li], m[:, li] + e[:, li],
                               color=color[name], alpha=0.25, lw=0)
            ax[1].axhline(0, color='0.6', lw=0.8)
            ax[1].axvline(0, color='0.6', lw=0.8, ls=':')
            ax[1].set_xlim(glo, ghi)
            ax[1].set_ylim(ylo - pad, yhi + pad)
            ax[1].set_xlabel('deg re. RF center')
            ax[1].set_ylabel('mean z')
            ax[1].set_title('spatial, at {:.0f} ms'.format(t[li] * 1000),
                            fontsize=14)
            ax[2].plot(t * 1000, m[pi], color=color[name], lw=2)
            ax[2].fill_between(t * 1000, m[pi] - e[pi], m[pi] + e[pi],
                               color=color[name], alpha=0.25, lw=0)
            ax[2].axhline(0, color='0.6', lw=0.8)
            ax[2].set_ylim(ylo - pad, yhi + pad)
            ax[2].set_xlabel('Lag (ms)')
            ax[2].set_ylabel('mean z')
            ax[2].set_title('temporal, at {:+.1f} deg'.format(g[pi]),
                            fontsize=14)
            fig.suptitle('{}  {}  {}  [{}]   {} flies, {} rois pooled'.format(
                name, k[0], k[1][:30], args.variant, out['n_fly'],
                out['n_roi']), fontsize=15)
            fig.savefig(os.path.join(
                args.output_directory,
                'group_MEAN_{}_{}_{}_{}.pdf'.format(name, k[0], k[1],
                                                    args.variant)))
            plt.close(fig)
            n_fig += 1

        # the comparison, when there is more than one group to compare
        if len(present) < 2:
            continue
        fig, ax = plt.subplots(1, 3, figsize=(17.5, 4.4), constrained_layout=True)
        for name in present:
            out = results[k, name]
            ks = coverage_bounds(out['n_roi_per_bin'])
            g, m, e, t = (out['grid_deg'][ks], out['mean'][ks],
                          out['sem'][ks], out['t_lag'])
            vis = (g >= glo) & (g <= ghi)
            mv = np.where(vis[:, None], np.abs(m), np.nan)
            li = int(np.nanargmax(np.nanmax(mv, axis=0)))
            pi = int(np.nanargmax(mv[:, li]))
            lab = '{} ({} flies, {} rois)'.format(name, out['n_fly'],
                                                  out['n_roi'])
            ax[0].plot(g, m[:, li], color=color[name], lw=2, label=lab)
            ax[0].fill_between(g, m[:, li] - e[:, li], m[:, li] + e[:, li],
                               color=color[name], alpha=0.22, lw=0)
            ax[1].plot(t * 1000, m[pi], color=color[name], lw=2, label=lab)
            ax[1].fill_between(t * 1000, m[pi] - e[pi], m[pi] + e[pi],
                               color=color[name], alpha=0.22, lw=0)
        ax[0].axhline(0, color='0.6', lw=0.8)
        ax[0].axvline(0, color='0.6', lw=0.8, ls=':')
        ax[0].set_xlim(glo, ghi)
        ax[0].set_ylim(ylo - pad, yhi + pad)
        ax[0].set_xlabel('deg re. RF center')
        ax[0].set_ylabel('mean z')
        ax[0].set_title('spatial', fontsize=14)
        ax[0].legend(fontsize=11)
        ax[1].axhline(0, color='0.6', lw=0.8)
        ax[1].set_ylim(ylo - pad, yhi + pad)
        ax[1].set_xlabel('Lag (ms)')
        ax[1].set_ylabel('mean z')
        ax[1].set_title('temporal', fontsize=14)

        # the honest denominator: one point per fly, shown individually
        rows = [r for r in per_fly_kin if r['ch'] == k[0] and r['geom'] == k[1]]
        for xi, name in enumerate(present):
            v = [r['peak_lag'] * 1000 for r in rows
                 if r['group'] == name and np.isfinite(r['peak_lag'])]
            if not v:
                continue
            ax[2].scatter(np.full(len(v), xi) + np.linspace(-.06, .06, len(v)),
                          v, s=70, color=color[name], edgecolors='k',
                          linewidths=0.5, zorder=3)
            ax[2].plot([xi - .2, xi + .2], [np.median(v)] * 2, color='k', lw=2)
        ax[2].set_xticks(range(len(present)))
        ax[2].set_xticklabels(present)
        ax[2].set_ylabel('peak lag (ms)')
        ax[2].set_title('one point per FLY\n(bar = median)', fontsize=14)
        fig.suptitle('{}  {}  [{}]   pooled means share axes; the fly panel is '
                     'the comparison unit'.format(k[0], k[1][:30], args.variant),
                     fontsize=15)
        fig.savefig(os.path.join(
            args.output_directory,
            'group_COMPARE_{}_{}_{}.pdf'.format(k[0], k[1], args.variant)))
        plt.close(fig)
        n_fig += 1
    print('\n{} figures written to {}'.format(n_fig, args.output_directory))


if __name__ == '__main__':
    main()
