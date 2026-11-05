"""
analyze_data_trf_contrast_switch.py

Contrast adaptation from the ContrastSwitchTernaryNoise protocol (labpack
JCS_protocol): full-field ternary noise whose contrast alternates between a low
and a high level at a fixed mean. Every block of the analysis lives in
visanalysis/util/contrast_switch.py; read its module docstring for the timing
convention, the pooled filter, and the Baccus & Meister normalisation.

Call structure (same as the other analyze_data scripts):
    python analyze_data_trf_contrast_switch.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --save_figs True \
        --tag raw \
        --output_directory "path/to/output"

Outputs, per contrast-switch series, in <output_directory>/<series>/ (or
<fly folder>/<tag>_contrast_switch/<series>/ when no output directory is given):
    contrast_switch.json   every number reported below, with 95% CIs
    contrast_switch.npz    filters, nonlinearities, gain dynamics, averages
    keep_<channel>.npy     the ROIs the gate kept
    *.png / *.pdf          figures

WHAT IT MEASURES, IN THE ORDER TO READ IT
-----------------------------------------
1. The gate: which ROIs have a high-contrast filter that beats the rolled
   null. Only these enter anything below.
2. Held-out r of the pooled filter in each contrast against the same null. A
   channel that fails has no filter and its gains mean nothing.
3. Steady-state gain per contrast (Baccus & Meister: own filter, var(F*s) =
   var(s)), the low/high ratio, and the slope over the generator range both
   contrasts share. Compare the ratio with 1 (no adaptation) and with the
   contrast ratio (complete adaptation).
4. The model-free gain ratio from the frozen noise, which uses no filter. If it
   disagrees with (3), the model is missing something.
5. Gain against time since the switch, relative to steady state.
6. The switch-triggered average of every channel: step, time constant, peak
   change. For a channel with no filter (cAMP so far) this is the readout.
7. Cross-fit: the gate run on one half of the stimulus and (3) measured on the
   other, both ways round, as a check that choosing ROIs on high-contrast data
   did not bias the ratio.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt

from visanalysis.analysis import imaging_data
from visanalysis.plugin import base as base_plugin
from visanalysis.util import contrast_switch as cs

try:
    from visanalysis.util.strf_plot import apply_figure_style
    apply_figure_style()
except Exception:
    pass
plt.ioff()

LOW, HIGH = '#3a7bd5', '#d0473c'


def _parse_map(s, cast=str):
    """'1:calcium,2:cAMP' -> {'1': 'calcium', '2': 'cAMP'}"""
    out = {}
    for item in str(s).split(','):
        if item.strip():
            k, v = item.split(':', 1)
            out[k.strip()] = cast(v.strip())
    return out


def _figure_main(res, arrs, title):
    ca = res['channels'][res['primary']]
    a = arrs[res['primary']]
    fig, ax = plt.subplots(2, 3, figsize=(15, 9))

    x = ax[0, 0]
    for q, col, nm in (('low', LOW, 'low'), ('high', HIGH, 'high')):
        h = a['h_' + q]
        x.plot(a['lags'], h / np.abs(h[a['lags'] < 0.6]).max(), color=col, lw=1.5, label=nm)
    x.axhline(0, color='0.8', lw=0.8)
    x.set_xlabel('lag (s)')
    x.set_title('{}: steady filters (peak-normalised)\nheld-out r low {:.3f} high {:.3f} (null max {:.3f}, {:.3f})'.format(
        res['primary'], ca['cv_r_low'], ca['cv_r_high'], ca['null_max_low'], ca['null_max_high']), fontsize=10)
    x.legend(frameon=False)

    x = ax[0, 1]
    for q, col, nm in (('low', LOW, 'low'), ('high', HIGH, 'high')):
        n = a['nl_' + q]
        x.errorbar(n[:, 0], n[:, 1], yerr=n[:, 2], fmt='o-', color=col, ms=4, lw=1.5, label=nm)
    lim = ca['shared_range']['range']
    x.axvspan(lim[0], lim[1], color='0.92', zorder=0)
    x.set_xlabel('generator (own filter, var = stimulus var)')
    x.set_ylabel('dF/F (each state at its own mean)')
    x.set_title('steady nonlinearity; grey = shared range\ngain ratio {:.2f} ({:.2f}-{:.2f}), shared range {:.2f} ({:.2f}-{:.2f})'.format(
        ca['gain']['ratio'], *ca['gain']['ci'], ca['shared_range']['ratio'], *ca['shared_range']['ci']), fontsize=10)
    x.legend(frameon=False)

    x = ax[0, 2]
    b = np.array(ca['gain_after_switch']['bins'])
    cen = 0.5 * (b[1:] + b[:-1])
    for nm, col, lab, off in (('to_low', LOW, 'after switch to low', 0.97), ('to_high', HIGH, 'after switch to high', 1.03)):
        g = ca['gain_after_switch'][nm]
        v, c_ = np.array(g['relative']), np.array(g['relative_ci'])
        x.plot(cen * off, v, 'o-', color=col, ms=4, lw=1.5, label=lab)
        x.fill_between(cen * off, c_[0], c_[1], color=col, alpha=0.15, lw=0)
    x.axhline(1, color='0.7', lw=0.8)
    x.axhline(0, color='0.85', lw=0.8)
    x.set_xscale('log')
    x.set_ylim(-1, 3)
    x.set_xlabel('time since switch (s)')
    x.set_ylabel('gain / steady gain')
    x.set_title('gain after a switch (95% CI)', fontsize=10)
    x.legend(frameon=False, fontsize=9)

    for j, ch in enumerate(list(res['channels'])[:2]):
        x = ax[1, j]
        s = arrs['sta_' + ch]
        sm = res['switch_triggered'][ch]
        for nm, col, lab in (('to_high', HIGH, 'to high'), ('to_low', LOW, 'to low')):
            x.plot(s['lag'], s[nm], color=col, lw=1.5, label='{} (step {:+.3f})'.format(lab, sm[nm]['step']))
            x.fill_between(s['lag'], s[nm] - s[nm + '_se'], s[nm] + s[nm + '_se'], color=col, alpha=0.2, lw=0)
        x.axvline(0, color='0.7', lw=0.8)
        x.set_xlabel('time from switch (s)')
        x.set_ylabel('dF/F')
        x.set_title('{}: switch-triggered average'.format(ch), fontsize=10)
        x.legend(frameon=False, fontsize=9)

    x = ax[1, 2]
    x.axis('off')
    mf = res.get('model_free')
    lines = ['{}: {} of {} ROIs kept by the gate'.format(res['primary'], res['gate']['n_keep'], res['gate']['n_roi']),
             'contrasts {:.2f} / {:.2f} (ratio {:.2f})'.format(*res['contrasts'], res['contrast_ratio']),
             'low/high (95% CI)  measured           vs non-adapting cell',
             '(own = each contrast filter, Baccus & Meister; hi filt = high filter in both)']
    sr = ca.get('static_reference')
    if ca['has_filter_both_states'] is False:
        lines.insert(1, '!! a filter fails its null: no ratio below is interpretable')
    for key, nm in (('gain', 'gain, own'), ('shared_range', 'shared, own'),
                    ('gain_high_filter', 'gain, hi filt'), ('shared_range_high_filter', 'shared, hi filt')):
        v = ca[key]
        lines.append('  {:<13} {:.2f} ({:.2f}-{:.2f})   {}'.format(
            nm, v['ratio'], *v['ci'], '{:.2f} ({:.2f}-{:.2f})'.format(v['ratio_vs_static'], *v['ci_vs_static'])
            if 'ratio_vs_static' in v else 'n/a'))
    if mf:
        lines.append('  {:<13} {:.2f} ({:.2f}-{:.2f})   {:.2f} ({:.2f}-{:.2f})'.format(
            'model-free', mf['gain_ratio'], *mf['gain_ratio_ci'], mf['gain_ratio_vs_static'], *mf['gain_ratio_ci_vs_static']))
    if sr:
        lines.append('a non-adapting cell reads {:.2f} / {:.2f} / {:.2f} / {:.2f} through these'.format(
            sr['gain'], sr['shared_range'], sr['gain_high_filter'], sr['shared_range_high_filter']))
        lines.append('  ({:.2f} exactly if linear; the model-free ratio is read against this)'.format(sr['truth_linear']))
    lines.append('filter trough low / high {:.0f} / {:.0f} ms'.format(ca['filter_low']['trough_lag'] * 1e3,
                                                                     ca['filter_high']['trough_lag'] * 1e3))
    lines.append('stretch low vs high {:.2f} ({:.2f}-{:.2f})'.format(ca['stretch_low_vs_high']['value'],
                                                                    *ca['stretch_low_vs_high']['ci']))
    for half in res.get('crossfit', []):
        g_ = half['gain']
        lines.append('cross-fit, measured on {:.0f}-{:.0f} s: gain vs non-adapting {}'.format(
            max(half['measure'][0], 0), min(half['measure'][1], 300), '{:.2f} ({:.2f}-{:.2f})'.format(
                g_['ratio_vs_static'], *g_['ci_vs_static']) if 'ratio_vs_static' in g_ else '{:.2f}'.format(g_['ratio'])))
    x.text(0, 1, '\n'.join(lines), va='top', fontsize=10, family='monospace', transform=x.transAxes)
    for x in ax.ravel():
        x.spines[['top', 'right']].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def _figure_gate(gate, title):
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    bins = np.linspace(min(np.nanmin(gate['r']), np.nanmin(gate['null'])), np.nanmax(gate['r']), 50)
    ax[0].hist(gate['null'], bins, density=True, color='0.7', label='null (rolled stimulus)')
    ax[0].hist(gate['r'][np.isfinite(gate['r'])], bins, density=True, histtype='step', color=HIGH, lw=1.5, label='ROIs')
    ax[0].set_xlabel('held-out r, high-contrast filter')
    ax[0].legend(frameon=False)
    ax[1].bar(np.arange(gate['plane_n'].size), gate['plane_n'], color='0.85', label='all ROIs')
    ax[1].bar(np.arange(gate['plane_k'].size), gate['plane_k'], color=HIGH, label='kept')
    ax[1].set_xlabel('z plane')
    ax[1].legend(frameon=False)
    for x in ax:
        x.spines[['top', 'right']].set_visible(False)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_file_directory', nargs='?')
    parser.add_argument('--rig', nargs='?', default='Bruker')
    parser.add_argument('--show_figs', nargs='?', default='False')
    parser.add_argument('--save_figs', nargs='?', default='True')
    parser.add_argument('--tag', nargs='?', default='raw')
    parser.add_argument('--experiment_file_name', nargs='?', default=None)
    parser.add_argument('--output_directory', nargs='?', default=None)
    parser.add_argument('--channels', nargs='?', default='1:calcium,2:cAMP',
                        help='functional channel number : name. The FIRST is the primary channel: '
                             'the gate runs on it and its kept ROIs are used for every channel.')
    parser.add_argument('--filter_length', nargs='?', default='1:3.0,2:5.0', help='s, per channel number')
    parser.add_argument('--lag_dt', nargs='?', default='1:0.02,2:0.1',
                        help='s, per channel number. Finer than the volume is resolvable because '
                             'the planes are sampled at different times.')
    parser.add_argument('--steady_start', nargs='?', default='',
                        help="s after a switch from which rows count as steady state; '' = each "
                             "channel's filter length, so no steady row's lag window crosses a switch")
    parser.add_argument('--first_dwell', nargs='?', type=int, default=2)
    parser.add_argument('--gate', nargs='?', default='high',
                        help="'high' (null test on each ROI's high-contrast filter) or a path to a .npy "
                             'of ROI indices to keep')
    parser.add_argument('--gate_q', nargs='?', type=float, default=0.05)
    parser.add_argument('--null_rolls', nargs='?', type=int, default=14)
    parser.add_argument('--n_boot', nargs='?', type=int, default=1000)
    parser.add_argument('--crossfit', nargs='?', default='True')
    parser.add_argument('--bleach', nargs='?', default='auto')
    parser.add_argument('--export_models', nargs='?', default='False',
                        help="True: write each channel's steady-state LN pair (SwitchedLN on contrast_state) "
                             "to the model library, for visanalysis.model / predict_ln.py")
    parser.add_argument('--model_library', nargs='?', default=None,
                        help='model library directory; default $VISANALYSIS_MODEL_LIBRARY, else '
                             '~/.visanalysis/config.yaml')
    args = parser.parse_args()

    efd = args.experiment_file_directory
    tag = args.tag
    if args.experiment_file_name:
        efn = args.experiment_file_name
    else:
        tagged = 'fly.hdf5' if tag == 'raw' else 'fly_' + tag + '.hdf5'
        efn = tagged if os.path.exists(os.path.join(efd, tagged)) else 'fly.hdf5'
    efp = os.path.join(efd, efn)
    if not os.path.exists(efp):
        print('No such file: ' + efp)
        sys.exit(1)
    print('reading ' + efp)
    if args.rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        plug = base_plugin.BasePlugin()

    chans = _parse_map(args.channels)
    flen = _parse_map(args.filter_length, float)
    ldt = _parse_map(args.lag_dt, float)
    primary = list(chans.values())[0]
    bleach_mode = {'auto': 'auto', 'on': True, 'off': False}[str(args.bleach).lower()]

    series_num = list(map(str, plug.getSeriesNumbers(efp)))
    ID = imaging_data.ImagingDataObject(efp, int(series_num[0]), quiet=True)
    prep = ID.getSubjectMetadata().get('prep', '')
    timing_channel_ind = 1 if prep == 'fly right optic lobe' else 0
    found = []
    for s in series_num:
        plug.updateImagingDataObject(efd, efn, int(s))
        if cs.is_contrast_switch(plug.ImagingDataObject.getRunParameters()):
            found.append(s)
    if not found:
        print('No ContrastSwitchTernaryNoise series found.')
        sys.exit(0)
    print('contrast-switch series: ' + ', '.join(found))

    for s in found:
        sn = 'series_{:03d}'.format(int(s))
        out_dir = (os.path.join(args.output_directory, sn) if args.output_directory
                   else os.path.join(efd, '{}_contrast_switch'.format(tag), sn))
        os.makedirs(out_dir, exist_ok=True)
        t_start = time.time()
        plug.updateImagingDataObject(efd, efn, int(s))
        ID = plug.ImagingDataObject
        # before any response call: epochs are cut on the photodiode
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6
        _, _, sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * sample_rate

        ctx = cs.build_context(ID, {nm: 'mask_ch' + n for n, nm in chans.items()}, bleach_mode)
        print('\n=== {} {}: {} trials, seeds {}{}, {} dropped frames, dwell {:g} s, contrasts {:.3f}/{:.3f}'.format(
            os.path.basename(efd), sn, ctx['n_ep'], sorted(set(ctx['seeds'])), ' (frozen)' if ctx['frozen'] else '',
            sum(ctx['n_drop']), ctx['dwell_s'], *ctx['realized']))
        for nm, dk in ctx['dark'].items():
            if dk:
                print('  WARNING: {} goes dark in trial(s) {} (population level < {:.0%} of typical); those '
                      'trials are dropped for {} only'.format(nm, [e + 1 for e in dk], cs.DARK_FRAC, nm))
        for nm, fd in ctx['faded'].items():
            if fd:
                print('  {}: {} ROI-trials dropped, baseline below {:.0%} of the channel median (ROIs {})'.format(
                    nm, len(fd), cs.FADE_FRAC, sorted({r for r, _ in fd})))
        n_of = {nm: n for n, nm in chans.items()}
        kw = dict(steady_start=float(args.steady_start) if args.steady_start else None)
        res = dict(fly=efd, series=sn, primary=primary, n_trials=ctx['n_ep'], seeds=sorted(set(ctx['seeds'])),
                   frozen=ctx['frozen'], dropped_frames=ctx['n_drop'], dwell_s=ctx['dwell_s'],
                   contrasts=list(ctx['realized']), contrast_ratio=ctx['contrast_ratio'],
                   stim_dur=ctx['stim_dur'].tolist(), dark_trials={k: [e + 1 for e in v] for k, v in ctx['dark'].items()},
                   faded_roi_trials={k: [[r, e + 1] for r, e in v] for k, v in ctx['faded'].items()},
                   channels={}, switch_triggered={},
                   settings=vars(args))
        arrs = {}

        # 1. gate
        R = ctx['resp'][primary]
        if args.gate == 'high':
            gate = cs.gate_rois(ctx, primary, flen[n_of[primary]], ldt[n_of[primary]], n_rolls=args.null_rolls,
                                q=args.gate_q, **kw)
            keep = gate['keep']
        else:
            keep = np.sort(np.load(args.gate).astype(int))
            gate = dict(keep=keep, p=np.full(R['n_roi'], np.nan), r=np.full(R['n_roi'], np.nan),
                        null=np.array([np.nan]), lam=np.nan)
            print('  gate: {} ROIs from {}'.format(keep.size, args.gate))
        np.save(os.path.join(out_dir, 'keep_{}.npy'.format(primary)), keep)
        gate['plane_n'] = np.bincount(R['plane'], minlength=R['frame_offsets'].size)
        gate['plane_k'] = np.bincount(R['plane'][keep], minlength=R['frame_offsets'].size)
        res['gate'] = dict(method=args.gate, q=args.gate_q, n_keep=int(keep.size), n_roi=int(R['n_roi']),
                           ridge=float(gate['lam']), per_plane_kept=gate['plane_k'].tolist(),
                           per_plane_all=gate['plane_n'].tolist())

        # 2-5. LN analysis per channel; 6. switch-triggered averages
        for n, nm in chans.items():
            out, a = cs.analyse_channel(ctx, nm, keep, flen[n], ldt[n], first_dwell=args.first_dwell,
                                        n_rolls=args.null_rolls, n_boot=args.n_boot, **kw)
            res['channels'][nm] = out
            arrs[nm] = a
            sta, sa = cs.switch_triggered(ctx, nm, keep, first_dwell=args.first_dwell)
            res['switch_triggered'][nm] = sta
            arrs['sta_' + nm] = sa
            if args.export_models == 'True':
                from visanalysis.model import export as mexport
                mp = mexport.contrast_switch_model(ctx, res, arrs, nm, efd, sn, library=args.model_library)
                print('  {}: LN model exported to {}'.format(nm, mp))
            print('  {} switch: step to high {:+.4f} ({:+.4f}, {:+.4f}), to low {:+.4f}; tau to high {:.2f} s, '
                  'to low {:.2f} s; peak change {:.4f} at {:.1f} s'.format(
                      nm, sta['to_high']['step'], *sta['to_high']['step_ci'], sta['to_low']['step'],
                      sta['to_high']['tau'], sta['to_low']['tau'], sta['peak_swing']['value'], sta['peak_swing']['at']))
        pc = res['channels'][primary]
        res['model_free'] = cs.model_free_gain(ctx, primary, keep, window=(pc['steady_start'], ctx['dwell_s']),
                                               first_dwell=args.first_dwell)
        if res['model_free']:
            mf = res['model_free']
            ref = pc.get('static_reference', {}).get('truth_linear', np.nan)
            mf['static_reference'] = ref
            mf['gain_ratio_vs_static'] = mf['gain_ratio'] / ref
            mf['gain_ratio_ci_vs_static'] = [x / ref for x in mf['gain_ratio_ci']]
            print('  model-free gain ratio {:.2f} ({:.2f}-{:.2f}); a non-adapting cell reads {:.2f}, so relative '
                  'to it {:.2f} ({:.2f}-{:.2f})'.format(mf['gain_ratio'], *mf['gain_ratio_ci'], ref,
                                                         mf['gain_ratio_vs_static'], *mf['gain_ratio_ci_vs_static']))

        # 7. cross-fit: select on one half of the stimulus, measure on the other
        res['crossfit'] = []
        if args.crossfit == 'True' and args.gate == 'high':
            half = float(np.median(ctx['stim_dur'])) / 2.0
            for sel, mea in (((-np.inf, half), (half, np.inf)), ((half, np.inf), (-np.inf, half))):
                g = cs.gate_rois(ctx, primary, flen[n_of[primary]], ldt[n_of[primary]], n_rolls=args.null_rolls,
                                 q=args.gate_q, pos_range=sel, **kw)
                o, _ = cs.analyse_channel(ctx, primary, g['keep'], flen[n_of[primary]], ldt[n_of[primary]],
                                          first_dwell=args.first_dwell, n_rolls=0, n_boot=args.n_boot,
                                          ridge=res['channels'][primary]['ridge'], pos_range=mea, **kw)
                res['crossfit'].append(dict(select=list(sel), measure=list(mea), n_keep=int(g['keep'].size),
                                            gain=o['gain'], shared_range=o['shared_range']))

        with open(os.path.join(out_dir, 'contrast_switch.json'), 'w') as f:
            json.dump(res, f, indent=1, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else float(v))
        flat = {'{}__{}'.format(k, kk): vv for k, v in arrs.items() for kk, vv in v.items()
                if isinstance(vv, np.ndarray)}
        np.savez(os.path.join(out_dir, 'contrast_switch.npz'), gate_r=gate['r'], gate_p=gate['p'],
                 gate_null=gate['null'], keep=keep, **flat)
        if args.save_figs == 'True':
            title = '{} {} {}'.format(os.path.basename(os.path.dirname(os.path.normpath(efd))),
                                      os.path.basename(os.path.normpath(efd)), sn)
            for name, fig in (('contrast_switch', _figure_main(res, arrs, title)),
                              ('contrast_switch_gate', _figure_gate(gate, title + ': ROI gate'))):
                fig.savefig(os.path.join(out_dir, name + '.png'), dpi=140)
                fig.savefig(os.path.join(out_dir, name + '.pdf'))
                if args.show_figs != 'True':
                    plt.close(fig)
        print('  saved to {} ({:.0f} s)'.format(out_dir, time.time() - t_start))
    if args.show_figs == 'True':
        plt.show()
