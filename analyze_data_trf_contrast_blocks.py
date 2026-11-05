"""
analyze_data_trf_contrast_blocks.py

Per-ROI TEMPORAL filters under switching contrast, from the
ContrastBlockUniformNoise protocol (labpack JCS_protocol): full-field uniform
noise whose range alternates between a low- and a high-contrast distribution
every ~5 s around one mean.

Call structure (same as analyze_data_trf_fullfield.py):
    python analyze_data_trf_contrast_blocks.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --show_figs False \
        --save_figs True \
        --tag final \
        --filter_length 2.0

Results saved to fly.hdf5 (or fly_<tag>.hdf5) under group
/TRF_CONTRAST_BLOCKS_{tag}, unless --output_directory is given, in which case
they go to <output_directory>/contrast_blocks_{tag}.hdf5 and the figures beside
it. Every block of the analysis lives in visanalysis/util/contrast_blocks.py;
read its module docstring for why the series is treated as one continuous
stimulus rather than cut into epochs.

WHAT IT MEASURES, IN THE ORDER TO READ IT
-----------------------------------------
1. cv_r of the steady-state filter in each contrast, with a circular-shift
   null. A filter that does not predict held-out data is not a filter, and
   nothing below means anything for an ROI that fails it.
2. Steady-state LN model per contrast, in the Baccus & Meister (2002)
   convention: each filter is scaled so var(filtered stimulus) = var(stimulus),
   so the FILTER carries only temporal processing (shape, kinetics) and the
   NONLINEARITY carries all the gain. gain = slope of the nonlinearity, dF/F
   per unit contrast. gain_ratio = gain(low) / gain(high); compare it with
   contrast_ratio = SD(high) / SD(low). normalization_index =
   log(gain_ratio) / log(contrast_ratio): 0 = no gain control, 1 = the response
   amplitude is independent of contrast.
3. The same nonlinearities over the generator range BOTH conditions cover. A
   static compressive nonlinearity would give them the same slope there;
   shared_gain_ratio tests exactly that.
4. Gain vs time since the switch, relative to steady state. A static
   nonlinearity would be adapted in the first bin; adaptation is not. This and
   (3) are the evidence that the gain change is adaptation.
5. Switch-triggered average: the contrast-envelope response. For a channel with
   no linear filter at this timescale (cAMP in Mi1, so far) it is the readout.
6. For a frozen-noise series (every epoch the same seed): the noise ceiling of
   the cycle-averaged response, and how much of it the steady-state filters,
   and the filters plus the gain dynamics, predict. Filters come from the fly's
   random-seed series, so this is an out-of-series test.
"""

import sys
import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import h5py

from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.util import nullcal
from visanalysis.util import contrast_blocks as cb
from visanalysis.util.fullfield_temporal import integration_index

try:
    from visanalysis.util.strf_plot import apply_figure_style
    apply_figure_style()
except Exception:
    pass
plt.ioff()

COND_COLORS = ['#3a7bd5', '#d0473c', '#6a9f3a', '#9b59b6']


def _peak(h):
    i = int(np.nanargmax(np.abs(h))) if np.any(np.isfinite(h)) else 0
    return i, h[i]


def _kinetics(h, dt):
    """Peak lag, signed peak, integration index, and time-to-half-decay."""
    i, pk = _peak(h)
    after = np.flatnonzero(np.abs(h[i:]) < 0.5 * abs(pk))
    t_half = (after[0] + i) * dt if after.size else np.nan
    return dict(peak_lag_s=i * dt, peak=pk, integration_index=integration_index(h),
                half_decay_s=t_half)


def _savefig(fig, path):
    """PDF for print, PNG beside it for viewing anywhere a PDF will not open."""
    fig.savefig(path)
    fig.savefig(os.path.splitext(path)[0] + '.png', dpi=130)


def _fmt(a):
    return '[' + ' '.join('{:.2f}'.format(x) for x in np.atleast_1d(a)) + ']'


def envelope_effect(lag, trials):
    """Contrast-envelope response per ROI from switch-aligned trials.

    effect = mean over the high block minus mean over the preceding low block,
    per switch, then mean and t across switches -- switches, not samples, are
    the independent unit. peak_s is the time of the extreme of the mean trace
    after the switch, and `late` its level over the last second, which together
    separate a sustained step (calcium in Mi1: the rectified mean rises and
    stays) from a transient (up and back).
    """
    if not trials.size:
        return None
    post = lag >= 0
    pre = lag < 0
    with np.errstate(invalid='ignore'):
        per_sw = (np.nanmean(trials[:, :, post], 2)
                  - np.nanmean(trials[:, :, pre], 2))
        m = np.nanmean(per_sw, 0)
        se = np.nanstd(per_sw, 0, ddof=1) / np.sqrt(per_sw.shape[0])
        mean_tr = np.nanmean(trials, 0)
    base = np.nanmean(mean_tr[:, pre], 1, keepdims=True)
    dev = mean_tr[:, post] - base
    pk = np.nanargmax(np.abs(dev), axis=1)
    late = np.nanmean(mean_tr[:, lag >= lag[-1] - 1.0], 1) - base[:, 0]
    return dict(effect=m, effect_se=se, t=m / se,
                peak_s=lag[post][pk], peak=dev[np.arange(len(pk)), pk],
                late=late)


def _print_envelope(env, n_sw, n_roi, label):
    print('    ENVELOPE ({}, {} switches, median over ROIs): high-minus-low '
          '{:+.4f} dF/F (t {:.1f}; {} of {} ROIs |t|>3); peak {:+.4f} at '
          '{:.2f} s after the switch; last 1 s {:+.4f} (close to peak = '
          'sustained, near 0 = transient)'.format(
              label, n_sw, np.median(env['effect']), np.median(env['t']),
              int(np.sum(np.abs(env['t']) > 3)), n_roi,
              np.median(env['peak']), np.median(env['peak_s']),
              np.median(env['late'])))


def _timing_kwargs(ID, args):
    """Lag grid and per-sample integration window for build_design.

    --sample_window 'frame' (default) is one plane's frame period from the
    acquisition metadata: a plane is scanned in that time, so that is how long
    each sample integrated the stimulus. 'volume' reproduces the old binned
    design. --lag_dt 0 means the volume period.
    """
    lag_dt = float(args.lag_dt) or None
    sw = str(args.sample_window).lower()
    if sw == 'volume':
        window = None
    elif sw == 'frame':
        window = float(ID.getAcquisitionMetadata('framePeriod'))
    else:
        window = float(sw)
    return dict(lag_dt=lag_dt, sample_window_s=window,
                penalty=str(args.penalty), lag_basis=str(args.lag_basis))


def analyze_random_series(ID, ch_name, args):
    """Everything computed from a random-seed series. Returns a result dict."""
    fl = float(args.filter_length)
    d = cb.build_design(ID, ch_name, fl, correct_bleach=args.bleach != 'off',
                        block_s=float(args.block_s), **_timing_kwargs(ID, args))
    n_cond = len(d['cond_keys'])
    ss_start = fl if args.steady_start is None else float(args.steady_start)
    ss = d['since_switch'] >= ss_start
    masks = [ss & (d['cond'] == q) for q in range(n_cond)]
    oi = d['onset_info']
    print('    {} rois, {} rows, {} lags of {:.1f} ms ({:.2f} s), volume period '
          '{:.1f} ms, sample window {:.1f} ms, {} z plane(s)'.format(
              d['n_roi'], d['Y'].shape[0], d['n_lags'], d['lag_dt'] * 1e3,
              d['filter_length_s'], d['dt'] * 1e3, d['window'] * 1e3,
              len(d['z_offsets'])))
    print('    onsets: {} detected for {} epochs; placed from neighbours: {}; '
          'clock residual (min/med/max) {} ms; dropped frames {}'.format(
              oi['n_detected'], oi['n_epochs'], list(oi['placed_epochs']),
              _fmt(oi['clock_resid_ms']), oi['n_dropped_frames']))
    if len(oi['placed_epochs']) > 0.05 * oi['n_epochs']:
        print('    WARNING: {} of {} epoch onsets had to be placed rather than '
              'detected. That is what the WRONG photodiode channel looks like '
              '(check the prep metadata); the stimulus timing below is then a '
              'unix-clock estimate, not a measurement.'.format(
                  len(oi['placed_epochs']), oi['n_epochs']))
    for q in range(n_cond):
        print('    condition {}: levels {} -> RMS contrast {:.3f} (displayed SD '
              '{:.3f}); steady-state rows (since switch >= {:.2f} s): {}'.format(
                  q, d['cond_keys'][q], d['rms_contrast'][q], d['stim_sd'][q],
                  ss_start, int(masks[q].sum())))

    # ── ridge, filters, cv_r ────────────────────────────────────────────────
    grid = None if args.ridge is None else np.array([args.ridge])
    if grid is None and d['penalty_kind'] == 'smooth':
        grid = np.logspace(0, 5, 11)
    ridge, grid, rscores = cb.choose_ridge(d, masks, grid=grid,
                                           n_fold=int(args.n_fold))
    h = np.array([cb.fit_filters(d, m, ridge)[0] for m in masks])
    b = np.array([cb.fit_filters(d, m, ridge)[1] for m in masks])
    cvr = np.array([cb.cv_r(d, m, ridge, int(args.n_fold)) for m in masks])
    print('    {} penalty {:.3g}; median cv_r by condition {}'.format(
        d['penalty_kind'], ridge, _fmt(np.nanmedian(cvr, axis=1))))

    # ── null: cycle-aligned shifts, prediction gate per condition ──────────
    keep_q = [np.arange(d['n_roi']) for _ in range(n_cond)]
    cv_p = np.full(cvr.shape, np.nan)
    null_thr = np.full(n_cond, np.nan)
    if int(args.null_shifts) > 0:
        _t = time.perf_counter()
        shifts = cb.cycle_shifts(d, int(args.null_shifts))
        for q in range(n_cond):
            nul = np.concatenate([cb.cv_r(d, masks[q], ridge, int(args.n_fold),
                                          shift=s) for s in shifts])
            nul = nul[np.isfinite(nul)]
            ok = np.isfinite(cvr[q])
            null_thr[q] = np.percentile(nul, 95)
            cv_p[q, ok] = nullcal.pvalues(cvr[q, ok], nul)[0]
            resolvable = int(np.ceil(1.0 / (0.05 * len(shifts))))
            if resolvable <= int(ok.sum()):
                sel, _ = nullcal.benjamini_hochberg(cv_p[q, ok], q=0.05)
                keep_q[q] = np.flatnonzero(ok)[sel] if len(sel) else np.array([], int)
            else:
                keep_q[q] = np.flatnonzero(ok & (cvr[q] > null_thr[q]))
        print('    null ({} cycle-aligned shifts, {:.1f} s): 95th pct cv_r {}; '
              'passing per condition {}'.format(
                  len(shifts), time.perf_counter() - _t, _fmt(null_thr),
                  [len(k) for k in keep_q]))
    keep = keep_q[0]
    for k in keep_q[1:]:
        keep = np.intersect1d(keep, k)
    print('    KEPT (prediction gate in every condition): {} of {} -> {}'.format(
        len(keep), d['n_roi'], list(keep)))
    sel = keep if len(keep) else np.arange(d['n_roi'])
    if not len(keep):
        print('    WARNING: no ROI passed in every condition; population numbers '
              'below are over ALL ROIs and are not a selected population.')

    # ── LN filters (Baccus & Meister normalisation) and kinetics ────────────
    # Scaled so var(F * s) = var(s): the filter keeps its shape, its gain moves
    # to the nonlinearity. amp is the RAW filter peak, kept for reference.
    h_ln = np.array([cb.normalize_filters(d, h[q], masks[q])[0]
                     for q in range(n_cond)])
    kin = [[_kinetics(h_ln[q, r], d['lag_dt']) for r in range(d['n_roi'])]
           for q in range(n_cond)]
    amp = np.array([[abs(_peak(h[q, r])[1]) for r in range(d['n_roi'])]
                    for q in range(n_cond)])
    lo, hi = 0, n_cond - 1
    contrast_ratio = d['stim_sd'][hi] / d['stim_sd'][lo]

    # ── gain after switch ───────────────────────────────────────────────────
    edges = np.asarray(args.gain_bins, dtype=float)
    gain, nrow, folds = cb.gain_after_switch(d, ss, ridge, edges,
                                             int(args.n_fold))
    # gain is ABSOLUTE here -- the slope of the nonlinearity, dF/F per unit
    # contrast -- because the folds are B&M-normalised. gain_rel divides by
    # the mean over bins lying wholly in the steady-state window:
    # out-of-fold slopes are not exactly 1 there (ridge shrinkage pushes them
    # up, estimation noise pulls them down), so "1" is set by the data.
    ss_bins = edges[:-1] >= ss_start - 1e-9
    ref = np.nanmean(gain[:, ss_bins, :], axis=1, keepdims=True)
    gain_rel = gain / ref
    centers = 0.5 * (edges[:-1] + edges[1:])
    exp_fit = np.array([cb.fit_exp_recovery(centers,
                                            np.nanmedian(gain_rel[q][:, sel], 1))
                        for q in range(n_cond)])

    # ── steady-state nonlinearities: where the gain lives ──────────────────
    nl = cb.ln_model(d, ss, folds)
    gain_ratio = nl['gain'][lo] / nl['gain'][hi]
    norm_index = np.log(gain_ratio) / np.log(contrast_ratio)
    shared_ratio = nl['shared_gain'][lo] / nl['shared_gain'][hi]
    full_ratio = gain_ratio

    # ── early vs late filters ──────────────────────────────────────────────
    windows = [(0.0, 0.5), (0.5, 1.0), (ss_start, 1e9)]
    wf = cb.window_filters(d, ridge, windows)

    # ── switch-triggered averages ──────────────────────────────────────────
    per = float(np.median(np.diff(d['onsets'])))
    sta_lag, sta_tr = cb.switch_triggered_average(
        d, hi, per, per, bin_s=float(args.sta_bin) or None)
    _, sta_tr_lo = cb.switch_triggered_average(
        d, lo, per, per, bin_s=float(args.sta_bin) or None)
    env = envelope_effect(sta_lag, sta_tr)
    if env is not None:
        _print_envelope(env, sta_tr.shape[0], d['n_roi'], 'random seeds')
    has_filter = bool(len(keep))
    if not has_filter:
        print('    no ROI has a filter that predicts held-out data in every '
              'condition, so gain, kinetics, nonlinearity and the gain time '
              'course are NOT reported for this channel: they would be '
              'descriptions of noise. The envelope does not depend on a '
              'filter and stands on its own.')
    else:
        print('    steady state (kept ROIs, median): peak lag {} s; |peak| {}; '
              'half-decay {} s'.format(
                  _fmt([np.nanmedian([kin[q][r]['peak_lag_s'] for r in sel])
                        for q in range(n_cond)]),
                  _fmt([np.nanmedian(amp[q, sel]) for q in range(n_cond)]),
                  _fmt([np.nanmedian([kin[q][r]['half_decay_s'] for r in sel])
                        for q in range(n_cond)])))
        # Sensitivity as Baccus & Meister define it (their Fig. 9C): the mean
        # slope dN/dg over the g values common to both conditions.
        print('    LN sensitivity (B&M: slope of N over the shared g range, '
              'dF/F per unit contrast): {}; whole-range slope {}'.format(
                  _fmt([np.nanmedian(nl['shared_gain'][q, sel])
                        for q in range(n_cond)]),
                  _fmt([np.nanmedian(nl['gain'][q, sel])
                        for q in range(n_cond)])))
        print('    gain ratio low/high {:.2f} (contrast ratio {:.2f}) -> '
              'normalization index {:.2f}'.format(
                  np.nanmedian(gain_ratio[sel]), contrast_ratio,
                  np.nanmedian(norm_index[sel])))
        print('    nonlinearity slope ratio low/high: whole range {:.2f}, SHARED '
              'input range {:.2f}  (static nonlinearity predicts 1 in the shared '
              'range)'.format(np.nanmedian(full_ratio[sel]),
                              np.nanmedian(shared_ratio[sel])))
        for q in range(n_cond):
            print('    gain after switch into cond {} (rel. steady state, median of '
                  'kept): {}  exp fit A={:+.2f} tau={:.2f} s'.format(
                      q, _fmt(np.nanmedian(gain_rel[q][:, sel], 1)),
                      exp_fit[q, 0], exp_fit[q, 1]))
        print('      bins (s): {}'.format(_fmt(edges)))

    return dict(design=d, ridge=ridge, ridge_grid=grid, ridge_scores=rscores,
                filters=h, intercepts=b, cv_r=cvr, cv_p=cv_p,
                cv_null_threshold=null_thr, keep=keep,
                keep_by_condition=keep_q, sel=sel, kinetics=kin, amp=amp,
                ln_filters=h_ln,
                gain_ratio=gain_ratio, contrast_ratio=contrast_ratio,
                normalization_index=norm_index, gain_bins=edges,
                gain=gain, gain_rel=gain_rel, gain_rows=nrow,
                exp_fit=exp_fit, nl=nl, shared_gain_ratio=shared_ratio,
                full_gain_ratio=full_ratio, windows=windows,
                window_filters=wf, sta_lag=sta_lag, sta_trials=sta_tr,
                sta_trials_to_low=sta_tr_lo,
                ss_start=ss_start, envelope=env, has_filter=has_filter)


def analyze_repeat_series(ID, ch_name, args, ref):
    """Noise ceiling and out-of-series prediction for a frozen series."""
    d = cb.build_design(ID, ch_name, float(args.filter_length),
                        correct_bleach=args.bleach != 'off',
                        block_s=float(args.block_s), **_timing_kwargs(ID, args))
    n_cond = len(d['cond_keys'])
    oi = d['onset_info']
    print('    frozen noise: {} epochs, onsets detected {}; placed {}'.format(
        oi['n_epochs'], oi['n_detected'], list(oi['placed_epochs'])))
    per = float(np.median(np.diff(d['onsets'])))
    sta_lag, sta_tr = cb.switch_triggered_average(
        d, n_cond - 1, per, per, bin_s=float(args.sta_bin) or None)
    out = dict(design=d, sta_lag=sta_lag, sta_trials=sta_tr)
    if ref is None:
        print('    no random-seed series in this fly to take filters from; '
              'reporting reliability only')
    if ref is not None and not ref.get('has_filter', True):
        print('    the random series has no validated filter on this channel; '
              'prediction skipped, reliability only')
        ref = None
    if ref is not None and ref['design']['n_roi'] != d['n_roi']:
        print('    ROI count differs from the random series ({} vs {}); '
              'prediction skipped'.format(d['n_roi'], ref['design']['n_roi']))
        ref = None
    tc, ymean, ceiling = cb.cycle_average(d, d['Y'], n_cond)
    out.update(t_cycle=tc, y_cycle=ymean, noise_ceiling=ceiling)
    if ref is not None:
        p_ss = cb.predict_series(d, ref['filters'], ref['intercepts'])
        p_dyn = cb.predict_series(d, ref['filters'], ref['intercepts'],
                                  gain=ref['gain_rel'], bin_edges=ref['gain_bins'])
        _, pm_ss, _ = cb.cycle_average(d, p_ss, n_cond)
        _, pm_dyn, _ = cb.cycle_average(d, p_dyn, n_cond)
        r_ss = cb._colcorr(pm_ss, ymean)
        r_dyn = cb._colcorr(pm_dyn, ymean)
        out.update(pred_cycle_ss=pm_ss, pred_cycle_dyn=pm_dyn, r_ss=r_ss,
                   r_dyn=r_dyn)
        sel = ref['sel']
        print('    cycle-average (kept ROIs, median): noise ceiling {:.3f}; '
              'r steady-state filters {:.3f}; r with gain dynamics {:.3f}'
              .format(np.nanmedian(ceiling[sel]), np.nanmedian(r_ss[sel]),
                      np.nanmedian(r_dyn[sel])))
    else:
        print('    cycle-average noise ceiling (median, all ROIs) {:.3f}'.format(
            np.nanmedian(ceiling)))
    env = envelope_effect(sta_lag, sta_tr)
    out['envelope'] = env
    if env is not None:
        _print_envelope(env, sta_tr.shape[0], d['n_roi'], 'frozen noise')
    return out


# ─────────────────────────────────────────────────────────────────────────────
# figures
# ─────────────────────────────────────────────────────────────────────────────

def _cond_label(res, q):
    return 'RMS c = {:.2f}'.format(res['design']['rms_contrast'][q])


def figure_summary(res, title, signal=''):
    d = res['design']
    n_cond = len(d['cond_keys'])
    sel = res['sel']
    lags = d['lags']
    fig, ax = plt.subplots(2, 3, figsize=(15, 8.5))

    a = ax[0, 0]
    for q in range(n_cond):
        m = np.nanmean(res['ln_filters'][q, sel], 0)
        sd = np.nanstd(res['ln_filters'][q, sel], 0) / np.sqrt(max(len(sel), 1))
        a.fill_between(lags, m - sd, m + sd, color=COND_COLORS[q], alpha=0.25,
                       lw=0)
        a.plot(lags, m, color=COND_COLORS[q], label=_cond_label(res, q))
    a.axhline(0, color='k', lw=0.5)
    a.set(xlabel='lag (s)', ylabel='filter (var(F*s) = var(s))',
          title='LN filters: temporal processing ({} ROIs)'.format(len(sel)))
    a.legend(frameon=False)

    a = ax[0, 1]
    nl = res['nl']
    for q in range(n_cond):
        x = np.nanmean(nl['x'][q][:, sel], 1)
        y = np.nanmean(nl['y'][q][:, sel], 1)
        a.plot(x, y, 'o-', color=COND_COLORS[q], ms=4,
               label='gain {:.3f}'.format(np.nanmedian(nl['gain'][q, sel])))
    lo, hi = np.nanmean(nl['shared_range'][:, sel], 1)
    a.axvspan(lo, hi, color='0.9', zorder=0)
    a.set(xlabel='generator (contrast units)', ylabel='dF/F',
          title='LN nonlinearities: gain; shared ratio {:.2f}'.format(
              np.nanmedian(res['shared_gain_ratio'][sel])))
    a.legend(frameon=False, fontsize=8)

    a = ax[0, 2]
    cen = 0.5 * (res['gain_bins'][:-1] + res['gain_bins'][1:])
    for q in range(n_cond):
        a.plot(cen, np.nanmedian(res['gain'][q][:, sel], 1), 'o-',
               color=COND_COLORS[q], lw=2, label='into ' + _cond_label(res, q))
        a.axhline(np.nanmedian(nl['gain'][q, sel]), color=COND_COLORS[q],
                  ls='--', lw=0.8)
    a.axvline(res['ss_start'], color='0.6', ls=':', lw=1)
    a.set(xlabel='time since switch (s)', ylabel='gain (dF/F per unit contrast)',
          title='gain after switch (dashed = steady state)')
    a.set_ylim(bottom=0)
    a.legend(frameon=False, fontsize=8)

    a = ax[1, 0]
    cen = 0.5 * (res['gain_bins'][:-1] + res['gain_bins'][1:])
    for q in range(n_cond):
        g = res['gain_rel'][q][:, sel]
        a.plot(cen, g, color=COND_COLORS[q], alpha=0.15, lw=0.7)
        a.plot(cen, np.nanmedian(g, 1), 'o-', color=COND_COLORS[q], lw=2,
               label='into {}; tau {:.2f} s'.format(_cond_label(res, q),
                                                    res['exp_fit'][q, 1]))
        A, tau = res['exp_fit'][q]
        tt = np.linspace(0, cen[-1], 200)
        if np.isfinite(tau):
            a.plot(tt, 1 + A * np.exp(-tt / tau), '--', color=COND_COLORS[q],
                   lw=1)
    a.axhline(1, color='k', lw=0.5)
    a.axvline(res['ss_start'], color='0.6', ls=':', lw=1)
    a.set(xlabel='time since switch (s)', ylabel='gain / steady-state gain',
          title='gain after switch (out-of-fold filters)')
    a.set_ylim(0, max(2.5, np.nanpercentile(res['gain_rel'][:, :, sel], 95)))
    a.legend(frameon=False, fontsize=8)

    a = ax[1, 1]
    wl = ['{:.1f}-{:.1f} s'.format(*w) if w[1] < 1e8 else 'steady'
          for w in res['windows']]
    ls = ['-', '--', ':']
    for q in range(n_cond):
        for w in range(len(res['windows'])):
            a.plot(lags, np.nanmean(res['window_filters'][q, w][sel], 0),
                   color=COND_COLORS[q], ls=ls[w % 3],
                   label='{} {}'.format(_cond_label(res, q), wl[w]))
    a.axhline(0, color='k', lw=0.5)
    a.set(xlabel='lag (s)', ylabel='dF/F per unit contrast',
          title='filters by time since switch')
    a.legend(frameon=False, fontsize=7)

    a = ax[1, 2]
    tr = res['sta_trials']
    if tr.size:
        m = np.nanmean(tr[:, sel, :], axis=(0, 1))
        s = np.nanstd(np.nanmean(tr[:, sel, :], 1), 0) / np.sqrt(tr.shape[0])
        a.axvspan(0, res['sta_lag'][-1], color=COND_COLORS[n_cond - 1],
                  alpha=0.08, lw=0)
        a.axvspan(res['sta_lag'][0], 0, color=COND_COLORS[0], alpha=0.08, lw=0)
        a.fill_between(res['sta_lag'], m - s, m + s, color='k', alpha=0.2, lw=0)
        a.plot(res['sta_lag'], m, 'k')
        a.axhline(0, color='k', lw=0.4)
    a.set(xlabel='time from switch to high contrast (s)', ylabel='dF/F',
          title='{}switch-triggered avg (n={})'.format(
              signal + ': ' if signal else '', tr.shape[0]))
    if not res.get('has_filter', True):
        for a in list(ax[0]) + [ax[1, 0], ax[1, 1]]:
            a.set_facecolor('0.93')
            a.text(0.5, 0.5, 'no ROI passes the\nprediction gate:\n'
                   'not interpretable', transform=a.transAxes, ha='center',
                   va='center', color='0.35', fontsize=11)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def figure_switch(entries, title):
    """Switch-triggered averages, every channel side by side, both directions.

    entries: list of (signal label, result dict) for ONE series. Each panel is
    the mean over that channel's ROIs (kept ROIs where a filter passed the
    gate, all ROIs otherwise), baseline = mean of the block before the switch,
    band = SEM across switches.
    """
    fig, ax = plt.subplots(1, len(entries), figsize=(6.5 * len(entries), 4.6),
                           squeeze=False)
    for a, (label, res) in zip(ax[0], entries):
        sel = res['sel']
        lag = res['sta_lag']
        for key, col, lab in (('sta_trials', COND_COLORS[1], 'low \u2192 high'),
                              ('sta_trials_to_low', COND_COLORS[0],
                               'high \u2192 low')):
            tr = res.get(key)
            if tr is None or not tr.size:
                continue
            with np.errstate(invalid='ignore'):
                per_sw = np.nanmean(tr[:, sel, :], 1)       # (n_switch, n_t)
                m = np.nanmean(per_sw, 0)
                n = np.sum(np.isfinite(per_sw), 0)
                se = np.nanstd(per_sw, 0) / np.sqrt(np.maximum(n, 1))
            base = np.nanmean(m[lag < 0])
            a.fill_between(lag, m - se - base, m + se - base, color=col,
                           alpha=0.25, lw=0)
            a.plot(lag, m - base, color=col, lw=1.2,
                   label='{} (n={})'.format(lab, tr.shape[0]))
        a.axvspan(lag[0], 0, color='0.93', lw=0, zorder=0)
        a.axhline(0, color='k', lw=0.5)
        a.axvline(0, color='k', lw=0.5, ls=':')
        a.set(xlabel='time from switch (s)',
              ylabel='dF/F re. pre-switch block mean',
              title='{} ({} ROIs{})'.format(
                  label, len(sel), '' if res.get('has_filter') else ', all'))
        a.legend(frameon=False, fontsize=8, loc='upper left')
    fig.suptitle(title)
    fig.tight_layout()
    return fig


def figure_rois(res, title):
    d = res['design']
    n_cond = len(d['cond_keys'])
    lags = d['lags']
    n = d['n_roi']
    nc = 6
    nr = int(np.ceil(n / nc))
    fig, ax = plt.subplots(nr, nc, figsize=(2.6 * nc, 2.0 * nr), squeeze=False)
    for r in range(nr * nc):
        a = ax.flat[r]
        if r >= n:
            a.axis('off')
            continue
        for q in range(n_cond):
            a.plot(lags, res['filters'][q, r], color=COND_COLORS[q], lw=1)
        a.axhline(0, color='k', lw=0.4)
        a.set_title('roi {}{}  r={}'.format(
            r, '*' if r in res['keep'] else '',
            '/'.join('{:.2f}'.format(x) for x in res['cv_r'][:, r])), fontsize=8)
        a.tick_params(labelsize=7)
    fig.suptitle(title + '   (* = passes prediction gate in every condition)')
    fig.tight_layout()
    return fig


def figure_repeat(rep, ref, title):
    fig, ax = plt.subplots(1, 2, figsize=(13, 4))
    sel = ref['sel'] if ref is not None else np.arange(rep['design']['n_roi'])
    a = ax[0]
    a.plot(rep['t_cycle'], np.nanmean(rep['y_cycle'][:, sel], 1), 'k',
           label='data (cycle mean)')
    if 'pred_cycle_ss' in rep:
        a.plot(rep['t_cycle'], np.nanmean(rep['pred_cycle_ss'][:, sel], 1),
               color='C1', label='steady-state filters')
        a.plot(rep['t_cycle'], np.nanmean(rep['pred_cycle_dyn'][:, sel], 1),
               color='C2', label='+ gain dynamics')
    a.set(xlabel='time in cycle (s)', ylabel='dF/F',
          title='frozen cycle, mean of {} ROIs'.format(len(sel)))
    a.legend(frameon=False, fontsize=8)
    a = ax[1]
    x = np.arange(len(sel))
    a.bar(x - 0.27, rep['noise_ceiling'][sel], 0.27, color='0.7',
          label='noise ceiling')
    if 'r_ss' in rep:
        a.bar(x, rep['r_ss'][sel], 0.27, color='C1', label='r steady-state')
        a.bar(x + 0.27, rep['r_dyn'][sel], 0.27, color='C2', label='r + dynamics')
    a.set_xticks(x)
    a.set_xticklabels([str(s) for s in sel], fontsize=7)
    a.set(xlabel='roi', ylabel='r with cycle mean')
    a.legend(frameon=False, fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# save
# ─────────────────────────────────────────────────────────────────────────────

def save_results(path, group_name, results, repeats, args):
    with h5py.File(path, 'a') as f:
        if group_name in f:
            del f[group_name]
        g = f.create_group(group_name)
        g.attrs['filter_length_s'] = float(args.filter_length)
        g.attrs['block_s'] = float(args.block_s)
        g.attrs['n_fold'] = int(args.n_fold)
        g.attrs['lag_dt'] = float(args.lag_dt)
        g.attrs['sample_window'] = str(args.sample_window)
        g.attrs['sta_bin'] = float(args.sta_bin)
        g.attrs['penalty'] = str(args.penalty)
        g.attrs['null_shifts'] = int(args.null_shifts)
        for (sn, ch), res in results.items():
            d = res['design']
            sg = g.require_group('{}/{}'.format(sn, ch))
            sg.attrs['ridge'] = res['ridge']
            sg.attrs['has_filter'] = bool(res['has_filter'])
            sg.attrs['steady_start_s'] = res['ss_start']
            sg.attrs['contrast_ratio'] = res['contrast_ratio']
            sg.attrs['placed_epochs'] = np.asarray(
                d['onset_info']['placed_epochs'], dtype=int)
            sg.attrs['condition_keys'] = str(d['cond_keys'])
            for k, v in dict(
                    filters=res['filters'], intercepts=res['intercepts'],
                    lags=d['lags'],
                    rms_contrast=d['rms_contrast'], stim_sd=d['stim_sd'],
                    cv_r=res['cv_r'], cv_p=res['cv_p'],
                    cv_null_threshold=res['cv_null_threshold'],
                    kept_rois=np.asarray(res['keep'], dtype=int),
                    amplitude=res['amp'], gain_ratio=res['gain_ratio'],
                    normalization_index=res['normalization_index'],
                    peak_lag_s=np.array([[k['peak_lag_s'] for k in kq]
                                         for kq in res['kinetics']]),
                    half_decay_s=np.array([[k['half_decay_s'] for k in kq]
                                           for kq in res['kinetics']]),
                    integration_index=np.array(
                        [[k['integration_index'] for k in kq]
                         for kq in res['kinetics']]),
                    gain_bins=res['gain_bins'], gain_after_switch=res['gain'],
                    gain_after_switch_rel=res['gain_rel'],
                    gain_rows=res['gain_rows'], gain_exp_fit=res['exp_fit'],
                    ln_filters=res['ln_filters'],
                    nl_x=res['nl']['x'], nl_y=res['nl']['y'],
                    ln_gain=res['nl']['gain'],
                    ln_shared_gain=res['nl']['shared_gain'],
                    ln_generator_sd=res['nl']['g_sd'],
                    nl_shared_range=res['nl']['shared_range'],
                    shared_gain_ratio=res['shared_gain_ratio'],
                    window_edges=np.array(res['windows']),
                    window_filters=res['window_filters'],
                    sta_lag=res['sta_lag'], sta_mean=np.nanmean(
                        res['sta_trials'], 0) if res['sta_trials'].size
                    else np.array([]),
                    envelope_effect=(res['envelope']['effect']
                                     if res['envelope'] else []),
                    envelope_t=res['envelope']['t'] if res['envelope'] else [],
                    envelope_peak_s=(res['envelope']['peak_s']
                                     if res['envelope'] else []),
                    ridge_grid=res['ridge_grid'],
                    ridge_cv_r=res['ridge_scores'],
                    roi_z_offset_s=d['z_offsets'][d['roi_group']],
                    epoch_onsets=d['onsets']).items():
                sg.create_dataset(k, data=np.asarray(v, dtype=float)
                                  if k != 'kept_rois' else v)
        for (sn, ch), rep in repeats.items():
            sg = g.require_group('{}/{}'.format(sn, ch))
            sg.attrs['frozen_noise'] = True
            for k in ('t_cycle', 'y_cycle', 'noise_ceiling', 'pred_cycle_ss',
                      'pred_cycle_dyn', 'r_ss', 'r_dyn', 'sta_lag'):
                if k in rep:
                    sg.create_dataset(k, data=np.asarray(rep[k], dtype=float))
            if rep['sta_trials'].size:
                sg.create_dataset('sta_mean', data=np.nanmean(rep['sta_trials'], 0))


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_file_directory', nargs='?')
    parser.add_argument('--rig', nargs='?', default='Bruker')
    parser.add_argument('--show_figs', nargs='?', default='False')
    parser.add_argument('--save_figs', nargs='?', default='False')
    parser.add_argument('--tag', nargs='?', default='raw')
    parser.add_argument('--experiment_file_name', nargs='?', default=None,
                        help='hdf5 to read. Default: fly_<tag>.hdf5 when it '
                             'exists, else fly.hdf5.')
    parser.add_argument('--output_directory', nargs='?', default=None,
                        help='write results to <dir>/contrast_blocks_<tag>.hdf5 '
                             'and figures to <dir>/<tag>_roi_figs instead of '
                             'into the fly folder')
    parser.add_argument('--channels', nargs='?', default='1,2')
    parser.add_argument('--filter_length', nargs='?', type=float, default=2.0,
                        help='seconds. Mi1 calcium filters are over by ~0.5 s; '
                             'the cost of a longer one is that the steady-state '
                             'window (since switch >= filter_length, by default) '
                             'shrinks inside a 5 s block.')
    parser.add_argument('--lag_dt', nargs='?', type=float, default=0.01,
                        help='filter lag spacing, s. Finer than the volume '
                             'period is resolvable because samples fall at '
                             'random phases of the 30 Hz stimulus; 0 = volume '
                             'period (the old design).')
    parser.add_argument('--sample_window', nargs='?', default='frame',
                        help="how long each sample integrated the stimulus: "
                             "'frame' (one plane's frame period, default), "
                             "'volume' (old design), or seconds")
    parser.add_argument('--penalty', nargs='?', default='auto',
                        help="filter regulariser: 'auto' (smooth on a sub-volume "
                             "lag grid, ridge otherwise), 'ridge', 'smooth'. "
                             "ridge on a fine grid leaves a 33 ms ripple where "
                             "the 30 Hz stimulus has no power.")
    parser.add_argument('--channel_names', nargs='?', default='',
                        help="labels for figure titles, e.g. "
                             "'1:calcium (GECO),2:cAMP (gFlamp2)'. What each "
                             "channel carries depends on the genotype, so it is "
                             "given here rather than assumed.")
    parser.add_argument('--lag_basis', nargs='?', default='point',
                        help="'point' (taps at each lag; for fine grids) or "
                             "'box' (stimulus averaged over each lag bin; for "
                             "coarse grids such as 0.04-0.05 s)")
    parser.add_argument('--sta_bin', nargs='?', type=float, default=0.02,
                        help='bin for jitter-binned switch-triggered averages, '
                             's; 0 = interpolate onto the volume grid')
    parser.add_argument('--steady_start', nargs='?', type=float, default=None,
                        help='time since switch from which a row counts as '
                             'steady state. Default = filter_length, the '
                             'earliest point at which the whole lag window lies '
                             'in the current block.')
    parser.add_argument('--gain_bins', nargs='?', default='0,0.25,0.5,0.75,1,'
                        '1.5,2,3,4,5.1',
                        help='edges (s since switch) for the gain time course')
    parser.add_argument('--block_s', nargs='?', type=float, default=60.0,
                        help='cross-validation chunk length, s')
    parser.add_argument('--n_fold', nargs='?', type=int, default=5)
    parser.add_argument('--ridge', nargs='?', type=float, default=None)
    parser.add_argument('--bleach', nargs='?', default='on',
                        help="'on' (default) or 'off'")
    parser.add_argument('--null_shifts', nargs='?', type=int, default=20)
    args = parser.parse_args()
    args.gain_bins = [float(x) for x in str(args.gain_bins).split(',')]
    channel_names = dict(kv.split(':', 1) for kv in
                         str(args.channel_names).split(',') if ':' in kv)

    experiment_file_directory = args.experiment_file_directory
    show_figs = args.show_figs == 'True'
    save_figs = args.save_figs == 'True'
    tag = args.tag
    if args.experiment_file_name:
        experiment_file_name = args.experiment_file_name
    else:
        _tagged = 'fly.hdf5' if tag == 'raw' else 'fly_' + tag + '.hdf5'
        experiment_file_name = (
            _tagged if os.path.exists(os.path.join(experiment_file_directory,
                                                   _tagged)) else 'fly.hdf5')
    experiment_file_path = os.path.join(experiment_file_directory,
                                        experiment_file_name)
    if not os.path.exists(experiment_file_path):
        print('No such file: ' + experiment_file_path)
        sys.exit(1)
    print('reading ' + experiment_file_path)

    if args.rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        plug = base_plugin.BasePlugin()

    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))
    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        int(series_num[0]), quiet=True)
    fly_metadata = ID.getSubjectMetadata()
    timing_channel_ind = (1 if fly_metadata.get('prep', '') ==
                          'fly right optic lobe' else 0)
    print('genotype {} x {}; notes: {}; timing_channel_ind {}'.format(
        fly_metadata.get('genotype_father'), fly_metadata.get('genotype_mother'),
        fly_metadata.get('notes'), timing_channel_ind))

    random_series, repeat_series = [], []
    for s in series_num:
        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name, int(s))
        ID = plug.ImagingDataObject
        if cb.is_contrast_blocks(ID.getRunParameters()):
            (repeat_series if cb.is_repeated(ID.getEpochParameters())
             else random_series).append(s)
    if not random_series and not repeat_series:
        print('No ContrastBlockUniformNoise series in this file.')
        sys.exit(0)
    print('random-seed series: {}; frozen-noise series: {}'.format(
        random_series, repeat_series))

    def load(s):
        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name, int(s))
        ID = plug.ImagingDataObject
        # before any timing call, or epochs are cut on the wrong photodiode
        ID.timing_channel_ind = timing_channel_ind
        return ID

    channels = [c.strip() for c in str(args.channels).split(',')]
    results, repeats = {}, {}
    for s in random_series:
        ID = load(s)
        for c in channels:
            print('\n== series {} ch{} (random seeds) =='.format(s, c))
            _t = time.perf_counter()
            try:
                results['sn' + s, 'ch' + c] = analyze_random_series(
                    ID, 'mask_ch' + c, args)
            except (KeyError, ValueError, AssertionError) as e:
                print('    skipped ({})'.format(e))
            print('    [{:.1f} s]'.format(time.perf_counter() - _t))
    for s in repeat_series:
        ID = load(s)
        for c in channels:
            print('\n== series {} ch{} (frozen noise) =='.format(s, c))
            # filters from the random series recorded CLOSEST in series order
            ref = None
            if random_series:
                near = min(random_series, key=lambda r: abs(int(r) - int(s)))
                ref = results.get(('sn' + near, 'ch' + c))
                if ref is not None:
                    print('    predicting from series {} filters'.format(near))
            try:
                repeats['sn' + s, 'ch' + c] = analyze_repeat_series(
                    ID, 'mask_ch' + c, args, ref)
                repeats['sn' + s, 'ch' + c]['ref'] = ref
            except (KeyError, ValueError, AssertionError) as e:
                print('    skipped ({})'.format(e))

    if args.output_directory:
        os.makedirs(args.output_directory, exist_ok=True)
        out_path = os.path.join(args.output_directory,
                                'contrast_blocks_{}.hdf5'.format(tag))
        figs_dir = os.path.join(args.output_directory, tag + '_roi_figs')
    else:
        out_path = experiment_file_path
        figs_dir = os.path.join(experiment_file_directory, tag + '_roi_figs')
    group_name = 'TRF_CONTRAST_BLOCKS_' + tag
    save_results(out_path, group_name, results, repeats, args)
    print('\nresults written to {} :/{}'.format(out_path, group_name))

    if save_figs or show_figs:
        os.makedirs(figs_dir, exist_ok=True)
        fly = os.path.basename(os.path.normpath(experiment_file_directory))
        for (sn, ch), res in results.items():
            sig = channel_names.get(ch[2:], '')
            ttl = '{} {} {}{}  contrast blocks, {} penalty {:.2g}'.format(
                fly, sn, ch, ' = ' + sig if sig else '',
                res['design']['penalty_kind'], res['ridge'])
            for name, fig in (('summary', figure_summary(res, ttl, sig)),
                              ('rois', figure_rois(res, ttl))):
                if save_figs:
                    _savefig(fig, os.path.join(figs_dir, 'trf_contrast_blocks_'
                                             '{}_{}_{}.pdf'.format(sn, ch, name)))
                if not show_figs:
                    plt.close(fig)
        for sn in sorted({k[0] for k in results}):
            entries = [(channel_names.get(ch[2:], ch), results[(s2, ch)])
                       for (s2, ch) in sorted(results) if s2 == sn]
            fig = figure_switch(entries, '{} {}  switch-triggered averages '
                                '({:g} ms jitter bins)'.format(
                                    fly, sn, 1e3 * float(args.sta_bin)
                                    if float(args.sta_bin) else float('nan')))
            if save_figs:
                _savefig(fig, os.path.join(figs_dir, 'trf_contrast_blocks_'
                                           '{}_switch.pdf'.format(sn)))
            if not show_figs:
                plt.close(fig)
        for (sn, ch), rep in repeats.items():
            fig = figure_repeat(rep, rep.get('ref'),
                                '{} {} {}{}  frozen contrast blocks'.format(
                                    fly, sn, ch,
                                    ' = ' + channel_names.get(ch[2:], '')
                                    if channel_names.get(ch[2:]) else ''))
            if save_figs:
                _savefig(fig, os.path.join(figs_dir, 'trf_contrast_blocks_'
                                         '{}_{}_frozen.pdf'.format(sn, ch)))
            if not show_figs:
                plt.close(fig)
        if save_figs:
            print('figures written to ' + figs_dir)
        if show_figs:
            plt.show()
