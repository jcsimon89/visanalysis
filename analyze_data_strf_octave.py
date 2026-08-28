"""
analyze_data_strf_octave.py

Per-ROI spatiotemporal linear filters (STRFs) from the multi-octave ternary
noise stimulus (labpack OctaveTernaryNoise), estimated in the CELL basis.

Call structure (same as analyze_data_strf_noizone.py):
    python analyze_data_strf_octave.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --show_figs False \
        --save_figs True \
        --tag final \
        --dff pre \
        --filter_length 5.0

Results saved to fly.hdf5 (or fly_final.hdf5) under group /STRF_OCTAVE_{tag}.

HOW THIS DIFFERS FROM THE NOIZONE PIPELINE
------------------------------------------
* No movie file. The stimulus is a pure function of the frame index, so the
  epoch metadata IS the stimulus. Reconstruction is a hash evaluation.

* Dropped frames need no repair of the stimulus array. The noizone path rebuilds
  stim_1d through reconstruct_displayed_frame_indices because its movie rows are
  positional. Here we simply ASK for the indices that actually rendered --
  reproduce_cells evaluates any index set exactly, gaps included.

* The regressor is cells on a sphere, not bars on a line. They are mutually
  independent with a common variance, so C_ss = 2p * I exactly: the plain STA is
  already ML up to a known constant. No whitening, no spatial ridge.

* The STA returns the filter convolved with the stimulus autocorrelation
  exp(-|lag|/tau), because the cells hold rather than refresh every frame. What
  that does, measured noiselessly, is narrower than it sounds:

      SPACE is untouched. For a filter separable within an octave,
      h(c,j) = a(c) k(j), the STA is a(c) (k*R)(L) -- R acts along the lag axis
      only, so the spatial profile comes through EXACTLY, r = 1.000000. Every
      spatial result (centroid, size, shape, the projection to visual space)
      therefore needs no correction at all.

      TIME is biased: the peak is late by 0.08-0.27 s depending on tau and the
      kernel, and the kernel is broadened by up to 2.2x.

  The bias is nearly constant within an octave -- 0.08 to 0.10 s at tau = 0.33
  across a 3.5x range of true peak -- so latency DIFFERENCES between ROIs are
  effectively unbiased and only the absolute number is off.

  Deconvolution (--deconvolve) is therefore OPT-IN, not the default. It is worth
  running only when you need absolute latency, and even then it earns its keep
  only at high SNR: measured against a known filter, it corrected the coarse
  layer from 1.88 s to 1.68 s against a true 1.60 s at 20 trials of 300 s, and
  was correctly refused everywhere else by the reliability guard.

  The case it cannot fix, and which no flag changes: a NON-separable filter
  inside one octave, centre and surround with different kinetics in the same
  layer. There the smear moves which lobe dominates the peak lag -- measured,
  peak 0.31 s -> 1.52 s and the spatial profile inverted to r = -0.60. Splitting
  centre and surround across octaves is what keeps you out of that regime, but
  the layers are overcomplete, so it is not guaranteed.

* The cell basis is overcomplete across octaves -- one 25 deg cell overlaps ~17
  of the 6 deg ones -- so the layers are NOT summed into a picture. Each octave
  independently measures the cell-averaged filter at its own scale.
"""

import sys
import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.util.octave_stim import (
    is_octave_ternary, resolve_octave_epoch_meta, load_octave_cells,
    octave_slices, cell_directions, cell_azel_deg, display_scale,
    effective_dof_frames,
)
from visanalysis.util.octave_strf import (
    compute_strf_cells_folds, peak_lag, independent_samples,
    deconvolve_octaves, zscore_by_octave, spherical_centroid,
    direction_to_azel_deg, coarse_assignment, fit_joint_cells,
)
from visanalysis.util import bleach, nullcal
from visanalysis.util.frame_timing import reconstruct_displayed_frame_indices
from visanalysis.util.strf_plot import _tighten
import h5py

plt.ioff()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_file_directory', nargs='?')
    parser.add_argument('--rig', nargs='?', default='Bruker')
    parser.add_argument('--show_figs', nargs='?', default='False')
    parser.add_argument('--save_figs', nargs='?', default='False')
    parser.add_argument('--tag', nargs='?', default='raw')
    parser.add_argument('--dff', nargs='?', default='pre')
    parser.add_argument('--null_shifts', nargs='?', type=int, default=4,
                        help='circular-shift null draws; 0 disables. Costs (n+1)x '
                             'the joint fit. See util.nullcal.')
    parser.add_argument('--bleach', nargs='?', default='auto',
                        help="bleach correction on raw F: 'auto' (default), 'on', 'off'")
    parser.add_argument('--filter_length', nargs='?', type=float, default=10.0,
                        help='STRF filter length in seconds')
    parser.add_argument('--monitor_hz', nargs='?', type=float, default=None,
                        help='override the display rate; default uses the ImagingDataObject')
    parser.add_argument('--lam', nargs='?', type=float, default=None,
                        help='deconvolution ridge. Omit to choose it by cross-validation, '
                             'which is the intended path -- a value that suits one indicator '
                             'will not suit another.')
    parser.add_argument('--lam_scale', nargs='?', type=float, default=1.0,
                        help='multiplies the cross-validated lambda. >1 is smoother.')
    parser.add_argument('--joint_ridge', nargs='?', type=float, default=1e-3,
                        help="ridge for the joint spatial solve, relative to "
                             "trace(X.T X)/n_fine. Conditioning only -- this is NOT "
                             "the temporal deconvolution, which stays off unless "
                             "--deconvolve is passed.")
    parser.add_argument('--deconvolve', action='store_true',
                        help='undo the temporal smear. OFF by default: it does nothing for '
                             'spatial results, which the smear leaves exact, and only '
                             'corrects absolute latency. Turn it on if you need that, and '
                             'check the per-octave decision it prints -- the reliability '
                             'guard refuses it wherever it would cost more than it buys.')
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    show_figs = args.show_figs == 'True'
    save_figs = args.save_figs == 'True'
    tag = args.tag
    dff = args.dff
    # acts on RAW F before dF/F; 'auto' skips epochs too short to separate
    # bleaching from signal. See visanalysis.util.bleach.
    bleach_mode = {'auto': 'auto', 'on': True, 'off': False}[str(args.bleach).lower()]
    null_shifts = int(args.null_shifts)
    filter_length_s = float(args.filter_length)
    monitor_hz_override = args.monitor_hz

    experiment_file_name = 'fly.hdf5' if tag == 'raw' else 'fly_' + tag + '.hdf5'
    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)
    if not os.path.exists(experiment_file_path):
        print('No such file: ' + experiment_file_path)
        sys.exit(1)

    if rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        plug = base_plugin.BasePlugin()

    response_set_name_prefix = 'mask_'
    func_channels_num = ['1', '2']

    # ── fly metadata and photodiode channel ──────────────────────────────────
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))
    ID = imaging_data.ImagingDataObject(experiment_file_path, int(series_num[0]), quiet=True)
    fly_metadata = ID.getSubjectMetadata()
    prep = fly_metadata.get('prep', '')
    timing_channel_ind = 1 if prep == 'fly right optic lobe' else 0
    print('fly_metadata: ' + repr(fly_metadata))
    print('timing_channel_ind: ' + repr(timing_channel_ind))

    # ── find the octave series ───────────────────────────────────────────────
    roi_data = {}
    run_parameters = {}
    octave_series_info = {}

    for current_series in series_num:
        sn = 'sn' + current_series
        plug.updateImagingDataObject(experiment_file_directory, experiment_file_name,
                                     int(current_series))
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6
        _, _, _sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * _sample_rate

        run_parameters[sn] = ID.getRunParameters()
        if not is_octave_ternary(run_parameters[sn]):
            continue

        metas = [resolve_octave_epoch_meta(ep) for ep in ID.getEpochParameters()]
        octave_series_info[sn] = {'epochs': metas}

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            _t0 = time.perf_counter()
            roi_data[sn, ch] = bleach.roi_responses(
                ID, response_set_name_prefix + ch,
                roi_prefix='aligned', dff=dff, correct_bleach=bleach_mode)
            print('    [timing] {} {}: getRoiResponses took {:.2f}s'.format(
                sn, ch, time.perf_counter() - _t0))

    if not octave_series_info:
        print('No octave ternary series identified. Check is_octave_ternary() and '
              'that process_data.py attached the epoch parameters for this series.')
        sys.exit(0)

    print('\nFound {} octave ternary series:'.format(len(octave_series_info)))
    for sn, info in octave_series_info.items():
        meta0 = info['epochs'][0][0]
        seeds = sorted({m[0]['seed'] for m in info['epochs']})
        fibs = sorted({m[0]['fibonacci_seed'] for m in info['epochs']})
        print('  {}: {} epochs  patches={} deg  taus={} s  cells={}  '
              '{} distinct seed(s), fibonacci_seed={}'.format(
                  sn, len(info['epochs']), tuple(meta0['patch_degrees']),
                  tuple(meta0['hold_taus']), sum(meta0['n_cells']), len(seeds), fibs))
        if len(fibs) > 1:
            print('    NOTE: the tessellation rotates within this series, so cells do '
                  'not mean the same thing across those epochs. They are grouped below.')

    # ── STRFs ────────────────────────────────────────────────────────────────
    strf_results = {}

    for sn, info in octave_series_info.items():
        pre_time = float(run_parameters[sn].get('pre_time', 0.0))

        plug.updateImagingDataObject(experiment_file_directory, experiment_file_name,
                                     int(sn.replace('sn', '')))
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6
        if monitor_hz_override is not None:
            ID.command_frame_rate = monitor_hz_override
        _, _, sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * sample_rate

        stimulus_timing = ID.getStimulusTiming()
        n_trials_timing = len(stimulus_timing['stimulus_start_times'])
        epoch_pre_times = np.array(ID.getEpochParameters('pre_time'))[:n_trials_timing]
        epoch_start_times = (stimulus_timing['stimulus_start_times'][:n_trials_timing]
                             - epoch_pre_times)
        drop_times_abs_s = stimulus_timing['dropped_frame_times'] / sample_rate
        frame_rate = ID.command_frame_rate
        raw_stim_dt = 1.0 / frame_rate

        print('  {}: {} dropped monitor frames (monitor_hz={})'.format(
            sn, len(drop_times_abs_s), frame_rate))

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            if (sn, ch) not in roi_data:
                continue

            epoch_response = roi_data[sn, ch]['epoch_response']
            time_vector_by_epoch = roi_data[sn, ch]['time_vector_by_epoch']
            n_roi = epoch_response.shape[0]
            roi_z_offsets = roi_data[sn, ch].get('roi_z_offsets', np.zeros(n_roi))

            # Period of the RESPONSE samples, which is the imaging rate and not the
            # display rate. n_valid counts imaging samples, so anything converting it
            # to a duration needs this and not raw_stim_dt -- the two differ by about
            # 12x and independent_samples is the difference between reporting a
            # recording as underdetermined and reporting it as ample.
            _dts = [np.median(np.diff(tv)) for tv in time_vector_by_epoch if len(tv) > 1]
            sample_dt = float(np.median(_dts)) if _dts else raw_stim_dt

            print('  {} {}: n_roi={}, n_trials={}, n_tp={}, imaging {:.2f} Hz'.format(
                sn, ch, n_roi, epoch_response.shape[1], epoch_response.shape[2],
                1.0 / sample_dt))

            # Group by geometry. Cells only mean the same thing across trials that
            # share a tessellation, so fibonacci_seed keys the group; the noise seed
            # varies freely within it and is exactly what we want varying.
            trials_by_geom = {}
            geom_meta = {}
            for trial_ind, (meta, _sampling) in enumerate(info['epochs']):
                if trial_ind >= epoch_response.shape[1]:
                    break
                tvec = time_vector_by_epoch[trial_ind]
                if len(tvec) == 0:
                    continue

                n_frames = int(np.ceil((tvec[-1] - pre_time) * frame_rate)) + 1
                if n_frames <= 0:
                    continue

                # Dropped frames: ask for the indices that actually rendered rather
                # than repairing an array afterwards. The display recomputes its frame
                # index from wall-clock time on every successful draw, so a dropped
                # frame is skipped permanently, not delayed -- and since every cell
                # value is a pure function of its index, evaluating the skipped
                # sequence is exact rather than an approximation.
                idx = np.arange(n_frames)
                if len(drop_times_abs_s) > 0 and trial_ind < len(epoch_start_times):
                    t0_abs = epoch_start_times[trial_ind] + epoch_pre_times[trial_ind]
                    t1_abs = tvec[-1] + epoch_start_times[trial_ind]
                    trial_drops = drop_times_abs_s[(drop_times_abs_s >= t0_abs)
                                                   & (drop_times_abs_s <= t1_abs)]
                    if len(trial_drops) > 0:
                        drop_inds = np.round((trial_drops - t0_abs) * frame_rate).astype(int)
                        actual = reconstruct_displayed_frame_indices(n_frames, drop_inds)
                        idx = np.where(actual == -1, np.arange(n_frames), actual)

                stim_cells = load_octave_cells(meta, n_frames, frame_indices=idx)

                key = int(meta['fibonacci_seed'])
                geom_meta[key] = meta
                grp = trials_by_geom.setdefault(key, [])
                for roi_ind in range(n_roi):
                    grp.append({
                        'roi_ind': roi_ind,
                        'response': epoch_response[roi_ind, trial_ind, :],
                        'time_vector': tvec + roi_z_offsets[roi_ind],
                        'stim_cells': stim_cells,
                        'stim_dt': raw_stim_dt,
                    })

            for fib, trials in trials_by_geom.items():
                meta = geom_meta[fib]
                n_cells = int(sum(meta['n_cells']))
                octs = octave_slices(meta)

                _t0 = time.perf_counter()
                strf_full, folds, t_lag, n_valid = compute_strf_cells_folds(
                    trials, n_roi, n_cells, filter_length_s, pre_time)
                print('    [timing] compute_strf_cells_folds ({} {} fib={}, {} '
                      'trial-roi entries) took {:.2f}s'.format(
                          sn, ch, fib, len(trials), time.perf_counter() - _t0))

                lag_dt = float(t_lag[1] - t_lag[0]) if len(t_lag) > 1 else raw_stim_dt

                if not args.deconvolve:
                    strf_dec, lam_used = strf_full, None
                    print('    not deconvolved (default). Spatial results are unaffected; '
                          'peak lags read late by roughly 0.1-0.3 s. Pass --deconvolve if '
                          'you need absolute latency.')
                else:
                    _t0 = time.perf_counter()
                    strf_dec, lam_used = deconvolve_octaves(
                        strf_full, octs, lag_dt, lam=args.lam,
                        lam_scale=args.lam_scale, cv_pair=folds)
                    print('    [timing] deconvolve_octaves took {:.2f}s'.format(
                        time.perf_counter() - _t0))
                    # Per octave, say what was decided and why. The guard can reject
                    # the deconvolution outright, and that has to be visible: an
                    # octave that fell back to the STA is still smeared by its hold,
                    # so its peak lag is late by roughly tau and its kernel is
                    # broadened. Silently returning either one would leave no way to
                    # tell which is on the page.
                    for _name, _sl, _d, _tau in octs:
                        _i = lam_used[_name]
                        _msg = '    {:>7}: {}'.format(_name, _i['how'])
                        if _i.get('lam') is not None:
                            _msg += '  lam={:.3g}'.format(_i['lam'])
                        if 'split_half_raw' in _i:
                            _msg += '  split-half {:.3f} -> {:.3f}'.format(
                                _i['split_half_raw'], _i['split_half_deconvolved'])
                        print(_msg)
                        if _i['how'].startswith('rejected'):
                            print('             -> kept the STA. It is still convolved '
                                  'with exp(-|t|/{:.2f}s), so read its peak lag as late.'
                                  .format(_tau))

                dirs = cell_directions(meta)
                azel = cell_azel_deg(meta)
                strf_z, z_info = zscore_by_octave(strf_dec, octs, lag_dt)

                # THE REPORTED FILTER. One filter on the fine cells, fitted from
                # every octave at once. There is no per-octave-then-merge path: the
                # two-step is the same estimator only in EXPECTATION, since it
                # replaces the Gram matrix with its expected value 2p*(I + B'B)
                # where the joint fit uses the one the stimulus actually realised.
                # On synthetic ground truth the joint fit wins wherever the two
                # differ -- 4.4% filter error against 17.0% at good SNR, 14.0
                # against 18.8 at moderate, 59 against 72 on a short recording --
                # and its closed form assumed a single coarse layer, so it was
                # never defined past two octaves anyway.
                fine_sl = octs[-1][1]
                joint_assign = [coarse_assignment(dirs[fine_sl], dirs[octs[j][1]])
                                for j in range(len(octs) - 1)]
                _t0 = time.perf_counter()
                strf_unified, _jlag, joint_info = fit_joint_cells(
                    trials, n_roi, list(meta['n_cells']), octs, dirs,
                    filter_length_s, pre_time, ridge=args.joint_ridge,
                    assignment=joint_assign)
                # Circular-shift null on the SAME estimator that produced
                # strf_unified -- a null from a different estimator would not
                # be the right reference. Only the per-ROI p and the threshold
                # are kept; no null arrays are stored.
                null_p, null_thr = None, None
                if null_shifts > 0:
                    _sdt = float(trials[0]['stim_dt'])

                    def _shifted(shift_s, _tr=trials, _sdt=_sdt):
                        k = int(round(shift_s / _sdt))
                        # roll ONCE per distinct stimulus array: it is shared
                        # across every ROI's trial dict, so rolling per trial
                        # would copy the same array hundreds of times
                        _rolled = {}
                        for t in _tr:
                            if id(t['stim_cells']) not in _rolled:
                                _rolled[id(t['stim_cells'])] = np.roll(
                                    t['stim_cells'], k, axis=0)
                        return fit_joint_cells(
                            [dict(t, stim_cells=_rolled[id(t['stim_cells'])])
                             for t in _tr],
                            n_roi, list(meta['n_cells']), octs, dirs,
                            filter_length_s, pre_time, ridge=args.joint_ridge,
                            assignment=joint_assign)[0]
                    _dur = trials[0]['stim_cells'].shape[0] * _sdt
                    _t1 = time.perf_counter()
                    _cal = nullcal.calibrate(
                        _shifted,
                        nullcal.shifts_for(_dur, filter_length_s, null_shifts))
                    _z = nullcal.peak_z(strf_unified)
                    null_p = np.full(n_roi, np.nan)
                    _ok = np.isfinite(_z)
                    null_p[_ok] = nullcal.pvalues(_z[_ok], _cal['null'])[0]
                    null_thr = _cal['threshold']
                    print('    [timing] null ({} shifts) took {:.2f}s; '
                          'threshold {:.2f}'.format(
                              null_shifts, time.perf_counter() - _t1, null_thr))

                print('    [timing] fit_joint_cells took {:.2f}s  '
                      '({} fine cells from {} octaves, ridge {:.0e}, {})'.format(
                          time.perf_counter() - _t0, strf_unified.shape[1], len(octs),
                          joint_info['ridge'],
                          'ANALYTIC Gram fallback' if joint_info['used_analytic']
                          else 'empirical Gram'))

                # Degrees of freedom differ by octave: a 1.00 s hold at 120 fps renews
                # every 120 frames against 40 for the 0.33 s one, so the same lag count
                # is three times less independent data. One z threshold applied to both
                # would be two different tests wearing the same number.
                dof = {name: independent_samples(n_valid, tau, sample_dt)
                       for name, _sl, _d, tau in octs}

                pk = np.array([peak_lag(strf_dec[r]) for r in range(n_roi)])

                # Centroid per ROI per octave, at that ROI's own peak lag. A weighted
                # mean of unit vectors renormalised to the sphere -- averaging az/el as
                # if Cartesian goes wrong at the edges of the lit band, which is exactly
                # where it would matter.
                centroids = {}
                n_kept_cells = {}
                for name, sl, _d, _tau in octs:
                    c = np.full((n_roi, 2), np.nan)
                    nk = np.zeros(n_roi, dtype=int)
                    for r in range(n_roi):
                        w = strf_dec[r, sl, pk[r]]
                        if not np.any(w):
                            continue
                        # The centroid of the DOMINANT LOBE, not of |w|. A centre
                        # surround filter has both signs at one lag, and |w| would
                        # average the two into a point between them that belongs to
                        # neither. Take the polarity of the strongest cell and keep
                        # only that lobe, which is also what the noizone workflow
                        # means by "cells above a fraction of the roi's own peak".
                        sign = np.sign(w[np.argmax(np.abs(w))])
                        lobe = np.clip(w * sign, 0.0, None)
                        d_hat, nk[r] = spherical_centroid(lobe, dirs[sl], frac=0.5)
                        if np.all(np.isfinite(d_hat)):
                            c[r] = direction_to_azel_deg(d_hat)
                    centroids[name] = c
                    n_kept_cells[name] = nk

                strf_results[sn, ch, fib] = {
                    'strf_raw': strf_full, 'strf': strf_dec, 'strf_z': strf_z,
                    't_lag': t_lag, 'n_valid': n_valid, 'lam': lam_used,
                    'meta': meta, 'octaves': octs, 'cell_azel_deg': azel,
                    'peak_lag_ind': pk, 'centroids': centroids, 'dof': dof,
                    'z_info': z_info, 'n_kept_cells': n_kept_cells,
                    'strf_unified': strf_unified, 'joint_info': joint_info,
                    # octs[-1], not octs[1] -- the fine layer is the LAST one, and
                    # those coincide only at two octaves.
                    'fine_azel_deg': azel[octs[-1][1]],
                    'display_scale': display_scale(meta),
                }
                print('    -> {} {} fib={}: STRF {} (n_roi, n_cells, n_lags) from {} '
                      'trial-roi entries'.format(sn, ch, fib, strf_dec.shape, len(trials)))
                for name, sl_, d, tau in octs:
                    n_cell_j = sl_.stop - sl_.start
                    indep = float(np.median(dof[name]))
                    print('       {:>7}: tau={:.2f}s  {:.0f} display frames per hold  '
                          '{:.0f} independent samples vs {} cells  noise sd={:.4g}'
                          .format(name, tau, effective_dof_frames(tau, lag_dt),
                                  indep, n_cell_j, z_info[name]['median_sd']))
                    if indep < n_cell_j:
                        print('                UNDERDETERMINED: fewer independent '
                              'stimulus configurations than cells, so the per-cell '
                              'estimate is poor however clean the recording is. '
                              'Needs about {:.0f}x more time.'
                              .format(n_cell_j / max(indep, 1e-9)))

    # ── save ─────────────────────────────────────────────────────────────────
    strf_group = 'STRF_OCTAVE_' + tag
    with h5py.File(experiment_file_path, 'a') as f:
        if strf_group in f:
            del f[strf_group]
        g = f.create_group(strf_group)
        g.attrs['filter_length_s'] = filter_length_s
        g.attrs['dff'] = dff
        g.attrs['deconvolved'] = bool(args.deconvolve)

        for (sn, ch, fib), res in strf_results.items():
            sg = g.require_group('{}/{}/fib_{}'.format(sn, ch, fib))
            for key in ('strf_raw', 'strf', 'strf_z'):
                sg.create_dataset(key, data=res[key], compression='gzip')
            # One filter on the fine cells, every octave folded in. This is the
            # thing to project to visual space; the per-octave arrays above are
            # kept because the unified fit cannot be undone from them.
            sg.create_dataset('strf_unified', data=res['strf_unified'],
                              compression='gzip')
            sg.create_dataset('unified_azel_deg', data=res['fine_azel_deg'])
            for k, v in res['joint_info'].items():
                sg.attrs['joint_' + k] = v
            sg.create_dataset('t_lag', data=res['t_lag'])
            sg.create_dataset('cell_azel_deg', data=res['cell_azel_deg'])
            sg.create_dataset('n_valid', data=res['n_valid'])
            sg.create_dataset('peak_lag_ind', data=res['peak_lag_ind'])
            sg.attrs['fibonacci_seed'] = fib
            sg.attrs['display_scale'] = res['display_scale']
            # Cell columns are only interpretable alongside the octave boundaries, so
            # the slices are recorded rather than left to be recomputed downstream.
            for name, sl, d, tau in res['octaves']:
                og = sg.require_group('octave_' + name)
                og.attrs['patch_deg'] = d
                og.attrs['hold_tau'] = tau
                og.attrs['cell_start'] = sl.start
                og.attrs['cell_stop'] = sl.stop
                og.create_dataset('centroid_azel_deg', data=res['centroids'][name])
                og.create_dataset('dof', data=res['dof'][name])
                og.attrs['median_noise_sd'] = res['z_info'][name]['median_sd']
                # How many cells actually defined each centroid -- a centroid backed
                # by two cells is a very different claim from one backed by forty.
                og.create_dataset('n_cells_in_centroid', data=res['n_kept_cells'][name])
                # Whether this octave was actually deconvolved, and on what evidence.
                # Recorded per octave because the guard decides per octave, so one
                # series can hold one deconvolved layer and one raw STA.
                if res['lam'] is not None:
                    _i = res['lam'][name]
                    og.attrs['deconv_how'] = _i['how']
                    og.attrs['deconv_applied'] = not _i['how'].startswith('rejected')
                    if _i.get('lam') is not None:
                        og.attrs['deconv_lam'] = float(_i['lam'])
                    if 'split_half_raw' in _i:
                        og.attrs['split_half_raw'] = _i['split_half_raw']
                        og.attrs['split_half_deconvolved'] = _i['split_half_deconvolved']
                else:
                    og.attrs['deconv_applied'] = False
                    og.attrs['deconv_how'] = 'not requested (default)'

    print('\nSaved STRFs to {} under /{}'.format(experiment_file_path, strf_group))

    # ── figures ──────────────────────────────────────────────────────────────
    if save_figs or show_figs:
        figs_dir = os.path.join(experiment_file_directory, tag + '_roi_figs')
        os.makedirs(figs_dir, exist_ok=True)

        for (sn, ch, fib), res in strf_results.items():
            strf, t_lag, octs = res['strf'], res['t_lag'], res['octaves']
            azel, pk = res['cell_azel_deg'], res['peak_lag_ind']
            n_roi = strf.shape[0]

            for r in range(n_roi):
                n_col = len(octs) + 1
                fig, axes = plt.subplots(2, n_col,
                                         figsize=(5.0 * n_col, 6.4), squeeze=False)
                lag_s = t_lag[pk[r]]
                for j, (name, sl, d, tau) in enumerate(octs):
                    w = strf[r, sl, pk[r]]
                    lim = np.abs(w).max() or 1.0

                    # One marker per cell, sized to the cell. Deliberately not
                    # interpolated onto a grid: the measurement IS per cell, and a
                    # smooth image would imply spatial detail the basis does not carry.
                    ax = axes[0][j]
                    ax.scatter(azel[sl, 0], azel[sl, 1], c=w, cmap='RdBu_r',
                               vmin=-lim, vmax=lim, s=(d * 2.2) ** 2,
                               edgecolors='none', alpha=0.85)
                    ax.set_title('{} cells @ lag {:.2f}s'.format(name, lag_s))
                    ax.set_xlabel('azimuth (deg)')
                    ax.set_ylabel('elevation (deg)')
                    ax.set_xlim(-69, 69)
                    ax.set_ylim(-51, 40)
                    ax.set_aspect('equal')
                    c = res['centroids'][name][r]
                    if np.all(np.isfinite(c)):
                        ax.plot(c[0], c[1], 'k+', markersize=12, markeredgewidth=1.6)

                    # The strongest cells' time courses, so the panel shows the filter
                    # rather than one arbitrary cell's trace.
                    ax = axes[1][j]
                    power = (strf[r, sl, :] ** 2).sum(axis=1)
                    for c_i in np.argsort(power)[-5:]:
                        ax.plot(t_lag, strf[r, sl, :][c_i], lw=1.0, alpha=0.75)
                    ax.axvline(lag_s, color='k', ls=':', lw=1.0)
                    ax.axhline(0, color='0.6', lw=0.8)
                    ax.set_xlabel('lag (s)')
                    ax.set_ylabel('filter (cell units)')
                    ax.set_title('{}: 5 strongest cells  tau={:.2f}s'.format(name, tau))

                # The unified filter, on the fine cells, at the same lag.
                mg = res['strf_unified'][r]
                m_azel = res['fine_azel_deg']
                lim = np.abs(mg[:, pk[r]]).max() or 1.0
                ax = axes[0][len(octs)]
                # octs[-1][2] is the FINE patch size; octs[1][2] happened to be the
                # same thing only while there were exactly two octaves.
                ax.scatter(m_azel[:, 0], m_azel[:, 1], c=mg[:, pk[r]], cmap='RdBu_r',
                           vmin=-lim, vmax=lim, s=(octs[-1][2] * 2.2) ** 2,
                           edgecolors='none', alpha=0.85)
                ax.set_title('JOINT @ lag {:.2f}s'.format(lag_s))
                ax.set_xlabel('azimuth (deg)')
                ax.set_ylabel('elevation (deg)')
                ax.set_xlim(-69, 69)
                ax.set_ylim(-51, 40)
                ax.set_aspect('equal')

                ax = axes[1][len(octs)]
                power = (mg ** 2).sum(axis=1)
                for c_i in np.argsort(power)[-5:]:
                    ax.plot(t_lag, mg[c_i], lw=1.0, alpha=0.75)
                ax.axvline(lag_s, color='k', ls=':', lw=1.0)
                ax.axhline(0, color='0.6', lw=0.8)
                ax.set_xlabel('lag (s)')
                ax.set_ylabel('filter (cell units)')
                ax.set_title('joint fit: 5 strongest cells')

                title = '{} {} fib={}  roi {}'.format(sn, ch, fib, r)
                if res['lam'] is None:
                    title += '  [STA, not deconvolved: peak lag reads late]'
                else:
                    _raw = [n for n, _s, _dd, _t in octs
                            if res['lam'][n]['how'].startswith('rejected')]
                    if _raw:
                        title += '  [{}: STA kept, still smeared]'.format(', '.join(_raw))
                fig.suptitle(title)
                _tighten(fig)
                if save_figs:
                    fig.savefig(os.path.join(
                        figs_dir,
                        'strf_octave_{}_{}_fib{}_roi{:03d}.pdf'.format(sn, ch, fib, r)))
                if not show_figs:
                    plt.close(fig)

        print('Figures written to ' + figs_dir)

    if show_figs:
        plt.show()
