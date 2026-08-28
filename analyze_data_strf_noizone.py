"""
analyze_data_strf.py

Computes per-ROI 1D spatiotemporal linear filters (STRFs) from flymax
stim_noizone bar-noise stimuli and combines orientation filters into 2D STRFs.

Call structure (same as analyze_data.py):
    python analyze_data_strf.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --show_figs False \
        --save_figs True \
        --tag raw \
        --dff pre \
        --filter_length 1.0

Results saved to fly.hdf5 (or fly_final.hdf5) under group /STRF_{tag}.
"""

import sys
import os
import argparse
import json
import pathlib
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.util.noise_stim import (
    is_noizone, resolve_noizone_epoch_noise_params,
    resolve_noizone_series_mesh_shape, load_bar_stim_1d,
)
from visanalysis.util import bleach, nullcal
from visanalysis.util.frame_timing import reconstruct_displayed_frame_indices
from visanalysis.util.strf_plot import (
    _pole_label, _panel_size, _grid_shape, _tighten,
)
import h5py

plt.ioff()

# point matplotlib at the ffmpeg binary bundled with this conda env, in case
# it's not on PATH when this script is invoked (e.g. via os.system from the
# shell wrapper without a fully-activated environment)
_ffmpeg_path = os.path.join(os.path.dirname(sys.executable), 'Library', 'bin', 'ffmpeg.exe')
if os.path.exists(_ffmpeg_path):
    plt.rcParams['animation.ffmpeg_path'] = _ffmpeg_path


# ────────────────────────────────────────────────────────────────────────────
# STRF computation
# ────────────────────────────────────────────────────────────────────────────

def compute_strf_1d(trials, n_roi, filter_length_s, pre_time=0.0):
    """
    Estimate 1D STRF by reverse correlation (response-triggered average),
    accumulated across trials that may each reference a different noise
    movie and/or stim_dt (e.g. different epochs within a series).

    Parameters
    ----------
    trials : list of dicts, one per (trial, roi) observation, each with:
        'roi_ind'     : int                          -- which ROI this is
                        (index into the strf output's first dimension)
        'response'    : (n_tp,)                       -- dF/F for this trial/roi
        'time_vector' : (n_tp,)                        -- actual imaging frame
                        times (s) for this ROI specifically -- already
                        includes any per-ROI timing correction (e.g. z-slice
                        offset within a volume); epoch-relative, stim onset
                        at t = pre_time
        'stim_1d'     : (n_stim_frames, n_positions)   -- mean-subtracted noise
                        for this trial's own movie
        'stim_dt'     : float                          -- this trial's stim
                        frame period (s)
        All entries must share the same n_positions (same orientation --
        group by orientation before calling this).
    n_roi          : int    -- total number of ROIs (size of output's first dim)
    filter_length_s: float  -- filter length to compute (s)
    pre_time       : float  -- baseline duration before stim onset (s)

    Returns
    -------
    strf     : (n_roi, n_positions, n_lags)  -- z-scored STA
    strf_raw : (n_roi, n_positions, n_lags)  -- raw STA (units: dF/F * stim)
    t_lag    : (n_lags,)  -- lag axis in seconds (0 = simultaneous, + = past)
    """
    n_positions = trials[0]['stim_1d'].shape[1]
    # stim_dt is the RAW display frame period (see the caller), so the lag axis
    # is deliberately oversampled: a noise image is held for noisetau ms, so
    # only filter_length_s / noisetau lags are independent. That is intentional
    # -- raw resolution is what makes the per-frame dropped-frame correction in
    # the caller exact, and the extra points interpolate rather than duplicate
    # because hold boundaries drift relative to the imaging clock.
    #
    # Two consequences to remember when reading the output:
    #   1. the kernel is smeared along tau by the stimulus autocorrelation,
    #      which for sample-and-hold is a TRIANGLE of half-width noisetau
    #      (hold = boxcar applied to a white sequence; its autocorrelation is
    #      boxcar (x) boxcar). It is symmetric, so it does NOT shift latency,
    #      and at typical noisetau it is small next to the indicator kernel.
    #   2. n_lags overstates the degrees of freedom by ~noisetau / lag_dt.
    #      Nothing here depends on that, but see peak_z_and_sign in
    #      average_strf_noizone.py, where it interacts with the z threshold.
    lag_dt = min(t['stim_dt'] for t in trials)
    n_lags = max(1, int(np.round(filter_length_s / lag_dt)))
    t_lag = np.arange(n_lags) * lag_dt

    strf_raw = np.zeros((n_roi, n_positions, n_lags), dtype=np.float64)
    n_valid_total = np.zeros(n_roi, dtype=np.int64)

    for trial in trials:
        roi_ind = trial['roi_ind']
        response, tvec = trial['response'], trial['time_vector']
        stim_1d, stim_dt = trial['stim_1d'], trial['stim_dt']
        if len(tvec) == 0:
            continue
        n_stim_frames = stim_1d.shape[0]
        n_tp = min(len(tvec), len(response))
        tvec_trial = tvec[:n_tp]
        resp_trial = response[:n_tp]

        # exclude response samples from the first filter_length_s of the noise
        # period -- they don't yet have a full filter_length_s of stimulus
        # history behind them, so every included sample gets a complete lag window
        time_valid = tvec_trial >= (pre_time + filter_length_s)

        # this trial's own frame 0 = its own onset (t = pre_time)
        stim_indices = np.round((tvec_trial - pre_time) / stim_dt).astype(int)

        for lag_idx, t in enumerate(t_lag):
            lag_frames = int(np.round(t / stim_dt))
            past_idx = stim_indices - lag_frames
            valid = time_valid & (past_idx >= 0) & (past_idx < n_stim_frames)
            if not valid.any():
                continue
            resp = resp_trial[valid].astype(np.float64)   # (n_valid,)
            stim = stim_1d[past_idx[valid], :]              # (n_valid, n_positions)
            strf_raw[roi_ind, :, lag_idx] += resp @ stim
            n_valid_total[roi_ind] += int(valid.sum())

    for roi_ind in range(n_roi):
        if n_valid_total[roi_ind] > 0:
            strf_raw[roi_ind] /= n_valid_total[roi_ind]

    strf = np.zeros_like(strf_raw)
    for roi in range(n_roi):
        s = strf_raw[roi].std()
        if s > 0:
            strf[roi] = strf_raw[roi] / s

    return strf, strf_raw, t_lag


def combine_strf_2d(strf_h, strf_v):
    """
    Combine two 1D STRFs from orthogonal orientations (separated by 90 deg)
    into a 2D STRF, assuming a separable (rank-1) receptive field.

    POLARITY
    --------
    A bare outer product cannot represent the receptive field's sign. For a
    separable field RF = A * f(p) * g(q), the marginals the 1D STRFs estimate
    are A*f*(int g) and A*g*(int f), so their product carries A SQUARED -- it
    comes out positive wherever the two marginals share a sign, regardless of
    whether that shared sign was + or -. An OFF-centre RF, whose 1D filters are
    both NEGATIVE at the centre, therefore appears as a strong POSITIVE peak in
    2D. That is a sign error, not a property of the data.

    The fix follows from the same algebra. With
        T = int int RF = the sum of either marginal
    we have outer(marginal_h, marginal_v) = RF * T exactly, so

        RF = outer(marginal_h, marginal_v) / T

    Dividing by T restores polarity AND relative scale, but is unstable
    wherever T ~ 0 -- a balanced centre-surround integrates to nothing -- which
    would blow up at weak lags. So only sign(T) is applied here, per
    (roi, lag): polarity is corrected, and the magnitude behaviour callers
    already scale their colour maps against is left alone.

    Both marginals estimate the same T up to a positive factor, so sign(T)
    computed from each must agree. Where it does not, separability is failing
    at that lag (or SNR is too low to tell) and the polarity is genuinely
    undetermined -- returned as `ambiguous` rather than silently resolved.

    Parameters
    ----------
    strf_h : (n_roi, n_pos_h, n_lags)
    strf_v : (n_roi, n_pos_v, n_lags)
      If n_lags differs between the two, truncates to the minimum.

    Returns
    -------
    strf_2d   : (n_roi, n_pos_h, n_pos_v, n_lags)
    polarity  : (n_roi, n_lags)      the sign applied, +1 / -1
    ambiguous : (n_roi, n_lags) bool the two marginals disagreed on sign(T)
    """
    n_lags = min(strf_h.shape[2], strf_v.shape[2])
    strf_h = strf_h[:, :, :n_lags]
    strf_v = strf_v[:, :, :n_lags]

    total_h = strf_h.sum(axis=1)              # (n_roi, n_lags) ~ T up to +ve scale
    total_v = strf_v.sum(axis=1)
    polarity = np.sign(total_h)
    polarity[polarity == 0] = 1.0
    ambiguous = np.sign(total_h) != np.sign(total_v)

    strf_2d = np.einsum('rpl,rql->rpql', strf_h, strf_v)
    return strf_2d * polarity[:, None, None, :], polarity, ambiguous


# ────────────────────────────────────────────────────────────────────────────
# Main
# ────────────────────────────────────────────────────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_file_directory', nargs='?')
    parser.add_argument('--rig', nargs='?', default='Bruker')
    parser.add_argument('--show_figs', nargs='?', default='False')
    parser.add_argument('--save_figs', nargs='?', default='False')
    parser.add_argument('--tag', nargs='?', default='raw')
    parser.add_argument('--null_shifts', nargs='?', type=int, default=4,
                        help='circular-shift null draws per group; 0 disables. '
                             'Costs (n+1)x the STRF computation. See util.nullcal.')
    parser.add_argument('--bleach', nargs='?', default='auto',
                        help="bleach correction on raw F: 'auto' (default, applied when the "
                             "epoch is long enough), 'on', or 'off'. See util.bleach.")
    parser.add_argument('--dff', nargs='?', default='pre')
    parser.add_argument('--filter_length', nargs='?', type=float, default=5.0,
                        help='STRF filter length in seconds')
    parser.add_argument('--monitor_hz', nargs='?', type=float, default=None,
                        help='Monitor refresh rate (Hz) for dropped frame detection. '
                             'Overrides ImagingDataObject default of 120 Hz. '
                             'Set to your actual monitor rate (e.g. 60).')
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    dff = args.dff
    # bleach correction acts on RAW F before dF/F -- see util.bleach for why it
    # must be multiplicative and why 'auto' skips short epochs rather than raising.
    bleach_mode = {'auto': 'auto', 'on': True, 'off': False}[str(args.bleach).lower()]
    null_shifts = int(args.null_shifts)
    filter_length_s = args.filter_length

    show_figs = args.show_figs == 'True'
    save_figs = args.save_figs == 'True'
    monitor_hz_override = args.monitor_hz  # None means use ImagingDataObject default (120 Hz)

    if args.tag in ('raw', 'final'):
        tag = args.tag
    else:
        raise NameError('--tag must be "raw" or "final"')

    experiment_file_name = 'fly.hdf5' if tag == 'raw' else 'fly_final.hdf5'
    json_file_name = 'fly.json'
    response_set_name_prefix = 'mask_'

    # ── fly metadata ─────────────────────────────────────────────────────────
    with open(pathlib.Path(experiment_file_directory, json_file_name), 'r') as f:
        fly_json = json.load(f)

    func_channels = (
        fly_json['functional_channel']
        .replace('[', '').replace(']', '').replace("'", '').split(',')
    )
    func_channels_num = [c.strip().split('_')[-1] for c in func_channels]
    struct_channel_num = [fly_json['structural_channel'].split('_')[-1]]

    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)
    print('experiment_file_path: ' + repr(experiment_file_path))

    # ── plugin ────────────────────────────────────────────────────────────────
    if rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
        print('****Bruker plugin****')
    elif rig == 'AODscope':
        from visanalysis.plugin import aodscope
        plug = aodscope.AodScopePlugin()
        print('****AODscope plugin****')
    else:
        plug = base_plugin.BasePlugin()
        print('****Unrecognized plugin name****')

    # ── series list and fly metadata ──────────────────────────────────────────
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))

    ID = imaging_data.ImagingDataObject(experiment_file_path, int(series_num[0]), quiet=False)
    fly_metadata = ID.getSubjectMetadata()
    print('fly_metadata: ' + repr(fly_metadata))

    prep = fly_metadata.get('prep', '')
    if prep == 'fly right optic lobe':
        timing_channel_ind = 1
    elif prep == 'fly left optic lobe':
        timing_channel_ind = 0
    else:
        timing_channel_ind = 0
        print('Could not determine photodiode channel from prep, defaulting to 0')
    print('timing_channel_ind: ' + repr(timing_channel_ind))

    # ── load noizone series only ────────────────────────────────────────────────
    # non-noise series are already handled by analyze_data.py in the shell workflow
    roi_data = {}
    run_parameters = {}
    acquisition_metadata = {}
    volume_frame_offsets = {}
    noise_series_info = {}   # sn -> {'mesh_shape', 'epochs': [{'bin_path','rtz','noisetau'}]}

    for current_series in series_num:
        current_series_int = int(current_series)
        sn = 'sn' + current_series

        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name,
                                     current_series_int)
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6  # photodiode trace threshold for up/down finding (0-1, normalized)
        _, _, _sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * _sample_rate  # widen tolerance to 1ms, scaled to sample_rate (datapoints)

        run_parameters[sn] = ID.getRunParameters()
        if not is_noizone(run_parameters[sn]):
            continue   # non-noise series already handled by analyze_data.py

        epoch_noise_params = [resolve_noizone_epoch_noise_params(ep) for ep in ID.getEpochParameters()]
        # stimpack writes protocol_parameters onto the series group
        # (experiment/data.py:107), so the display mesh -- which the bar-position
        # correction depends on -- is read once here rather than assumed or
        # re-read per epoch. Raises if the series ever lists more than one.
        series_mesh_shape = resolve_noizone_series_mesh_shape(run_parameters[sn])

        # Screen dimensions and bar geometry are no longer read from epoch
        # attrs -- each movie's master is authoritative and is loaded per
        # bin_path below. rtz/noisetau/bin_path can legitimately vary per epoch
        # (e.g. interleaved orientations), so those stay a per-epoch list.
        noise_series_info[sn] = {
            'mesh_shape': series_mesh_shape,
            'epochs': [
                {'bin_path': ep['bin_path'], 'rtz': float(ep['rtz']),
                 'noisetau': float(ep['noisetau'])}
                for ep in epoch_noise_params
            ],
        }

        acquisition_metadata[sn] = ID.getAcquisitionMetadata()
        volume_frame_offsets[sn] = ID.getVolumeFrameOffsets()

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            response_set_name = response_set_name_prefix + ch
            _t0 = time.perf_counter()
            roi_data[sn, ch] = bleach.roi_responses(ID, response_set_name,
                                                    roi_prefix='aligned',
                                                    dff=dff,
                                                    correct_bleach=bleach_mode)
            print(f'    [timing] {sn} {ch}: getRoiResponses took {time.perf_counter() - _t0:.2f}s')

    if not noise_series_info:
        print('No noizone series identified. '
              'Check is_noizone() and run_parameters key names above.')
        sys.exit(0)

    print(f'\nFound {len(noise_series_info)} noizone series:')
    for sn, info in noise_series_info.items():
        n_epochs = len(info['epochs'])
        unique_rtz = sorted({e['rtz'] for e in info['epochs']})
        unique_noisetau = sorted({e['noisetau'] for e in info['epochs']})
        unique_files = sorted({pathlib.Path(e['bin_path']).name for e in info['epochs']})
        print(f'  {sn}: {n_epochs} epochs  rtz values={unique_rtz}  '
              f'noisetau values={unique_noisetau} ms  '
              f'{len(unique_files)} unique bin file(s)  '
              f'display mesh={info["mesh_shape"]}')

    # ── compute 1D STRFs ──────────────────────────────────────────────────────
    # (n_roi, n_positions, n_lags) per (sn, ch, rtz, noisetau). Epochs within a
    # series can reference different bin files/rtz/noisetau, so trials are
    # grouped by (exact angle, noisetau) -- n_positions differs by angle, and
    # different update rates are kept as separate filters rather than blended.
    strf_results = {}
    stim_1d_cache = {}   # (bin_path, mesh_shape) -> (stim_1d, bar_pos_deg, meta)

    def get_cached_stim_1d(bin_path, mesh_shape):
        # The exact per-bar noise sequence, from the movie's noiz1d master
        # (_master.json + _values.bin + _barindex.bin). Each bar is one
        # regressor, so there is no projection or binning step and nothing is
        # approximated. Bar width, phase, pole, pose and timebase all come from
        # the master, so nothing is re-derived from filenames or params.json.
        # Values are expanded to raw display-frame resolution on load, so the
        # dropped-frame correction below still applies at the true raw-frame
        # level. Cached per bin_path: one movie, one orientation.
        key = (bin_path, mesh_shape)
        if key not in stim_1d_cache:
            _t0 = time.perf_counter()
            stim_1d, bar_pos_deg, meta = load_bar_stim_1d(
                bin_path, mesh_shape=mesh_shape)
            print(f'    [timing] load_bar_stim_1d took {time.perf_counter() - _t0:.2f}s '
                  f'-> stim_1d shape={stim_1d.shape}, '
                  f'bars {meta["bar_index"].min()}-{meta["bar_index"].max()} @ {meta["bar_width_deg"]:g} deg, '
                  f'position axis {bar_pos_deg[0]:+.1f} to {bar_pos_deg[-1]:+.1f} deg')
            _shift = meta['clip_shift_deg'] + meta['mesh_shift_deg']
            print(f'    [geometry] mesh {meta["mesh_shape"]}: bar position correction '
                  f'max {np.abs(_shift).max():.2f} deg, median {np.median(np.abs(_shift)):.2f} deg '
                  f'(clip max {np.abs(meta["clip_shift_deg"]).max():.2f}, '
                  f'mesh max {np.abs(meta["mesh_shift_deg"]).max():.2f})')
            stim_1d_cache[key] = (stim_1d - stim_1d.mean(axis=0), bar_pos_deg, meta)
        return stim_1d_cache[key]

    for sn, info in noise_series_info.items():
        mesh_shape = info['mesh_shape']   # all bar geometry comes from the master
        epochs = info['epochs']

        pre_time = float(run_parameters[sn].get('pre_time', 0.0))

        # ── dropped frame correction ─────────────────────────────────────────
        # Re-load ID for this series to access stimulus timing
        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name,
                                     int(sn.replace('sn', '')))
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6  # photodiode trace threshold for up/down finding (0-1, normalized)
        if monitor_hz_override is not None:
            ID.command_frame_rate = monitor_hz_override

        _, _, sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * sample_rate  # widen tolerance to 1ms, scaled to sample_rate (datapoints)

        stimulus_timing = ID.getStimulusTiming()

        n_trials_timing = len(stimulus_timing['stimulus_start_times'])
        epoch_pre_times = np.array(ID.getEpochParameters('pre_time'))[:n_trials_timing]
        epoch_start_times = stimulus_timing['stimulus_start_times'][:n_trials_timing] - epoch_pre_times

        n_dropped = len(stimulus_timing['dropped_frame_times'])
        print(f'  {sn}: {n_dropped} dropped monitor frames detected '
              f'(monitor_hz={ID.command_frame_rate})')

        _durs_ms = stimulus_timing['all_frame_durations'] / sample_rate * 1000.0
        _ideal_ms = stimulus_timing['ideal_frame_len'] / sample_rate * 1000.0
        _slop_ms = np.asarray(ID.frame_slop).item() / sample_rate * 1000.0
        if _durs_ms.size > 0:
            _pct = np.percentile(_durs_ms, [1, 5, 25, 50, 75, 95, 99])
            print(f'    measured inter-frame intervals (ms): '
                  f'1%={_pct[0]:.3f} 5%={_pct[1]:.3f} 25%={_pct[2]:.3f} 50%={_pct[3]:.3f} '
                  f'75%={_pct[4]:.3f} 95%={_pct[5]:.3f} 99%={_pct[6]:.3f}  '
                  f'| ideal={_ideal_ms:.3f}ms +/- {_slop_ms:.3f}ms tolerance  '
                  f'(n={_durs_ms.size})')

        for _label, _key in [('up-to-up', 'up_to_up_intervals'), ('down-to-down', 'down_to_down_intervals')]:
            _same_type_ms = stimulus_timing[_key] / sample_rate * 1000.0
            if _same_type_ms.size > 0:
                _pct2 = np.percentile(_same_type_ms, [1, 5, 25, 50, 75, 95, 99])
                print(f'    {_label} intervals (ms): '
                      f'1%={_pct2[0]:.3f} 5%={_pct2[1]:.3f} 25%={_pct2[2]:.3f} 50%={_pct2[3]:.3f} '
                      f'75%={_pct2[4]:.3f} 95%={_pct2[5]:.3f} 99%={_pct2[6]:.3f}  '
                      f'(n={_same_type_ms.size})')

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            if (sn, ch) not in roi_data:
                continue

            epoch_response        = roi_data[sn, ch]['epoch_response']        # (n_roi, n_trials, n_tp)
            time_vector_by_epoch  = roi_data[sn, ch]['time_vector_by_epoch']  # list of n_trials arrays

            n_roi = epoch_response.shape[0]
            # per-ROI timing correction for z-slice offset within a volume
            # (zeros if not volumetric / not available -- see getRoiResponses)
            roi_z_offsets = roi_data[sn, ch].get('roi_z_offsets', np.zeros(n_roi))

            print(f'  {sn} {ch}: n_roi={n_roi}, '
                  f'n_trials={epoch_response.shape[1]}, n_tp={epoch_response.shape[2]}')

            # dropped-frame times, in seconds -- used below to correct the
            # stimulus content (not response timing; see
            # reconstruct_displayed_frame_indices) for whichever trials they
            # actually fall within.
            drop_times_abs_s = stimulus_timing['dropped_frame_times'] / sample_rate
            frame_rate = ID.command_frame_rate  # true raw display frame rate
            # stim_1d (from get_cached_stim_1d) is at native raw display-frame
            # resolution -- every bin-file row is one real display frame, not
            # one noise-image update (a noise image is simply held across
            # many consecutive identical rows). Indexing stim_1d must use this
            # raw frame period, not noisetau, so per-raw-frame dropped-frame
            # correction lines up with it exactly (see below).
            raw_stim_dt = 1.0 / frame_rate

            # flatten to one dict per (trial, roi) observation, grouped by
            # (exact rtz, noisetau) -- n_positions differs by angle, and
            # different update rates are kept as separate filters rather than
            # blended together. Each observation's own time_vector already
            # has that ROI's z-slice offset baked in, so compute_strf_1d
            # doesn't need to know about ROI timing at all.
            trials_by_group = {}   # (rtz, noisetau) -> list of (trial, roi) dicts
            group_position_info = {}   # (rtz, noisetau) -> bar_pos_deg
            for trial_ind, ep in enumerate(epochs):
                if trial_ind >= epoch_response.shape[1]:
                    break
                stim_1d, bar_pos_deg, _bar_meta = get_cached_stim_1d(
                    ep['bin_path'], mesh_shape)
                group_position_info[ep['rtz'], ep['noisetau']] = bar_pos_deg
                tvec = time_vector_by_epoch[trial_ind]

                # correct stim_1d for this trial's dropped display frames, at
                # native raw-frame resolution: flymax recomputes its
                # displayed frame index fresh from wall-clock time on every
                # successful draw, so a dropped frame is permanently skipped
                # (not delayed). Since stim_1d is already raw-frame indexed,
                # each dropped raw frame is corrected individually here --
                # no collapsing to noise-image granularity, so a drop doesn't
                # get treated as if it corrupted the whole surrounding
                # noise-image window when only one raw frame within it did.
                if trial_ind < len(epoch_start_times) and len(drop_times_abs_s) > 0 and len(tvec) > 0:
                    noise_start_abs = epoch_start_times[trial_ind] + epoch_pre_times[trial_ind]
                    noise_end_abs = tvec[-1] + epoch_start_times[trial_ind]
                    trial_drops_s = drop_times_abs_s[
                        (drop_times_abs_s >= noise_start_abs) & (drop_times_abs_s <= noise_end_abs)
                    ]
                    if len(trial_drops_s) > 0:
                        duration = noise_end_abs - noise_start_abs
                        n_raw_frames_trial = min(
                            stim_1d.shape[0],
                            int(np.ceil(duration * frame_rate)) + 1,
                        )
                        drop_raw_inds = np.round(
                            (trial_drops_s - noise_start_abs) * frame_rate
                        ).astype(int)
                        actual_raw_idx = reconstruct_displayed_frame_indices(n_raw_frames_trial, drop_raw_inds)
                        actual_raw_idx = np.where(
                            actual_raw_idx == -1, np.arange(n_raw_frames_trial), actual_raw_idx,
                        )
                        stim_1d = stim_1d.copy()
                        stim_1d[:n_raw_frames_trial] = stim_1d[actual_raw_idx]

                group = trials_by_group.setdefault((ep['rtz'], ep['noisetau']), [])
                for roi_ind in range(n_roi):
                    group.append({
                        'roi_ind': roi_ind,
                        'response': epoch_response[roi_ind, trial_ind, :],
                        'time_vector': tvec + roi_z_offsets[roi_ind],
                        'stim_1d': stim_1d,
                        'stim_dt': raw_stim_dt,
                    })

            for (rtz, noisetau), trials in trials_by_group.items():
                _t0 = time.perf_counter()
                strf, strf_raw, t_lag = compute_strf_1d(trials, n_roi, filter_length_s, pre_time)
                print(f'    [timing] compute_strf_1d ({sn} {ch} rtz={rtz} noisetau={noisetau}ms, '
                      f'{len(trials)} trial-roi entries) took {time.perf_counter() - _t0:.2f}s')
                # Circular-shift null. Rolling the stimulus destroys the
                # pairing while keeping its autocorrelation and the response's
                # own structure, so peak |z| under it is the right reference
                # for "does this ROI respond at all". Costs (n+1)x this
                # computation; what gets STORED is only the per-ROI p-value and
                # the threshold, so q stays adjustable downstream and no null
                # arrays go into the h5.
                null_p, null_thr = None, None
                if null_shifts > 0:
                    # from the trials, not the loop variable: raw_stim_dt is
                    # bound inside the epoch loop above and would otherwise
                    # carry whatever the LAST epoch set
                    _sdt = float(trials[0]['stim_dt'])

                    def _shifted(shift_s, _tr=trials, _sdt=_sdt):
                        k = int(round(shift_s / _sdt))
                        # roll ONCE per distinct stimulus array: it is shared
                        # across every ROI's trial dict, so rolling per trial
                        # would copy the same array hundreds of times
                        _rolled = {}
                        for t in _tr:
                            if id(t['stim_1d']) not in _rolled:
                                _rolled[id(t['stim_1d'])] = np.roll(
                                    t['stim_1d'], k, axis=0)
                        return compute_strf_1d(
                            [dict(t, stim_1d=_rolled[id(t['stim_1d'])])
                             for t in _tr],
                            n_roi, filter_length_s, pre_time)[0]
                    _dur = trials[0]['stim_1d'].shape[0] * _sdt
                    _t1 = time.perf_counter()
                    _cal = nullcal.calibrate(
                        _shifted,
                        nullcal.shifts_for(_dur, filter_length_s, null_shifts))
                    _z = nullcal.peak_z(strf)
                    null_p = np.full(n_roi, np.nan)
                    _ok = np.isfinite(_z)
                    null_p[_ok] = nullcal.pvalues(_z[_ok], _cal['null'])[0]
                    null_thr = _cal['threshold']
                    print('    [timing] null ({} shifts) took {:.2f}s; '
                          'threshold {:.2f}'.format(
                              null_shifts, time.perf_counter() - _t1, null_thr))

                strf_results[sn, ch, rtz, noisetau] = {
                    'strf': strf,
                    'strf_raw': strf_raw,
                    't_lag': t_lag,
                    'null_p': null_p,
                    'null_threshold': null_thr,
                    'null_shifts': null_shifts,
                    'rtz': rtz,
                    'noisetau': noisetau,
                    'bar_pos_deg': group_position_info[rtz, noisetau],
                }
                print(f'    -> {sn} {ch} rtz={rtz} noisetau={noisetau}ms: '
                      f'STRF shape {strf.shape}  (n_roi, n_positions, n_lags)  '
                      f'from {len(trials)} trials')

    # ── combine orthogonal angle pairs into 2D STRFs (outer product) ──────────
    # (n_roi, n_pos_h, n_pos_v, n_lags) per (sn_h, sn_v, ch, noisetau, h_angle,
    # v_angle). Assumes noise angles come in pairs separated by 90 deg (mod
    # 180); pairing is across (channel, noisetau) only -- NOT restricted to a
    # single series -- since a fly's two orthogonal orientations may each be
    # presented in their own separate series (e.g. sn2=rtz0, sn3=rtz90) rather
    # than interleaved within one series. Each pair is combined independently
    # via a separable (rank-1) outer product. Angles without an orthogonal
    # partner are skipped with a warning.
    strf_2d_results = {}

    channel_noisetau = {(c, nt) for (s, c, o, nt) in strf_results}
    for ch, noisetau in channel_noisetau:
        group_results = {}   # rtz -> (res, series it came from)
        for (s, c, o, nt), res in strf_results.items():
            if c != ch or nt != noisetau:
                continue
            if o in group_results:
                raise ValueError(
                    f'{ch} noisetau={noisetau}ms: rtz={o} found in more than one series '
                    f'({group_results[o][1]} and {s}) -- ambiguous for 2D pairing.'
                )
            group_results[o] = (res, s)

        available = sorted(group_results)
        used = set()
        for theta in available:
            if theta in used:
                continue
            partner = next(
                (r for r in available if r not in used and r != theta
                 and abs((r - theta - 90) % 180) < 1e-6),
                None,
            )
            if partner is None:
                _, sn_theta = group_results[theta]
                print(f'\n{sn_theta} {ch} noisetau={noisetau}ms: rtz={theta} has no orthogonal '
                      f'partner (90 deg away) among {available} -- skipping')
                used.add(theta)
                continue
            used.add(theta)
            used.add(partner)

            res_h, sn_h = group_results[theta]
            res_v, sn_v = group_results[partner]
            t_lag = res_h['t_lag']
            strf_2d, polarity, ambiguous = combine_strf_2d(res_h['strf'], res_v['strf'])
            n_lags_2d = strf_2d.shape[3]
            strf_2d_results[sn_h, sn_v, ch, noisetau, theta, partner] = {
                'strf_2d': strf_2d, 't_lag': t_lag[:n_lags_2d],
                'bar_pos_h_deg': res_h['bar_pos_deg'],
                'bar_pos_v_deg': res_v['bar_pos_deg'],
                'polarity': polarity, 'polarity_ambiguous': ambiguous,
            }
            n_off = int((polarity < 0).sum())
            print(f'    [polarity] {n_off}/{polarity.size} (roi, lag) entries are '
                  f'OFF-sign and were inverted; {int(ambiguous.sum())} ambiguous '
                  f'(the two marginals disagreed on sign)')
            print(f'\n{sn_h}(rtz={theta}) + {sn_v}(rtz={partner})  {ch} noisetau={noisetau}ms: '
                  f'2D STRF shape {strf_2d.shape}  (n_roi, n_pos_h, n_pos_v, n_lags)')

    # ── save to HDF5 ───────────────────────────────────────────────────────────
    strf_group = f'STRF_{tag}'
    with h5py.File(experiment_file_path, 'a') as f:
        if strf_group in f:
            del f[strf_group]
        g = f.create_group(strf_group)
        g.attrs['filter_length_s'] = filter_length_s
        g.attrs['dff'] = dff

        for (sn, ch, orientation, noisetau), res in strf_results.items():
            sg = g.require_group(f'{sn}/{ch}/rtz_{orientation}/noisetau_{noisetau}')
            sg.create_dataset('strf',     data=res['strf'],     compression='gzip')
            sg.create_dataset('strf_raw', data=res['strf_raw'], compression='gzip')
            sg.create_dataset('t_lag',    data=res['t_lag'])
            sg.create_dataset('bar_pos_deg', data=res['bar_pos_deg'])
            sg.attrs['rtz'] = res['rtz']
            sg.attrs['noisetau'] = res['noisetau']

        for (sn_h, sn_v, ch, noisetau, theta, partner), res in strf_2d_results.items():
            g2 = g.require_group(f'combined_2d/{sn_h}_{sn_v}/{ch}/noisetau_{noisetau}/rtz_{theta}_{partner}')
            g2.create_dataset('strf_2d', data=res['strf_2d'], compression='gzip')
            g2.create_dataset('t_lag',   data=res['t_lag'])
            g2.create_dataset('bar_pos_h_deg', data=res['bar_pos_h_deg'])
            g2.create_dataset('bar_pos_v_deg', data=res['bar_pos_v_deg'])
            g2.create_dataset('polarity', data=res['polarity'])
            g2.create_dataset('polarity_ambiguous', data=res['polarity_ambiguous'])
            g2.attrs['rtz_h'] = theta
            g2.attrs['rtz_v'] = partner

    print(f'\nSaved STRFs to {experiment_file_path} under /{strf_group}')

    # ── figures ────────────────────────────────────────────────────────────────
    figs_dir = os.path.join(experiment_file_directory, f'{tag}_roi_figs')
    os.makedirs(figs_dir, exist_ok=True)
    fig_format = '.pdf'
    movie_fps = 30

    # extract n_roi from first available entry
    first_key = next(iter(roi_data))
    n_roi = len(roi_data[first_key]['roi_response'])

    # figure 1: 1D STRF -- one figure per ROI, per series/channel/orientation/noisetau
    for (sn, ch, orientation, noisetau), res in strf_results.items():
        strf   = res['strf']      # (n_roi, n_positions, n_lags)
        t_lag  = res['t_lag']
        rot    = res['rtz']
        bar_pos_deg = res['bar_pos_deg']   # (n_positions,), ascending, matches strf's position axis
        ori_label = f'rot{int(rot)}'
        vmax = np.abs(strf).max() or 1.0

        for roi_ind in range(n_roi):
            fig_name = f'strf_1d_{ori_label}_noisetau{int(noisetau)}_{ch}_{sn}_{tag}_roi_{roi_ind}_'
            fh, ax = plt.subplots(figsize=(4, 4), constrained_layout=True)
            im = ax.imshow(
                strf[roi_ind],           # (n_positions, n_lags)
                aspect='auto',
                origin='lower',
                extent=[t_lag[0] * 1000, t_lag[-1] * 1000, bar_pos_deg[0], bar_pos_deg[-1]],
                cmap='RdBu_r',
                vmin=-vmax, vmax=vmax,
            )
            ax.set_xlabel('Lag (ms)')
            ax.set_ylabel(f'{_pole_label(rot)}  (deg, 90-psi)')
            plt.colorbar(im, ax=ax, label='z-score')
            plt.suptitle(f'1D STRF  {ori_label}  noisetau={noisetau}ms  {ch}  {sn}  roi {roi_ind}  {tag}')
            if save_figs:
                plt.savefig(os.path.join(figs_dir, fig_name + fig_format),
                            dpi=200, transparent=True)
            if show_figs:
                plt.show()
            plt.close()

    # figure 2: 2D STRF at peak lag -- one figure per ROI, per series pair/channel/noisetau/angle pair
    for (sn_h, sn_v, ch, noisetau, theta, partner), res in strf_2d_results.items():
        strf_2d = res['strf_2d']   # (n_roi, n_pos_h, n_pos_v, n_lags)
        t_lag   = res['t_lag']
        n_roi_2d = strf_2d.shape[0]
        pos_h_deg = res['bar_pos_h_deg']
        pos_v_deg = res['bar_pos_v_deg']
        h_label = _pole_label(theta)
        v_label = _pole_label(partner)
        extent_2d = [pos_v_deg[0], pos_v_deg[-1], pos_h_deg[0], pos_h_deg[-1]]

        peak_lag_idx = int(np.abs(strf_2d).mean(axis=(0, 1, 2)).argmax())
        peak_lag_ms  = t_lag[peak_lag_idx] * 1000

        for roi_ind in range(n_roi_2d):
            fig_name = f'strf_2d_peak_lag_noisetau{int(noisetau)}_rtz{int(theta)}_{int(partner)}_{sn_h}_{sn_v}_{ch}_{tag}_roi_{roi_ind}_'
            # single color scale across all lags for this ROI (not just this
            # panel's lag), so peak-lag/subset/movie figures are comparable
            vmax = np.abs(strf_2d[roi_ind]).max() or 1.0
            fh, ax = plt.subplots(figsize=_panel_size(extent_2d, 4.2),
                                  constrained_layout=True)
            # rows = h-group bar position, cols = v-group bar position, both
            # in true visual-angle degrees -- aspect='equal' with a
            # degrees-based extent keeps equal spacing = equal true visual
            # angle on both axes (rather than relying on bin counts happening
            # to be proportional to each axis's true FOV)
            im = ax.imshow(
                strf_2d[roi_ind, :, :, peak_lag_idx],
                aspect='equal',
                origin='lower',
                extent=extent_2d,
                cmap='RdBu_r',
                vmin=-vmax, vmax=vmax,
            )
            ax.set_xlabel(f'{v_label}  (deg, 90-psi)')
            ax.set_ylabel(f'{h_label}  (deg, 90-psi)')
            plt.colorbar(im, ax=ax)
            plt.suptitle(f'2D STRF at peak lag ({peak_lag_ms:.0f} ms)  {sn_h}+{sn_v}  {ch}  '
                         f'noisetau={noisetau}ms  rtz=({theta},{partner})  roi {roi_ind}  {tag}')
            if save_figs:
                plt.savefig(os.path.join(figs_dir, fig_name + fig_format),
                            dpi=200, transparent=True)
            if show_figs:
                plt.show()
            plt.close()

    # figure 3: 2D STRF across a subset of lags -- per ROI, per series/channel/noisetau
    # showing every lag gets unreadable for long filter_length/fine noisetau
    # (e.g. 5s / 50ms = 100 lags), so a fixed number of evenly-spaced lags are
    # shown instead, regardless of the total lag count.
    n_display_lags = 24
    for (sn_h, sn_v, ch, noisetau, theta, partner), res in strf_2d_results.items():
        strf_2d = res['strf_2d']   # (n_roi, n_pos_h, n_pos_v, n_lags)
        t_lag   = res['t_lag']
        n_roi_2d, _, _, n_lags = strf_2d.shape
        pos_h_deg = res['bar_pos_h_deg']
        pos_v_deg = res['bar_pos_v_deg']
        extent_2d = [pos_v_deg[0], pos_v_deg[-1], pos_h_deg[0], pos_h_deg[-1]]
        lag_idxs = np.unique(np.linspace(0, n_lags - 1, min(n_display_lags, n_lags)).round().astype(int))

        for roi_ind in range(n_roi_2d):
            fig_name = f'strf_2d_lags_noisetau{int(noisetau)}_rtz{int(theta)}_{int(partner)}_{sn_h}_{sn_v}_{ch}_{tag}_roi_{roi_ind}_'
            # single color scale across all lags for this ROI, not per panel
            vmax = np.abs(strf_2d[roi_ind]).max() or 1.0
            # Grid and panel size both follow the data aspect, so panels sit
            # flush instead of each image being letterboxed inside a 4x2 box.
            n_rows, n_cols = _grid_shape(len(lag_idxs), extent_2d)
            pw, ph = _panel_size(extent_2d, 1.7)
            fh, axes = plt.subplots(n_rows, n_cols,
                                    figsize=(pw * n_cols, ph * n_rows + 0.35),
                                    constrained_layout=True)
            _tighten(fh)
            axes_flat = np.array(axes).flatten()
            for panel_ind, lag_idx in enumerate(lag_idxs):
                ax = axes_flat[panel_ind]
                ax.imshow(
                    strf_2d[roi_ind, :, :, lag_idx],
                    aspect='equal',
                    origin='lower',
                    extent=extent_2d,
                    cmap='RdBu_r',
                    vmin=-vmax, vmax=vmax,
                )
                ax.text(0.03, 0.97, f'{t_lag[lag_idx]*1000:.0f} ms',
                        transform=ax.transAxes, ha='left', va='top', fontsize=7,
                        bbox=dict(fc='white', ec='none', alpha=0.75, pad=1))
                ax.set_xticks([])
                ax.set_yticks([])
                for _s in ax.spines.values():
                    _s.set_visible(False)
            for ax in axes_flat[len(lag_idxs):]:
                ax.set_visible(False)
            plt.suptitle(f'2D STRF ({len(lag_idxs)} of {n_lags} lags)  {sn_h}+{sn_v}  {ch}  '
                         f'noisetau={noisetau}ms  rtz=({theta},{partner})  roi {roi_ind}  {tag}')
            if save_figs:
                plt.savefig(os.path.join(figs_dir, fig_name + fig_format),
                            dpi=150, transparent=True)
            if show_figs:
                plt.show()
            plt.close()

    # figure 4: 2D STRF across every lag -- per ROI, per series/channel/noisetau
    # (same as figure 3 but unsubsampled). Commented out for now: with n_lags
    # now ~600 (raw display-frame resolution) this grid gets huge and mostly
    # unusable as a static image, and figure 6's movie shows the same thing
    # far more usefully.
    # for (sn_h, sn_v, ch, noisetau, theta, partner), res in strf_2d_results.items():
    #     strf_2d = res['strf_2d']   # (n_roi, n_pos_h, n_pos_v, n_lags)
    #     t_lag   = res['t_lag']
    #     n_roi_2d, _, _, n_lags = strf_2d.shape
    #
    #     for roi_ind in range(n_roi_2d):
    #         fig_name = f'strf_2d_all_lags_noisetau{int(noisetau)}_rtz{int(theta)}_{int(partner)}_{sn_h}_{sn_v}_{ch}_{tag}_roi_{roi_ind}_'
    #         n_cols = min(n_lags, 8)
    #         n_rows = int(np.ceil(n_lags / n_cols))
    #         fh, axes = plt.subplots(n_rows, n_cols,
    #                                 figsize=(2 * n_cols, 2 * n_rows),
    #                                 constrained_layout=True)
    #         axes_flat = np.array(axes).flatten()
    #         for lag_idx in range(n_lags):
    #             ax = axes_flat[lag_idx]
    #             ax.imshow(
    #                 strf_2d[roi_ind, :, :, lag_idx],
    #                 aspect='auto',
    #                 origin='lower',
    #                 cmap='RdBu_r',
    #                 vmin=-1, vmax=1,
    #             )
    #             ax.set_title(f'{t_lag[lag_idx]*1000:.0f} ms', fontsize=8)
    #             ax.set_xticks([])
    #             ax.set_yticks([])
    #         for ax in axes_flat[n_lags:]:
    #             ax.set_visible(False)
    #         plt.suptitle(f'2D STRF all lags  {sn_h}+{sn_v}  {ch}  noisetau={noisetau}ms  '
    #                      f'rtz=({theta},{partner})  roi {roi_ind}  {tag}')
    #         if save_figs:
    #             plt.savefig(os.path.join(figs_dir, fig_name + fig_format),
    #                         dpi=150, transparent=True)
    #         if show_figs:
    #             plt.show()
    #         plt.close()

    # figure 5: movie of the 1D STRF across every lag. Commented out for now
    # -- only the 2D movie (figure 6) is wanted.
    # for (sn, ch, orientation, noisetau), res in strf_results.items():
    #     strf  = res['strf']   # (n_roi, n_positions, n_lags)
    #     t_lag = res['t_lag']
    #     rot   = res['rtz']
    #     n_roi_1d, n_positions, n_lags = strf.shape
    #     vmax = np.abs(strf).max() or 1.0
    #
    #     if not save_figs:
    #         continue
    #
    #     for roi_ind in range(n_roi_1d):
    #         fig_name = f'strf_1d_movie_rot{int(rot)}_noisetau{int(noisetau)}_{ch}_{sn}_{tag}_roi_{roi_ind}_'
    #         fh, ax = plt.subplots(figsize=(5, 4), constrained_layout=True)
    #         positions = np.arange(n_positions)
    #         bars = ax.bar(positions, strf[roi_ind, :, 0], color='k')
    #         ax.set_ylim(-vmax, vmax)
    #         ax.set_xlabel('Position')
    #         ax.set_ylabel('STRF (z-scored)')
    #         title = ax.set_title('')
    #
    #         def _update_1d(lag_idx, _strf=strf, _roi_ind=roi_ind, _t_lag=t_lag,
    #                         _bars=bars, _title=title):
    #             values = _strf[_roi_ind, :, lag_idx]
    #             for bar, v in zip(_bars, values):
    #                 bar.set_height(v)
    #             _title.set_text(f'{sn}  {ch}  roi {_roi_ind}  lag={_t_lag[lag_idx]*1000:.1f} ms')
    #             return list(_bars) + [_title]
    #
    #         ani = animation.FuncAnimation(fh, _update_1d, frames=n_lags, blit=False)
    #         movie_path = os.path.join(figs_dir, fig_name + '.mp4')
    #         ani.save(movie_path, writer='ffmpeg', fps=movie_fps, dpi=150)
    #         plt.close(fh)

    # figure 6: movie of the 2D STRF across every lag -- one video per ROI,
    # per series pair/channel/noisetau/angle pair (separate files, not one
    # combined grid, so playback is readable at full size per ROI)
    for (sn_h, sn_v, ch, noisetau, theta, partner), res in strf_2d_results.items():
        strf_2d = res['strf_2d']   # (n_roi, n_pos_h, n_pos_v, n_lags)
        t_lag   = res['t_lag']
        n_roi_2d, _, _, n_lags = strf_2d.shape
        pos_h_deg = res['bar_pos_h_deg']
        pos_v_deg = res['bar_pos_v_deg']
        h_label = _pole_label(theta)
        v_label = _pole_label(partner)
        extent_2d = [pos_v_deg[0], pos_v_deg[-1], pos_h_deg[0], pos_h_deg[-1]]

        if not save_figs:
            continue

        for roi_ind in range(n_roi_2d):
            fig_name = f'strf_2d_movie_noisetau{int(noisetau)}_rtz{int(theta)}_{int(partner)}_{sn_h}_{sn_v}_{ch}_{tag}_roi_{roi_ind}_'
            # single color scale across all lags, so a weak lag doesn't get
            # rendered as if it were as strong as the peak lag
            vmax = np.abs(strf_2d[roi_ind]).max() or 1.0
            fh, ax = plt.subplots(figsize=_panel_size(extent_2d, 4.2),
                                  constrained_layout=True)
            im = ax.imshow(
                strf_2d[roi_ind, :, :, 0],
                aspect='equal',
                origin='lower',
                extent=extent_2d,
                cmap='RdBu_r',
                vmin=-vmax, vmax=vmax,
            )
            ax.set_xlabel(f'{v_label}  (deg, 90-psi)')
            ax.set_ylabel(f'{h_label}  (deg, 90-psi)')
            title = ax.set_title('')

            def _update(lag_idx, _strf_2d=strf_2d, _roi_ind=roi_ind, _t_lag=t_lag,
                        _im=im, _title=title):
                _im.set_data(_strf_2d[_roi_ind, :, :, lag_idx])
                _title.set_text(f'{sn_h}+{sn_v}  {ch}  roi {_roi_ind}  '
                                 f'lag={_t_lag[lag_idx]*1000:.1f} ms')
                return _im, _title

            ani = animation.FuncAnimation(fh, _update, frames=n_lags, blit=False)
            movie_path = os.path.join(figs_dir, fig_name + '.mp4')
            ani.save(movie_path, writer='ffmpeg', fps=movie_fps, dpi=150)
            plt.close(fh)
