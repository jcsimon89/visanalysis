"""
analyze_data_trf_fullfield.py

Per-ROI TEMPORAL filters from the full-field ternary noise stimulus
(labpack FullFieldTernaryNoise).

Call structure (same as analyze_data_strf_octave.py):
    python analyze_data_trf_fullfield.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --show_figs False \
        --save_figs True \
        --tag final \
        --filter_length 10.0

Results saved to fly.hdf5 (or fly_final.hdf5) under group /TRF_FULLFIELD_{tag}.

HOW THIS DIFFERS FROM THE OCTAVE STRF PIPELINE
----------------------------------------------
* No space. The whole field takes one value, so every ROI sees the identical
  signal and the answer is one filter per ROI rather than a filter per cell per
  ROI. That is the entire point of the stimulus: the per-cell SNR division is
  what makes a spatiotemporal STRF need a long recording, and if you only want
  the temporal filter you should not pay it.

* REGRESSION, not a spike-triggered average. The octave path takes an STA and
  then has to undo the stimulus autocorrelation, because holding rather than
  refreshing convolves the answer with exp(-|lag|/tau). Here the filter is
  solved from the normal equations, so the empirical stimulus autocorrelation is
  divided out as part of the fit. There is nothing to deconvolve and no
  --deconvolve flag: the peak lag is already the peak lag.

* The ridge is therefore the ONLY regularisation, and it is cross-validated
  rather than assumed. It is chosen once per series and applied to every ROI --
  per-ROI selection would let each ROI pick the smoothing that flatters it, and
  the differences between ROIs are what the experiment is for.

* No movie, and no dropped-frame repair. Every value is a pure function of
  (seed, band, frame), so the stimulus is regenerated from the recorded
  `fullfield_stimulus` metadata. What the metadata cannot supply is WHICH frames
  rendered; this assumes the nominal mapping frame = round(t * fps) and reports
  the photodiode's dropped-frame count so a drift can be seen rather than
  silently absorbed.

WHAT TO READ FIRST, AND WHY IN THAT ORDER
-----------------------------------------
cv_r -- the correlation between held-out response and its prediction -- comes
before everything else. A filter that does not predict held-out data is not a
filter, and its peak time, gain and shape are then descriptions of noise. The
circular-shift null is the second gate and answers a different question: not
"how well does it predict" but "is this distinguishable from no coupling at
all". Both are per ROI, and both are saved.

MULTI-BAND SERIES. The regressor is the DISPLAYED contrast, which is what the
cell saw, so a multi-band series yields one filter for the composite. That is
correct but it is not the multi-band experiment: several bands exist to ask what
MODULATES the filter, which means conditioning the fast band's estimate on the
slow band's state. That analysis is not implemented here. --band fits one band
alone as a diagnostic.
"""

import sys
import os
import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
from visanalysis.util.strf_plot import apply_figure_style
import h5py

from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.util import bleach, nullcal
from visanalysis.util import fullfield_temporal as ff

plt.ioff()
apply_figure_style()


def _summarise(h, lags, dt):
    """The per-ROI scalars, all derived from the filter itself."""
    n_roi = h.shape[0]
    peak_ind = np.array([int(np.nanargmax(np.abs(x))) if np.any(np.isfinite(x))
                         else 0 for x in h])
    out = dict(
        peak_lag_s=lags[peak_ind],
        peak_signed=np.array([h[r, peak_ind[r]] for r in range(n_roi)]),
        # DC gain: dF/F per unit SUSTAINED contrast. Signed on purpose -- an ON
        # and an OFF cell differ by its sign, and |.| would hide that.
        dc_gain=np.nansum(h, axis=1),
        integration_index=np.array([ff.integration_index(x) for x in h]),
    )
    # Does the filter decay inside its window? This, not the window length, is
    # the condition under which the low frequencies mean anything.
    rz = [bleach.returns_to_zero(x, dt) for x in h]
    out['returns_to_zero'] = np.array([bool(a) for a, _ in rz])
    out['tail_fraction'] = np.array([b for _, b in rz])
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_file_directory', nargs='?')
    parser.add_argument('--rig', nargs='?', default='Bruker')
    parser.add_argument('--show_figs', nargs='?', default='False')
    parser.add_argument('--save_figs', nargs='?', default='False')
    parser.add_argument('--tag', nargs='?', default='raw')
    parser.add_argument('--experiment_file_name', nargs='?', default=None,
                        help="hdf5 to read. Default: fly_<tag>.hdf5 when that "
                             "file exists (a recording whose rois were selected "
                             "from a search stimulus), else fly.hdf5 -- so --tag "
                             "can be a plain label on flies that never had one.")
    parser.add_argument('--filter_length', nargs='?', type=float, default=10.0,
                        help='filter length in seconds. Must cover the slowest '
                             'kernel you expect: cAMP indicators still carry '
                             'power past 3 s.')
    parser.add_argument('--block_s', nargs='?', type=float, default=60.0,
                        help='cross-validation block length in seconds. Folds '
                             'are contiguous blocks rather than epochs, because '
                             'the protocol ships one trial and a fold per epoch '
                             'would leave nothing to hold out.')
    parser.add_argument('--n_fold', nargs='?', type=int, default=6)
    parser.add_argument('--ridge', nargs='?', type=float, default=None,
                        help='fix the ridge instead of cross-validating it. '
                             'Leave unset -- the right value depends on the hold '
                             'time and the indicator, so a number that suits one '
                             'recording will not suit the next.')
    parser.add_argument('--band', nargs='?', type=int, default=None,
                        help='fit against ONE band of a multi-band stimulus '
                             'instead of the displayed composite. Diagnostic '
                             'only; the composite is what the cell saw.')
    parser.add_argument('--bleach', nargs='?', default='auto',
                        help="bleach correction on raw F: 'auto' (default), "
                             "'on', 'off'")
    parser.add_argument('--peak_window', nargs='?', type=float,
                        default=nullcal.DEFAULT_PEAK_WINDOW_S,
                        help='lag window, in seconds from lag 0, searched for the '
                             'peak |z| -- on the rois and on every null draw alike; '
                             '0 searches the whole filter. See util.nullcal.')
    parser.add_argument('--null_shifts', nargs='?', type=int, default=20,
                        help='circular-shift null draws; 0 disables. Costs '
                             '(n+1)x the fit. See util.nullcal. Default 20 '
                             'because min_detectable_k ~= 1/(q*n_shifts), so '
                             'below 20 a single responsive ROI cannot be '
                             'resolved at q=0.05 however strong it is -- and '
                             'between 5 and 7 shifts FDR returns an EMPTY set '
                             'on real data for that reason, which reads as '
                             '"nothing is real" when it means "too few '
                             'shifts".')
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    show_figs = args.show_figs == 'True'
    save_figs = args.save_figs == 'True'
    tag = args.tag
    bleach_mode = {'auto': 'auto', 'on': True, 'off': False}[str(args.bleach).lower()]
    filter_length_s = float(args.filter_length)
    null_shifts = int(args.null_shifts)
    peak_window_s = float(args.peak_window)

    # Cross-validation scores only rows whose whole lag window lies inside their
    # own block, so a block no longer than the filter leaves nothing to score
    # and the ridge would be chosen from an empty set. Caught here rather than
    # after the fits, because it is a configuration error and it is silent.
    if float(args.block_s) <= filter_length_s:
        print('block_s ({:g} s) must exceed filter_length ({:g} s), or no row '
              'can be scored: every row in a block has a lag window reaching '
              'into the previous one.'.format(args.block_s, filter_length_s))
        sys.exit(1)
    if float(args.block_s) < 2 * filter_length_s:
        print('WARNING: block_s ({:g} s) is less than twice filter_length '
              '({:g} s), so most of each block is spent on the guard and the '
              'cross-validation is scored on very little.'
              .format(args.block_s, filter_length_s))

    # --tag is a LABEL for the output -- the hdf5 group and the figure folder --
    # not a filename. It used to be both: 'raw' meant fly.hdf5 and anything else
    # meant fly_<tag>.hdf5, the filtered copy select_final_rois.py writes after
    # picking responsive rois from a search stimulus. That step is being retired,
    # because responsive rois are now identified from the noise itself against a
    # calibrated null, which needs no separate stimulus and does not require the
    # receptive field to sit near the middle of the screen.
    #
    # Recordings that DO have a search stim still have their fly_<tag>.hdf5, and
    # must keep working, so the tagged file is used WHEN IT EXISTS and fly.hdf5
    # otherwise. That leaves --tag free to be a plain label like 'mi1_1d' or
    # 'lag60' on a fly that never had one, which is also how two analyses can sit
    # side by side in one file instead of overwriting each other.
    if getattr(args, 'experiment_file_name', None):
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
    print('reading ' + experiment_file_name)

    if rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        plug = base_plugin.BasePlugin()

    response_set_name_prefix = 'mask_'
    func_channels_num = ['1', '2']

    # ── fly metadata and photodiode channel ──────────────────────────────────
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))
    ID = imaging_data.ImagingDataObject(experiment_file_path,
                                        int(series_num[0]), quiet=True)
    fly_metadata = ID.getSubjectMetadata()
    prep = fly_metadata.get('prep', '')
    timing_channel_ind = 1 if prep == 'fly right optic lobe' else 0
    print('fly_metadata: ' + repr(fly_metadata))
    print('timing_channel_ind: ' + repr(timing_channel_ind))

    # ── find the full-field series ───────────────────────────────────────────
    fullfield_series = []
    run_parameters = {}
    for current_series in series_num:
        sn = 'sn' + current_series
        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name, int(current_series))
        ID = plug.ImagingDataObject
        run_parameters[sn] = ID.getRunParameters()
        if ff.is_fullfield_ternary(run_parameters[sn]):
            fullfield_series.append(current_series)

    if not fullfield_series:
        print('No full-field ternary series identified. Check '
              'is_fullfield_ternary() and that process_data.py attached the '
              'epoch parameters for this series.')
        sys.exit(0)

    print('\nFound {} full-field ternary series: {}'.format(
        len(fullfield_series), ', '.join(fullfield_series)))

    # ── filters ──────────────────────────────────────────────────────────────
    results = {}

    for current_series in fullfield_series:
        sn = 'sn' + current_series
        plug.updateImagingDataObject(experiment_file_directory,
                                     experiment_file_name, int(current_series))
        ID = plug.ImagingDataObject
        # BEFORE any response call: getRoiResponses splits epochs on the
        # photodiode, so the timing channel has to be set first or every epoch
        # is cut on the wrong trace.
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6
        _, _, sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * sample_rate

        rp = run_parameters[sn]
        duration_s = (float(rp.get('pre_time', 0.0))
                      + float(rp.get('stim_time', 0.0))
                      + float(rp.get('tail_time', 0.0)))

        # Dropped display frames are NOT repaired. Every stimulus value is a
        # pure function of its frame index, so a drop shifts the time base of
        # everything after it by 1/fps -- it does not corrupt a value. Reported
        # rather than corrected, because only the display machine knows which
        # index was skipped; if the cumulative shift approaches an imaging
        # sample the lag axis is what pays for it.
        try:
            st = ID.getStimulusTiming()
            n_drop = len(st['dropped_frame_times'])
        except Exception as e:
            n_drop = None
            print('  {}: could not read stimulus timing ({})'.format(sn, e))
        if n_drop:
            print('  {}: {} dropped monitor frames -> up to {:.3f} s of '
                  'cumulative time-base shift'.format(
                      sn, n_drop, n_drop / float(rp.get('fps', 120.0))))

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            _t0 = time.perf_counter()
            try:
                design = ff.build_roi_design(
                    ID, response_set_name_prefix + ch,
                    filter_length_s, roi_prefix='aligned', band=args.band,
                    correct_bleach=bleach_mode, block_s=float(args.block_s))
            # AssertionError: getRoiResponses asserts the ROI set exists, so a
            # fly with one functional channel would otherwise crash here
            # instead of skipping the channel it does not have.
            except (KeyError, ValueError, OSError, AssertionError) as e:
                print('  {} {}: skipped ({})'.format(sn, ch, e))
                continue
            print('    [timing] build_roi_design took {:.2f}s'.format(
                time.perf_counter() - _t0))

            meta0 = design['metas'][0]
            lags = np.arange(design['n_lags']) * design['dt']
            print('  {} {}: n_roi={}, {} epochs, {} rows, {} lags, imaging '
                  '{:.2f} Hz, {} z plane(s), {} block(s)'.format(
                      sn, ch, design['n_roi'], design['n_ep'],
                      design['Y'].shape[0], design['n_lags'],
                      1.0 / design['dt'], len(design['z_offsets']),
                      len(np.unique(design['fold']))))
            print('       taus={} s  clip_k={}  p={:.6f}  fps={:g}  seeds={}{}'
                  .format(tuple(meta0['hold_taus']), meta0['clip_k'],
                          meta0['p'], meta0['fps'], design['seeds'],
                          '  [repeat_stimulus]' if design['repeated'] else ''))
            if design['n_bands'] > 1 and args.band is None:
                print('       NOTE: {} bands. The filter below is for the '
                      'DISPLAYED composite, which is what the cell saw. The '
                      'modulation analysis -- the fast band conditioned on the '
                      'slow band\'s state -- is not implemented here.'
                      .format(design['n_bands']))
            if design['n_rows_dropped']:
                print('       {} rows dropped for non-finite responses'
                      .format(design['n_rows_dropped']))

            # ── ridge ────────────────────────────────────────────────────────
            _t0 = time.perf_counter()
            grid = None if args.ridge is None else np.array([args.ridge])
            grid, scores = ff.ridge_sweep(design, grid=grid,
                                          n_fold=int(args.n_fold))
            med = np.nanmedian(scores, axis=1)
            cv_failed = not np.any(np.isfinite(med))
            # Middle of the grid, NOT grid[0], when nothing could be scored.
            # The low end is the least regularised fit, which is the worst thing
            # to fall back on precisely when there is no evidence to justify it.
            best_i = (len(grid) // 2 if cv_failed else int(np.nanargmax(med)))
            ridge = float(grid[best_i])
            cv_r = scores[best_i]
            print('    [timing] ridge_sweep took {:.2f}s'.format(
                time.perf_counter() - _t0))
            if cv_failed:
                print('       WARNING: cross-validation produced no score for '
                      'any ROI, so the ridge fell back to the middle of the '
                      'grid ({:.3g}) and is NOT chosen from this data. Usual '
                      'cause: too few blocks ({} here) -- shorten --block_s or '
                      'record longer. Treat every filter below as unvalidated.'
                      .format(ridge, len(np.unique(design['fold']))))
            elif int(np.sum(np.isfinite(cv_r))) < design['n_roi']:
                print('       {} of {} ROIs have no cv_r (flat or all-NaN '
                      'response); they are still fitted and still tested '
                      'against the null.'.format(
                          int(np.sum(~np.isfinite(cv_r))), design['n_roi']))
            if cv_failed:
                pass                                   # already reported above
            elif args.ridge is not None:
                print('       ridge {:.3g} (fixed), median cv_r {:.3f}'
                      .format(ridge, med[best_i]))
            else:
                print('       ridge {:.3g} (cross-validated), median cv_r '
                      '{:.3f}'.format(ridge, med[best_i]))
                if best_i in (0, len(grid) - 1):
                    print('       WARNING: the ridge landed on the edge of its '
                          'grid, so it was not really chosen. Widen the grid, '
                          'or read this as "the data cannot constrain it".')

            # ── filters ──────────────────────────────────────────────────────
            h = ff.fit_rois(design, ridge)
            split_half = ff.split_half_rois(design, ridge)
            summary = _summarise(h, lags, design['dt'])

            # ── null: TWO gates, one set of shifts ───────────────────────────
            # Same estimator, same ridge, stimulus rolled. A null from a
            # different estimator would not be the right reference.
            #
            # SHAPE gate: is the filter's peak distinguishable from one the
            # estimator manufactures by chance? PREDICTION gate: does it predict
            # held-out blocks better than chance? These disagree, and the
            # disagreement is not academic -- on fly_014 ch1, ROIs 0 and 1
            # passed the shape gate at every shift count with cv_r of 0.024 and
            # 0.026, and went into the population mean and the -3 dB band
            # alongside two ROIs at 0.27 and 0.33. A filter that does not
            # predict is not a filter, so inclusion now requires both.
            #
            # The prediction gate is NULL-CALIBRATED rather than a fixed cv_r
            # floor, on purpose. cAMP filters are slower and weaker than the
            # calcium filters measured in the same fly -- here ch2 tops out near
            # cv_r 0.06 where ch1 reaches 0.33 -- so any floor tuned on calcium
            # silently discards the cAMP channel, and the comparison between
            # those two channels is the experiment. Chance is measured per
            # channel instead of assumed.
            null_p = np.full(design['n_roi'], np.nan)
            cv_p = np.full(design['n_roi'], np.nan)
            null_thr, keep, null_info = np.nan, np.array([], int), {}
            keep_shape = keep_cv = np.array([], int)
            cv_null_thr = np.nan
            cv_method = 'not run'
            # Lags searched for the peak, rois and null alike; the noise sd
            # still comes from the late lags of the full axis. See
            # nullcal.DEFAULT_PEAK_WINDOW_S.
            _pl = nullcal.peak_lags_for(peak_window_s, design['dt'],
                                        design['n_lags'])
            if null_shifts > 0:
                _t0 = time.perf_counter()
                try:
                    # One set of shifts for both gates, so they are calibrated
                    # against the same surrogate data.
                    shifts = nullcal.shifts_for(duration_s, filter_length_s,
                                                null_shifts)
                    cal = nullcal.calibrate(
                        lambda s: ff.fit_rois(design, ridge, shift_s=s), shifts,
                        peak_lags=_pl)
                    z = nullcal.peak_z(h, peak_lags=_pl)
                    ok = np.isfinite(z)
                    null_p[ok] = nullcal.pvalues(z[ok], cal['null'])[0]
                    null_thr = cal['threshold']
                    keep_shape, _signs, null_info = nullcal.select_rois(h, cal)

                    # Prediction null: cv_r at the CHOSEN ridge on each rolled
                    # stimulus. Pooled over ROIs and shifts exactly as the shape
                    # null is, so the two carry the same resolution limit.
                    cv_null = np.concatenate([
                        ff.ridge_sweep(design, grid=np.array([ridge]),
                                       n_fold=int(args.n_fold), shift_s=s)[1][0]
                        for s in shifts])
                    cv_null = cv_null[np.isfinite(cv_null)]
                    ok_cv = np.isfinite(cv_r)
                    if cv_null.size and ok_cv.any():
                        cv_null_thr = float(np.percentile(cv_null, 95))
                        cv_p[ok_cv] = nullcal.pvalues(cv_r[ok_cv], cv_null)[0]
                        # Same rule select_rois applies to the shape gate: FDR
                        # only when the shift count can resolve at least one
                        # ROI, the percentile otherwise. Without this, too few
                        # shifts put the p floor above q, BH rejects nothing,
                        # and the gate empties silently -- which is what an
                        # empty keep and an all-ROI 'mean' came from on fly_014
                        # when the shell wrapper still passed 4 shifts.
                        cv_resolvable = int(np.ceil(1.0 / (0.05 * len(shifts))))
                        if cv_resolvable <= int(ok_cv.sum()):
                            sel_cv, _cut = nullcal.benjamini_hochberg(
                                cv_p[ok_cv], q=0.05)
                            keep_cv = (np.flatnonzero(ok_cv)[sel_cv]
                                       if len(sel_cv) else np.array([], int))
                            cv_method = 'fdr'
                        else:
                            keep_cv = np.flatnonzero(ok_cv & (cv_r > cv_null_thr))
                            cv_method = ('percentile (min_detectable_k={} > '
                                         'n_roi={})'.format(cv_resolvable,
                                                            int(ok_cv.sum())))
                    keep = np.intersect1d(keep_shape, keep_cv)

                    print('    [timing] null ({} shifts, both gates) took '
                          '{:.2f}s'.format(null_shifts,
                                           time.perf_counter() - _t0))
                    print('       shape gate: {} of {} ({}); null 95th '
                          'percentile z = {:.2f}; reported threshold {:.2f} '
                          '({})'.format(
                              len(keep_shape), design['n_roi'],
                              null_info.get('method_used', ''), null_thr,
                              null_info.get('threshold', np.nan),
                              null_info.get('threshold_kind', 'unspecified')))
                    print('       prediction gate: {} of {} ({}; null 95th '
                          'percentile cv_r = {:.3f} from {} draws)'.format(
                              len(keep_cv), design['n_roi'], cv_method,
                              cv_null_thr, cv_null.size))
                    print('       KEPT (both gates): {} of {} -> {}'.format(
                        len(keep), design['n_roi'], list(keep)))
                    if len(keep_shape) and not len(keep):
                        print('       NOTE: every ROI that passed the shape '
                              'gate failed the prediction gate. Filters with a '
                              'convincing peak that predict nothing are what '
                              'this combination exists to catch.')
                    if 'reason' in null_info:
                        print('       fallback: ' + str(null_info['reason']))
                    if 'warning' in null_info:
                        print('       ' + null_info['warning'])
                except ValueError as e:
                    print('    null skipped: {}'.format(e))

            # ── spectrum ─────────────────────────────────────────────────────
            # Of the mean over SURVIVING ROIs, unaligned. Do not sign-align each
            # ROI to its own peak first: that makes noise add constructively
            # there, so the average acquires a peak whether or not one is real.
            sel = keep if len(keep) else np.arange(design['n_roi'])
            if not len(keep):
                print('       WARNING: no ROI passed both gates, so the '
                      'population mean and the -3 dB band below are over ALL '
                      '{} ROIs and are NOT a selected population.'
                      .format(design['n_roi']))
            h_mean = np.nanmean(h[sel], axis=0)
            # Saved, not only plotted: this is the per-fly population filter
            # a cross-fly average is built from.
            h_sem = (np.nanstd(h[sel], axis=0) / np.sqrt(len(sel))
                     if len(sel) > 1 else np.full(h.shape[1], np.nan))
            # Not sign-aligning is right, but it means a population that
            # disagrees in sign averages toward zero. nullcal's rule is to
            # split such a population and report the groups; at minimum it
            # must be visible, because the mean then describes neither group.
            kept_signs = np.sign(summary['peak_signed'][keep]).astype(int)
            mixed = bool(len(np.unique(kept_signs)) > 1)
            if mixed:
                print('       WARNING: kept ROIs disagree in the sign of '
                      'their peak {}, so the mean filter partly CANCELS: its '
                      'amplitude, shape and -3 dB band describe neither '
                      'group. Split the population, or find out which '
                      'component sets each ROI\'s peak.'.format(
                          {int(r): int(v) for r, v in zip(keep, kept_signs)}))
            f_ax, mag = ff.filter_spectrum(h_mean, design['dt'])
            f_lo, f_hi, f_pk = ff.band_3db(f_ax, mag)
            print('       population -3 dB band {:.3f} - {:.3f} Hz '
                  '(peak {:.3f} Hz) from {} ROI(s)'.format(
                      f_lo, f_hi, f_pk, len(sel)))
            if design['bleach_info'].get('applied'):
                print('       bleach correction high-passes the response: it is '
                      '0.99 at 0.02 Hz but 0.71 at 0.01 Hz, so read the low edge '
                      'as a floor near 0.01 Hz rather than as a measurement.')
            n_flat = int(np.sum(~summary['returns_to_zero']))
            if n_flat:
                print('       {} of {} filters have not returned to zero by '
                      '{:.0f} s. Their low frequencies are not resolved -- '
                      'raise --filter_length.'.format(
                          n_flat, design['n_roi'], filter_length_s))

            # ── trial by trial: repeated stimuli only ────────────────────────
            # With the same stimulus on every trial, differences between
            # per-trial filters are the cell's own -- drift, run-down, state --
            # not a different stimulus sampling the filter differently. So this
            # runs when a stimulus repeats and is skipped when none does.
            trial_rel = ff.trial_reliability(design, ridge, h)
            if trial_rel is None:
                print('       no repeated stimulus: trial-by-trial analysis '
                      'skipped')
            else:
                show = keep if len(keep) else np.arange(design['n_roi'])
                fmt = lambda a: '[' + ' '.join('{:.2f}'.format(x) for x in a) + ']'
                for grp in trial_rel:
                    print('       repeated stimulus, seed {}: {} trials (epochs '
                          '{})'.format(grp['seed'], len(grp['epochs']),
                                       list(grp['epochs'])))
                    for r in show:
                        nc = grp['noise_ceiling'][r]
                        print('         roi {}: filter trial-vs-rest r {}  gain {}'
                              .format(r, fmt(grp['filter_loo_r'][:, r]),
                                      fmt(grp['trial_gain'][:, r])))
                        print('                response repeatability {:.3f} +- '
                              '{:.3f}; noise ceiling {}; cv_r / ceiling {}'.format(
                                  grp['response_pair_r'][r],
                                  grp['response_pair_r_se'][r],
                                  'unresolved' if not np.isfinite(nc)
                                  else '{:.3f}'.format(nc),
                                  '-' if not np.isfinite(nc)
                                  else '{:.2f}'.format(cv_r[r] / nc)))

            results[sn, ch] = dict(
                filter=h, lags=lags, cv_r=cv_r, split_half_r=split_half,
                ridge=ridge, ridge_grid=grid, ridge_scores=scores,
                null_p=null_p, null_threshold=null_thr, cv_p=cv_p,
                cv_null_threshold=cv_null_thr, keep_shape=keep_shape,
                keep_cv=keep_cv, mean_filter=h_mean, mean_filter_sem=h_sem,
                mean_rois=sel, mean_is_selected=bool(len(keep)),
                kept_signs=kept_signs, mixed_polarity=mixed,
                trial_rel=trial_rel,
                null_peak_z=nullcal.peak_z(h, peak_lags=_pl), keep=keep,
                null_info=null_info,
                design=design, meta=meta0, freqs=f_ax, mag_mean=mag,
                band_3db=(f_lo, f_hi, f_pk), duration_s=duration_s,
                n_dropped_frames=(-1 if n_drop is None else int(n_drop)),
                **summary)

    if not results:
        print('\nNothing fitted.')
        sys.exit(0)

    # ── save ─────────────────────────────────────────────────────────────────
    group_name = 'TRF_FULLFIELD_' + tag
    with h5py.File(experiment_file_path, 'a') as f:
        if group_name in f:
            del f[group_name]
        g = f.create_group(group_name)
        g.attrs['filter_length_s'] = filter_length_s
        g.attrs['block_s'] = float(args.block_s)
        g.attrs['band'] = 'composite' if args.band is None else int(args.band)
        g.attrs['bleach'] = str(args.bleach)
        g.attrs['null_shifts'] = null_shifts

        for (sn, ch), res in results.items():
            d = res['design']
            sg = g.require_group('{}/{}'.format(sn, ch))
            sg.create_dataset('filter', data=res['filter'], compression='gzip')
            sg.create_dataset('lags', data=res['lags'])
            for key in ('cv_r', 'split_half_r', 'peak_lag_s', 'peak_signed',
                        'dc_gain', 'integration_index', 'null_p', 'cv_p',
                        'null_peak_z', 'returns_to_zero', 'tail_fraction'):
                sg.create_dataset(key, data=res[key])
            # kept_rois is the INTERSECTION; the per-gate sets are kept beside
            # it so a reader can see which gate excluded what.
            sg.create_dataset('kept_rois', data=np.asarray(res['keep'], dtype=int))
            sg.create_dataset('kept_rois_shape_gate',
                              data=np.asarray(res['keep_shape'], dtype=int))
            sg.create_dataset('kept_rois_prediction_gate',
                              data=np.asarray(res['keep_cv'], dtype=int))
            # Per-fly population filter. mean_filter_is_selected is False when
            # nothing passed both gates and the mean fell back to ALL ROIs --
            # never pool such a fly into a cross-fly average unflagged.
            sg.create_dataset('mean_filter', data=res['mean_filter'])
            sg.create_dataset('mean_filter_sem', data=res['mean_filter_sem'])
            sg.create_dataset('mean_filter_rois',
                              data=np.asarray(res['mean_rois'], dtype=int))
            sg.create_dataset('kept_peak_sign',
                              data=np.asarray(res['kept_signs'], dtype=int))
            sg.attrs['mean_filter_is_selected'] = res['mean_is_selected']
            sg.attrs['kept_mixed_polarity'] = res['mixed_polarity']
            if res['trial_rel'] is None:
                sg.attrs['trial_reliability'] = 'skipped: no repeated stimulus'
            else:
                sg.attrs['trial_reliability'] = '{} repeat group(s)'.format(
                    len(res['trial_rel']))
                for gi, grp in enumerate(res['trial_rel']):
                    tg = sg.create_group('trial_reliability/group_{}'.format(gi))
                    tg.attrs['seed'] = grp['seed']
                    tg.attrs['definitions'] = (
                        'filter_loo_r: each trial filter vs mean of the others, '
                        'over support_lags; trial_gain: LS scale onto the full '
                        'filter; response_pair_r: mean corr between single-trial '
                        'responses (model-free); noise_ceiling: '
                        'sqrt(response_pair_r), NaN unless > 2 SE')
                    tg.create_dataset('trial_filters', data=grp['trial_filters'],
                                      compression='gzip')
                    for key in ('epochs', 'support_lags', 'filter_loo_r',
                                'trial_gain', 'response_pair_r',
                                'response_pair_r_se', 'response_loo_r',
                                'noise_ceiling'):
                        tg.create_dataset(key, data=np.asarray(grp[key]))
            sg.create_dataset('roi_z_offset_s',
                              data=d['z_offsets'][d['roi_group']])
            sg.create_dataset('ridge_grid', data=res['ridge_grid'])
            sg.create_dataset('ridge_cv_r', data=res['ridge_scores'],
                              compression='gzip')
            sg.create_dataset('freqs', data=res['freqs'])
            sg.create_dataset('mag_mean', data=res['mag_mean'])

            # What the filter IS, so a reader does not have to guess its units
            # or which normalisation produced it.
            sg.attrs['filter_units'] = ('dF/F per unit displayed contrast, per '
                                        'imaging sample; divide by dt for a '
                                        'continuous kernel, sum for DC gain')
            sg.attrs['dff'] = ('per-ROI, per-epoch median baseline, after '
                                'bleach correction when bleach_applied')
            sg.attrs['ridge'] = res['ridge']
            sg.attrs['ridge_is_fixed'] = args.ridge is not None
            sg.attrs['dt'] = d['dt']
            sg.attrs['peak_window_s'] = (float(peak_window_s)
                                         if peak_window_s > 0 else 0.0)
            sg.attrs['n_epochs'] = d['n_ep']
            sg.attrs['n_rows'] = int(d['Y'].shape[0])
            sg.attrs['n_rows_dropped'] = d['n_rows_dropped']
            sg.attrs['n_blocks'] = int(len(np.unique(d['fold'])))
            sg.attrs['repeat_stimulus'] = bool(d['repeated'])
            sg.attrs['seeds'] = np.asarray(d['seeds'], dtype=np.int64)
            sg.attrs['hold_taus'] = np.asarray(res['meta']['hold_taus'],
                                               dtype=float)
            sg.attrs['clip_k'] = int(res['meta']['clip_k'])
            sg.attrs['p'] = float(res['meta']['p'])
            sg.attrs['fps'] = float(res['meta']['fps'])
            sg.attrs['n_bands'] = int(d['n_bands'])
            sg.attrs['null_threshold'] = res['null_threshold']
            sg.attrs['cv_null_threshold'] = res['cv_null_threshold']
            sg.attrs['inclusion_rule'] = ('shape gate (peak z vs circular-shift '
                                          'null) AND prediction gate (cv_r vs '
                                          'the same null)')
            sg.attrs['band_3db_lo_hz'] = res['band_3db'][0]
            sg.attrs['band_3db_hi_hz'] = res['band_3db'][1]
            sg.attrs['band_3db_peak_hz'] = res['band_3db'][2]
            sg.attrs['n_dropped_frames'] = res['n_dropped_frames']
            sg.attrs['bleach_applied'] = bool(
                d['bleach_info'].get('applied', False))
            # 'null_select_', not 'null_': select_rois reports its own
            # 'threshold' (the weakest ROI it kept), which under a shared prefix
            # silently overwrites the null threshold written above.
            for k, v in res['null_info'].items():
                sg.attrs['null_select_' + k] = v

    print('\nSaved temporal filters to {} under /{}'.format(
        experiment_file_path, group_name))

    # ── summary table ────────────────────────────────────────────────────────
    print('\n{:<10} {:>5} {:>7} {:>8} {:>9} {:>9} {:>9}'.format(
        'series/ch', 'n_roi', 'kept', 'cv_r', 'peak lag', 'DC gain', 'band Hz'))
    for (sn, ch), res in results.items():
        sel = res['keep'] if len(res['keep']) else np.arange(len(res['cv_r']))
        print('{:<10} {:>5} {:>7} {:>8.3f} {:>8.2f}s {:>9.3f} {:>4.2f}-{:.2f}'
              .format('{}/{}'.format(sn, ch), len(res['cv_r']),
                      len(res['keep']), float(np.nanmedian(res['cv_r'][sel])),
                      float(np.nanmedian(res['peak_lag_s'][sel])),
                      float(np.nanmedian(res['dc_gain'][sel])),
                      res['band_3db'][0], res['band_3db'][1]))

    # ── figures ──────────────────────────────────────────────────────────────
    if save_figs or show_figs:
        figs_dir = os.path.join(experiment_file_directory, tag + '_roi_figs')
        os.makedirs(figs_dir, exist_ok=True)

        for (sn, ch), res in results.items():
            h, lags, cv_r = res['filter'], res['lags'], res['cv_r']
            keep = res['keep']
            n_roi = h.shape[0]
            in_keep = np.zeros(n_roi, dtype=bool)
            in_keep[keep] = True

            fig, axes = plt.subplots(2, 3, figsize=(15, 7.5))

            ax = axes[0][0]
            for r in range(n_roi):
                ax.plot(lags, h[r], lw=0.7, alpha=0.35,
                        color='C0' if in_keep[r] else '0.75')
            ax.axhline(0, color='0.6', lw=0.8)
            ax.set_xlabel('lag (s)')
            ax.set_ylabel('filter (dF/F per unit contrast)')
            ax.set_title('all {} ROIs (grey = below null)'.format(n_roi))

            ax = axes[0][1]
            sel = keep if len(keep) else np.arange(n_roi)
            m = np.nanmean(h[sel], axis=0)
            sd = np.nanstd(h[sel], axis=0) / max(np.sqrt(len(sel)), 1)
            ax.fill_between(lags, m - sd, m + sd, alpha=0.3)
            ax.plot(lags, m, lw=1.6)
            ax.axhline(0, color='0.6', lw=0.8)
            ax.set_xlabel('lag (s)')
            ax.set_ylabel('filter')
            ax.set_title('mean +- sem over {} ROI(s), NOT sign-aligned{}'
                         .format(len(sel),
                                 '\nMIXED POLARITY: mean partly cancels'
                                 if res['mixed_polarity'] else
                                 ('' if res['mean_is_selected'] else
                                  '\nNOTHING KEPT: mean over ALL ROIs')))

            ax = axes[0][2]
            ax.semilogx(res['freqs'][1:], res['mag_mean'][1:], lw=1.4)
            lo, hi, pk = res['band_3db']
            if np.isfinite(lo):
                ax.axvspan(max(lo, res['freqs'][1]), hi, alpha=0.15,
                           color='C1')
            ax.set_xlabel('frequency (Hz)')
            ax.set_ylabel('|H(f)|')
            ax.set_title('-3 dB band {:.3f} - {:.3f} Hz'.format(lo, hi))

            ax = axes[1][0]
            # Both gates really are the two axes of this panel, so draw both
            # cutoffs: an ROI is kept only in the upper-right quadrant.
            ax.scatter(res['null_peak_z'], cv_r, s=18,
                       c=np.where(in_keep, 'C0', '0.7'))
            if np.isfinite(res['null_threshold']):
                ax.axvline(res['null_threshold'], color='k', ls=':', lw=1.0)
            if np.isfinite(res['cv_null_threshold']):
                ax.axhline(res['cv_null_threshold'], color='k', ls=':', lw=1.0)
            # Both dotted lines are the 95th percentile of the null for
            # that axis. On the x axis that is a reference, not the decision:
            # under fdr the shape gate cuts on p, not on z.
            ax.set_xlabel('peak |z|  (dotted = null 95th pct, reference)')
            ax.set_ylabel('cross-validated r  (dotted = null 95th pct)')
            ax.set_title('the two gates: kept = upper right')

            ax = axes[1][1]
            ax.scatter(res['peak_lag_s'][sel], res['dc_gain'][sel], s=18)
            ax.axhline(0, color='0.6', lw=0.8)
            ax.set_xlabel('peak lag (s)')
            ax.set_ylabel('DC gain (signed)')
            ax.set_title('kinetics vs polarity')

            ax = axes[1][2]
            med = np.nanmedian(res['ridge_scores'], axis=1)
            ax.semilogx(res['ridge_grid'], med, marker='o', ms=3)
            ax.axvline(res['ridge'], color='k', ls=':', lw=1.0)
            ax.set_xlabel('ridge')
            ax.set_ylabel('median cv_r')
            ax.set_title('ridge chosen by cross-validation')

            fig.suptitle('{} {} {}  |  {} epochs, {} blocks, ridge {:.3g}'
                         .format(sn, ch,
                                 os.path.basename(
                                     os.path.normpath(experiment_file_directory)),
                                 res['design']['n_ep'],
                                 len(np.unique(res['design']['fold'])),
                                 res['ridge']))
            fig.tight_layout()
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'trf_fullfield_{}_{}_summary.pdf'.format(sn, ch)))
            if not show_figs:
                plt.close(fig)

            # Per-ROI small multiples. One page rather than one file per ROI:
            # a temporal filter is one trace, and the comparison between ROIs is
            # the thing worth looking at.
            n_col = 6
            n_row = int(np.ceil(n_roi / n_col))
            fig, axes = plt.subplots(n_row, n_col, squeeze=False,
                                     figsize=(2.4 * n_col, 1.9 * n_row),
                                     sharex=True)
            for r in range(n_row * n_col):
                ax = axes[r // n_col][r % n_col]
                if r >= n_roi:
                    ax.axis('off')
                    continue
                ax.plot(lags, h[r], lw=1.0,
                        color='C0' if in_keep[r] else '0.6')
                ax.axhline(0, color='0.7', lw=0.6)
                ax.set_title('roi {}  r={:.2f}{}'.format(
                    r, cv_r[r], '' if in_keep[r] else '  (ns)'), fontsize=14)
                ax.tick_params(labelsize=7)
                if r // n_col == n_row - 1:
                    ax.set_xlabel('lag (s)', fontsize=14)
            fig.suptitle('{} {}: per-ROI temporal filters'.format(sn, ch))
            fig.tight_layout()
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'trf_fullfield_{}_{}_rois.pdf'.format(sn, ch)))
            if not show_figs:
                plt.close(fig)

            # Trial by trial, largest repeat group. Rows are ROIs (kept ones,
            # or all when none were kept); trials coloured in recording order so
            # drift or run-down reads as a colour gradient.
            if res['trial_rel'] is not None:
                grp = res['trial_rel'][0]
                show = list(keep if len(keep) else range(n_roi))[:8]
                n_tr = len(grp['epochs'])
                cols = plt.cm.viridis(np.linspace(0, 0.9, n_tr))
                fig, axes = plt.subplots(len(show), 3, squeeze=False,
                                         figsize=(13, 2.4 * len(show)))
                for row, r in enumerate(show):
                    ax = axes[row][0]
                    xmax = min(lags[-1], max(1.0, 1.5 * grp['support_lags'][r]
                                             * res['design']['dt']))
                    for j in range(n_tr):
                        ax.plot(lags, grp['trial_filters'][j, r], color=cols[j],
                                lw=0.9, label='trial {}'.format(j))
                    ax.plot(lags, h[r], 'k', lw=1.8, label='all trials')
                    ax.axhline(0, color='0.7', lw=0.6)
                    ax.set_xlim(0, xmax)
                    nc = grp['noise_ceiling'][r]
                    ax.set_title('roi {}  cv_r {:.2f}  ceiling {}'.format(
                        r, cv_r[r], 'unresolved' if not np.isfinite(nc)
                        else '{:.2f}'.format(nc)), fontsize=14)
                    if row == 0:
                        ax.legend(fontsize=14, ncol=2)
                    ax = axes[row][1]
                    ax.plot(range(n_tr), grp['filter_loo_r'][:, r], 'o-')
                    ax.set_ylim(min(-0.1, np.nanmin(grp['filter_loo_r'][:, r])
                                    - 0.05), 1.02)
                    ax.set_ylabel('trial vs rest r', fontsize=14)
                    ax = axes[row][2]
                    ax.plot(range(n_tr), grp['trial_gain'][:, r], 'o-')
                    ax.axhline(1, color='0.6', lw=0.8, ls=':')
                    ax.set_ylabel('gain vs all-trial filter', fontsize=14)
                    for c in (1, 2):
                        axes[row][c].set_xticks(range(n_tr))
                for c, lab in enumerate(('lag (s)', 'trial (recording order)',
                                         'trial (recording order)')):
                    axes[-1][c].set_xlabel(lab)
                fig.suptitle('{} {}: per-trial filters, repeated stimulus (seed '
                             '{}, {} trials)'.format(sn, ch, grp['seed'], n_tr))
                fig.tight_layout()
                if save_figs:
                    fig.savefig(os.path.join(
                        figs_dir, 'trf_fullfield_{}_{}_trials.pdf'.format(sn, ch)))
                if not show_figs:
                    plt.close(fig)

        print('Figures written to ' + figs_dir)

    if show_figs:
        plt.show()
