"""
analyze_data_strf_one_d.py

Per-ROI 1D spatiotemporal linear filters (STRFs) from the single-layer
bar/wedge ternary noise stimulus (labpack SphericalOneDTernaryNoise),
estimated in the CELL basis.

Call structure (same as analyze_data_strf_octave.py):
    python analyze_data_strf_one_d.py \
        --experiment_file_directory "path/to/fly_folder" \
        --rig Bruker \
        --show_figs False \
        --save_figs True \
        --tag final \
        --dff pre \
        --filter_length 10.0

Results saved to fly.hdf5 (or fly_final.hdf5) under group /STRF_ONE_D_{tag}.

HOW THIS DIFFERS FROM THE FLYMAX NOIZONE PIPELINE
--------------------------------------------------
* No movie file. SphericalOneDTernaryNoise is not a flymax pre-rendered .bin
  movie -- it never sets `movie_bin_name` -- so `is_noizone`/
  `resolve_noizone_epoch_noise_params` in noise_stim.py do not, and cannot,
  recognize these epochs. The stimulus is a pure function of the frame index,
  so the epoch's own metadata IS the stimulus, and reconstruction
  (`one_d_stim.load_one_d_cells`) is a hash evaluation via labpack, not a
  memmap read.

* Dropped frames need no repair of the stimulus array, for the same reason:
  ask for the indices that actually rendered and evaluate them exactly, gaps
  included, rather than reconstructing a positional movie row.

* The hold is exponential (mean `hold_tau`), not the fixed sample-and-hold
  `noisetau` a pre-rendered movie has, so the smearing it puts on the STA is
  not the small triangle the noizone pipeline could ignore -- see
  `--deconvolve` below and `one_d_strf.py` / `octave_strf.deconvolve_octaves`,
  which this script reuses directly: a cell basis with one common variance and
  one exponential hold is the same computation whether the cells are HEALPix
  patches or bars, and duplicating that numerics a second time would only be
  a second way for it to be subtly wrong.

HOW THIS DIFFERS FROM analyze_data_strf_octave.py
---------------------------------------------------
* ONE layer. Bars/wedges are disjoint rings or sectors, not nested scales, so
  the cell basis is NOT overcomplete and there is no merge step and no
  `fit_joint_cells` -- the reverse-correlation STA already IS the reported
  per-cell filter.

* Position is 1D, not a direction on the sphere. A cell is a whole ring or
  sector, so there is no (az, el) and no `spherical_centroid`; position
  summaries use `one_d_strf.centroid_1d` along the stimulus's own axis
  (beta for bars, alpha for wedges -- alpha wraps at 360).

* Visibility is read from the epoch, not recomputed: the protocol already
  records exactly which cells the projector could ever draw a texel from
  (`one_d_visible_cells`, from the real screen geometry at trial time), so
  this script uses that directly rather than re-deriving it from the patch
  size the way `octave_stim.visible_cells` has to.

WHICH ROIS ARE KEPT
-------------------
The circular-shift null (util.nullcal) is run on ONE channel, by default the
structural channel from fly.json -- the calcium indicator in this lab's
convention -- and the kept set it produces is handed to every other channel
of the same family. The other channels are not tested at all. The kept set
answers "which cells responded to the stimulus", and that is read from the
indicator that reports a response reliably; a slow or weak second indicator
(cAMP) is then measured on exactly those cells instead of being asked to
re-earn each one against its own noise. Each channel still reads its own
polarity off its own filter.

A family is a geometry AND a hold_tau (see _geom_key), so the null runs
separately per tau and a tau-0.04 kept set never leaks into a tau-0.4 one:
different taus are different experiments. --select_channel each restores the
old rule, every channel tested on itself; the choice is recorded with the
saved STRFs and a file selected the other way is recomputed, not reused.

POLARITY AND THE REGISTRATION LAG
---------------------------------
A kept roi's sign is the sign of its filter's FIRST lobe, not of its largest
value (`peak_lag_ind` is that lobe's peak for kept rois, the energy peak for
the rest). The extremum is the detection statistic the null calibrates, but
a biphasic calcium filter's rebound can outgrow its first lobe once the
stimulus hold has smeared both: on fly_022 at hold_tau 0.4 the extremum rule
called a third of the kept Mi1 rois negative and contradicted itself between
the two pole families -- independent data on the same cells -- 8 times in 61.
Polarity is read by projecting each roi's pooled time course onto the
family's sign-aligned template and orienting by the template's own first
lobe (`nullcal.first_lobe_polarity`, which carries the measurements behind
each choice; 61 of 61 agreements on the same fly). The extremum signs are
saved alongside as `kept_signs_extremum`, with `polarity_confidence` (|cos|
to the template) and `first_lobe_lag_ind` per kept roi; the rule is the
group attr `polarity_rule` and is part of the reuse check.

The CENTROID of a kept roi is not read at any single lag. Its spatial profile
is the filter projected onto the same family template across lags -- each
cell's time course times the template, summed over the template's support --
which uses every lag that carries signal and has no lag to choose. A single
slice fails either way: at the first-lobe peak when the rebound dominates
(the smeared first lobe is a sliver), at the energy peak when a weak roi's
maximum is a noise lag. See the comment at the centroid block; recorded as
the group attr `registration_rule`, also part of the reuse check.
"""

import sys
import os
import argparse
import hashlib
import io
import json
import pathlib
import time
import numpy as np
import matplotlib.pyplot as plt
from visanalysis.plugin import base as base_plugin
from visanalysis.analysis import imaging_data
from visanalysis.analysis.imaging_data import roi_bool_masks
from visanalysis.util.one_d_stim import (
    is_one_d_ternary, resolve_one_d_epoch_meta, epoch_visible_cells,
    load_one_d_cells, cell_positions_deg, cell_widths_deg, runt_cell_ids,
    pole_label, cell_offsets_deg, screen_axis_label, screen_center_coord,
    azel_from_pole_coords, direction_from_two_families, azel_of,
    screen_azel_of, pole_frame,
)
from visanalysis.util.one_d_strf import (
    one_layer, centroid_1d, compute_strf_cells, compute_strf_cells_folds,
    peak_lag, independent_samples, deconvolve_octaves, zscore_from_noise_cells,
    trial_reliability,
)
from visanalysis.util import bleach, nullcal, separability
from visanalysis.util.frame_timing import reconstruct_displayed_frame_indices
from visanalysis.util.strf_plot import _tighten, apply_figure_style
from visanalysis.util.one_d_group import (
    register_average, polarity_groups, coverage_bounds,
    radial_from_projections, radial_to_image)
import h5py

plt.ioff()

apply_figure_style()

#: How kept rois get their sign and their first-lobe lag. Recorded with the
#: saved STRFs; a file written under another rule is recomputed, not reused.
#: See nullcal, section POLARITY, and the module docstring above.
POLARITY_RULE = 'first_lobe'
#: Where a kept roi's spatial profile is read for its centroid: the filter
#: projected onto the family template across lags (see the module docstring).
#: Recorded with the saved STRFs likewise.
REGISTRATION_RULE = 'matched_profile'


def parse_per_channel(spec, channels, default, name):
    """'10' -> every channel; 'ch1:3,ch2:25' -> per channel.

    Same spelling as average_strf_noizone.parse_thresholds, so one convention
    covers both. A channel missing from a per-channel spec is an error rather
    than a silent default: the whole point of writing one is that the channels
    differ, so a typo that quietly hands a channel the wrong value is the
    failure worth refusing.
    """
    if spec is None:
        return {c: float(default) for c in channels}
    spec = str(spec).strip()
    if ':' not in spec and '=' not in spec:
        return {c: float(spec) for c in channels}
    out = {}
    for part in spec.replace(';', ',').split(','):
        if not part.strip():
            continue
        k, v = part.replace('=', ':').split(':')
        k = k.strip()
        out[k if k.startswith('ch') else 'ch' + k] = float(v)
    missing = [c for c in channels if c not in out]
    if missing:
        raise ValueError('{} gives no value for {}; channels present are {}'
                         .format(name, missing, channels))
    return out


def _geom_key(meta):
    """A hashable key for 'these trials share one tessellation and one hold'.

    Cells (and the single deconvolution layer built from them) only mean the
    same thing across trials that share every field the geometry and the
    per-cell distribution depend on. Unlike the octave stimulus's single
    `rotation_seed`, this stimulus has no one field that captures that -- any
    of mode, cell_width_deg, phase_deg, pole_az_deg, pole_el_deg, hold_tau or
    p can be swept independently -- so the key is all of them together. Seed
    is deliberately excluded: a new seed is a new noise DRAW of the same
    layer, exactly what should vary freely inside one group.
    """
    return (str(meta['mode']), round(float(meta['cell_width_deg']), 6),
            round(float(meta['phase_deg']), 6), round(float(meta['pole_az_deg']), 6),
            round(float(meta['pole_el_deg']), 6), round(float(meta['hold_tau']), 6),
            round(float(meta['p']), 6))


def _geom_label(key):
    """A stable name for one geometry group, derived from its CONTENT.

    Deliberately not a positional index. `trials_by_geom` is rebuilt inside the
    per-channel loop and a channel can drop a trial the other keeps (the empty
    time-vector check), so an enumerate() index would let sn/ch1/geom_0 and
    sn/ch2/geom_0 name different tessellations under the same hdf5 path. The
    octave twin gets this property for free by keying on rotation_seed; this
    stimulus has no single such field, so the label is built from all of them.

    Readable prefix so a figure folder can be scanned by eye, plus a short
    digest of the whole key so two geometries cannot collide even where the
    prefix formats them identically. hold_tau is in the prefix because the
    protocol can sweep it trial by trial, and two families that differ only
    in timescale should be tellable apart by name, not by digest. p stays in
    the digest alone. The digest is over the whole key and did not change
    when tau was added to the prefix, so a group an earlier run saved as
    `bars_w5_ph0_az0_el0_<digest>` is the one now named
    `bars_w5_ph0_az0_el0_tau0.45_<digest>`.
    """
    mode, width, phase, pole_az, pole_el, tau, _p = key
    digest = hashlib.sha1(repr(key).encode('utf-8')).hexdigest()[:6]
    return '{}_w{:g}_ph{:g}_az{:g}_el{:g}_tau{:g}_{}'.format(
        mode, width, phase, pole_az, pole_el, tau, digest)


class _Tee:
    """Mirror a stream to a log file.

    Flushed on every write rather than at close, so a run that dies partway
    through still leaves a complete log up to the point of failure -- which is
    the case a log is most useful for. stderr is teed as well as stdout, so a
    traceback lands in the file instead of only in a terminal that may be gone
    by the time anyone looks.
    """

    def __init__(self, stream, fh):
        self._stream, self._fh = stream, fh

    def write(self, s):
        self._stream.write(s)
        self._fh.write(s)
        self._fh.flush()
        return len(s)

    def flush(self):
        self._stream.flush()
        self._fh.flush()

    def isatty(self):
        return getattr(self._stream, 'isatty', lambda: False)()


#: One log per fly folder, shared by every analysis that writes one. The file
#: name stays generic on purpose -- what distinguishes one run from another goes
#: in the banner, not in the name, so the history of a fly reads in one place
#: and in order.
ANALYSIS_LOG_NAME = 'analysis_log.txt'

#: Figures from THIS analysis go in their own folder, not the one
#: analyze_data.py fills with mask overlays and search-stim responses.
#: Mixing them is how a fly folder ends up with 2600 files of which the
#: ones you want are a minority -- fly_008's raw_roi_figs did exactly
#: that. The tag still prefixes it, so raw and final stay separate.
FIGURE_SUBDIR = 'noise_figs'


def start_log(directory, analysis, tag, argv):
    """Tee stdout and stderr into the fly's running analysis log.

    APPENDED, never overwritten, and named generically: re-running is normal
    here -- to redraw figures, or with --force_reanalysis after something
    changed -- and the useful question afterwards is what was run against this
    fly and in what order. A file per run answers that only if you already know
    which file to open. Each entry opens with a banner naming the analysis, the
    time and the exact command, which is what the name would otherwise carry.

    Returns the log path.
    """
    path = os.path.join(directory, ANALYSIS_LOG_NAME)
    try:
        fh = io.open(path, 'a', encoding='utf-8')
    except (IOError, OSError) as e:
        # The log is a record, not the job. A fly folder that is read-only to
        # whoever is running -- another lab member's data, a share mounted
        # without write access -- must not become a traceback before any
        # analysis has been attempted.
        print('(no analysis log: {}; continuing without one)'.format(e))
        return None
    sys.stdout = _Tee(sys.stdout, fh)
    sys.stderr = _Tee(sys.stderr, fh)
    print('')
    print('=' * 78)
    print('analysis : ' + analysis)
    print('run at   : ' + time.strftime('%Y-%m-%d %H:%M:%S'))
    print('tag      : ' + str(tag))
    print('command  : ' + ' '.join(argv))
    print('=' * 78)
    return path


def pole_of(meta):
    """The pole a geometry is built on, rounded so two runs of the same pole
    group together. A bars and a wedges family about ONE pole are a complete
    spherical coordinate, whatever that pole is -- see
    one_d_stim.azel_from_pole_coords -- so this, not the mode, is what decides
    which geometries can be paired into a position."""
    return (round(float(meta['pole_az_deg']), 6),
            round(float(meta['pole_el_deg']), 6))


def _full_meta_from_epochs(f, sn, saved_meta):
    """The complete one_d_stimulus record of a saved family, read from the
    epoch records of its series, for files whose family group predates the
    `one_d_stimulus_json` attr. Matched on the geometry key, so a series that
    swept tau or pole hands back the right record for each family. None when
    no epoch matches."""
    want = _geom_key(saved_meta)
    series = 'series_{:03d}'.format(int(str(sn).replace('sn', '')))
    subjects = f.get('Subjects')
    if subjects is None:
        return None
    for sub in subjects:
        runs = subjects[sub].get('epoch_runs')
        if runs is None or series not in runs or 'epochs' not in runs[series]:
            continue
        for ek in runs[series]['epochs']:
            attrs = dict(runs[series]['epochs'][ek].attrs)
            try:
                meta, _s = resolve_one_d_epoch_meta(attrs)
            except (KeyError, ValueError):
                continue
            if _geom_key(meta) == want:
                return meta
    return None


def load_saved_strfs(path, group, filter_length_spec, filter_length_by_ch,
                     dff, deconvolve, lag_dt, null_method, select_channel,
                     peak_window_s=nullcal.DEFAULT_PEAK_WINDOW_S):
    """Read back a previous run's STRFs, when they describe the SAME analysis.

    Re-running to refresh figures should not recompute every STA. But saved
    results from a different filter length, dF/F mode or deconvolution setting
    are a different analysis wearing the same group name, so compatibility is
    checked against the attrs that run recorded rather than assumed -- silently
    redrawing stale filters under new parameters is the failure worth spending
    a few lines to prevent.

    Returns (results, reused, reason). Anything that does not line up gives
    ({}, False, why), which puts the caller back on the full path.
    """
    if not os.path.exists(path):
        return {}, False, 'no experiment file yet - computing'
    with h5py.File(path, 'r') as f:
        if group not in f:
            return {}, False, 'no saved /{} - computing'.format(group)
        g = f[group]
        # FILTER LENGTH: compared as the RESOLVED number per channel, with a
        # tolerance, not as the spec string. The spelling is not the analysis:
        # '10', '10.0' and 'ch1=10,ch2=10' are one filter length, and a string
        # comparison refused to reuse across them -- which cost a full
        # recompute, null test included, for nothing. What IS the analysis is
        # the number each channel ran with, so that is what is stored (as a
        # JSON dict, written below) and what is compared, for every channel
        # this run wants. A file from before the dict existed falls back to
        # the string, exactly as before, so nothing already saved is reused
        # under settings that cannot be checked.
        want_len = {str(k): float(v) for k, v in dict(filter_length_by_ch).items()}
        saved_len_raw = g.attrs.get('filter_length_by_ch', None)
        if saved_len_raw is not None:
            try:
                saved_len = {str(k): float(v) for k, v in json.loads(
                    saved_len_raw.decode('utf-8') if isinstance(saved_len_raw, bytes)
                    else str(saved_len_raw)).items()}
            except (ValueError, AttributeError):
                saved_len = None
        else:
            saved_len = None
        if saved_len is not None:
            _bad = {c: (saved_len.get(c), v) for c, v in want_len.items()
                    if c not in saved_len
                    or not np.isclose(saved_len[c], v, rtol=1e-9, atol=1e-9)}
            if _bad:
                return ({}, False,
                        'saved /{} has filter_length {} per channel, this run '
                        'wants {} - recomputing'.format(
                            group,
                            {c: saved_len.get(c) for c in sorted(want_len)},
                            {c: want_len[c] for c in sorted(want_len)}))
            len_ok = True
            len_saved_desc = str({c: saved_len[c] for c in sorted(want_len)})
        else:
            len_ok = (str(g.attrs.get('filter_length_spec', '<none>'))
                      == str(filter_length_spec))
            len_saved_desc = str(g.attrs.get('filter_length_spec', '<none>'))
        saved = (len_ok,
                 str(g.attrs.get('dff', '')),
                 bool(g.attrs.get('deconvolved', False)),
                 float(g.attrs.get('lag_dt', -1.0)))
        want = (True, str(dff), bool(deconvolve), float(lag_dt))
        # The kept set is part of what is saved, so the rule that made it is
        # part of the analysis identity. A file written before the attr
        # existed was selected by FDR (the only rule there was) and reads as
        # such, so switching the default to the percentile recomputes it
        # rather than redrawing an FDR selection under the new name.
        saved_null = str(g.attrs.get('null_method', 'fdr'))
        if saved_null != str(null_method):
            return ({}, False,
                    'saved /{} selected rois by {}, this run wants {} - '
                    'recomputing'.format(group, saved_null, null_method))
        # Likewise WHICH CHANNEL the kept set was decided on. Before the attr
        # existed every channel was tested on itself, so a file without it
        # reads as 'each'; a run that inherits from the structural channel
        # then recomputes rather than redrawing per-channel sets as inherited.
        _ss = g.attrs.get('roi_selection_channel', 'each')
        saved_sel = _ss.decode('utf-8') if isinstance(_ss, bytes) else str(_ss)
        want_sel = 'each' if select_channel is None else str(select_channel)
        if saved_sel != want_sel:
            return ({}, False,
                    'saved /{} decided kept rois on {}, this run wants {} - '
                    'recomputing'.format(group, saved_sel, want_sel))
        # And the POLARITY rule. A file from before it existed read signs at
        # the |z| extremum and registered every roi at its energy peak; those
        # are the kept_signs, peak_lag_ind and centroids a reuse would redraw.
        _pr = g.attrs.get('polarity_rule', 'extremum')
        saved_pol = _pr.decode('utf-8') if isinstance(_pr, bytes) else str(_pr)
        if saved_pol != POLARITY_RULE:
            return ({}, False,
                    'saved /{} read polarity by the {} rule, this run uses {} - '
                    'recomputing'.format(group, saved_pol, POLARITY_RULE))
        # And the LAG WINDOW the peak was searched over, rois and null alike.
        # A file from before the attr existed searched the whole axis, which
        # is what 0 means, so it reads as 0 and recomputes under any window.
        saved_win = float(g.attrs.get('peak_window_s', 0.0))
        want_win = float(peak_window_s) if float(peak_window_s) > 0 else 0.0
        if not np.isclose(saved_win, want_win, rtol=1e-9, atol=1e-9):
            return ({}, False,
                    'saved /{} searched the peak over {}, this run wants {} - '
                    'recomputing'.format(
                        group,
                        '{:g} s of lag'.format(saved_win) if saved_win > 0
                        else 'the whole filter',
                        '{:g} s of lag'.format(want_win) if want_win > 0
                        else 'the whole filter'))
        _rr = g.attrs.get('registration_rule', 'energy_peak')
        saved_reg = _rr.decode('utf-8') if isinstance(_rr, bytes) else str(_rr)
        if saved_reg != REGISTRATION_RULE:
            return ({}, False,
                    'saved /{} registered rois by {}, this run uses {} - '
                    'recomputing'.format(group, saved_reg, REGISTRATION_RULE))
        # lag_dt compared with a tolerance, the rest exactly: it is a float a
        # user may retype rather than a label, and 1/60 written back and read
        # again should not count as a different analysis.
        if (saved[:3] != want[:3]
                or not np.isclose(saved[3], want[3], rtol=1e-9, atol=0.0)):
            return ({}, False,
                    'saved /{} is filter_length={} dff={} deconvolved={} '
                    'lag_dt={:.5f}s, this run wants {} / {} / {} / {:.5f}s - '
                    'recomputing'.format(
                        group, len_saved_desc, saved[1], saved[2], saved[3],
                        filter_length_spec, want[1], want[2], want[3]))
        out = {}
        try:
          for sn in g:
            for ch in g[sn]:
                for gname in g[sn][ch]:
                    sg = g[sn][ch][gname]
                    glabel = (gname[len('geom_'):] if gname.startswith('geom_')
                              else gname)
                    meta = {k: sg.attrs[k] for k in
                            ('mode', 'cell_width_deg', 'phase_deg', 'pole_az_deg',
                             'pole_el_deg', 'hold_tau', 'p', 'seed')
                            if k in sg.attrs}
                    if 'mode' in meta:
                        meta['mode'] = (meta['mode'].decode('utf-8')
                                        if isinstance(meta['mode'], bytes)
                                        else str(meta['mode']))
                    # The FULL record, not the eight labelling attrs: the
                    # geometry helpers need the screen frame (frame_az_deg /
                    # frame_el_deg) to place anything on a pitched screen.
                    # From the attr written since 2026-11-17 when present;
                    # otherwise from the series' own epoch records, which have
                    # carried it all along; failing both, the KeyError below
                    # sends the file to a recompute rather than a wrong map.
                    _mj = sg.attrs.get('one_d_stimulus_json', None)
                    if _mj is not None:
                        _full = json.loads(_mj.decode('utf-8')
                                           if isinstance(_mj, bytes) else str(_mj))
                    else:
                        _full = _full_meta_from_epochs(f, sn, meta)
                        if _full is None:
                            raise KeyError(
                                'no one_d_stimulus_json on /{}/{}/{}/{} and no '
                                'matching epoch record; cannot recover the '
                                'screen frame'.format(group, sn, ch, gname))
                    _full['mode'] = str(_full.get('mode', meta.get('mode', '')))
                    meta = _full
                    res = {
                        'strf': sg['strf'][()], 'strf_z': sg['strf_z'][()],
                        't_lag': sg['t_lag'][()],
                        'positions_deg': sg['positions_deg'][()],
                        'widths_deg': sg['widths_deg'][()],
                        'screen_offset_deg': sg['screen_offset_deg'][()],
                        'center_screen_deg': sg['center_screen_deg'][()],
                        'jitter_deg': (lambda v: np.nan if v < 0 else v)(
                            float(sg.attrs.get('registration_jitter_deg', -1.0))),
                        'axis_label': str(sg.attrs.get('screen_axis_label', '')),
                        'peak_lag_ind': sg['peak_lag_ind'][()],
                        'visible_cells': sg['visible_cells'][()],
                        'center_deg': sg['center_deg'][()],
                        'meta': meta,
                        'deconv_applied': bool(sg.attrs.get('deconv_applied', False)),
                        'deconv_how': str(sg.attrs.get('deconv_how', 'unknown')),
                        'keep_rois': (sg['kept_rois'][()] if 'kept_rois' in sg
                                      else None),
                        # These three are read by the group section and the rf
                        # center map. They were saved but not loaded, and every
                        # consumer reaches them through .get with a default, so
                        # a reuse run produced plausible-but-wrong output in
                        # silence: every roi treated as positive polarity, and
                        # the patch center taken as (0, 0) instead of its real
                        # position, which shifts the whole rf map.
                        'kept_signs': (sg['kept_signs'][()] if 'kept_signs' in sg
                                       else None),
                        'polarity_confidence': (sg['polarity_confidence'][()]
                                                if 'polarity_confidence' in sg
                                                else None),
                        'center_az_deg': float(sg.attrs.get('center_az_deg', 0.0)),
                        'center_el_deg': float(sg.attrs.get('center_el_deg', 0.0)),
                        'trial_rel': None,
                        'sep': None,
                        'sep_window_lags': int(sg['strf'].shape[-1]),
                    }
                    if 'separability' in sg:
                        _sg = sg['separability']
                        res['sep'] = {k: _sg[k][()] for k in _sg}
                        res['sep']['n_tested'] = int(_sg.attrs.get('n_tested', 0))
                        res['sep_window_lags'] = int(_sg.attrs.get(
                            'window_lags', res['sep_window_lags']))
                    if 'trial_reliability' in sg:
                        grps = []
                        for gk in sorted(sg['trial_reliability']):
                            tg = sg['trial_reliability'][gk]
                            d = {k: tg[k][()] for k in tg}
                            d['seed'] = int(tg.attrs.get('seed', -1))
                            grps.append(d)
                        res['trial_rel'] = grps or None
                    out[sn, ch, glabel] = res
        except KeyError as e:
            # A group written by an older version of this script: it passed the
            # settings check above but does not carry every field the figures now
            # read. Recompute rather than half-load -- the alternative is a crash
            # partway through, or figures drawn from a mix of two formats.
            return ({}, False,
                    'saved /{} predates a field this version needs ({}) - '
                    'recomputing'.format(group, e))
        if not out:
            return {}, False, 'saved /{} is empty - computing'.format(group)
        return (out, True,
                'reusing saved STRFs from /{} ({} geometry group(s)) - pass '
                '--force_reanalysis to recompute'.format(group, len(out)))


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
    parser.add_argument('--dff', nargs='?', default='pre')
    parser.add_argument('--null_shifts', nargs='?', type=int, default=20,
                        help='circular-shift null draws; 0 disables. Costs (n+1)x '
                             'the STA computation. 20 is the shared default of every '
                             'noise analysis: the p floor of the pooled null is '
                             '1/(n_roi*n_shifts+1), and at 10 shifts BH could not '
                             'start on a family with only a handful of responsive '
                             'rois (fly_022 ch1 tau 0.04 kept 5 of ~10 real). Past '
                             '~20 the shifts start sampling overlapping stretches of '
                             'one recording, so the draws stop being independent. '
                             'See util.nullcal.')
    parser.add_argument('--peak_window', nargs='?', type=float,
                        default=nullcal.DEFAULT_PEAK_WINDOW_S,
                        help='lag window, in seconds from lag 0, searched for the '
                             'peak |z| -- on the rois and on every null draw alike; '
                             '0 searches the whole filter. See util.nullcal.')
    parser.add_argument('--separability', nargs='?', default='False',
                        help="'True' tests every kept roi for space-time separability "
                             "across the two trial folds: the cross-validated rank-1 "
                             "fraction and a second-component test against the shift "
                             "null. Raw STA, visible cells, the peak window. 'False' "
                             "(default) skips the per-roi test: separability is a "
                             "question for GROUP data, where the second component is "
                             "pooled across rois and flies, not for one roi of one fly, "
                             "whose per-roi call is noise-limited at the SNR these "
                             "recordings reach. Nothing per fold is computed or saved "
                             "unless the test runs. See util.separability.")
    parser.add_argument('--bleach', nargs='?', default='auto',
                        help="bleach correction on raw F: 'auto' (default), 'on', 'off'")
    parser.add_argument('--filter_length', nargs='?', default='10.0',
                        help="STRF filter length in seconds. One number for "
                             "every channel, or per channel as "
                             "'ch1:3,ch2:25' -- a slow indicator needs a long "
                             "filter while a fast one only pays for it, in "
                             "compute and in a lag axis that is mostly empty.")
    parser.add_argument('--lag_dt', nargs='?', type=float, default=1.0 / 60.0,
                        help='lag axis step in SECONDS. Default 1/60 s, i.e. half '
                             'the display rate. One lag per display frame is about '
                             "14 bins per response sample at this rig's imaging "
                             'rate, which costs compute without resolving '
                             'anything; measured, 1/60 s matches the 1/120 s '
                             'filter to r = 0.9998 for half the cost. A value '
                             'finer than the display frame period is raised to it.')
    parser.add_argument('--monitor_hz', nargs='?', type=float, default=None,
                        help='override the display rate; default uses the ImagingDataObject')
    parser.add_argument('--lam', nargs='?', type=float, default=None,
                        help='deconvolution ridge. Omit to choose it by cross-validation, '
                             'which is the intended path -- a value that suits one indicator '
                             'will not suit another.')
    parser.add_argument('--lam_scale', nargs='?', type=float, default=1.0,
                        help='multiplies the cross-validated lambda. >1 is smoother.')
    parser.add_argument('--majority_sign_only', action='store_true',
                        help='also average only the rois whose filter carries '
                             'the channel-majority polarity. Defensible when the '
                             "cell type's sign is known independently -- Mi1's "
                             'calcium response is ON -- because the criterion is '
                             'then prior knowledge rather than the data being '
                             'averaged. It is NOT the same as sign-aligning, '
                             'which flips rois to match and makes noise add '
                             'constructively at every peak (see util.nullcal); '
                             'this only drops dissenters. Note the majority is '
                             'still identified empirically here, so on a channel '
                             'whose expected sign is not established the argument '
                             'is weaker. Produced beside the full mean, not '
                             'instead of it.')
    parser.add_argument('--edge_margins', nargs='?', default='',
                        help="comma-separated margins in degrees, e.g. "
                             "'0,10,20'. For each one, an extra registered mean "
                             "is made from the rois at least that far inside the "
                             "OBSERVED extent of the population's receptive-field "
                             "centers, in azimuth and in elevation. The observed "
                             "extent, not the lit region: the bowl's illuminated "
                             "area is what the projector reached, while the roi "
                             "centers are the only evidence of what the fly "
                             "actually resolved, and the two need not agree. Run "
                             "several margins and read them as a series -- if the "
                             "kinetics do not move, edge rois are not affecting "
                             "them and the full population can be kept on "
                             "evidence rather than assumption.")
    parser.add_argument('--mean_span_deg', nargs='?', type=float, default=75.0,
                        help='half-width, in degrees, of the registered-mean '
                             'figures. The registered grid runs as far as the '
                             'widest roi reached, which on a wedges family is '
                             'most of a hemisphere and leaves the receptive '
                             'field a few pixels across. Axes are clamped to '
                             'this or to the data, whichever is tighter, and '
                             'the z limits are measured INSIDE it so zooming '
                             'magnifies the filter instead of leaving it flat '
                             'against an excursion off-screen.')
    parser.add_argument('--n_display_lags', nargs='?', type=int, default=12,
                        help='how many lags the 2D lag grid shows, spread over '
                             'the filter support. A single frame at the peak is '
                             'not enough when center and surround peak at '
                             'different times, which is exactly when the 2D '
                             'picture is worth having.')
    parser.add_argument('--radial_symmetry', action='store_true',
                        help='also recover a radial profile h(r) from the '
                             'registered means, by forward-modelling the strip '
                             'projections and solving one h against EVERY family '
                             'at once. The reported residual is what falsifies '
                             'the assumption -- a single family cannot.')
    parser.add_argument('--figure_rois', nargs='?', default='kept',
                        help="which rois get per-roi figures: 'kept' (default, "
                             "those that passed the null), 'all', 'none' to "
                             "write no per-roi figures at all, or an integer N "
                             "for the N strongest kept. 'all' on a 289-roi fly "
                             "with two geometries and two channels is 2312 "
                             "files, most of them noise in a channel where "
                             "nothing passed; 'none' is what you want when only "
                             "the population figures have changed, which is "
                             "every rerun that only moves a display option.")
    parser.add_argument('--localize_channel', nargs='?', default='structural',
                        help="channel whose centroids place the receptive-field "
                             "map: 'structural' (default, read from fly.json), "
                             "'each' for one map per channel, or a channel "
                             "number such as '1'.")
    parser.add_argument('--null_method', nargs='?', default='fdr',
                        choices=['fdr', 'percentile', 'auto'],
                        help="how the circular-shift null becomes a kept set. "
                             "'fdr' (default): Benjamini-Hochberg at q=0.05 over the "
                             "rois, i.e. about one in twenty KEPT rois may be noise. "
                             "It is not a family-wise test: where there is signal it "
                             "keeps most of what p<0.05 keeps (fly_022 ch1: 78 of 90), "
                             "and where there is none it keeps nothing. Its one "
                             "weakness is the p floor 1/(n_null+1) from few shifts, "
                             "which costs real rois in a family with only a handful; "
                             "raise --null_shifts there. 'percentile': peak z above "
                             "the 95th percentile of the pooled null, p < 0.05 per roi "
                             "with NO correction -- lets ~5%% of ALL rois through as "
                             "noise, about 11 per 213-roi family, whether or not the "
                             "family has signal. 'auto': fdr where the shift count can "
                             "resolve it, percentile otherwise. The choice is recorded "
                             "with the saved STRFs; a file selected by another rule is "
                             "recomputed, not reused.")
    parser.add_argument('--select_channel', nargs='?', default='structural',
                        help="which channel's null decides the kept rois for ALL "
                             "channels of a family: 'structural' (default, read from "
                             "fly.json; calcium by lab convention), a channel number "
                             "such as '1', or 'each' to test every channel on itself "
                             "(the old rule). Other channels run no null; they take "
                             "the kept set and read their own polarity off their own "
                             "filter. Families are per hold_tau, so each tau is "
                             "selected on its own. Recorded with the saved STRFs; a "
                             "file selected the other way is recomputed, not reused.")
    parser.add_argument('--force_reanalysis', action='store_true',
                        help='recompute even when this file already holds STRFs for '
                             'this tag, overwriting them. Without it a run that finds '
                             'compatible saved results reuses them and only redraws '
                             'the figures.')
    parser.add_argument('--deconvolve', action='store_true',
                        help='undo the exponential-hold temporal smear. OFF by default: it '
                             'does nothing to the spatial (positional) result for a separable '
                             'filter, and only corrects absolute latency. Turn it on if you '
                             'need that, and check the printed decision -- the reliability '
                             'guard refuses it wherever it would cost more than it buys.')
    args = parser.parse_args()

    experiment_file_directory = args.experiment_file_directory
    rig = args.rig
    show_figs = args.show_figs == 'True'
    save_figs = args.save_figs == 'True'
    tag = args.tag
    dff = args.dff
    bleach_mode = {'auto': 'auto', 'on': True, 'off': False}[str(args.bleach).lower()]
    null_shifts = int(args.null_shifts)
    null_method = str(args.null_method)
    peak_window_s = float(args.peak_window)
    do_sep = str(args.separability) == 'True'
    filter_length_spec = str(args.filter_length)
    lag_dt_s = float(args.lag_dt)
    monitor_hz_override = args.monitor_hz

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
    experiment_file_path = os.path.join(experiment_file_directory, experiment_file_name)
    if not os.path.exists(experiment_file_path):
        print('No such file: ' + experiment_file_path)
        sys.exit(1)
    # Opened before anything is decided, so the log holds the WHOLE run: which
    # file was read, the reuse decision, the per-geometry summaries and any
    # traceback. It used to open after the reuse check, which left the one line
    # saying why a run recomputed -- or did not -- in the terminal only, and
    # that line is the first thing worth knowing when rereading a log.
    log_path = start_log(experiment_file_directory,
                         '1D ternary STRF (analyze_data_strf_one_d.py)',
                         tag, sys.argv)
    print('reading ' + experiment_file_name)

    strf_group = 'STRF_ONE_D_' + tag
    if rig == 'Bruker':
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        plug = base_plugin.BasePlugin()

    response_set_name_prefix = 'mask_'

    # Channels come from fly.json, as process_data.py and analyze_data.py read
    # them, rather than being assumed to be 1 and 2: a fly with one functional
    # channel would otherwise have a second one looked up and quietly skipped.
    # The STRUCTURAL channel is what places the receptive-field map by default --
    # it is the one set to calcium in this lab's convention, and localising every
    # channel on itself lets a noisy indicator put its own noise on the map.
    struct_channel_num = None
    try:
        with open(pathlib.Path(experiment_file_directory, 'fly.json'), 'r') as _f:
            _fly_json = json.load(_f)
        struct_channel_num = str(_fly_json['structural_channel']).split('_')[-1]
        _fc = _fly_json['functional_channel'].replace('[', '').replace(']', '')
        func_channels_num = [c.split('_')[-1].strip()
                             for c in _fc.replace("'", '').split(',') if c.strip()]
        print('structural channel: ch{}   functional: {}'.format(
            struct_channel_num, ['ch' + c for c in func_channels_num]))
    except (IOError, OSError, KeyError, ValueError) as _e:
        func_channels_num = ['1', '2']
        print('could not read fly.json ({}); assuming functional channels 1 and 2 '
              'and no structural channel'.format(_e))

    # Resolved once, after the channel list is settled either way -- which is
    # also why the settings banner below cannot be printed any earlier.
    filter_length_by_ch = parse_per_channel(
        filter_length_spec, ['ch' + c for c in func_channels_num], 10.0,
        '--filter_length')

    # Which channel's null decides the kept set (None: each on itself). Needs
    # fly.json for 'structural', so it too is resolved only here. A channel
    # that is not functional is an error, not a fallback: a typo must not
    # silently change what every other channel is measured on.
    _sel = str(args.select_channel).lower()
    if _sel == 'structural' and struct_channel_num is not None:
        select_ch = 'ch' + str(struct_channel_num)
    elif _sel in ('each', 'all'):
        select_ch = None
    elif _sel == 'structural':
        select_ch = None
        print('no structural channel in fly.json; every channel selects its '
              'own rois (--select_channel each)')
    else:
        select_ch = _sel if _sel.startswith('ch') else 'ch' + _sel
    if select_ch is not None and select_ch not in ['ch' + c for c in func_channels_num]:
        print('--select_channel {}: not a functional channel of this fly ({})'
              .format(args.select_channel, ['ch' + c for c in func_channels_num]))
        sys.exit(1)

    # Decided here rather than at the top, because the reuse check compares
    # the filter length per channel and the channel list has only now been
    # read from fly.json. Nothing between the top and here depends on it.
    strf_results, reused = {}, False
    if args.force_reanalysis:
        print('--force_reanalysis: recomputing and overwriting /' + strf_group)
    else:
        strf_results, reused, _why = load_saved_strfs(
            experiment_file_path, strf_group, filter_length_spec,
            filter_length_by_ch, dff, args.deconvolve, float(args.lag_dt),
            null_method, select_ch, peak_window_s)
        print(_why)

    print('settings : filter_length={}  dff={}  deconvolve={}  '
          'null_shifts={}  null_method={}  roi_selection={}  bleach={}  rig={}  '
          'polarity={}  registration={}'.format(
              {k: '{:g}s'.format(v) for k, v in sorted(filter_length_by_ch.items())},
              dff, args.deconvolve, null_shifts, null_method,
              ('each channel on itself' if select_ch is None
               else 'on {}, inherited by the rest'.format(select_ch)),
              args.bleach, rig, POLARITY_RULE, REGISTRATION_RULE))
    print('           lag_dt {:.5f}s ({:g} Hz)'.format(
        lag_dt_s, (1.0 / lag_dt_s) if lag_dt_s > 0 else float('nan')))
    print('           peak window {}'.format(
        '{:g}s of lag, rois and null alike'.format(peak_window_s)
        if peak_window_s > 0 else 'the whole filter'))

    # ── fly metadata and photodiode channel ──────────────────────────────────
    series_num = list(map(str, plug.getSeriesNumbers(experiment_file_path)))
    ID = imaging_data.ImagingDataObject(experiment_file_path, int(series_num[0]), quiet=True)
    fly_metadata = ID.getSubjectMetadata()
    prep = fly_metadata.get('prep', '')
    timing_channel_ind = 1 if prep == 'fly right optic lobe' else 0
    print('fly_metadata: ' + repr(fly_metadata))
    print('timing_channel_ind: ' + repr(timing_channel_ind))

    # ── find the 1D ternary series ───────────────────────────────────────────
    roi_data = {}
    run_parameters = {}
    one_d_series_info = {}

    # Skipped wholesale when reusing: this loop is where roi responses are read,
    # which is the expensive part of the run.
    for current_series in ([] if reused else series_num):
        sn = 'sn' + current_series
        plug.updateImagingDataObject(experiment_file_directory, experiment_file_name,
                                     int(current_series))
        ID = plug.ImagingDataObject
        ID.timing_channel_ind = timing_channel_ind
        ID.threshold = 0.6
        _, _, _sample_rate = ID.getVoltageData()
        ID.frame_slop = 0.001 * _sample_rate

        run_parameters[sn] = ID.getRunParameters()
        if not is_one_d_ternary(run_parameters[sn]):
            continue

        epochs = []
        for ep in ID.getEpochParameters():
            meta, sampling = resolve_one_d_epoch_meta(ep)
            visible = epoch_visible_cells(ep, meta, sampling)
            epochs.append((meta, sampling, visible))
        one_d_series_info[sn] = {'epochs': epochs}

        for current_channel in func_channels_num:
            ch = 'ch' + current_channel
            _t0 = time.perf_counter()
            roi_data[sn, ch] = bleach.roi_responses(
                ID, response_set_name_prefix + ch,
                roi_prefix='aligned', dff=dff, correct_bleach=bleach_mode)
            print('    [timing] {} {}: getRoiResponses took {:.2f}s'.format(
                sn, ch, time.perf_counter() - _t0))

    if not one_d_series_info and not reused:
        print('No 1D ternary series identified. Check is_one_d_ternary() and '
              'that process_data.py attached the epoch parameters for this series.')
        sys.exit(0)

    if one_d_series_info:
        print('\nFound {} 1D ternary series:'.format(len(one_d_series_info)))
    for sn, info in one_d_series_info.items():
        seeds = sorted({m[0]['seed'] for m in info['epochs']})

        # Every field but the seed can be swept trial by trial -- the pole
        # routinely is, and hold_tau can be -- so each is reported as the set of
        # values the series actually showed. Reading them off the first epoch
        # presented one trial's tau as the series' tau.
        def _vals(field, fmt):
            vs = sorted({float(m[0][field]) for m in info['epochs']})
            return '/'.join(fmt.format(v) for v in vs)
        modes = '/'.join(sorted({str(m[0]['mode']) for m in info['epochs']}))
        print('  {}: {} epochs  mode={}  cell_width={}deg  phase={}deg  '
              'pole=(az {}, el {})  hold_tau={}s  p={}  n_cells={}  '
              '{} distinct seed(s)'.format(
                  sn, len(info['epochs']), modes, _vals('cell_width_deg', '{:.3g}'),
                  _vals('phase_deg', '{:.3g}'), _vals('pole_az_deg', '{:.3g}'),
                  _vals('pole_el_deg', '{:.3g}'), _vals('hold_tau', '{:.3g}'),
                  _vals('p', '{:.5f}'), _vals('n_cells', '{:g}'), len(seeds)))

    # ── STRFs ────────────────────────────────────────────────────────────────

    for sn, info in one_d_series_info.items():
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

        # The selection channel goes first, so its kept set already exists
        # for each family when the other channels reach for it below.
        _ch_order = sorted(func_channels_num, key=lambda c: 'ch' + c != select_ch)
        for current_channel in _ch_order:
            ch = 'ch' + current_channel
            if (sn, ch) not in roi_data:
                continue

            epoch_response = roi_data[sn, ch]['epoch_response']
            time_vector_by_epoch = roi_data[sn, ch]['time_vector_by_epoch']
            n_roi = epoch_response.shape[0]
            roi_z_offsets = roi_data[sn, ch].get('roi_z_offsets', np.zeros(n_roi))

            # This channel's own filter length: a slow indicator needs a long
            # lag axis and a fast one is only slowed down by it.
            filter_length_s = filter_length_by_ch[ch]
            _dts = [np.median(np.diff(tv)) for tv in time_vector_by_epoch if len(tv) > 1]
            sample_dt = float(np.median(_dts)) if _dts else raw_stim_dt

            print('  {} {}: n_roi={}, n_trials={}, n_tp={}, imaging {:.2f} Hz, '
                  'filter_length {:g}s'.format(
                      sn, ch, n_roi, epoch_response.shape[1],
                      epoch_response.shape[2], 1.0 / sample_dt, filter_length_s))

            # Group by geometry+hold: cells (and the one deconvolution layer built
            # from them) only mean the same thing across trials that share those.
            trials_by_geom = {}
            geom_meta = {}
            geom_sampling = {}
            geom_visible_sets = {}
            for trial_ind, (meta, _sampling, visible) in enumerate(info['epochs']):
                if trial_ind >= epoch_response.shape[1]:
                    break
                tvec = time_vector_by_epoch[trial_ind]
                if len(tvec) == 0:
                    continue

                n_frames = int(np.ceil((tvec[-1] - pre_time) * frame_rate)) + 1
                if n_frames <= 0:
                    continue

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

                stim_cells = load_one_d_cells(meta, n_frames, frame_indices=idx)

                key = _geom_key(meta)
                geom_meta[key] = meta
                # The patch geometry is what turns the stimulus's pole-anchored
                # coordinate into a screen-referenced one.
                geom_sampling.setdefault(key, _sampling)
                # Visible-cell set depends on the PATCH (sampling), not just the
                # tessellation (meta), so in principle it can differ between trials
                # grouped by _geom_key alone -- e.g. if 'center' were swept. Collect
                # every trial's set and union them below rather than assuming they
                # agree: a cell the patch showed in even one grouped trial carried
                # real driven signal in that trial and belongs in the visible set,
                # while one shown in none of them is pure STA noise regardless.
                geom_visible_sets.setdefault(key, []).append(frozenset(int(c) for c in visible))
                grp = trials_by_geom.setdefault(key, [])
                for roi_ind in range(n_roi):
                    grp.append({
                        'roi_ind': roi_ind,
                        # trial_ind and seed identify the TRIAL and the stimulus
                        # it showed. Repeats are detected from the seed (see
                        # one_d_strf.repeat_groups), and trial_ind additionally
                        # makes compute_strf_cells_folds split folds on a stated
                        # trial rather than on an array id.
                        'trial_ind': trial_ind,
                        'seed': int(meta['seed']),
                        'response': epoch_response[roi_ind, trial_ind, :],
                        'time_vector': tvec + roi_z_offsets[roi_ind],
                        'stim_cells': stim_cells,
                        'stim_dt': raw_stim_dt,
                    })

            for key, trials in trials_by_geom.items():
                glabel = _geom_label(key)
                meta = geom_meta[key]
                distinct_sets = geom_visible_sets[key]
                if len(set(distinct_sets)) > 1:
                    print('    NOTE: visible-cell set differs across trials grouped '
                          'under geom={} (e.g. a swept patch center) -- using the '
                          'union of {} distinct sets.'.format(glabel,
                                                             len(set(distinct_sets))))
                visible = np.array(sorted(set().union(*distinct_sets)), dtype=np.int64)
                n_cells = int(meta['n_cells'])
                layer = one_layer(meta)
                sampling = geom_sampling[key]
                positions = cell_positions_deg(meta)
                # Where each cell sits RELATIVE TO THE MIDDLE OF THE DISPLAY.
                # `positions` is anchored to the pole and, for wedges, runs
                # phase..phase+360, so neighbouring cells can sit at opposite
                # ends of it. Anchoring to the patch center puts the branch cut
                # behind the fly where nothing is lit, which makes the visible
                # cells contiguous by construction, and makes 0 mean the middle
                # of the screen. Both are recorded: `positions` is the
                # stimulus's own coordinate and matches cell_bounds.
                offsets, center_defined = cell_offsets_deg(meta, sampling)
                widths = cell_widths_deg(meta)
                runts = runt_cell_ids(meta)
                vis_mask = np.zeros(n_cells, dtype=bool)
                vis_mask[visible] = True

                _t0 = time.perf_counter()
                strf_full, folds, t_lag, n_valid = compute_strf_cells_folds(
                    trials, n_roi, n_cells, filter_length_s, pre_time, lag_dt_s)
                print('    [timing] compute_strf_cells_folds ({} {} geom={}, {} '
                      'trial-roi entries) took {:.2f}s'.format(
                          sn, ch, glabel, len(trials), time.perf_counter() - _t0))

                lag_dt = float(t_lag[1] - t_lag[0]) if len(t_lag) > 1 else raw_stim_dt

                if not args.deconvolve:
                    strf_dec, lam_used = strf_full, None
                    print('    not deconvolved (default). Positional results are '
                          'unaffected for a separable filter; peak lags read late by '
                          'roughly the hold_tau. Pass --deconvolve for absolute latency.')
                else:
                    _t0 = time.perf_counter()
                    strf_dec, lam_used = deconvolve_octaves(
                        strf_full, layer, lag_dt, lam=args.lam,
                        lam_scale=args.lam_scale, cv_pair=folds)
                    print('    [timing] deconvolve took {:.2f}s'.format(
                        time.perf_counter() - _t0))
                    _i = lam_used[layer[0][0]]
                    _msg = '    {}'.format(_i['how'])
                    if _i.get('lam') is not None:
                        _msg += '  lam={:.3g}'.format(_i['lam'])
                    if 'split_half_raw' in _i:
                        _msg += '  split-half {:.3f} -> {:.3f}'.format(
                            _i['split_half_raw'], _i['split_half_deconvolved'])
                    print(_msg)
                    # lam, not the 'rejected' prefix: deconvolve_octaves has THREE
                    # pass-through exits -- the reliability guard, 'cv unusable (no
                    # signal)' (few driven cells, which this stimulus reaches easily
                    # with only a couple of tens of them) and 'none (returned
                    # un-deconvolved)' when there are no folds to cross-validate on.
                    # All three leave the STA untouched and all three leave lam None.
                    if _i.get('lam') is None:
                        print('        -> kept the STA. It is still convolved with '
                              'exp(-|t|/{:.2f}s), so read its peak lag as late.'
                              .format(meta['hold_tau']))

                strf_z, z_info = zscore_from_noise_cells(strf_dec, vis_mask)
                if z_info.get('n_sd_fallback'):
                    print('    WARNING: {} of {} rois had no usable noise '
                          'estimate ({}); their z is raw filter units, not z'
                          .format(z_info['n_sd_fallback'], n_roi, z_info['how']))

                # Significance: a circular-shift null on the SAME estimator being
                # reported, i.e. the un-deconvolved STA -- --deconvolve changes
                # latency, not whether a cell counts as signal. See util.nullcal.
                # Restricted to the VISIBLE cells, on BOTH sides. The reported
                # statistic is a max, so every never-displayed column is one more
                # pure-noise comparison it maxes over -- here roughly half of them.
                # Including them on both sides stays VALID (the null inherits the
                # same columns) but needlessly raises the threshold; dropping them
                # from both is the same test over fewer comparisons, which is
                # simply more powerful. What would be wrong is restricting one
                # side only.
                null_p, null_threshold = None, None
                keep_rois, kept_signs, null_select_info = None, None, {}
                mixed_polarity = False
                # Bound here, not only inside the try below: keep_rois is set
                # before it, so a failure in between would leave the saved
                # 'minority_fraction' reading a name that does not exist.
                minority_frac = 0.0
                # Inherit, or test. Every channel but the selection channel
                # takes the kept set that channel's null produced for THIS
                # family (same geometry, same hold_tau) and runs no null of
                # its own; its signs are still its own. See the module
                # docstring, WHICH ROIS ARE KEPT.
                # Lags searched for the peak, on the rois and on every null
                # draw alike; the noise sd still comes from the late lags of
                # the full axis. See nullcal.DEFAULT_PEAK_WINDOW_S.
                _sep_null, _sep_va, _sep_vb = None, None, None
                _pl = nullcal.peak_lags_for(
                    peak_window_s,
                    float(t_lag[1] - t_lag[0]) if len(t_lag) > 1 else lag_dt_s,
                    strf_full.shape[-1])
                _inherit_from = (select_ch if select_ch is not None and ch != select_ch
                                 else None)
                if null_shifts > 0 and _inherit_from is not None:
                    _src = strf_results.get((sn, _inherit_from, glabel))
                    null_select_info = {'method_used': 'inherited',
                                        'from_channel': _inherit_from, 'n_roi': n_roi}
                    if _src is None or _src.get('keep_rois') is None:
                        null_select_info['warning'] = (
                            'no kept set on {} for this family (its null did not '
                            'run); nothing inherited'.format(_inherit_from))
                        print('    ' + null_select_info['warning'])
                    else:
                        keep_rois = np.asarray(_src['keep_rois'], dtype=int)
                        kept_signs = np.asarray(nullcal.peak_signs(
                            strf_full[:, vis_mask, :], keep_rois, peak_lags=_pl),
                            dtype=int)
                        null_select_info['from_method'] = str(
                            _src.get('null_select_info', {}).get('method_used', '?'))
                        print('    kept {} of {} rois: the set {} selected for this '
                              'family ({}); no null run on {}'.format(
                                  len(keep_rois), n_roi, _inherit_from,
                                  null_select_info['from_method'], ch))
                elif null_shifts > 0:
                    _sdt = float(trials[0]['stim_dt'])

                    def _shifted(shift_s, _tr=trials, _sdt=_sdt, _nc=n_cells):
                        k = int(round(shift_s / _sdt))
                        _rolled = {}
                        for t in _tr:
                            if id(t['stim_cells']) not in _rolled:
                                _rolled[id(t['stim_cells'])] = np.roll(
                                    t['stim_cells'], k, axis=0)
                        return compute_strf_cells(
                            [dict(t, stim_cells=_rolled[id(t['stim_cells'])])
                             for t in _tr],
                            n_roi, _nc, filter_length_s, pre_time,
                            lag_dt_s)[0][:, vis_mask, :]

                    _dur = trials[0]['stim_cells'].shape[0] * _sdt
                    _t1 = time.perf_counter()
                    try:
                        # The separability test rides on the same shifts:
                        # its statistic is evaluated on every shifted filter,
                        # with the test vectors (each fold's second spatial
                        # profile) computed once here. See util.separability.
                        _sep_extra = None
                        if do_sep and folds is not None:
                            _sep_va, _sep_vb = separability.test_vectors(
                                folds[0][:, vis_mask, :], folds[1][:, vis_mask, :],
                                strf_full[:, vis_mask, :], _pl)
                            _sep_extra = separability.extra_stats(
                                _sep_va, _sep_vb, peak_lags=_pl)
                        _cal = nullcal.calibrate(
                            _shifted, nullcal.shifts_for(_dur, filter_length_s, null_shifts),
                            peak_lags=_pl, extra_stats=_sep_extra)
                        _sep_null = _cal.get('null_extra') or None
                        _vis_strf = strf_full[:, vis_mask, :]
                        _z = nullcal.peak_z(_vis_strf, peak_lags=_pl)
                        null_p = np.full(n_roi, np.nan)
                        _ok = np.isfinite(_z)
                        null_p[_ok] = nullcal.pvalues(_z[_ok], _cal['null'])[0]
                        null_threshold = _cal['threshold']
                        print('    [timing] null ({} shifts) took {:.2f}s; '
                              'threshold {:.2f}; peak searched over {} of {} lags'.format(
                                  null_shifts, time.perf_counter() - _t1, null_threshold,
                                  strf_full.shape[-1] if _pl is None else _pl,
                                  strf_full.shape[-1]))

                        # Turn the null into an actual kept set. Default is FDR
                        # at q=0.05, the rule analyze_data_trf_fullfield.py uses:
                        # about one in twenty KEPT rois may be noise. It was
                        # briefly replaced (2026-11-16) by the uncorrected
                        # percentile on the belief that FDR was family-wise and
                        # over-strict; on fly_022 it is not (78 of 90 kept where
                        # there is signal, 0 where there is none), while the
                        # percentile admitted ~11 noise rois per family, so it
                        # went back the same day. See --null_method.
                        keep_rois, kept_signs, null_select_info = nullcal.select_rois(
                            _vis_strf, _cal, method=null_method)
                        kept_signs = np.asarray(kept_signs, dtype=int)
                        _how_used = null_select_info.get('method_used', '?')
                        print('    kept {} of {} rois by {}'.format(
                            len(keep_rois), n_roi,
                            ('z > 95th percentile of the null, p < 0.05 per roi, '
                             'uncorrected' if _how_used == 'percentile'
                             else 'FDR (q=0.05)' if _how_used == 'fdr'
                             else _how_used)))
                        if 'warning' in null_select_info:
                            print('        {}'.format(null_select_info['warning']))
                    except ValueError as e:
                        print('    null calibration skipped: {}'.format(e))

                # POLARITY AND THE REGISTRATION LAG, from the FIRST lobe
                # rather than the extremum -- see nullcal, section POLARITY,
                # for the rule and the measurements behind it. kept_signs is
                # replaced; the extremum signs travel alongside for comparison.
                # Read from strf_z: strf_dec, the array the centroid is taken
                # from, over each roi's own noise sd -- so the lag it returns
                # indexes the same axis, and each roi weighs in the family
                # template by its SNR (see first_lobe_polarity).
                polarity = None
                if keep_rois is not None and len(keep_rois):
                    polarity = nullcal.first_lobe_polarity(
                        strf_z[:, vis_mask, :], keep_rois)
                    _ext = np.asarray(polarity['extremum_signs'], dtype=int)
                    kept_signs = np.asarray(polarity['signs'], dtype=int)
                    _tl = polarity['template_lobe']
                    print('    polarity by {} over {} kept rois: {} flipped vs the '
                          'extremum rule ({} pos / {} neg there)'.format(
                              polarity['rule_used'], len(keep_rois),
                              polarity['n_flipped_vs_extremum'],
                              int((_ext > 0).sum()), int((_ext < 0).sum())))
                    if _tl is not None:
                        _c = np.asarray(polarity['confidence'], dtype=float)
                        _src = np.asarray(polarity['lag_source']).astype(str)
                        print('        template first lobe peaks {:.0f} ms (run '
                              '{:.0f}-{:.0f} ms); |cos| to template median {:.2f}, '
                              '{} rois below 0.3; first-lobe lag from own lobe '
                              'for {} rois, from the template for {}'.format(
                                  t_lag[_tl[2]] * 1000, t_lag[_tl[0]] * 1000,
                                  t_lag[min(_tl[1], len(t_lag) - 1)] * 1000,
                                  float(np.nanmedian(_c)), int(np.sum(_c < 0.3)),
                                  int(np.sum(_src == 'own_lobe')),
                                  int(np.sum(_src == 'template'))))
                    else:
                        print('        too few kept rois for a template; each '
                              'roi by its own first lobe')

                if kept_signs is not None:
                    # Checked, never used to flip traces before averaging:
                    # aligning each ROI to its own peak makes noise add
                    # constructively there, so the mean grows a peak even when
                    # nothing carries signal. See util.nullcal.
                    #
                    # The MINORITY FRACTION, not 'are both signs present'. A
                    # population with two dissenters out of 243 is uniform for
                    # every practical purpose, and a warning that fires on it
                    # is one an operator learns to ignore. On an inherited set
                    # these are THIS channel's signs over the other channel's
                    # cells, which is exactly the comparison wanted.
                    _npos = int((kept_signs > 0).sum())
                    _nneg = int((kept_signs < 0).sum())
                    minority_frac = (min(_npos, _nneg) / float(len(kept_signs))
                                     if len(kept_signs) else 0.0)
                    mixed_polarity = bool(minority_frac > 0.1)
                    print('        polarity: {} positive, {} negative '
                          '({:.1%} minority)'.format(
                              _npos, _nneg, minority_frac))
                    if mixed_polarity:
                        print('        NOTE: genuinely MIXED in polarity; '
                              'split them rather than averaging together.')

                # SPACE-TIME SEPARABILITY of the kept rois, judged across the
                # two trial folds (util.separability): raw STA, visible cells,
                # the null's lag window. The second-component test needs the
                # null's extra draws, so a channel that inherited its kept set
                # and ran no null of its own gets the fractions only.
                sep, fold_blocks = None, None
                _n_win = int(_pl) if _pl is not None else int(strf_full.shape[-1])
                # The per-fold blocks are computed and saved only WITH the
                # test: a split-half analysis works from half the data, and
                # whether to pay that is decided at group-analysis time, from
                # the full STAs saved here, not ahead of it for every fly.
                if do_sep and folds is not None:
                    fold_blocks = np.stack([
                        separability.window(folds[0][:, vis_mask, :], _n_win),
                        separability.window(folds[1][:, vis_mask, :], _n_win)]
                    ).astype(np.float32)
                if (do_sep and folds is not None and keep_rois is not None
                        and len(keep_rois)):
                    _t1 = time.perf_counter()
                    sep = separability.analyze(
                        folds[0][:, vis_mask, :], folds[1][:, vis_mask, :],
                        strf_full[:, vis_mask, :], keep_rois, _pl,
                        v_a=_sep_va, v_b=_sep_vb, null=_sep_null)
                    _f = sep['frac'][keep_rois]
                    _sz = sep['second_z'][keep_rois]
                    _ins = sep['inseparable'][keep_rois]
                    _verdict = ('second-component test not run (no null on this '
                                'channel)' if _sep_null is None else
                                '{} of {} inseparable (second-component peak z '
                                'passes BH q=0.05; null 95th pct {:.2f})'.format(
                                    int((_ins == 1).sum()), len(keep_rois),
                                    float(np.percentile(_sep_null['sep_second'], 95))
                                    if np.asarray(_sep_null.get('sep_second', [])).size
                                    else float('nan')))
                    print('    separability over {} kept rois: rank-1 fraction median '
                          '{:.2f} (IQR {:.2f}-{:.2f}), rank-2 increment median {:.3f}; '
                          'second-component peak z median {:.1f} (90th {:.1f}); {}  '
                          '[{:.2f}s]'.format(
                              len(keep_rois), float(np.nanmedian(_f[:, 0])),
                              float(np.nanpercentile(_f[:, 0], 25)),
                              float(np.nanpercentile(_f[:, 0], 75)),
                              float(np.nanmedian(_f[:, 1] - _f[:, 0])),
                              float(np.nanmedian(_sz)), float(np.nanpercentile(_sz, 90)),
                              _verdict, time.perf_counter() - _t1))
                elif do_sep and folds is None:
                    print('    separability skipped: a single trial, no folds')

                dof = independent_samples(n_valid, float(meta['hold_tau']), sample_dt)
                n_indep = float(np.median(dof))
                n_vis = int(vis_mask.sum())
                print('    -> {} {} geom={}: STRF {} (n_roi, n_cells, n_lags) from {} '
                      'trial-roi entries; tau={:.2f}s  {} of {} cells visible  '
                      '{:.0f} independent samples  noise sd={:.4g}'.format(
                          sn, ch, glabel, strf_dec.shape, len(trials), meta['hold_tau'],
                          n_vis, n_cells, n_indep, z_info['median_sd']))
                if n_indep < n_vis:
                    print('       UNDERDETERMINED: fewer independent stimulus '
                          'configurations than visible cells, so the per-cell estimate '
                          'is poor however clean the recording is. Needs about {:.0f}x '
                          'more time.'.format(n_vis / max(n_indep, 1e-9)))

                # Visible cells only. The never-displayed columns are pure STA
                # noise at the same per-cell variance as the real ones, and here
                # they can outnumber the driven cells -- bars about a zenith pole
                # span 180 deg while the patch covers about 91. Everything
                # positional is read AT this lag, so a lag set by noise corrupts
                # center_deg, spread_deg and the saved peak_lag_ind together.
                # Repeated trials, if the protocol ran with repeat_stimulus=True.
                # Detected rather than declared: nothing in the run parameters has
                # to be trusted, since the seed actually recorded per trial says
                # whether two trials really showed the same frames.
                _t0 = time.perf_counter()
                trial_rel = trial_reliability(
                    trials, n_roi, n_cells, strf_full, filter_length_s,
                    pre_time, visible_mask=vis_mask, lag_dt=lag_dt_s)
                if trial_rel is None:
                    print('       no repeated stimulus: trial-by-trial '
                          'reliability skipped')
                else:
                    print('       [timing] trial_reliability took {:.2f}s'.format(
                        time.perf_counter() - _t0))
                    # Over the KEPT rois, the way the full-field twin's
                    # `show = keep if len(keep) else ...` already does it. The
                    # unrestricted median mixes in the rois the null just
                    # rejected -- 46 of 289 on this fly's ch1, which pulled the
                    # reported filter_loo_r from 0.794 down to 0.757 -- and on a
                    # channel where nothing passed it reports the reliability of
                    # pure noise in the same field as a measurement. The
                    # population is named in the line so the number can never be
                    # read against the wrong denominator.
                    _rel_kept = keep_rois is not None and len(keep_rois) > 0
                    _rel_sel = (np.asarray(keep_rois, dtype=int) if _rel_kept
                                else np.arange(n_roi))
                    _rel_lbl = ('{} kept rois'.format(len(_rel_sel)) if _rel_kept
                                else 'all {} rois -- none passed the null'
                                     .format(n_roi))
                    for _g in trial_rel:
                        _nc = _g['noise_ceiling'][_rel_sel]
                        print('       repeated stimulus, seed {}: {} trials '
                              '(epochs {})'.format(_g['seed'], len(_g['trials']),
                                                   list(_g['trials'])))
                        print('         over {}:'.format(_rel_lbl))
                        print('         response_pair_r median {:.3f}   '
                              'filter_loo_r median {:.3f}   '
                              'noise ceiling median {:.3f} ({} of {} '
                              'resolved)'.format(
                                  float(np.nanmedian(
                                      _g['response_pair_r'][_rel_sel])),
                                  float(np.nanmedian(
                                      _g['filter_loo_r'][:, _rel_sel])),
                                  float(np.nanmedian(_nc)),
                                  int(np.isfinite(_nc).sum()), len(_rel_sel)))

                pk = np.array([peak_lag(strf_dec[r][vis_mask]) for r in range(n_roi)])
                # For kept rois peak_lag_ind is the FIRST-LOBE lag (polarity
                # block above); the energy peak stays for the rest. It is a
                # description of the filter's timing -- the centroid below is
                # not read at this one lag but across all of them.
                pk_energy = pk.copy()
                if polarity is not None:
                    pk[np.asarray(keep_rois, dtype=int)] = np.asarray(
                        polarity['lag_ind'], dtype=int)

                # THE SPATIAL PROFILE THAT IS CENTROIDED: the filter projected
                # onto the family template across lags, not one lag's slice.
                # For each cell, its own time course times the template, summed
                # over the template's support. Where a cell's time course has
                # the family's shape every lag adds with the same sign -- the
                # rebound reinforces the first lobe instead of competing with
                # it -- while noise adds in quadrature, so this is the
                # highest-SNR profile a separable filter allows and there is no
                # lag to choose. Why not a slice: at the first-lobe peak a roi
                # whose rebound dominates has only a sliver left after the
                # hold's smear (fly_022, tau 0.4: z 2-4 at 33 ms against 6-9 at
                # the rebound) and its half-data centroids landed 70-110 deg
                # apart; at the energy peak a weak roi can land on a noise lag
                # anywhere on the axis. The weights are unit norm, so the
                # profile keeps one slice's nominal noise sd and centroid_1d's
                # noise cut applies unchanged; they are the same for every
                # cell, so a roi whose kinetics differ from the family's loses
                # efficiency, not position. Assumes only that each cell's time
                # course is a scaled copy of the template; a surround with
                # other kinetics comes out attenuated, which the centroid --
                # the dominant-sign lobe -- never used anyway. A kept roi in a
                # family too small for a template uses its own pooled trace
                # over its own first lobe; a roi the null did not keep has
                # neither and uses its energy-peak slice.
                _pol_idx = ({int(_r): _i for _i, _r in enumerate(keep_rois)}
                            if polarity is not None else {})
                _tw = (nullcal.template_lag_weights(polarity['template'])
                       if polarity is not None else None)

                def _lag_weights(r):
                    i = _pol_idx.get(int(r))
                    if i is not None and _tw is not None:
                        return _tw
                    lw = np.zeros(strf_dec.shape[-1])
                    if i is not None:
                        a, b = polarity['runs'][i]
                        seg = np.abs(polarity['traces'][i, a:b])
                        nrm = float(np.sqrt(np.sum(seg ** 2)))
                        if nrm > 0:
                            lw[a:b] = seg / nrm
                            return lw
                    lw[pk_energy[r]] = 1.0
                    return lw

                def _centroid_profile(block, lw):
                    # the dominant-sign lobe of the profile, the convention
                    # centroid_1d documents; None when there is nothing to read
                    w = np.where(vis_mask, block @ lw, 0.0)
                    if not np.any(w):
                        return None
                    sign = np.sign(w[np.argmax(np.abs(w))])
                    return np.clip(w * sign, 0.0, None)

                _noise_cells = ~vis_mask
                _noise_ratio = []

                def _profile_noise(block, lw, nominal):
                    # THE NOISE OF THE PROFILE, not of one filter value. The
                    # nominal sd is one cell at one lag. A profile summed over
                    # lags with unit-norm weights keeps that sd only if the
                    # noise were independent across lags, and it is not: the
                    # stimulus hold puts its autocorrelation on every STA
                    # column, and the template's lobes are long runs of
                    # same-sign weights, so neighbouring lags add their shared
                    # noise in step. Measured on fly_022 the profile's real sd
                    # was 4.4-4.9x nominal at hold_tau 0.4 and about 2x at 0.04,
                    # and two thirds of pure-noise cells cleared a '2 sd' cut
                    # at 0.4 -- the cut was admitting any noise cell that
                    # touched the receptive field. So measure it: the same
                    # weights applied to the never-displayed cells give profile
                    # values that are noise with exactly these statistics, no
                    # model of the correlation needed. A one-lag profile IS a
                    # filter value and keeps the nominal sd, which is estimated
                    # from far more samples; so does a family with too few
                    # never-displayed cells to estimate anything from.
                    if (nominal is None or np.count_nonzero(lw) <= 1
                            or int(_noise_cells.sum()) < 3):
                        return nominal
                    v = block[_noise_cells] @ lw
                    v = v[np.isfinite(v)]
                    if v.size < 3:
                        return nominal
                    sd = float(np.sqrt(np.mean(v ** 2)))
                    return sd if sd > 0 else nominal

                def _fold_jitter(cf):
                    # sd of the difference between two half-data centroids is
                    # sqrt(2) x a half's error; the full-data estimate uses
                    # twice the data, so its own sd is that difference sd over 2.
                    d = cf[0] - cf[1]
                    if str(meta['mode']) == 'wedges':
                        d = ((d + 180.0) % 360.0) - 180.0
                    d = d[np.asarray(_use, dtype=int)]
                    d = d[np.isfinite(d)]
                    return (float(np.std(d) / 2.0) if d.size > 1 else np.nan), int(d.size)

                _use = (keep_rois if keep_rois is not None and len(keep_rois)
                        else np.arange(n_roi))
                centers = np.full(n_roi, np.nan)
                centers_screen = np.full(n_roi, np.nan)
                centers_fold = [np.full(n_roi, np.nan), np.full(n_roi, np.nan)]
                spreads = np.full(n_roi, np.nan)
                n_kept = np.zeros(n_roi, dtype=int)
                for r in range(n_roi):
                    lw = _lag_weights(r)
                    lobe = _centroid_profile(strf_dec[r], lw)
                    if lobe is None:
                        continue
                    # Noise-relative inclusion, using THIS roi's own noise: a
                    # peak-relative cut snaps the estimate to a cell center
                    # whenever the lobe is narrow. See centroid_1d. The noise
                    # is that of the profile being centroided (_profile_noise).
                    _nsd_nom = z_info['sd'][r] if 'sd' in z_info else None
                    _nsd = _profile_noise(strf_dec[r], lw, _nsd_nom)
                    if (_nsd_nom and _nsd is not None and _nsd != _nsd_nom
                            and int(r) in _pol_idx):
                        _noise_ratio.append(_nsd / _nsd_nom)
                    centers[r], spreads[r], n_kept[r] = centroid_1d(
                        lobe, positions, meta['mode'], cell_widths_deg=widths,
                        noise_sd=_nsd)
                    # Computed in the screen frame rather than converted from the
                    # pole frame: for wedges the circular mean has to be taken in
                    # the frame it will be reported in, or a lobe spanning the
                    # pole-frame seam averages to the wrong place before the
                    # conversion can fix it.
                    centers_screen[r], _sp, _nk = centroid_1d(
                        lobe, offsets, meta['mode'], cell_widths_deg=widths,
                        noise_sd=_nsd)
                    # REGISTRATION JITTER. The mean filter is the true one
                    # convolved with the distribution of registration errors, so
                    # a noisy centroid makes the population average look wider
                    # than any cell in it. Centroid each independent half of the
                    # trials at the SAME lag -- the question is spatial, not
                    # temporal -- and the spread of the difference measures that
                    # error without assuming anything about it.
                    if folds is not None:
                        for _fi, _fa in enumerate(folds):
                            _lf = _centroid_profile(_fa[r], lw)
                            if _lf is None:
                                continue
                            # each half's own profile noise: half the data
                            # is noisier, and its cut should say so
                            centers_fold[_fi][r], _s2, _n2 = centroid_1d(
                                _lf, offsets, meta['mode'], cell_widths_deg=widths,
                                noise_sd=_profile_noise(_fa[r], lw, _nsd))

                jitter_deg, _n_pairs = _fold_jitter(centers_fold)
                if _noise_ratio:
                    print('    matched-profile noise is {:.2f}x the single-value sd '
                          '(median over {} kept rois, from the never-displayed '
                          'cells); the centroid cut uses it'.format(
                              float(np.median(_noise_ratio)), len(_noise_ratio)))
                if np.isfinite(jitter_deg):
                    print('    registration jitter {:.2f} deg (from {} split-half '
                          'centroid pairs, matched profile)'.format(
                              jitter_deg, _n_pairs))

                strf_results[sn, ch, glabel] = {
                    'strf_raw': strf_full, 'strf': strf_dec, 'strf_z': strf_z,
                    't_lag': t_lag, 'n_valid': n_valid, 'lam': lam_used,
                    'meta': meta, 'positions_deg': positions, 'widths_deg': widths,
                    'filter_length_s': float(filter_length_s),
                    'screen_offset_deg': offsets,
                    'center_screen_deg': centers_screen,
                    'center_fold_deg': np.vstack(centers_fold),
                    'jitter_deg': jitter_deg,
                    'center_defined': bool(center_defined),
                    'axis_label': screen_axis_label(meta, sampling),
                    'center_az_deg': float(sampling.get('center_az_deg', 0.0)),
                    'center_el_deg': float(sampling.get('center_el_deg', 0.0)),
                    'runt_cells': runts, 'visible_cells': visible,
                    'peak_lag_ind': pk, 'center_deg': centers, 'spread_deg': spreads,
                    'n_kept_cells': n_kept, 'dof': dof, 'z_info': z_info,
                    'deconv_applied': bool(
                        lam_used is not None
                        and lam_used[layer[0][0]].get('lam') is not None),
                    'deconv_how': ('not requested (default)' if lam_used is None
                                   else lam_used[layer[0][0]]['how']),
                    'null_p': null_p, 'null_threshold': null_threshold,
                    'keep_rois': keep_rois, 'kept_signs': kept_signs,
                    'polarity': polarity,
                    'polarity_confidence': (polarity['confidence']
                                            if polarity is not None else None),
                    'null_select_info': null_select_info,
                    'mixed_polarity': mixed_polarity,
                    'minority_fraction': (minority_frac
                                          if keep_rois is not None else 0.0),
                    'trial_rel': trial_rel,
                    'sep': sep, 'fold_blocks': fold_blocks,
                    'sep_window_lags': _n_win,
                }

    # ── save ─────────────────────────────────────────────────────────────────
    # Nothing new to write when the results came out of this very file.
    if not reused:
        with h5py.File(experiment_file_path, 'a') as f:
            if strf_group in f:
                del f[strf_group]
            g = f.create_group(strf_group)
            g.attrs['filter_length_spec'] = filter_length_spec
            # the numbers the reuse check actually compares; the spec above is
            # kept for the log and for files read by older code
            g.attrs['filter_length_by_ch'] = json.dumps(
                {k: float(v) for k, v in sorted(filter_length_by_ch.items())})
            g.attrs['lag_dt'] = float(lag_dt_s)
            g.attrs['dff'] = dff
            g.attrs['deconvolved'] = bool(args.deconvolve)
            g.attrs['null_method'] = null_method
            g.attrs['roi_selection_channel'] = ('each' if select_ch is None
                                                else select_ch)
            g.attrs['polarity_rule'] = POLARITY_RULE
            g.attrs['registration_rule'] = REGISTRATION_RULE
            g.attrs['peak_window_s'] = (float(peak_window_s)
                                        if peak_window_s > 0 else 0.0)

            for (sn, ch, glabel), res in strf_results.items():
                sg = g.require_group('{}/{}/geom_{}'.format(sn, ch, glabel))
                # strf_raw is the un-deconvolved STA. Without --deconvolve
                # it IS strf -- the same object, not merely equal -- so writing
                # both puts two identical gzipped copies of a
                # (n_roi, n_cells, n_lags) array in the file, which on this
                # rig's first fly is about 100 MB of pure duplicate. Write it
                # only when it differs, and say so where it would have been.
                for key in ('strf', 'strf_z'):
                    sg.create_dataset(key, data=res[key], compression='gzip')
                if res['strf_raw'] is res['strf']:
                    sg.attrs['strf_raw'] = ('omitted: identical to strf '
                                            '(not deconvolved)')
                else:
                    sg.create_dataset('strf_raw', data=res['strf_raw'],
                                      compression='gzip')
                sg.create_dataset('t_lag', data=res['t_lag'])
                sg.create_dataset('n_valid', data=res['n_valid'])
                sg.create_dataset('peak_lag_ind', data=res['peak_lag_ind'])
                sg['peak_lag_ind'].attrs['rule'] = (
                    'kept rois: peak of the FIRST lobe (first_lobe_lag_ind); '
                    'other rois: lag of greatest summed squared filter over '
                    'the visible cells. A timing descriptor only: the centroid '
                    'of a kept roi is read from the filter projected onto the '
                    'family template across lags (group attr registration_rule)')
                sg.create_dataset('positions_deg', data=res['positions_deg'])
                sg.create_dataset('screen_offset_deg', data=res['screen_offset_deg'])
                sg.create_dataset('center_screen_deg', data=res['center_screen_deg'])
                sg.create_dataset('center_fold_deg', data=res['center_fold_deg'])
                sg.attrs['registration_jitter_deg'] = (
                    res['jitter_deg'] if np.isfinite(res['jitter_deg']) else -1.0)
                sg.attrs['screen_frame_defined'] = res['center_defined']
                sg.attrs['screen_axis_label'] = res['axis_label']
                sg.attrs['center_az_deg'] = res['center_az_deg']
                sg.attrs['center_el_deg'] = res['center_el_deg']
                sg.create_dataset('widths_deg', data=res['widths_deg'])
                sg.create_dataset('runt_cells', data=res['runt_cells'])
                sg.create_dataset('visible_cells', data=res['visible_cells'])
                sg.create_dataset('center_deg', data=res['center_deg'])
                sg.create_dataset('spread_deg', data=res['spread_deg'])
                sg.create_dataset('n_cells_in_centroid', data=res['n_kept_cells'])
                sg.create_dataset('dof', data=res['dof'])
                # Separability (util.separability), kept rois only, NaN
                # elsewhere; cells are the visible ones in visible_cells order
                # and lags the leading window_lags of t_lag.
                _sep = res.get('sep')
                if _sep is not None:
                    _g = sg.create_group('separability')
                    for _k in ('frac', 'cv_corr', 'second_z', 'second_p', 'second_z_dir',
                               'second_p_dir', 'inseparable', 'second_trace',
                               'test_vector', 'energy_frac', 'profile', 'time_course'):
                        _g.create_dataset(_k, data=np.asarray(_sep[_k]))
                    _g.attrs['bh_cutoff'] = float(_sep['bh_cutoff'])
                    _g.attrs['n_tested'] = int(_sep['n_tested'])
                    _g.attrs['window_lags'] = int(res['sep_window_lags'])
                    _g.attrs['cells'] = 'visible cells, visible_cells order'
                    _g.attrs['frac'] = ('cross-validated <A_k, B>/<A, B> for rank k = 1, 2, '
                                        'averaged over both fold directions; raw STA; '
                                        'biased low below peak z ~7')
                    _g.attrs['second'] = (
                        'second-component test: each fold\'s second spatial profile '
                        'with the full block\'s rank-1 profile projected out '
                        '(test_vector, rows A and B) applied to the OTHER fold, the '
                        'two traces summed sign-aligned (second_trace); second_z = '
                        'its peak |z| (late-lag noise, peak window), second_p against '
                        'the shift null of the same statistic; inseparable = 1 where '
                        'BH q=0.05 passes, 0 tested and not, -1 untested. '
                        'second_z_dir / second_p_dir are the two directions alone.')
                    _g.attrs['null'] = (
                        'shift null shared with roi selection'
                        if np.isfinite(np.asarray(_sep['second_p'])).any()
                        else 'none: no null ran on this channel (inherited kept set)')
                if res.get('fold_blocks') is not None:
                    sg.create_dataset('fold_blocks', data=res['fold_blocks'],
                                      compression='gzip')
                    sg['fold_blocks'].attrs['layout'] = (
                        '(fold, roi, visible cell, window lag): raw STA of each trial '
                        'fold over the separability window, float32')
                meta = res['meta']
                for k in ('mode', 'cell_width_deg', 'phase_deg', 'pole_az_deg',
                         'pole_el_deg', 'hold_tau', 'p', 'seed'):
                    sg.attrs[k] = meta[k]
                # The WHOLE stimulus record as well. The eight attrs above are
                # what the figures label things by, but the geometry helpers
                # read more -- in particular frame_az_deg / frame_el_deg, the
                # screen's pitch, which pole_frame uses to turn a recorded pole
                # into a fly-frame direction. A reuse run that rebuilt meta
                # from the eight attrs alone had no frame and silently drew
                # the rf-centre map with the screen taken as level: on this
                # rig that is 38 degrees of pitch, and the map was wrong in
                # every reuse run while right in every fresh one.
                sg.attrs['one_d_stimulus_json'] = json.dumps(
                    {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                     for k, v in dict(meta).items()}, default=float)
                sg.attrs['filter_length_s'] = res['filter_length_s']
                sg.attrs['median_noise_sd'] = res['z_info']['median_sd']
                # Which pool the sd came from; the two are not comparable.
                sg.attrs['noise_pool'] = res['z_info']['how']
                if res['null_p'] is not None:
                    sg.create_dataset('null_p', data=res['null_p'])
                    sg.attrs['null_threshold'] = res['null_threshold']
                if res['keep_rois'] is not None:
                    # Same names analyze_data_trf_fullfield.py writes.
                    sg.create_dataset('kept_rois',
                                      data=np.asarray(res['keep_rois'], dtype=int))
                    sg.create_dataset('kept_signs',
                                      data=np.asarray(res['kept_signs'], dtype=int))
                    _pol = res.get('polarity')
                    if _pol is not None:
                        # Positional within kept_rois, like kept_signs. The
                        # extremum signs are what this pipeline reported before
                        # 2026-11-17, kept so a reader can see what changed.
                        sg.create_dataset('kept_signs_extremum', data=np.asarray(
                            _pol['extremum_signs'], dtype=int))
                        sg.create_dataset('first_lobe_lag_ind', data=np.asarray(
                            _pol['lag_ind'], dtype=int))
                        sg.create_dataset('polarity_confidence', data=np.asarray(
                            _pol['confidence'], dtype=float))
                        sg.create_dataset('first_lobe_amp_ratio', data=np.asarray(
                            _pol['amp_ratio'], dtype=float))
                        sg.create_dataset('first_lobe_lag_source', data=np.asarray(
                            [str(x) for x in _pol['lag_source']], dtype='S'))
                        sg.create_dataset('polarity_template', data=np.asarray(
                            _pol['template'], dtype=float))
                        sg.attrs['polarity_rule'] = str(_pol['rule_used'])
                        sg.attrs['polarity_template_lobe'] = (
                            np.asarray(_pol['template_lobe'], dtype=int)
                            if _pol['template_lobe'] is not None
                            else np.array([-1, -1, -1]))
                        sg.attrs['polarity_n_flipped_vs_extremum'] = int(
                            _pol['n_flipped_vs_extremum'])
                        sg.attrs['polarity_smooth_samples'] = int(_pol['smooth'])
                        sg.attrs['polarity_z_min'] = float(_pol['z_min'])
                        sg.attrs['polarity_definitions'] = (
                            'kept_signs: sign of the FIRST lobe, by projection '
                            'of each roi\'s pooled time course (nullcal.lobe_trace) '
                            'onto the family\'s sign-aligned template, oriented by '
                            'the template\'s first lobe; kept_signs_extremum: sign '
                            'at the |z| extremum (the pre-2026-11-17 rule); '
                            'polarity_confidence: |cos| between trace and template; '
                            'first_lobe_lag_ind: lag index of the first-lobe peak '
                            '(= peak_lag_ind for kept rois); first_lobe_amp_ratio: '
                            '|trace| at that lag over max |trace|, < 1 when a later '
                            'lobe is bigger; first_lobe_lag_source: own_lobe, '
                            'template or energy_peak; polarity_template: the '
                            'smoothed, late-lag-zeroed template (n_lags); '
                            'polarity_template_lobe: (start, stop, peak) lag indices '
                            'of its first lobe, -1 if none')
                    sg.attrs['kept_mixed_polarity'] = res['mixed_polarity']
                    sg.attrs['kept_minority_fraction'] = res['minority_fraction']
                    sg.attrs['selection_basis'] = (
                        'visible cells only; kept set inherited from {}, no null '
                        'run on this channel'.format(
                            res['null_select_info']['from_channel'])
                        if res['null_select_info'].get('method_used') == 'inherited'
                        else 'visible cells only')
                    # 'null_select_', not 'null_': select_rois reports its own
                    # 'threshold' (the weakest ROI it kept), which under a shared
                    # prefix would silently overwrite the null threshold above.
                    for _k, _v in res['null_select_info'].items():
                        sg.attrs['null_select_' + _k] = _v
                # Trial-by-trial reliability, same group layout and the same measure
                # names analyze_data_trf_fullfield.py writes, so one reader serves both.
                if res['trial_rel'] is None:
                    sg.attrs['trial_reliability'] = 'skipped: no repeated stimulus'
                else:
                    sg.attrs['trial_reliability'] = '{} repeat group(s)'.format(
                        len(res['trial_rel']))
                    for _gi, _grp in enumerate(res['trial_rel']):
                        tg = sg.create_group('trial_reliability/group_{}'.format(_gi))
                        tg.attrs['seed'] = _grp['seed']
                        tg.attrs['definitions'] = (
                            'filter_loo_r: each trial STRF vs mean of the others, over '
                            'the visible cells and support_lags; trial_gain: LS scale '
                            'onto the full STRF; response_pair_r: mean corr between '
                            'single-trial responses (model-free); noise_ceiling: '
                            'sqrt(response_pair_r), NaN unless > 2 SE')
                        for key in ('trials', 'support_lags', 'filter_loo_r',
                                    'trial_gain', 'response_pair_r',
                                    'response_pair_r_se', 'response_loo_r',
                                    'noise_ceiling', 'full_time_course'):
                            tg.create_dataset(key, data=_grp[key])
                        tg.create_dataset('trial_time_courses',
                                          data=_grp['trial_time_courses'],
                                          compression='gzip')
                sg.attrs['deconv_how'] = res['deconv_how']
                sg.attrs['deconv_applied'] = res['deconv_applied']
                if res['lam'] is not None:
                    _i = res['lam'][next(iter(res['lam']))]
                    if _i.get('lam') is not None:
                        sg.attrs['deconv_lam'] = float(_i['lam'])
                    if 'split_half_raw' in _i:
                        sg.attrs['split_half_raw'] = _i['split_half_raw']
                        sg.attrs['split_half_deconvolved'] = _i['split_half_deconvolved']

        print('\nSaved STRFs to {} under /{}'.format(experiment_file_path, strf_group))


    # ── registered mean across rois: what does this cell type's filter look like ──
    # Each roi is shifted onto its OWN receptive-field center and the population
    # averaged, which is the only way a mean means anything: averaging unshifted
    # filters of cells that look at different places gives the average of the
    # POPULATION's coverage, not of its filter.
    #
    # Derived from strf_results, so it costs nothing on a reuse run and is
    # rebuilt whenever the figures are refreshed. z units, never sign-aligned --
    # see one_d_group for both reasons.
    # ── the channel's majority polarity, if asked ────────────────────────
    # Pooled over every family of the channel rather than decided per family:
    # the sign belongs to the cell type, not to the orientation it was measured
    # with, and an roi reading one way in one family and the other way in the
    # other is telling you something about that roi rather than about Mi1.
    majority_sign = {}
    if args.majority_sign_only:
        _pool = {}
        for (sn, ch, glabel), res in sorted(strf_results.items()):
            _s = res.get('kept_signs')
            if _s is not None and len(_s):
                _pool.setdefault(ch, []).append(np.asarray(_s, int))
        for ch, _arrs in sorted(_pool.items()):
            _a = np.concatenate(_arrs)
            _np_, _nn = int((_a > 0).sum()), int((_a < 0).sum())
            majority_sign[ch] = 1 if _np_ >= _nn else -1
            print('  {}: majority polarity {} ({} positive, {} negative over all '
                  'families)'.format(ch, '+' if majority_sign[ch] > 0 else '-',
                                     _np_, _nn))

    # ── rois held back from the edge of the observed extent ──────────────
    # ONE window, derived from a single reference channel and applied to every
    # channel unchanged. Deriving it per channel would put each in a different
    # physical region -- the channels keep different roi sets -- and a
    # difference in region would then read as a difference in kinetics, which
    # is exactly the comparison this is meant to protect.
    #
    # The extent uses the 2nd/98th percentile rather than min/max: the raw
    # extremes are two rois, so one mislocalised cell stretches the span and
    # drags the whole window with it.
    # WHICH CHANNEL'S CENTRES place an roi. The selection channel -- the
    # structural, calcium one by default -- defines where a cell's receptive
    # field is. A channel that inherited its kept set from it is registered,
    # localised and trimmed at THAT channel's centroid of the same roi, never
    # at its own: its own centroid is a measurement of a weaker indicator, and
    # for an roi whose filter is mostly noise it is a random position, so
    # registering there smears the population mean spatially and scatters the
    # map. Same cell, same receptive field; the position is not a property of
    # the indicator. Falls back to the channel's own centres only when the
    # selection channel has no result for the family.
    def _centre_res(sn, ch, glabel, res):
        if select_ch is None or ch == select_ch:
            return res, ch
        _r = strf_results.get((sn, select_ch, glabel))
        if _r is None or len(_r.get('center_deg', [])) != len(res['center_deg']):
            return res, ch
        return _r, select_ch

    edge_sets = {}
    _margins = [float(x) for x in str(args.edge_margins).replace(';', ',').split(',')
                if x.strip() != '']
    if _margins:
        _pos, _byc = {}, {}
        for (sn, ch, glabel), res in sorted(strf_results.items()):
            if str(res['meta']['mode']) in ('bars', 'wedges'):
                # the result whose CENTRES localise this channel's rois
                _byc.setdefault(ch, []).append(
                    (glabel, _centre_res(sn, ch, glabel, res)[0]))
        for ch, fams in sorted(_byc.items()):
            if len(fams) < 2:
                continue
            (gl_a, ra), (gl_b, rb) = fams[0], fams[1]
            _caz = 0.5 * (float(ra.get('center_az_deg', 0.0))
                          + float(rb.get('center_az_deg', 0.0)))
            _cel = 0.5 * (float(ra.get('center_el_deg', 0.0))
                          + float(rb.get('center_el_deg', 0.0)))
            _pref = np.array([
                np.cos(np.radians(_cel)) * np.cos(np.radians(_caz)),
                np.cos(np.radians(_cel)) * np.sin(np.radians(_caz)),
                np.sin(np.radians(_cel))])
            _n = min(len(ra['center_deg']), len(rb['center_deg']))
            _d = direction_from_two_families(ra['meta'], ra['center_deg'][:_n],
                                             rb['meta'], rb['center_deg'][:_n],
                                             prefer=_pref)
            # Relative to the patch centre by ROTATION, not az/el subtraction:
            # see one_d_stim.screen_azel_of.
            _az, _el = screen_azel_of(_d, _caz, _cel)
            _ka, _kb = ra.get('keep_rois'), rb.get('keep_rois')
            _e = (np.arange(_n) if _ka is None or _kb is None
                  else np.intersect1d(np.asarray(_ka, int), np.asarray(_kb, int)))
            _e = _e[_e < _n]
            _e = _e[np.isfinite(_az[_e]) & np.isfinite(_el[_e])]
            if _e.size >= 8:
                _pos[ch] = (_e, _az, _el)
        if _pos:
            _ref = max(_pos, key=lambda c: _pos[c][0].size)
            _e, _az, _el = _pos[_ref]
            _abox = np.percentile(_az[_e], [2, 98])
            _ebox = np.percentile(_el[_e], [2, 98])
            print('  edge margins: window from {} ({} localised rois); observed '
                  'extent az {:+.1f}..{:+.1f}, el {:+.1f}..{:+.1f} deg'.format(
                      _ref, _e.size, _abox[0], _abox[1], _ebox[0], _ebox[1]))
            for _m in _margins:
                # MARGIN 0 MEANS NO TRIM. The percentile extent is the reference
                # the other margins are measured from, not a filter in its own
                # right, so applying it at zero would silently drop the outer 2%
                # per side and make the baseline something other than 'no
                # trimming'. Kept separate so the series isolates one thing per
                # step: all -> m0 is the cost of requiring a 2D position (which
                # needs the roi to pass the null in BOTH families), and m0 -> m10
                # -> m20 is the cost of trimming.
                _lo_a, _hi_a = _abox[0] + _m, _abox[1] - _m
                _lo_e, _hi_e = _ebox[0] + _m, _ebox[1] - _m
                if _m > 0 and (_hi_a <= _lo_a or _hi_e <= _lo_e):
                    print('    margin {:g} deg: window is empty; skipped'.format(_m))
                    continue
                _per_ch = {}
                for ch, (_ec, _azc, _elc) in sorted(_pos.items()):
                    _sel = _ec if _m <= 0 else _ec[
                        (_azc[_ec] >= _lo_a) & (_azc[_ec] <= _hi_a)
                        & (_elc[_ec] >= _lo_e) & (_elc[_ec] <= _hi_e)]
                    _per_ch[ch] = set(int(r) for r in _sel)
                    print('    margin {:>4g} deg  {}: {:>4} of {:<4} rois '
                          '({:.0%})'.format(_m, ch, len(_sel), _ec.size,
                                            len(_sel) / max(_ec.size, 1)))
                edge_sets[_m] = _per_ch

    avg_results = {}
    for (sn, ch, glabel), res in sorted(strf_results.items()):
        keep = res.get('keep_rois')
        if keep is None or not len(keep):
            print('  {} {} {}: no kept rois, no registered mean'.format(
                sn, ch, glabel))
            continue
        keep = np.asarray(keep, dtype=int)
        vis = np.asarray(res['visible_cells'], dtype=int)
        signs = np.asarray(res.get('kept_signs', np.ones(len(keep))), dtype=int)
        groups, minority = polarity_groups(signs)
        _cres, _cch = _centre_res(sn, ch, glabel, res)
        if _cch != ch:
            print('  {} {} {}: registered at the {} centroid of each roi'.format(
                sn, ch, glabel, _cch))
        # Parent first, then its trimmed subsets, so each can be pinned
        # to the parent's bin width. auto_grid_step is data-driven and half as
        # many rois earns coarser bins -- correct on its own terms, but it puts
        # the two means on different grids, and comparing them is the whole
        # reason the subset exists.
        # The polarity cut is applied to EVERY positional variant, not only to
        # the full population: 'central' and 'correct sign' are independent
        # claims about an roi, and stacking them is the only way to see whether
        # they are removing the same cells or different ones. Built in two
        # passes -- positional first, then sign on each -- so the naming reads
        # as what was done, in order: all_m10_maj is the 10 deg margin with the
        # dissenters dropped.
        _base = [(gname, sel, None) for gname, sel in groups.items()]
        for _mg in sorted(edge_sets):
            _cen = edge_sets[_mg].get(ch)
            if not _cen:
                continue
            for gname, sel in list(groups.items()):
                _m = np.array([int(r) in _cen for r in keep[sel]], dtype=bool)
                if _m.sum() >= 3:
                    _base.append(('{}_m{:g}'.format(gname, _mg), sel[_m], gname))

        _variants = list(_base)
        if ch in majority_sign:
            _want = majority_sign[ch]
            _sg = res.get('kept_signs')
            if _sg is not None and len(_sg) == len(keep):
                _sg = np.asarray(_sg, int)
                _same = []
                for _nm, _sel, _par in _base:
                    # Only on the pooled group. 'positive' and 'negative' are
                    # already polarity-pure, so the cut is either a no-op or
                    # removes everything -- and the everything case fell below
                    # the 3-roi floor and vanished with no message, leaving a
                    # fly that happened to split missing a variant its siblings
                    # had. The split already answers what _maj asks.
                    if (_par if _par is not None else _nm).split('_m')[0] != 'all':
                        continue
                    _m = _sg[_sel] == _want
                    # every variant is pinned to the BASE group's bin width, so
                    # all of them land on one grid and can be read against each
                    # other rather than each choosing its own by occupancy
                    if _m.sum() >= 3:
                        _variants.append(('{}_maj'.format(_nm), _sel[_m],
                                          _par if _par is not None else _nm))
                        # Written even when it removes nothing. A family with no
                        # dissenters would otherwise be the one family missing an
                        # _maj file, and anything reading a fixed set of variant
                        # names -- a cross-fly step, a figure panel -- would have
                        # to special-case it. Asking for the cut and getting the
                        # same rois back is a result; a missing file is a gap.
                        if _m.sum() == _sel.size:
                            _same.append(_nm)
                if _same:
                    print('    {} {} {}: no dissenters, so _maj is identical to '
                          '{}'.format(sn, ch, glabel, ', '.join(_same)))
        _step_of = {}
        for gname, sel, _parent in _variants:
            rois = keep[sel]
            if rois.size < 3:
                continue
            out = register_average(
                res['strf_z'][np.ix_(rois, vis)],
                np.asarray(res['screen_offset_deg'])[vis],
                np.asarray(_cres['center_screen_deg'])[rois],
                float(res['meta']['cell_width_deg']),
                mode=str(res['meta']['mode']),
                grid_step=_step_of.get(_parent),
                # each cell's OWN width: a wedge seam cell can be a partial one
                cell_widths_deg=np.asarray(res['widths_deg'], dtype=float)[vis])
            _step_of[gname] = out['grid_step_deg']
            out.update(jitter_deg=float(_cres.get('jitter_deg', np.nan)),
                       centres_from=_cch,
                       polarity=gname, minority_fraction=float(minority),
                       rois=rois, axis_label=res.get('axis_label', ''),
                       t_lag=res['t_lag'],
                       cell_width_deg=float(res['meta']['cell_width_deg']))
            avg_results[sn, ch, glabel, gname] = out
            _j = out['jitter_deg']
            print('  {} {} {} [{}]: registered mean over {} rois, '
                  '{:g} deg bins, {} bins{}'.format(
                      sn, ch, glabel, gname, out['n_roi'], out['grid_step_deg'],
                      int(np.sum(np.any(np.isfinite(out['mean']), axis=1))),
                      '' if not np.isfinite(_j)
                      else ', registration jitter {:.2f} deg'.format(_j)))

    # ── radial profile, if asked: one h(r) against every family at once ──────
    radial_results = {}
    if args.radial_symmetry:
        # The DOMINANT polarity group of each family, not a hardcoded 'positive'.
        # Which sign dominates is a property of the cell type and of where the
        # filter peaks -- a biphasic kernel read at its peak lag can be mostly
        # negative -- so naming one is a coin flip that silently drops the
        # population when it lands wrong.
        #
        # Families are pooled per CONDITION (hold_tau, p), not per channel
        # alone: a tau sweep gives two families about each pole, and one h(r)
        # fit over all of them would be asked to explain two timescales at
        # once. With a single condition the key is the bare channel name, so
        # the hdf5 path and figure name are what they always were; with
        # several it carries the tau.
        _cond_of = {}
        for (sn, ch, glabel), res in strf_results.items():
            _cond_of[glabel] = (round(float(res['meta']['hold_tau']), 6),
                                round(float(res['meta']['p']), 6))
        _n_cond = len(set(_cond_of.values()))
        by_ch = {}
        for ch in sorted({k[1] for k in avg_results}):
            fams = {}
            for (sn, c, glabel, gname), out in avg_results.items():
                if c != ch:
                    continue
                best = fams.get(glabel)
                if best is None or out['n_roi'] > best[1]['n_roi']:
                    fams[glabel] = (glabel, out, gname)
            for glabel in sorted(fams):
                _tau = _cond_of[glabel][0]
                _key = ch if _n_cond == 1 else '{}_tau{:g}'.format(ch, _tau)
                by_ch.setdefault(_key, []).append(fams[glabel])
        for ch, fams in sorted(by_ch.items()):
            if len(fams) < 2:
                print('  {}: radial symmetry needs two families, have {} -- '
                      'a single projection cannot even test it'.format(ch, len(fams)))
                continue
            # Two families whose DOMINANT polarity groups differ are not two
            # views of one population -- they are different rois. The gain
            # cannot paper over it (it is clamped non-negative), so this would
            # surface as a residual near 1; saying so here names the cause.
            _pols = {f[2] for f in fams}
            if len(_pols) > 1:
                print('  {}: families disagree on dominant polarity ({}) -- '
                      'they describe different roi populations, so the radial '
                      'fit below is not meaningful'.format(
                          ch, ', '.join('{}={}'.format(f[0][:16], f[2])
                                        for f in fams)))
            n_lags = min(f[1]['mean'].shape[1] for f in fams)
            profs = [f[1]['mean'][:, :n_lags] for f in fams]
            grids = [f[1]['grid_deg'] for f in fams]
            width = float(fams[0][1]['cell_width_deg'])
            # smooth_width_deg: the registered means spread each cell over its
            # width, so the forward model must too, or h(r) comes out a box
            # too wide. See register_average / projection_matrix.
            rad = radial_from_projections(profs, grids, width,
                                          smooth_width_deg=width)
            rad.update(families=[f[0] for f in fams], t_lag=fams[0][1]['t_lag'][:n_lags])
            radial_results[ch] = rad
            med = float(np.nanmedian(rad['residual']))
            print('  {}: radial fit over {} families, median residual {:.3f}{}'
                  .format(ch, len(fams), med,
                          '  -- LARGE: the field is not radially symmetric'
                          if med > 0.15 else ''))
            # The per-family amplitudes are fitted, not assumed equal, because
            # two roi populations have no common scale. The values are worth
            # reading: far from 1 means the two families disagree about how
            # loud the same cell type is, which is a recording or selection
            # difference rather than anything about the receptive field.
            _g = np.asarray(rad['gain'], dtype=float)
            print('     family gains (re. {}): {}'.format(
                rad['families'][0][:22],
                '  '.join('{}={:.2f}'.format(f[:16], v)
                          for f, v in zip(rad['families'], _g))))
            if np.nanmin(_g) <= 0.0:
                print('     NOTE: a family was clamped to zero gain -- it could '
                      'not be reconciled with the others at any positive scale, '
                      'which usually means opposite polarity, not elongation')
            elif np.nanmax(_g) / np.nanmin(_g) > 3.0:
                print('     NOTE: more than 3x apart -- check the kept-roi sets '
                      'before reading the residual as elongation')

    # ── save both ───────────────────────────────────────────────────────────
    avg_group = 'STRF_ONE_D_AVG_' + tag
    if avg_results:
        with h5py.File(experiment_file_path, 'a') as f:
            if avg_group in f:
                del f[avg_group]
            g = f.create_group(avg_group)
            g.attrs['registered_on'] = (
                'each roi at its center_screen_deg on the channel named by '
                'centers_from: the selection (structural) channel for a channel '
                'that inherited its kept set, otherwise its own')
            g.attrs['registration'] = (
                'each cell spread over its own width, overlap-weighted per bin '
                '(one_d_group.register_average); n_roi_per_bin is summed weight')
            g.attrs['units'] = 'z (per-roi noise)'
            g.attrs['replicate'] = 'roi'
            for (sn, ch, glabel, gname), out in avg_results.items():
                sg = g.require_group('{}/{}/geom_{}/{}'.format(sn, ch, glabel, gname))
                for k in ('grid_deg', 'mean', 'sem', 'n_roi_per_bin', 'rois'):
                    sg.create_dataset(k, data=out[k])
                sg.create_dataset('t_lag', data=out['t_lag'])
                sg.attrs['n_roi'] = out['n_roi']
                sg.attrs['centers_from'] = str(out.get('centres_from', ch))
                sg.attrs['grid_step_deg'] = out['grid_step_deg']
                sg.attrs['minority_fraction'] = out['minority_fraction']
                sg.attrs['registration_jitter_deg'] = (
                    out['jitter_deg'] if np.isfinite(out['jitter_deg']) else -1.0)
                sg.attrs['jitter_note'] = (
                    'the mean is broadened by this in quadrature: an apparent '
                    'width w reads as sqrt(true^2 + jitter^2 + cell^2/12), the '
                    'last term from spreading each cell over its width')
                sg.attrs['axis_label'] = out['axis_label']
            for ch, rad in radial_results.items():
                rg = g.require_group('{}/radial'.format(ch))
                rg.create_dataset('r_centers_deg', data=rad['r_centers_deg'])
                rg.create_dataset('h', data=rad['h'])
                rg.create_dataset('residual', data=np.atleast_1d(rad['residual']))
                rg.create_dataset('residual_per_family',
                                  data=np.atleast_2d(rad['residual_per_family']))
                rg.create_dataset('family_gain',
                                  data=np.atleast_1d(rad['gain']))
                rg.attrs['families'] = ', '.join(rad['families'])
                rg.attrs['note'] = ('h(r) constrained by every family at once; '
                                    'residual is the test of radial symmetry')
        print('\nSaved registered means to {} under /{}'.format(
            experiment_file_path, avg_group))

    # ── figures ──────────────────────────────────────────────────────────────
    if save_figs or show_figs:
        figs_dir = os.path.join(experiment_file_directory,
                                tag + '_' + FIGURE_SUBDIR)
        os.makedirs(figs_dir, exist_ok=True)
        # Per-roi files one level down. A few hundred kept rois give roughly
        # 1500 of them against a dozen population figures, and the population
        # ones -- the registered mean, the radial fit, the rf map, the repeat
        # summary -- are what gets read first. Mixing them buries the answer
        # in the evidence.
        rois_dir = os.path.join(figs_dir, 'rois')
        os.makedirs(rois_dir, exist_ok=True)

        # ── where the kept rois sit in the volume ──────────────────────────
        # One figure per family, one panel per slice: the mean image under
        # every roi of the aligned set, the kept ones filled in colour and the
        # rejected ones a thin outline -- a grey FILL vanishes on the bright
        # neuropil the rois sit on -- so the question 'is the responsive
        # population a layer, a column, or scattered noise?' is answered by
        # looking rather than by cross-referencing roi indices. Colour is
        # POLARITY, majority blue and minority red, the same convention as the
        # rf-centre map below, so a kept roi of the dissenting sign can be
        # placed too. Masks are read from the file here rather than kept from
        # the analysis pass, because a reused run never made that pass.
        _masks_cache = {}

        def _roi_masks(sn, ch):
            if (sn, ch) not in _masks_cache:
                _rd = roi_data.get((sn, ch))
                if _rd is None or 'roi_mask' not in _rd:
                    plug.updateImagingDataObject(
                        experiment_file_directory, experiment_file_name,
                        int(sn.replace('sn', '')))
                    _rd = plug.ImagingDataObject.getRoiMasks(
                        response_set_name_prefix + ch, roi_prefix='aligned')
                _masks_cache[(sn, ch)] = (_rd['roi_mask'], _rd['roi_image'])
            return _masks_cache[(sn, ch)]

        _slice_cache = {}

        def _roi_slices(sn, ch, n_roi):
            """The slice each roi mostly lives in (most voxels), or None."""
            if (sn, ch) not in _slice_cache:
                try:
                    _raw, _img = _roi_masks(sn, ch)
                    _b = roi_bool_masks(_raw, n_roi)
                    _slice_cache[(sn, ch)] = np.array(
                        [int(np.argmax(m.sum(axis=(0, 1)))) for m in _b])
                except (KeyError, AssertionError, ValueError) as _e:
                    print('  {} {}: roi slices unavailable ({})'.format(sn, ch, _e))
                    _slice_cache[(sn, ch)] = None
            return _slice_cache[(sn, ch)]

        for (sn, ch, glabel), res in sorted(strf_results.items()):
            _keep = res.get('keep_rois')
            if _keep is None:
                continue                      # no null ran: nothing to show
            _n_roi = res['strf'].shape[0]
            try:
                _raw_mask, _img = _roi_masks(sn, ch)
                _bool = roi_bool_masks(_raw_mask, _n_roi)
            except (KeyError, AssertionError, ValueError) as _e:
                print('  {} {} {}: no kept-roi map ({})'.format(sn, ch, glabel, _e))
                continue
            _img = np.asarray(_img, dtype=float)
            if _img.ndim == 2:
                _img = _img[:, :, np.newaxis]
            _n_z = _bool[0].shape[-1]
            if _img.shape[-1] != _n_z:
                print('  {} {} {}: no kept-roi map (image has {} slices, masks '
                      '{})'.format(sn, ch, glabel, _img.shape[-1], _n_z))
                continue
            _keep = np.asarray(_keep, dtype=int)
            _signs = res.get('kept_signs')
            _signs = (np.asarray(_signs, dtype=int) if _signs is not None
                      else np.ones(len(_keep), dtype=int))
            _maj = 1 if (_signs > 0).sum() >= (_signs < 0).sum() else -1
            _col = {1: (0.12, 0.47, 0.71), -1: (0.84, 0.15, 0.16)}
            _kept_set = {int(r): int(s) for r, s in zip(_keep, _signs)}
            # Contrast from the whole volume, not per slice, so a dim slice
            # reads as dim instead of being stretched to look like the others.
            _lo, _hi = np.nanpercentile(_img, [1.0, 99.5])
            _hi = _hi if _hi > _lo else _lo + 1.0

            # Panel shape follows the image: this rig's slices are several
            # times taller than wide, and a square panel for each leaves
            # most of the page blank. Width is floored so the title fits.
            _aspect = _img.shape[1] / float(_img.shape[0])     # imshow: cols/rows
            _ph = 4.2
            _pw = max(2.4, _ph * _aspect)
            _ncol = min(_n_z, 5)
            _nrow = int(np.ceil(_n_z / float(_ncol)))
            fig, axes = plt.subplots(_nrow, _ncol, squeeze=False,
                                     figsize=(_pw * _ncol, _ph * _nrow + 0.9),
                                     constrained_layout=True)
            for _z in range(_n_z):
                ax = axes[_z // _ncol, _z % _ncol]
                ax.imshow(np.clip((_img[:, :, _z] - _lo) / (_hi - _lo), 0, 1),
                          cmap='gray', vmin=0, vmax=1, interpolation='nearest')
                _ov = np.zeros(_img.shape[:2] + (4,))
                _rej = np.zeros(_img.shape[:2], dtype=bool)
                _n_here, _n_kept_here = 0, 0
                for _r in range(_n_roi):
                    _m = _bool[_r][:, :, _z]
                    if not _m.any():
                        continue
                    _n_here += 1
                    if _r in _kept_set:
                        _n_kept_here += 1
                        _ov[_m] = _col[_kept_set[_r] * _maj] + (0.65,)
                    else:
                        _rej |= _m
                ax.imshow(_ov, interpolation='nearest')
                if _rej.any():
                    # One contour over the union: neighbouring rejected rois
                    # share edges, and per-roi outlines would draw them twice.
                    ax.contour(_rej.astype(float), levels=[0.5],
                               colors=['#e0e0e0'], linewidths=0.35, alpha=0.9)
                ax.set_title('slice {}: {} of {} rois kept'.format(
                    _z, _n_kept_here, _n_here), fontsize=12)
                ax.set_axis_off()
            for _z in range(_n_z, _nrow * _ncol):
                axes[_z // _ncol, _z % _ncol].set_axis_off()
            fig.suptitle('{} {} {}\nkept {} of {} rois  (blue: majority sign, '
                         'red: minority sign, outlined: rejected by the null)'.format(
                             sn, ch, glabel, len(_keep), _n_roi), fontsize=13)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_KEPT_MAP_{}_{}_{}.pdf'.format(
                        sn, ch, glabel)))
            if not show_figs:
                plt.close(fig)

        # ── SEPARABILITY, per family: a summary and a grid of example rois ──
        # See util.separability. The fraction is biased low at low SNR, so
        # the summary draws it against each roi's peak z; the yes/no is the
        # residual test. The grid overlays every driven cell's own time
        # course (normalised, oriented by the rank-1 profile) on the rank-1
        # time course: a separable roi collapses onto one curve.
        for (sn, ch, glabel), res in strf_results.items():
            _sep = res.get('sep')
            _keep = res.get('keep_rois')
            if _sep is None or _keep is None or not len(_keep):
                continue
            _keep = np.asarray(_keep, dtype=int)
            _vis = np.asarray(res['visible_cells'], dtype=int)
            _tl = np.asarray(res['t_lag'])
            _nw = int(res.get('sep_window_lags', len(_tl)))
            _frac = np.asarray(_sep['frac'], dtype=float)
            _f1 = _frac[_keep, 0]
            _inc = _frac[_keep, 1] - _f1
            _rz = np.asarray(_sep['second_z'], dtype=float)[_keep]
            _flag = np.asarray(_sep['inseparable'])[_keep]
            _pz = np.nanmax(np.abs(np.asarray(res['strf_z'])[_keep][:, _vis, :_nw]),
                            axis=(1, 2))
            _tested = bool(np.isfinite(np.asarray(_sep['second_p'], dtype=float)).any())
            _ins = _flag == 1
            _col = np.where(_ins, '#d62728', '#1f77b4')

            fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6), constrained_layout=True)
            ax = axes[0]
            _bins = np.linspace(0, 1.1, 23)
            ax.hist([_f1[~_ins], _f1[_ins]], bins=_bins, stacked=True,
                    color=['#1f77b4', '#d62728'],
                    label=['separable / untested', 'inseparable (second-component test)'])
            ax.set_xlabel('cross-validated rank-1 fraction')
            ax.set_ylabel('kept rois')
            ax.legend(fontsize=10, loc='upper left')
            ax.set_title('one profile x one time course accounts for', fontsize=12)
            ax = axes[1]
            ax.scatter(_pz, _f1, s=26, c=_col, edgecolors='k', linewidths=0.4)
            ax.set_xlabel('roi peak |z| (window)')
            ax.set_ylabel('rank-1 fraction')
            ax.set_ylim(min(0.0, float(np.nanmin(_f1)) if np.isfinite(_f1).any() else 0.0),
                        1.15)
            ax.axhline(1.0, color='0.6', lw=0.8, ls=':')
            ax.set_title('fraction is biased low at low z: read against z', fontsize=12)
            ax = axes[2]
            ax.scatter(_rz, _inc, s=26, c=_col, edgecolors='k', linewidths=0.4)
            ax.set_xlabel('second-component peak |z| (cross-fold, summed)')
            ax.set_ylabel('rank-2 increment of the fraction')
            ax.axhline(0.0, color='0.6', lw=0.8, ls=':')
            ax.axvline(float(np.nanpercentile(_rz, 50)), color='0.85', lw=0.6)
            ax.set_title('second component: energy vs significance', fontsize=12)
            fig.suptitle('{} {} {}\nseparability over {} kept rois: {}'.format(
                sn, ch, glabel, len(_keep),
                ('{} inseparable by the second-component test (BH q=0.05)'.format(
                    int(_ins.sum())) if _tested
                 else 'second-component test not run on this channel (no null)')),
                fontsize=13)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_SEPARABILITY_{}_{}_{}.pdf'.format(
                        sn, ch, glabel)))
            if not show_figs:
                plt.close(fig)

            # Example rois: every inseparable one (strongest residual first,
            # up to 8) plus the highest-z separable ones (up to 4).
            _ord_ins = _keep[_ins][np.argsort(-_rz[_ins])][:8]
            _ord_sep = _keep[~_ins][np.argsort(-_pz[~_ins])][:4]
            _show = list(_ord_ins) + list(_ord_sep)
            if not _show:
                continue
            _ncol = 3 if len(_show) > 4 else min(len(_show), 4)
            _nrow = int(np.ceil(len(_show) / float(_ncol)))
            fig, axes = plt.subplots(_nrow, _ncol, squeeze=False,
                                     figsize=(4.6 * _ncol, 3.4 * _nrow + 0.8),
                                     constrained_layout=True)
            _off = np.asarray(res['screen_offset_deg'], dtype=float)[_vis]
            _cmap = plt.get_cmap('viridis')
            _lo, _hi = float(np.nanmin(_off)), float(np.nanmax(_off))
            _strf = np.asarray(res.get('strf_raw', res['strf']))
            _prof = np.asarray(_sep['profile'], dtype=float)
            _tc = np.asarray(_sep['time_course'], dtype=float)
            _ef = np.asarray(_sep['energy_frac'], dtype=float)
            _str = np.asarray(_sep['second_trace'], dtype=float)
            for _i, _r in enumerate(_show):
                ax = axes[_i // _ncol, _i % _ncol]
                _blk = _strf[_r][_vis, :_nw]
                _zb = np.asarray(res['strf_z'])[_r][_vis, :_nw]
                _cells = np.flatnonzero(np.nanmax(np.abs(_zb), axis=1) > 2.5)
                for _c in _cells:
                    _y = _blk[_c] * (np.sign(_prof[0, _r, _c]) or 1.0)
                    _m = float(np.nanmax(np.abs(_y)))
                    if not _m > 0:
                        continue
                    _t = 0.5 if _hi <= _lo else (_off[_c] - _lo) / (_hi - _lo)
                    ax.plot(_tl[:_nw], _y / _m, color=_cmap(_t), lw=0.9, alpha=0.85)
                _t1 = _tc[0, _r]
                _m1 = float(np.nanmax(np.abs(_t1))) or 1.0
                ax.plot(_tl[:_nw], _t1 / _m1, color='k', lw=2.2, label='rank-1 time course')
                _t2 = _str[_r]
                if np.isfinite(_t2).any() and np.nanmax(np.abs(_t2)) > 0:
                    ax.plot(_tl[:_nw], _t2 / float(np.nanmax(np.abs(_t2))), color='0.45',
                            lw=1.4, ls='--',
                            label='tested 2nd-component trace (cross-fold)')
                ax.axhline(0, color='0.7', lw=0.6)
                _k = int(np.flatnonzero(_keep == _r)[0])
                ax.set_title('roi {}  z {:.1f}  frac {:.2f} (+{:.2f} rank 2)\n'
                             '2nd-component z {:.1f}  {}'.format(
                                 _r, _pz[_k], _f1[_k], _inc[_k], _rz[_k],
                                 'INSEPARABLE' if _ins[_k] else
                                 ('separable' if _tested else 'untested')),
                             fontsize=10,
                             color='#d62728' if _ins[_k] else 'k')
                ax.set_xlabel('lag (s)')
                if _i % _ncol == 0:
                    ax.set_ylabel('normalised filter')
                if _i == 0:
                    ax.legend(fontsize=8, loc='upper right')
            for _i in range(len(_show), _nrow * _ncol):
                axes[_i // _ncol, _i % _ncol].set_axis_off()
            fig.suptitle('{} {} {}: driven cells (|z| > 2.5, colour = position along the '
                         'axis, oriented by the rank-1 profile) over the rank-1 time '
                         'course'.format(sn, ch, glabel), fontsize=12)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_SEP_ROIS_{}_{}_{}.pdf'.format(
                        sn, ch, glabel)))
            if not show_figs:
                plt.close(fig)

        for (sn, ch, glabel), res in strf_results.items():
            strf, t_lag = res['strf'], res['t_lag']
            # Plotted in the SCREEN frame throughout: 0 is the middle of the
            # display, and for wedges the visible cells are contiguous in it.
            positions, pk = res['screen_offset_deg'], res['peak_lag_ind']
            # The centre marked on every per-roi panel is the one the
            # population analysis registers at: the selection channel's for
            # an inherited channel (see _centre_res), so the figure shows the
            # position the roi is actually averaged at, not a second one.
            centers_plot = _centre_res(sn, ch, glabel, res)[0]['center_screen_deg']
            vis_mask = np.zeros(strf.shape[1], dtype=bool)
            vis_mask[res['visible_cells']] = True
            meta = res['meta']
            n_roi = strf.shape[0]
            axis_label = res.get('axis_label') or pole_label(meta)

            # WHICH ROIS GET A FIGURE. Drawing every roi means drawing the ones
            # the null rejected, which in a channel where nothing passed is a
            # folder of noise -- 578 files on this rig's first real fly. Ranked
            # by peak |filter| over the visible cells so a truncated list is the
            # strongest rois rather than the lowest-numbered ones.
            _keep_fig = res.get('keep_rois')
            _spec = str(args.figure_rois).lower()
            if _spec in ('none', 'no', '0'):
                # Checked before anything else so it always wins. 'none' used to
                # fall through the int() below into the except branch, which set
                # _n_want to the whole ranked list -- i.e. it drew EVERY kept roi,
                # the exact opposite of what it says. An integer 0 worked, but
                # only by accident of arithmetic.
                _roi_list = []
                _why_rois = 'none (--figure_rois {})'.format(_spec)
            elif _spec == 'all':
                _roi_list = list(range(n_roi))
                _why_rois = 'all rois'
            elif _keep_fig is None:
                _roi_list = list(range(n_roi))
                _why_rois = 'all rois (no null selection ran)'
            elif not len(_keep_fig):
                # The null RAN and kept nothing. That is an answer, not an
                # absent selection, and the two were conflated: an empty kept
                # set fell through to 'draw everything', so the channel where
                # nothing passed produced the most figures. On this fly's ch2,
                # 289 rois x 3 figures x 2 geometries = 1734 files of the exact
                # noise the null had just rejected. Draw none; --figure_rois
                # all is there for looking anyway.
                _roi_list = []
                _why_rois = 'none -- the null kept no rois in this channel'
            else:
                _strength = np.array([np.abs(strf[r][vis_mask]).max()
                                      for r in _keep_fig])
                _ranked = [int(r) for r in
                           np.asarray(_keep_fig)[np.argsort(_strength)[::-1]]]
                try:
                    _n_want = int(_spec)
                except ValueError:
                    _n_want = len(_ranked)
                _roi_list = _ranked[:max(_n_want, 0)]
                _why_rois = '{} of {} kept rois, strongest first'.format(
                    len(_roi_list), len(_ranked))
            print('  {} {} {}: {} per-roi figures ({})'.format(
                sn, ch, glabel, 2 * len(_roi_list), _why_rois))

            # ── the STRF itself: angle on y, lag on x, as the noizone analysis
            # draws it. z-scored, so the colour scale means the same thing in
            # every panel, and shared across rois of one geometry so they are
            # directly comparable. Only the cells the projector could show are
            # drawn -- a band of never-displayed cells is a band of pure noise,
            # and on this stimulus it can be half the range.
            # Rows go in SCREEN order, not cell-index order. For wedges the two
            # differ whenever the lit arc crosses the pole frame's seam, which is
            # the common case -- meridians about the zenith at phase 0 give
            # visible cells [0..4, 19..23]. Sorting by offset makes them one run.
            _vis_idx = np.flatnonzero(vis_mask)
            _vis_idx = _vis_idx[np.argsort(positions[_vis_idx])]
            _steps = np.diff(positions[_vis_idx]) if _vis_idx.size > 1 else np.array([1.0])
            _contiguous = bool(np.allclose(_steps, _steps[0], atol=1e-6))
            # THE LIT ARC, shared by both figures so the map's y axis and the
            # profile's x axis are the same range and can be read together.
            # Cells are ordered by screen offset, which for wedges is not cell
            # index order: a family straddling the pole frame's seam is one run
            # here but two in index. Unequal steps (a runt cell inside the arc)
            # would put uneven angular spacing on an even image axis, so that
            # case falls back to every cell and says so in the title.
            _rows = (_vis_idx if (_vis_idx.size and _contiguous)
                     else np.argsort(positions))
            _pos = positions[_rows]
            _wid = res['widths_deg'][_rows]
            # Edge to edge, not center to center. imshow spreads N rows evenly
            # across the extent, so passing cell CENTRES squeezes the axis by
            # (N-1)/N and misplaces every row but the middle one -- 6.4 deg at
            # the ends for 7 fifteen-degree cells, which is most of a cell.
            _lo = float(_pos[0] - _wid[0] / 2.0)
            _hi = float(_pos[-1] + _wid[-1] / 2.0)
            _n_lit = int(vis_mask.sum())
            if _vis_idx.size:
                _img_all = res['strf_z'][:, _rows, :]
                _m = float(np.nanmax(np.abs(_img_all)))
                _vmax = _m if np.isfinite(_m) and _m > 0 else 1.0
                for r in _roi_list:
                    fig, ax = plt.subplots(figsize=(5.8, 5.0), constrained_layout=True)
                    im = ax.imshow(
                        _img_all[r], aspect='auto', origin='lower',
                        extent=[t_lag[0] * 1000, t_lag[-1] * 1000, _lo, _hi],
                        cmap='RdBu_r', vmin=-_vmax, vmax=_vmax)
                    ax.axvline(t_lag[pk[r]] * 1000, color='k', ls=':', lw=1.0)
                    if np.isfinite(centers_plot[r]):
                        ax.axhline(centers_plot[r], color='k', ls='--', lw=0.8)
                    ax.set_xlabel('Lag (ms)')
                    ax.set_ylabel(axis_label)
                    fig.colorbar(im, ax=ax, label='z-score')
                    _t = '1D STRF  {} {} {}  roi {}'.format(sn, ch, glabel, r)
                    if not _contiguous:
                        _t += '  [ALL cells: unequal cell steps]'
                    fig.suptitle(_t, fontsize=15)
                    if save_figs:
                        fig.savefig(os.path.join(
                            rois_dir,
                            'strf_one_d_map_{}_{}_{}_roi{:03d}.pdf'.format(
                                sn, ch, glabel, r)))
                    if not show_figs:
                        plt.close(fig)

            for r in _roi_list:
                fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.7),
                                         constrained_layout=True)
                lag_s = t_lag[pk[r]]

                w = strf[r, :, pk[r]]
                # Symmetric about zero and scaled to the VISIBLE cells, so the
                # filter's polarity is readable and an off-screen noise excursion
                # cannot set the scale. Chosen explicitly rather than with `or`,
                # since np.abs(...).max() is NaN-truthy and would slip past it.
                _peak = np.nanmax(np.abs(w[_rows])) if _rows.size else np.nan
                lim = float(_peak) if np.isfinite(_peak) and _peak > 0 else 1.0
                ax = axes[0]
                # Cropped to the lit arc, on the same axis range as the map.
                # The cells outside it are never displayed, so their coefficients
                # are noise; the title keeps the count so the crop is visible.
                ax.plot(_pos, w[_rows], '-o', ms=3, color='0.15', lw=1.0)
                ax.set_xlim(_lo, _hi)
                ax.axhline(0, color='0.6', lw=0.8)
                if np.isfinite(centers_plot[r]):
                    ax.axvline(centers_plot[r], color='r', ls=':', lw=1.2)
                ax.set_ylim(-1.05 * lim, 1.05 * lim)
                ax.set_xlabel(axis_label)
                ax.set_ylabel('filter (cell units)')
                ax.set_title('{} cells @ lag {:.2f}s  ({} of {} lit)'.format(
                    meta['mode'], lag_s, _n_lit, vis_mask.size), fontsize=15)

                ax = axes[1]
                power = (strf[r] ** 2).sum(axis=1)
                power = np.where(vis_mask, power, -1.0)
                _cmax = 0.0
                for c_i in np.argsort(power)[-5:]:
                    if power[c_i] < 0:
                        continue
                    ax.plot(t_lag, strf[r, c_i], lw=1.0, alpha=0.75)
                    _cmax = max(_cmax, float(np.nanmax(np.abs(strf[r, c_i]))))
                # The trace the polarity and the first-lobe lag were read
                # from: the receptive-field cells pooled with non-negative
                # weights (nullcal.lobe_trace), scaled to the strongest cell so
                # it shares the axis. The dotted line marks its first lobe.
                _tr, _ = nullcal.lobe_trace(strf[r][vis_mask])
                _tmax = float(np.nanmax(np.abs(_tr))) if np.any(_tr) else 0.0
                if _tmax > 0 and _cmax > 0:
                    ax.plot(t_lag, _tr * (_cmax / _tmax), 'k', lw=1.6, alpha=0.9,
                            label='pooled trace (scaled)')
                    ax.legend(fontsize=11, loc='upper right')
                ax.axvline(lag_s, color='k', ls=':', lw=1.0)
                ax.axhline(0, color='0.6', lw=0.8)
                ax.set_xlabel('lag (s)')
                ax.set_ylabel('filter (cell units)')
                _ttl = '5 strongest visible cells  tau={:.2f}s'.format(meta['hold_tau'])
                _pc = res.get('polarity_confidence')
                if _pc is not None and res.get('keep_rois') is not None:
                    _pcm = {int(k): float(v) for k, v in zip(res['keep_rois'], _pc)}
                    if r in _pcm and np.isfinite(_pcm[r]):
                        _ttl += '  polarity |cos| {:.2f}'.format(_pcm[r])
                ax.set_title(_ttl, fontsize=15)

                title = '{} {} geom={}  roi {}'.format(sn, ch, glabel, r)
                if not res['deconv_applied']:
                    title += '  [{}: peak lag reads late]'.format(res['deconv_how'])
                fig.suptitle(title, fontsize=15)
                # Not the default crush: _tighten's 0.01 h_pad leaves the axes
                # titles sitting on top of the suptitle.
                _tighten(fig, w_pad=0.02, h_pad=0.06, wspace=0.06, hspace=0.05)
                if save_figs:
                    fig.savefig(os.path.join(
                        rois_dir,
                        'strf_one_d_{}_{}_{}_roi{:03d}.pdf'.format(sn, ch, glabel, r)))
                if not show_figs:
                    plt.close(fig)

            # ── repeatability, when the stimulus repeated. Same three columns
            # analyze_data_trf_fullfield.py draws: per-trial filters coloured in
            # recording order so drift reads as a colour gradient, then each
            # trial against the rest, then its gain on the all-trial filter.
            if res.get('trial_rel'):
                _grp = res['trial_rel'][0]
                _ntr = len(_grp['trials'])
                _lim = min(n_roi, _grp['trial_time_courses'].shape[1])
                _cols = plt.cm.viridis(np.linspace(0, 0.9, _ntr))

                # ── per roi, one file each ──────────────────────────────────
                # Same roi list as the other per-roi figures, so --figure_rois
                # governs all of them together. Three columns, as the full-field
                # pipeline draws them: the per-trial filters coloured in
                # recording order so drift reads as a colour gradient, then each
                # trial against the mean of the others, then its gain on the
                # all-trial filter.
                for r in [int(_x) for _x in _roi_list if 0 <= int(_x) < _lim]:
                    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.0),
                                             constrained_layout=True)
                    ax = axes[0]
                    for j in range(_ntr):
                        ax.plot(t_lag, _grp['trial_time_courses'][j, r],
                                color=_cols[j], lw=0.9,
                                label='trial {}'.format(int(_grp['trials'][j])))
                    ax.plot(t_lag, _grp['full_time_course'][r], 'k', lw=1.8,
                            label='all trials')
                    ax.axhline(0, color='0.7', lw=0.6)
                    ax.set_xlabel('lag (s)')
                    ax.set_ylabel('filter (cell units)')
                    ax.legend(fontsize=14, ncol=2)
                    ax = axes[1]
                    ax.plot(range(_ntr), _grp['filter_loo_r'][:, r], 'o-')
                    # _ylo, not _lo: _lo is the lit arc's lower bound, set
                    # once per geometry and read by the map extent and the
                    # profile's xlim above. Reusing the name worked only
                    # because it is recomputed next iteration.
                    _ylo = float(np.nanmin(_grp['filter_loo_r'][:, r]))
                    ax.set_ylim(min(-0.1, _ylo - 0.05), 1.02)
                    ax.set_xticks(range(_ntr))
                    ax.set_xlabel('trial (recording order)')
                    ax.set_ylabel('trial vs rest r', fontsize=14)
                    ax = axes[2]
                    ax.plot(range(_ntr), _grp['trial_gain'][:, r], 'o-')
                    ax.axhline(1, color='0.6', lw=0.8, ls=':')
                    ax.set_xticks(range(_ntr))
                    ax.set_xlabel('trial (recording order)')
                    ax.set_ylabel('gain vs all-trial filter', fontsize=14)
                    _nc = _grp['noise_ceiling'][r]
                    fig.suptitle(
                        '{} {} {}  roi {}   repeated stimulus seed {}, {} trials'
                        '   pair r {:.2f}   ceiling {}'.format(
                            sn, ch, glabel, r, _grp['seed'], _ntr,
                            _grp['response_pair_r'][r],
                            'unresolved' if not np.isfinite(_nc)
                            else '{:.2f}'.format(_nc)), fontsize=15)
                    if save_figs:
                        fig.savefig(os.path.join(
                            rois_dir,
                            'strf_one_d_repeat_{}_{}_{}_roi{:03d}.pdf'.format(
                                sn, ch, glabel, r)))
                    if not show_figs:
                        plt.close(fig)

                # ── population summary over every kept roi ──────────────────
                # The per-roi figures do not scale to 243 rois -- nobody reads
                # 243 of them -- so the distribution is the thing that answers
                # "was this recording reliable" for the population at once.
                _pop = res.get('keep_rois')
                _pop = ([int(_x) for _x in _pop if 0 <= int(_x) < _lim]
                        if _pop is not None and len(_pop) else list(range(_lim)))
                if _pop:
                    _pair = _grp['response_pair_r'][_pop]
                    _ceil = _grp['noise_ceiling'][_pop]
                    _loo = np.nanmedian(_grp['filter_loo_r'][:, _pop], axis=0)
                    _gain = _grp['trial_gain'][:, _pop]
                    fig, axes = plt.subplots(1, 4, figsize=(18, 4.1),
                                             constrained_layout=True)
                    ax = axes[0]
                    ax.hist(_pair[np.isfinite(_pair)], bins=30, color='#1f77b4')
                    ax.axvline(float(np.nanmedian(_pair)), color='k', ls='--', lw=1)
                    ax.set_xlabel('response_pair_r')
                    ax.set_ylabel('rois')
                    ax.set_title('single-trial repeatability\nmedian {:.3f}'.format(
                        float(np.nanmedian(_pair))), fontsize=14)
                    ax = axes[1]
                    _fin = _ceil[np.isfinite(_ceil)]
                    if _fin.size:
                        ax.hist(_fin, bins=30, color='#2ca02c')
                        ax.axvline(float(np.median(_fin)), color='k', ls='--', lw=1)
                    ax.set_xlabel('noise_ceiling')
                    ax.set_title('noise ceiling\n{} of {} resolved'.format(
                        int(np.isfinite(_ceil).sum()), len(_pop)), fontsize=14)
                    ax = axes[2]
                    ax.hist(_loo[np.isfinite(_loo)], bins=30, color='#ff7f0e')
                    ax.axvline(float(np.nanmedian(_loo)), color='k', ls='--', lw=1)
                    ax.set_xlabel('median filter_loo_r')
                    ax.set_title('filter consistency\nmedian {:.3f}'.format(
                        float(np.nanmedian(_loo))), fontsize=14)
                    ax = axes[3]
                    _m = np.nanmedian(_gain, axis=1)
                    _q1 = np.nanpercentile(_gain, 25, axis=1)
                    _q3 = np.nanpercentile(_gain, 75, axis=1)
                    ax.fill_between(range(_ntr), _q1, _q3, alpha=0.25,
                                    color='#9467bd', lw=0)
                    ax.plot(range(_ntr), _m, 'o-', color='#9467bd')
                    ax.axhline(1, color='0.6', lw=0.8, ls=':')
                    ax.set_xticks(range(_ntr))
                    ax.set_xlabel('trial (recording order)')
                    ax.set_ylabel('gain vs all-trial filter', fontsize=14)
                    # A trend here is run-down or drift across the recording,
                    # which a per-roi panel shows one roi at a time and a
                    # population median makes obvious.
                    ax.set_title('gain across trials\nmedian and IQR over rois',
                                 fontsize=14)
                    _pop_lbl = ('{} kept rois'.format(len(_pop))
                                if _keep_fig is not None and len(_keep_fig)
                                else 'all {} rois -- none passed the null'
                                     .format(len(_pop)))
                    fig.suptitle('{} {} {}: repeat reliability over {} '
                                 '(seed {}, {} trials)'.format(
                                     sn, ch, glabel, _pop_lbl,
                                     _grp['seed'], _ntr), fontsize=15)
                    if save_figs:
                        fig.savefig(os.path.join(
                            figs_dir,
                            'strf_one_d_repeat_summary_{}_{}_{}.pdf'.format(
                                sn, ch, glabel)))
                    if not show_figs:
                        plt.close(fig)

        # ── receptive-field centers on the screen ────────────────────────
        # One point per roi: two 1D families combined into a position, both
        # already relative to the patch center so they share an origin without
        # any registration.
        #
        # Pairs ACROSS SERIES as well as within one: rois are drawn once per fly
        # and the same masks are applied to every series, so roi index means the
        # same cell whether the two families were interleaved in one series or
        # run as two. Channels are kept apart, since a roi can carry signal in
        # one indicator and not the other.
        #
        # Included only when a roi passed the null in BOTH families. That is a
        # conjunction, which for a COMPARISON of two orientations would bias the
        # sample toward radially symmetric cells -- but a position is not a
        # comparison: without both coordinates there is no point to place.
        _loc = str(args.localize_channel).lower()
        if _loc == 'structural' and struct_channel_num is not None:
            _loc_ch = 'ch' + str(struct_channel_num)
        elif _loc in ('each', 'all'):
            _loc_ch = None
        elif _loc == 'structural':
            _loc_ch = None
            print('  no structural channel in fly.json; mapping every channel')
        else:
            _loc_ch = _loc if _loc.startswith('ch') else 'ch' + _loc
        if _loc_ch is not None:
            print('  receptive-field map localised on {}'.format(_loc_ch))

        # A map is made from two families of the SAME MODE about DIFFERENT
        # poles, within ONE condition (hold_tau, p). Each family pins the
        # direction to a curve and two curves meet in a point, so geometrically
        # bars+wedges about one pole would localise too (see
        # one_d_stim.direction_from_two_families) -- but the two modes probe
        # the field with different cell shapes and are not combined here, by
        # jcsimon's call (2026-11-13). What this rig runs is two WEDGE families
        # about orthogonal poles, or two BARS families likewise.
        _fams = {}
        for (sn, ch, glabel), res in strf_results.items():
            if _loc_ch is not None and ch != _loc_ch:
                continue
            if str(res['meta']['mode']) in ('bars', 'wedges'):
                _fams.setdefault(ch, []).append((sn, glabel, res))

        for ch, fams in sorted(_fams.items()):
            if len(fams) < 2:
                print('  {}: only {} family -- a position needs two'.format(
                    ch, len(fams)))
                continue
            for _i in range(len(fams)):
                for _j in range(_i + 1, len(fams)):
                    (sn_a, gl_a, res_a) = fams[_i]
                    (sn_e, gl_e, res_e) = fams[_j]
                    # Same mode, same condition, different axis -- each test
                    # names its reason, since a run that sweeps pole and tau
                    # together produces every kind of rejected pair.
                    _ma, _me = res_a['meta'], res_e['meta']
                    if str(_ma['mode']) != str(_me['mode']):
                        print('    skip {} x {}: bars and wedges are not '
                              'combined into one map'.format(gl_a, gl_e))
                        continue
                    if (abs(float(_ma['hold_tau']) - float(_me['hold_tau'])) > 1e-6
                            or abs(float(_ma['p']) - float(_me['p'])) > 1e-6):
                        print('    skip {} x {}: different hold_tau / p, one '
                              'map per condition'.format(gl_a, gl_e))
                        continue
                    # The same pole or its antipode gives parallel curves (two
                    # bar rings, or two wedge fans meeting only at the pole):
                    # direction_from_two_families returns NaN for every roi,
                    # and the pair would print and save an empty map.
                    _axis_dot = abs(float(np.dot(pole_frame(_ma)[0],
                                                 pole_frame(_me)[0])))
                    if _axis_dot > np.cos(np.radians(1.0)):
                        print('    skip {} x {}: same axis, parallel families '
                              'pin no point'.format(gl_a, gl_e))
                        continue
                    # Both families', not just the second's. The center is
                    # what the plotted coordinates are relative to and what
                    # breaks the two-root ambiguity, so taking it from one
                    # family of the pair silently mis-centers the whole map if
                    # the patch was ever moved between them.
                    _az_a = float(res_a.get('center_az_deg', 0.0))
                    _el_a = float(res_a.get('center_el_deg', 0.0))
                    _az_e = float(res_e.get('center_az_deg', 0.0))
                    _el_e = float(res_e.get('center_el_deg', 0.0))
                    _daz = ((_az_a - _az_e + 180.0) % 360.0) - 180.0
                    if abs(_daz) > 1.0 or abs(_el_a - _el_e) > 1.0:
                        print('    NOTE: {} and {} were shown about different '
                              'patch centers ((az {:g}, el {:g}) vs '
                              '(az {:g}, el {:g})); the map is drawn relative '
                              'to their midpoint'.format(
                                  gl_a, gl_e, _az_a, _el_a, _az_e, _el_e))
                    _c_az = _az_e + _daz / 2.0      # circular midpoint
                    _c_el = 0.5 * (_el_a + _el_e)
                    # The patch center breaks the two-root ambiguity: the lit
                    # region is a cap about it and the other root is generally
                    # behind the fly.
                    _pref = np.array([
                        np.cos(np.radians(_c_el)) * np.cos(np.radians(_c_az)),
                        np.cos(np.radians(_c_el)) * np.sin(np.radians(_c_az)),
                        np.sin(np.radians(_c_el))])
                    n_common = min(len(res_a['center_deg']),
                                   len(res_e['center_deg']))
                    # POLE-FRAME centroids: combining two families into a point
                    # has to go through the frame they were both binned in.
                    _d = direction_from_two_families(
                        res_a['meta'], res_a['center_deg'][:n_common],
                        res_e['meta'], res_e['center_deg'][:n_common],
                        prefer=_pref)
                    az_c, el_c = screen_azel_of(_d, _c_az, _c_el)
                    # sign lookups, keyed by roi index: kept_signs is
                    # positional within kept_rois, not indexed by roi
                    _sa, _se = {}, {}
                    for _r, _dst in ((res_a, _sa), (res_e, _se)):
                        _kk, _ss = _r.get('keep_rois'), _r.get('kept_signs')
                        if _kk is not None and _ss is not None:
                            _dst.update(zip(np.asarray(_kk, int).tolist(),
                                            np.asarray(_ss, int).tolist()))
                    keep_a = res_a.get('keep_rois')
                    keep_e = res_e.get('keep_rois')
                    if keep_a is None or keep_e is None:
                        sel = np.arange(n_common)
                        how = 'all rois (no null selection available)'
                    else:
                        sel = np.intersect1d(np.asarray(keep_a, dtype=int),
                                             np.asarray(keep_e, dtype=int))
                        sel = sel[sel < n_common]
                        how = 'passed the null in both'
                    ok = sel[np.isfinite(az_c[sel]) & np.isfinite(el_c[sel])]
                    print('  {}: rf center map, {} of {} rois ({})'.format(
                        ch, len(ok), n_common, how))
                    if not len(ok):
                        continue

                    fig, ax = plt.subplots(figsize=(7.2, 5.8),
                                           constrained_layout=True)
                    # COLOURED BY POLARITY, and the minority named per
                    # CHANNEL rather than by a fixed sign: the cAMP indicator's
                    # majority is negative and calcium's is positive, so a rule
                    # written around one sign would mark every cAMP roi as
                    # anomalous. What is worth seeing is each population's own
                    # dissenters, and where they sit -- an roi whose receptive
                    # field CENTER falls outside what the fly can see, but whose
                    # antagonistic surround does not, is driven only through the
                    # surround and reads with the opposite sign. If that is what
                    # these are, they lie at the rim of the cloud.
                    _sgn = {}
                    for _r in ok:
                        _v = [d.get(int(_r)) for d in (_sa, _se)]
                        _v = [x for x in _v if x is not None]
                        _sgn[int(_r)] = (None if len(set(_v)) != 1 else _v[0])
                    _vals = [v for v in _sgn.values() if v is not None]
                    _maj = 1 if sum(1 for v in _vals if v > 0) >= \
                        sum(1 for v in _vals if v < 0) else -1
                    _cat = {r: ('disagree' if v is None else
                                ('majority' if v == _maj else 'minority'))
                            for r, v in _sgn.items()}
                    _style = {'majority': ('#1f77b4', 'o', 42),
                              'minority': ('#d62728', 'o', 52),
                              'disagree': ('0.55', 's', 42)}
                    for _k, (_col, _mk, _sz) in _style.items():
                        _m = np.array([_cat[int(r)] == _k for r in ok], dtype=bool)
                        if not _m.any():
                            continue
                        ax.scatter(az_c[ok][_m], el_c[ok][_m], s=_sz, c=_col,
                                   marker=_mk, edgecolors='k', linewidths=0.5,
                                   zorder=3, label='{} (n={}, sign {})'.format(
                                       _k, int(_m.sum()),
                                       '+' if _k == 'majority' and _maj > 0
                                       else ('-' if _k == 'majority'
                                             else ('-' if _maj > 0 else '+'))
                                       if _k != 'disagree' else '?'))
                    ax.legend(fontsize=12, loc='best')

                    # ECCENTRICITY BY SIGN. The surround-driven explanation above
                    # is a prediction about position, so it is worth testing
                    # rather than asserting: measured from the population's own
                    # median direction, the minority should sit further out.
                    _md = np.median(_d[ok], axis=0)
                    _md = _md / np.linalg.norm(_md)
                    _ecc = np.degrees(np.arccos(np.clip(_d[ok] @ _md, -1.0, 1.0)))
                    _lbl = np.array([_cat[int(r)] for r in ok])
                    _parts = []
                    for _k in ('majority', 'minority', 'disagree'):
                        _e = _ecc[_lbl == _k]
                        if _e.size:
                            _parts.append('{} n={} median ecc {:.1f} deg '
                                          '(90th {:.1f})'.format(
                                              _k, _e.size, float(np.median(_e)),
                                              float(np.percentile(_e, 90))))
                    print('    ' + ';  '.join(_parts))
                    _mj, _mn = _ecc[_lbl == 'majority'], _ecc[_lbl == 'minority']
                    if _mj.size >= 3 and _mn.size >= 3:
                        # rank-based, so a handful of rois and a skewed spread
                        # cannot carry it: the fraction of majority-minority
                        # pairs in which the minority roi is the further out.
                        _auc = float(np.mean(_mn[:, None] > _mj[None, :]))
                        print('        minority is further from the population '
                              'center in {:.0%} of pairs (0.5 = no relation)'
                              .format(_auc))
                    ax.axhline(0, color='0.6', lw=0.8, ls=':')
                    ax.axvline(0, color='0.6', lw=0.8, ls=':')
                    ax.plot(0, 0, 'k+', ms=11, mew=1.6, zorder=4)
                    ax.set_xlabel('Azimuth re. patch center (deg)')
                    ax.set_ylabel('Elevation re. patch center (deg)')
                    ax.set_aspect('equal', adjustable='datalim')
                    ax.set_title('receptive-field centers  {}\n{}\n'
                                 '{} {}\n  x  {} {}'.format(
                                     ch, how, sn_a, gl_a, sn_e, gl_e),
                                 fontsize=14)
                    if save_figs:
                        fig.savefig(os.path.join(
                            figs_dir,
                            'strf_one_d_rf_centers_{}_{}_{}.pdf'.format(
                                ch, gl_a, gl_e)))
                    if not show_figs:
                        plt.close(fig)

                    # ── the same centres, BY SLICE ─────────────────────────
                    # Which imaging plane each roi lives in, against where its
                    # receptive field sits. A retinotopic layer should show the
                    # planes as ordered bands across the map; a cloud with the
                    # planes mixed through it says the stack cuts across the
                    # retinotopy, or that the centres are noise. One panel
                    # with every slice coloured, then one per slice with the
                    # others greyed, since ten colours on one axis are hard to
                    # read apart. Minority-sign rois keep a red edge so the
                    # two figures can be read against each other.
                    _zs = _roi_slices(sn_a, ch, len(az_c))
                    if _zs is not None and len(ok):
                        _zok = _zs[ok]
                        _planes = sorted(set(int(z) for z in _zok))
                        _cmap = plt.get_cmap('viridis', max(len(_planes), 2))
                        _zcol = {z: _cmap(i) for i, z in enumerate(_planes)}
                        _edge = np.array(['#d62728' if _cat[int(r)] == 'minority'
                                          else 'k' for r in ok])
                        _n_pan = 1 + len(_planes)
                        _ncol = min(_n_pan, 4)
                        _nrow = int(np.ceil(_n_pan / float(_ncol)))
                        fig, axes = plt.subplots(
                            _nrow, _ncol, squeeze=False,
                            figsize=(4.6 * _ncol, 4.3 * _nrow + 1.2),
                            constrained_layout=True)
                        _xl = (float(np.nanmin(az_c[ok])) - 5, float(np.nanmax(az_c[ok])) + 5)
                        _yl = (float(np.nanmin(el_c[ok])) - 5, float(np.nanmax(el_c[ok])) + 5)
                        for _pi in range(_nrow * _ncol):
                            ax = axes[_pi // _ncol, _pi % _ncol]
                            if _pi >= _n_pan:
                                ax.set_axis_off()
                                continue
                            ax.axhline(0, color='0.6', lw=0.8, ls=':')
                            ax.axvline(0, color='0.6', lw=0.8, ls=':')
                            ax.plot(0, 0, 'k+', ms=10, mew=1.4, zorder=4)
                            if _pi == 0:
                                for z in _planes:
                                    _m = _zok == z
                                    ax.scatter(az_c[ok][_m], el_c[ok][_m], s=40,
                                               color=_zcol[z], edgecolors=_edge[_m],
                                               linewidths=0.7, zorder=3,
                                               label='slice {} (n={})'.format(
                                                   z, int(_m.sum())))
                                ax.legend(fontsize=9, loc='best', ncol=2)
                                ax.set_title('all slices', fontsize=12)
                            else:
                                z = _planes[_pi - 1]
                                _m = _zok == z
                                ax.scatter(az_c[ok][~_m], el_c[ok][~_m], s=22,
                                           color='0.82', edgecolors='none', zorder=2)
                                ax.scatter(az_c[ok][_m], el_c[ok][_m], s=44,
                                           color=_zcol[z], edgecolors=_edge[_m],
                                           linewidths=0.7, zorder=3)
                                ax.set_title('slice {}  (n={})'.format(
                                    z, int(_m.sum())), fontsize=12)
                            ax.set_xlim(*_xl)
                            ax.set_ylim(*_yl)
                            ax.set_aspect('equal', adjustable='datalim')
                            ax.set_xlabel('Azimuth re. patch center (deg)', fontsize=9)
                            ax.set_ylabel('Elevation re. patch center (deg)', fontsize=9)
                        fig.suptitle('receptive-field centers by slice  {}  ({})\n'
                                     '{} {}  x  {} {}\n'
                                     'colour: imaging slice; red edge: minority sign'
                                     .format(ch, how, sn_a, gl_a, sn_e, gl_e),
                                     fontsize=12)
                        if save_figs:
                            fig.savefig(os.path.join(
                                figs_dir,
                                'strf_one_d_rf_centers_by_slice_{}_{}_{}.pdf'.format(
                                    ch, gl_a, gl_e)))
                        if not show_figs:
                            plt.close(fig)

        # ── the registered mean: the answer to "what does this filter look like" ──
        # COMMON AXES PER CHANNEL. Each geometry's registered grid ends where
        # its own coverage does, and each mean was being scaled to its own
        # peak, so two panels of the same cell type came out at different
        # degrees-per-inch and different z-per-colour. The eye reads that as a
        # difference in the filter. Every mean of one channel is in the same
        # units -- z against each roi's own noise, degrees from its own centre
        # -- so they can and should share limits. Per CHANNEL, not globally: a
        # calcium and a cAMP indicator have no common scale, and forcing one
        # would flatten whichever is smaller.
        _span = float(args.mean_span_deg)
        _mean_ax = {}
        for (_sn, _ch, _gl, _gn), _o in avg_results.items():
            _ks = coverage_bounds(_o['n_roi_per_bin'], grid_deg=_o['grid_deg'])
            if _ks.size < 3:
                continue
            _g, _m, _e = _o['grid_deg'][_ks], _o['mean'][_ks], _o['sem'][_ks]
            # Clamped BEFORE anything is measured, so the z limits and the
            # colour scale describe what is actually drawn.
            _in = np.abs(_g) <= _span
            if int(_in.sum()) < 3:
                continue
            _g, _m, _e = _g[_in], _m[_in], _e[_in]
            _amp = float(np.nanmax(np.abs(_m)))
            _cur = _mean_ax.get(_ch)
            _new = (float(_g[0]), float(_g[-1]),
                    float(np.nanmin(_m - _e)), float(np.nanmax(_m + _e)),
                    _amp if np.isfinite(_amp) and _amp > 0 else 1.0)
            _mean_ax[_ch] = _new if _cur is None else (
                min(_cur[0], _new[0]), max(_cur[1], _new[1]),
                min(_cur[2], _new[2]), max(_cur[3], _new[3]),
                max(_cur[4], _new[4]))

        for (sn, ch, glabel, gname), out in sorted(avg_results.items()):
            grid, mean, sem = out['grid_deg'], out['mean'], out['sem']
            t_lag = out['t_lag']
            ks = coverage_bounds(out['n_roi_per_bin'], grid_deg=grid)
            if ks.size < 3:
                continue
            # ks is contiguous (coverage_bounds guarantees it), which the
            # imshow extent below relies on: rows must be evenly spaced in
            # degrees for a uniform extent to place them where they belong.
            g, m, e = grid[ks], mean[ks], sem[ks]
            _glo, _ghi, _ylo, _yhi, vmax = _mean_ax[ch]
            _pad = 0.05 * (_yhi - _ylo) if _yhi > _ylo else 1.0

            fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.4),
                                     constrained_layout=True)
            ax = axes[0]
            im = ax.imshow(m, aspect='auto', origin='lower', cmap='RdBu_r',
                           vmin=-vmax, vmax=vmax,
                           extent=[t_lag[0] * 1000, t_lag[-1] * 1000,
                                   g[0] - out['grid_step_deg'] / 2,
                                   g[-1] + out['grid_step_deg'] / 2])
            ax.set_ylim(_glo, _ghi)
            ax.axhline(0, color='k', lw=0.8, ls=':')
            ax.set_xlabel('Lag (ms)')
            ax.set_ylabel('deg re. each roi\'s own RF center')
            fig.colorbar(im, ax=ax, label='z')

            # each marginal at the mean's OWN peak, with the population SEM.
            # Restricted to the drawn range: a peak outside it would pick the
            # lag and the row for a feature the figure does not show.
            _vis = (g >= _glo - 1e-9) & (g <= _ghi + 1e-9)
            _mv = np.where(_vis[:, None], np.abs(m), np.nan)
            if not np.any(np.isfinite(_mv)):
                _mv = np.abs(m)
            li = int(np.nanargmax(np.nanmax(_mv, axis=0)))
            pi = int(np.nanargmax(_mv[:, li]))
            ax = axes[1]
            ax.plot(g, m[:, li], 'k', lw=2)
            ax.fill_between(g, m[:, li] - e[:, li], m[:, li] + e[:, li],
                            color='0.6', alpha=0.4, lw=0)
            ax.axhline(0, color='0.6', lw=0.8)
            ax.axvline(0, color='0.6', lw=0.8, ls=':')
            ax.set_xlim(_glo, _ghi)
            ax.set_ylim(_ylo - _pad, _yhi + _pad)
            ax.set_xlabel('deg re. RF center')
            ax.set_ylabel('mean z')
            ax.set_title('spatial, at {:.0f} ms'.format(t_lag[li] * 1000),
                         fontsize=14)
            ax = axes[2]
            ax.plot(t_lag * 1000, m[pi], 'k', lw=2)
            ax.fill_between(t_lag * 1000, m[pi] - e[pi], m[pi] + e[pi],
                            color='0.6', alpha=0.4, lw=0)
            ax.axhline(0, color='0.6', lw=0.8)
            ax.set_ylim(_ylo - _pad, _yhi + _pad)
            ax.set_xlabel('Lag (ms)')
            ax.set_ylabel('mean z')
            ax.set_title('temporal, at {:+.1f} deg'.format(g[pi]), fontsize=14)

            fig.suptitle('registered mean  {} {} {} [{}]   n = {} rois, '
                         '{:g} deg bins   (axes shared across {}, '
                         '+-{:g} deg)'.format(
                             sn, ch, glabel, gname, out['n_roi'],
                             out['grid_step_deg'], ch, _span),
                         fontsize=15)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_MEAN_{}_{}_{}_{}.pdf'.format(
                        sn, ch, glabel, gname)))
            if not show_figs:
                plt.close(fig)

        # ── the radial profile, and what it costs to assume it ──────────────
        for ch, rad in sorted(radial_results.items()):
            r, h, t_lag = rad['r_centers_deg'], rad['h'], rad['t_lag']
            resid = np.atleast_1d(rad['residual'])
            _m = float(np.nanmax(np.abs(h)))
            vmax = _m if np.isfinite(_m) and _m > 0 else 1.0
            li = int(np.nanargmax(np.nanmax(np.abs(h), axis=0)))

            fig, axes = plt.subplots(1, 3, figsize=(17.5, 4.4),
                                     constrained_layout=True)
            ax = axes[0]
            im = ax.imshow(h, aspect='auto', origin='lower', cmap='RdBu_r',
                           vmin=-vmax, vmax=vmax,
                           extent=[t_lag[0] * 1000, t_lag[-1] * 1000,
                                   r[0], r[-1]])
            ax.set_xlabel('Lag (ms)')
            ax.set_ylabel('radius from RF center (deg)')
            fig.colorbar(im, ax=ax, label='h(r)')
            ax.set_title('radial profile', fontsize=14)

            ax = axes[1]
            axpix, img = radial_to_image(r, h[:, li])
            im = ax.imshow(img, origin='lower', cmap='RdBu_r',
                           vmin=-vmax, vmax=vmax,
                           extent=[axpix[0], axpix[-1], axpix[0], axpix[-1]])
            ax.set_xlabel('deg')
            ax.set_ylabel('deg')
            ax.set_aspect('equal')
            fig.colorbar(im, ax=ax, label='h(r)')
            # a rendering of h(r), not a measurement of 2D structure
            ax.set_title('peak lag {:.0f} ms -- one frame only\n'
                         '(rotation of h, not 2D data)'
                         .format(t_lag[li] * 1000), fontsize=13)

            ax = axes[2]
            ax.plot(t_lag * 1000, resid, 'k', lw=2, label='all families')
            per = np.atleast_2d(rad['residual_per_family'])
            for i, fam in enumerate(rad['families']):
                ax.plot(t_lag * 1000, per[i], lw=1.2, alpha=0.8, label=fam[:22])
            ax.axhline(0.15, color='r', ls=':', lw=1.2)
            ax.set_xlabel('Lag (ms)')
            ax.set_ylabel('relative misfit')
            ax.set_title('does one h explain BOTH families?\n'
                         'red line: elongation becomes visible', fontsize=13)
            ax.legend(fontsize=11)

            fig.suptitle('radial reconstruction  {}   median residual {:.3f}'
                         '   family gains {}'
                         .format(ch, float(np.nanmedian(resid)),
                                 np.array2string(
                                     np.asarray(rad['gain'], dtype=float),
                                     precision=2)), fontsize=15)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_RADIAL_{}.pdf'.format(ch)))
            if not show_figs:
                plt.close(fig)

            # ── the same reconstruction across lags ─────────────────────────
            # One frame at the peak shows whichever lobe happens to be largest
            # there. Center and surround peak at different times in this
            # preparation, so the lag series is where the opponent structure
            # actually appears -- and a shared colour scale is what makes the
            # panels comparable rather than each one renormalised to itself.
            amp = np.nanmax(np.abs(h), axis=0)
            strong = np.flatnonzero(amp >= 0.15 * np.nanmax(amp))
            lo, hi = (int(strong[0]), int(strong[-1])) if strong.size else (0, len(t_lag) - 1)
            sel = np.unique(np.linspace(lo, hi, max(int(args.n_display_lags), 2)
                                        ).astype(int))
            n_col = int(np.ceil(np.sqrt(len(sel))))
            n_row = int(np.ceil(len(sel) / n_col))
            fig, axes = plt.subplots(n_row, n_col, squeeze=False,
                                     figsize=(3.0 * n_col, 3.2 * n_row),
                                     constrained_layout=True)
            for a in axes.ravel():
                a.set_axis_off()
            im = None
            for k, li2 in enumerate(sel):
                a = axes[k // n_col][k % n_col]
                a.set_axis_on()
                axpix, img2 = radial_to_image(r, h[:, li2])
                im = a.imshow(img2, origin='lower', cmap='RdBu_r',
                              vmin=-vmax, vmax=vmax,
                              extent=[axpix[0], axpix[-1], axpix[0], axpix[-1]])
                a.set_aspect('equal')
                a.set_title('{:.0f} ms'.format(t_lag[li2] * 1000), fontsize=13)
                a.set_xticks([])
                a.set_yticks([])
            if im is not None:
                fig.colorbar(im, ax=axes.ravel().tolist(), label='h(r)',
                             shrink=0.7)
            fig.suptitle('radial reconstruction across lags  {}   '
                         'median residual {:.3f}'.format(
                             ch, float(np.nanmedian(resid))), fontsize=15)
            if save_figs:
                fig.savefig(os.path.join(
                    figs_dir, 'strf_one_d_RADIAL_LAGS_{}.pdf'.format(ch)))
            if not show_figs:
                plt.close(fig)

        print('Figures written to ' + figs_dir)
        _n_pop = len([_f for _f in os.listdir(figs_dir)
                      if os.path.isfile(os.path.join(figs_dir, _f))])
        print('  {} population figures, {} per-roi figures in rois/'.format(
            _n_pop, len(os.listdir(rois_dir))))

    if log_path:
        print('\nappended to analysis log: ' + log_path)

    if show_figs:
        plt.show()
