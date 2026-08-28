"""
octave_stim.py

Stimulus-side helpers for the multi-octave ternary noise stimulus, the
analogue of noise_stim.py for the bar-noise ("noizone") protocols.

WHAT IS DIFFERENT FROM noizone
------------------------------
* There is no movie file. The stimulus is a pure function of the frame index,
  so the epoch's two metadata dicts ARE the stimulus. Reconstruction is a
  hash evaluation, not a file read, and it is exact for any set of indices --
  including a set with gaps where frames dropped.

* The regressor is a set of CELLS on the sphere, not bars on a line. Cells are
  mutually independent within and across octaves, so in this basis the spatial
  covariance is exactly (2p/clip_k^2) * I -- a scalar multiple of the identity,
  not merely diagonal. A plain STA is therefore already the ML estimate up to a
  known constant: no whitening, no spatial ridge, no C_ss to invert.

* The stimulus is TEMPORALLY correlated, per octave, with an exponential
  autocorrelation exp(-|lag|/tau_j). The STA returns the filter convolved with
  that, and unlike the noizone case the smearing is NOT small next to the
  indicator kernel -- see hold_transfer() and the deconvolution step in
  analyze_data_strf_octave.py. This is the one place regularisation is needed.

* The cell basis is overcomplete across octaves (one 25 deg cell overlaps ~17
  of the 6 deg cells), so the recovered weights are not a picture. Each octave
  independently measures the cell-averaged filter at its own scale:
  w_jc = (integral of k over cell jc) / clip_k. They are combined by fitting
  through both footprints, never by summing the layers.
"""

import ast
import json
import pathlib
import re
import sys

import numpy as np

# The generator lives in labpack. Imported lazily so this module can be
# imported (for is_octave_ternary etc.) on a machine without labpack on the
# path -- only the reconstruction calls actually need it.
_gen = None
_tess = None


def _find_labpack():
    """Put a sibling clandinin_labpack checkout on the path if labpack is absent.

    The analysis environment has visanalysis but not labpack, and unlike the
    noizone pipeline this one genuinely needs the generator: there is no movie
    file to read, so reconstruction IS a call into labpack. Rather than make the
    analysis env depend on the stimulus package, look next to this repo, which is
    where both checkouts live.

    Silent when labpack is already importable, and silent when there is no
    sibling either -- _generator() then raises with something actionable.
    """
    import importlib.util
    if importlib.util.find_spec('labpack') is not None:
        return
    # Anchored on the module actually imported below, not on labpack/__init__.py:
    # labpack is a PEP 420 namespace package and has no __init__.py, so looking
    # for one finds nothing and the fallback silently never fires.
    here = pathlib.Path(__file__).resolve()
    needed = pathlib.Path('labpack', 'visual_stim', 'clandinin', 'octave_ternary.py')
    for parent in here.parents:
        candidate = parent / 'clandinin_labpack'
        if (candidate / needed).exists():
            sys.path.insert(0, str(candidate))
            return


def _generator():
    global _gen, _tess
    if _gen is None:
        _find_labpack()
        try:
            from labpack.visual_stim.clandinin import octave_ternary, sphere_tessellation
        except ImportError as e:
            raise ImportError(
                'the octave ternary stimulus is regenerated from its metadata rather '
                'than read from a movie file, so reconstructing it needs labpack, and '
                'labpack is not importable here. Either pip install -e the '
                'clandinin_labpack checkout into this environment, or place it beside '
                'the visanalysis checkout so it can be found automatically.') from e
        _gen, _tess = octave_ternary, sphere_tessellation
    return _gen, _tess


# ─────────────────────────────────────────────────────────────────────────────
# series / epoch identification
# ─────────────────────────────────────────────────────────────────────────────

def is_octave_ternary(run_parameters):
    """True if a series ran the octave ternary protocol.

    Keyed on protocol_ID rather than on a filename, because unlike the movie
    protocols there is no file to name.
    """
    pid = str(run_parameters.get('protocol_ID', ''))
    return 'OctaveTernary' in pid


#: numpy scalar reprs left inside a stringified dict, e.g. np.float64(0.5)
_NP_SCALAR = re.compile(r'\bnp\.\w+\(([^()]*)\)')


def _parse_meta_attr(v, name):
    """Read one metadata dict back out of an hdf5 attribute.

    stimpack stores these by stringifying the dict, so what comes back is a
    Python repr and not JSON -- single-quoted keys, True rather than true. That
    rules out json.loads and leaves ast.literal_eval.

    literal_eval alone is not quite enough either. Any numpy scalar that reached
    as_metadata reprs as "np.float64(0.707...)", which is a call expression and
    therefore a SyntaxError to literal_eval. The generator no longer emits those
    -- as_metadata now casts everything to plain Python types -- but files
    recorded before that fix carry them in `weights`, and those recordings are
    not reproducible by any other route. Unwrapping the call is exact, since the
    argument is the literal value, so accepting them costs nothing.
    """
    if isinstance(v, dict):
        return v
    if isinstance(v, bytes):
        v = v.decode('utf-8')
    if not isinstance(v, str):
        raise TypeError('{} is a {}, expected a dict or a stringified one'
                        .format(name, type(v).__name__))
    for attempt in (lambda s: json.loads(s),
                    lambda s: ast.literal_eval(s),
                    lambda s: ast.literal_eval(_NP_SCALAR.sub(r'\1', s))):
        try:
            return attempt(v)
        except (ValueError, SyntaxError):
            continue
    raise ValueError(
        '{} could not be parsed. It is neither JSON nor a Python literal, even '
        'after unwrapping numpy scalar reprs. First 200 characters:\n{}'
        .format(name, v[:200]))


def resolve_octave_epoch_meta(epoch_params):
    """Return (meta, sampling) for one epoch.

    Both dicts are written by the protocol under the names it actually uses --
    `octave_stimulus` and `octave_sampling`, matching
    JCS_protocol.get_trial_parameters. Nothing is defaulted: a default that
    changed later would silently rewrite the past, so a missing field is an
    error rather than something to paper over.
    """
    meta = epoch_params.get('octave_stimulus')
    sampling = epoch_params.get('octave_sampling')
    if meta is None or sampling is None:
        present = sorted(str(k) for k in epoch_params if 'octave' in str(k).lower())
        raise KeyError(
            'epoch is missing octave_stimulus / octave_sampling. Without them the '
            'stimulus cannot be reconstructed -- check that process_data.py '
            'attached the epoch parameters for this series. Octave-ish keys that '
            'ARE present: {}'.format(present or 'none'))
    meta = _parse_meta_attr(meta, 'octave_stimulus')
    sampling = _parse_meta_attr(sampling, 'octave_sampling')
    if meta.get('kind') != 'octave_ternary_sphere':
        raise ValueError('unexpected stimulus kind {!r}'.format(meta.get('kind')))
    return meta, sampling


# ─────────────────────────────────────────────────────────────────────────────
# reconstruction
# ─────────────────────────────────────────────────────────────────────────────

def load_octave_cells(meta, n_frames, frame_indices=None):
    """Exact cell values at raw display-frame resolution.

    Returns (n_frames, n_cells_total) float32, mean-subtracted.

    Rows are RAW display frames, one per frame the projector drew, matching the
    convention noise_stim.load_bar_stim_1d uses. That is what lets the
    dropped-frame correction downstream operate per raw frame.

    frame_indices lets a caller ask for an arbitrary set instead -- scattered,
    out of order, or with gaps -- which is exact here rather than approximate,
    because every value is a pure function of its index.
    """
    gen, _ = _generator()
    idx = np.arange(int(n_frames)) if frame_indices is None else np.asarray(frame_indices)
    cells = gen.OctaveTernaryStimulus.reproduce_cells(meta, idx)   # (n_cells, T) int8
    out = np.ascontiguousarray(cells.T, dtype=np.float32)          # (T, n_cells)
    # Zero-mean by construction, but subtract the realised mean so the STA is
    # not biased by a finite-sample offset in a short trial.
    out -= out.mean(axis=0, keepdims=True)
    return out


def octave_slices(meta):
    """[(name, slice, patch_deg, tau), ...] -- which columns belong to which octave.

    Everything downstream that treats the two layers differently (deconvolution,
    effective dof, cell footprint) indexes through this rather than hard-coding
    a split.
    """
    n_cells = list(meta['n_cells'])
    patch = list(meta['patch_degrees'])
    taus = list(meta['hold_taus'])
    off = np.cumsum([0] + n_cells)
    return [(f'{d:g}deg', slice(int(off[j]), int(off[j + 1])), float(d), float(t))
            for j, (d, t) in enumerate(zip(patch, taus))]


def cell_directions(meta):
    """(n_cells_total, 3) unit vectors -- where each cell sits on the sphere.

    Needed for the centroid, which must be a weighted mean of unit vectors
    renormalised to the sphere. Averaging az/el as if they were Cartesian goes
    wrong near the edges of the lit band, which is exactly where it matters.
    """
    _, tess = _generator()
    fib_seed = int(meta['fibonacci_seed'])
    return np.vstack([tess.fibonacci_directions(int(n), fib_seed)
                      for n in meta['n_cells']])


def cell_azel_deg(meta):
    """(n_cells_total, 2) azimuth, elevation in degrees. For plotting only."""
    d = cell_directions(meta)
    az = np.degrees(np.arctan2(d[:, 1], d[:, 0]))
    el = np.degrees(np.arcsin(np.clip(d[:, 2], -1, 1)))
    return np.column_stack([az, el])


def cell_variance(meta):
    """Per-cell variance of the regressor, 2p.

    UNITS. reproduce_cells returns the RAW ternary values in {-1, 0, +1}, not
    the display-scaled ones -- the 1/clip_k that converts a cell value to
    display contrast is a constant that folds into the recovered filter's
    units, so it is left out of the regressor. The variance here is therefore
    2p = 1/3, not 2p/clip_k^2.

    What matters for the analysis is that it is identical for every cell in
    every octave, so C_ss = 2p * I -- a scalar multiple of the identity, not
    merely diagonal. That is why a plain STA is already the ML estimate up to a
    known constant, with no whitening and no spatial ridge.

    To express a recovered filter in units of display contrast rather than cell
    value, multiply by clip_k.
    """
    return 2.0 * float(meta['p'])


def display_scale(meta):
    """clip_k -- multiply a cell-basis filter by this to get contrast units."""
    return float(meta['clip_k'])


# ─────────────────────────────────────────────────────────────────────────────
# temporal: the one place a correction is genuinely needed
# ─────────────────────────────────────────────────────────────────────────────

def hold_autocorr(tau, lag_dt, n_lags):
    """exp(-|lag|/tau) sampled on the lag axis, normalised to 1 at zero lag.

    The STA returns the filter convolved with this. For the noizone stimulus
    the equivalent was a triangle of half-width noisetau and was small enough
    to ignore; here it is not, because tau is comparable to the kernels being
    measured.
    """
    lags = np.arange(n_lags) * float(lag_dt)
    return np.exp(-lags / float(tau))


def hold_transfer(tau, lag_dt, n_lags):
    """Frequency-domain version: the Lorentzian the deconvolution divides out.

    Returns (freqs_hz, transfer) with transfer normalised to 1 at DC. The
    amplification a naive deconvolution would apply is 1/transfer, which grows
    as f^2 -- 5x at 1 Hz and 28x at 2.5 Hz for tau = 0.33 s. That is why the
    deconvolution is regularised rather than done by division.
    """
    f = np.fft.rfftfreq(int(n_lags), d=float(lag_dt))
    return f, 1.0 / (1.0 + (2.0 * np.pi * f * float(tau)) ** 2)


def effective_dof_frames(tau, lag_dt):
    """Frames per independent sample for one octave.

    n_lags overstates the degrees of freedom by roughly this factor, and it
    differs between octaves -- 120 frames for a 1.0 s hold against 40 for a
    0.33 s one. A single z threshold applied to both treats three times as much
    independent data as if it were the same amount.
    """
    return float(tau) / float(lag_dt)


def visible_cells(meta, sampling):
    """(n_cells,) bool -- which cells fall inside the displayed patch.

    The tessellation covers the whole sphere but the patch covers about a
    quarter of it, so roughly three quarters of the cells are generated and
    never shown. Their STA is pure noise, and including them does two kinds of
    damage: it triples the multiple-comparisons burden behind any peak
    statistic, and it lets a centroid be dragged to a cell on the back of the
    sphere. Measured on the first real recording, 325 of 1212 cells were in
    view, and the strongest cell for several ROIs sat at azimuths beyond 160
    degrees -- entirely outside the lit region.

    Membership is tested in patch coordinates rather than by az/el bounds, so
    it stays correct when the patch is not centred on the origin.
    """
    d = cell_directions(meta)
    az0 = np.radians(float(sampling['center_az_deg']))
    el0 = np.radians(float(sampling['center_el_deg']))
    # rotate the patch centre onto +x, the same convention the display uses
    ca, sa = np.cos(-az0), np.sin(-az0)
    Rz = np.array([[ca, -sa, 0.0], [sa, ca, 0.0], [0.0, 0.0, 1.0]])
    ce, se = np.cos(-el0), np.sin(-el0)
    Ry = np.array([[ce, 0.0, se], [0.0, 1.0, 0.0], [-se, 0.0, ce]])
    p = d @ Rz.T @ Ry.T
    az = np.degrees(np.arctan2(p[:, 1], p[:, 0]))
    el = np.degrees(np.arcsin(np.clip(p[:, 2], -1, 1)))
    return ((np.abs(az) <= float(sampling['width_deg']) / 2) &
            (np.abs(el) <= float(sampling['height_deg']) / 2) &
            (p[:, 0] > 0))
