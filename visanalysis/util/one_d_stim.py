"""
one_d_stim.py

Stimulus-side helpers for the 1D ternary noise stimulus (labpack
SphericalOneDTernaryNoise -- bars or wedges about a pole), the 1D analogue of
octave_stim.py for the bar-noise ("noizone") protocols.

WHAT IS DIFFERENT FROM noizone (flymax NoizOne bar movies)
------------------------------------------------------------
* There is no movie file. The stimulus is a pure function of the frame index,
  so the epoch's two metadata dicts ARE the stimulus. Reconstruction is a hash
  evaluation, not a file read, and it is exact for any set of indices --
  including a set with gaps where frames dropped.

* The regressor is a set of CELLS (bars or wedges), not positions read off a
  memmap'd movie. Cells are mutually independent, so the spatial covariance in
  the cell basis is exactly 2p * I -- a scalar multiple of the identity, not
  merely diagonal. A plain STA is therefore already the ML estimate up to a
  known constant: no whitening, no spatial ridge, no C_ss to invert.

* The hold is EXPONENTIAL (mean hold_tau), not a fixed sample-and-hold period.
  The STA returns the filter convolved with exp(-|lag|/hold_tau), and unlike
  the noizone case that smearing is not small next to the indicator kernel --
  see one_d_strf.py and octave_strf.deconvolve_octaves, which this module
  feeds.

WHAT IS DIFFERENT FROM octave_stim.py
--------------------------------------
* ONE layer, one hold_tau. There is no octave_slices, no coarse/fine merge and
  no fit_joint_cells here: bars/wedges are disjoint rings or sectors, not
  nested scales, so the cell basis is NOT overcomplete and the plain STA is
  already the final per-cell estimate.

* A cell is a RING or a SECTOR, not a point on the sphere -- it covers a whole
  range of directions, so there is no (az, el) to give it and no
  cell_directions()/spherical centroid. Position is one coordinate along the
  stimulus's own axis: beta (angle from the pole, bars) or alpha (azimuth
  about the pole, wedges, wraps at 360). See `cell_positions_deg`.

* The protocol records exactly which cells the projector could ever draw a
  texel from -- `one_d_visible_cells`, computed at trial time from the real
  screen geometry -- so `epoch_visible_cells` reads that back rather than
  recomputing it here. `recompute_visible_cells` is kept for old recordings
  that predate the field, or as a cross-check.
"""

import ast
import json
import pathlib
import re
import sys

import numpy as np

# The generator lives in labpack. Imported lazily so this module can be
# imported (for is_one_d_ternary etc.) on a machine without labpack on the
# path -- only the reconstruction calls actually need it.
_gen = None
_tess = None


def _find_labpack():
    """Put a sibling clandinin_labpack checkout on the path if labpack is absent.

    Same trick as octave_stim._find_labpack, and for the same reason: there is
    no movie file, so reconstruction IS a call into labpack.
    """
    import importlib.util
    if importlib.util.find_spec('labpack') is not None:
        return
    # Anchored on the module actually imported below, not on labpack/__init__.py:
    # labpack is a PEP 420 namespace package and has no __init__.py.
    here = pathlib.Path(__file__).resolve()
    needed = pathlib.Path('labpack', 'visual_stim', 'clandinin', 'one_d_ternary.py')
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
            from labpack.visual_stim.clandinin import one_d_ternary, one_d_tessellation
        except ImportError as e:
            raise ImportError(
                'the 1D ternary stimulus is regenerated from its metadata rather '
                'than read from a movie file, so reconstructing it needs labpack, and '
                'labpack is not importable here. Either pip install -e the '
                'clandinin_labpack checkout into this environment, or place it beside '
                'the visanalysis checkout so it can be found automatically.') from e
        _gen, _tess = one_d_ternary, one_d_tessellation
    return _gen, _tess


def pole_frame(meta):
    """(pole, u, v, reference) of this family, as unit vectors IN THE FLY FRAME.

    The stimulus records its pole in a frame whose forward direction is
    `frame_az_deg`/`frame_el_deg` -- the rig's screen_center, since
    one_d_ternary/3 -- so a recorded pole of (0, 90) means "up as the screen
    sees it", which on a mount pitched by 38 degrees is fly (0, 52). The wedge
    seam reference is fixed in that frame too. A record without the keys (/2
    and earlier) is the identity frame, exactly as from_metadata treats it.

    Everything that turns a binned coordinate back into a direction has to go
    through this. Reading the pole as fly-frame on a tilted mount puts every
    localisation in the screen's coordinates while the caller believes they
    are the fly's, and the root-selection preference -- which IS fly-frame --
    then picks the antipode for anything high on the screen.
    """
    _gen, tess = _generator()
    R = tess.frame_rotation(float(meta.get('frame_az_deg', 0.0)),
                            float(meta.get('frame_el_deg', 0.0)))
    pole_f = tess.pole_vector(float(meta['pole_az_deg']), float(meta['pole_el_deg']))
    u_f, v_f, ref = tess.wedge_reference_frame(pole_f)
    return R @ pole_f, R @ u_f, R @ v_f, ref


# ─────────────────────────────────────────────────────────────────────────────
# series / epoch identification
# ─────────────────────────────────────────────────────────────────────────────

def is_one_d_ternary(run_parameters):
    """True if a series ran the 1D ternary (bar/wedge noise) protocol.

    Keyed on protocol_ID -- stimpack sets it to the protocol class's own name,
    SphericalOneDTernaryNoise -- rather than on a filename, because unlike the
    movie protocols there is no file to name.
    """
    pid = str(run_parameters.get('protocol_ID', ''))
    return 'OneDTernaryNoise' in pid


#: numpy scalar reprs left inside a stringified dict, e.g. np.float64(0.5)
_NP_SCALAR = re.compile(r'\bnp\.\w+\(([^()]*)\)')


def _parse_meta_attr(v, name):
    """Read one metadata dict back out of an hdf5 attribute.

    Same parsing problem octave_stim._parse_meta_attr solves and the same
    fix: stimpack stores these by stringifying the dict, so what comes back is
    a Python repr (literal_eval), not JSON, and a stray numpy scalar repr
    (a call expression) has to be unwrapped before literal_eval can read it.
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


def resolve_one_d_epoch_meta(epoch_params):
    """Return (meta, sampling) for one epoch.

    Both dicts are written by the protocol under the names it actually uses --
    `one_d_stimulus` and `one_d_sampling`, matching
    JCS_protocol.SphericalOneDTernaryNoise.get_trial_parameters. Nothing is
    defaulted: a default that changed later would silently rewrite the past,
    so a missing field is an error rather than something to paper over.
    """
    meta = epoch_params.get('one_d_stimulus')
    sampling = epoch_params.get('one_d_sampling')
    if meta is None or sampling is None:
        present = sorted(str(k) for k in epoch_params if 'one_d' in str(k).lower())
        raise KeyError(
            'epoch is missing one_d_stimulus / one_d_sampling. Without them the '
            'stimulus cannot be reconstructed -- check that process_data.py '
            'attached the epoch parameters for this series. one_d-ish keys that '
            'ARE present: {}'.format(present or 'none'))
    meta = _parse_meta_attr(meta, 'one_d_stimulus')
    sampling = _parse_meta_attr(sampling, 'one_d_sampling')
    if meta.get('kind') != 'one_d_ternary_sphere':
        raise ValueError('unexpected stimulus kind {!r}'.format(meta.get('kind')))
    return meta, sampling


def epoch_visible_cells(epoch_params, meta=None, sampling=None):
    """(n_visible,) int, ascending -- cell ids the texture could ever show.

    Prefers the value the protocol recorded at trial time
    (`one_d_visible_cells`), which is exact because it was computed from the
    real screen geometry rather than re-derived here. Falls back to
    `recompute_visible_cells` (needs `meta` and `sampling`) only for a
    recording made before that field existed.

    This is the TEXTURE-visible set, not the lit region: the bowl's circular
    aperture cuts some of these cells further, and no fixed set names that --
    compare `one_d_stim.cell_solid_angles(meta)` against the area actually
    displayed for that.
    """
    raw = epoch_params.get('one_d_visible_cells')
    if raw is not None:
        return np.asarray(sorted(int(c) for c in raw), dtype=np.int64)
    if meta is None or sampling is None:
        raise KeyError(
            "epoch has no 'one_d_visible_cells' (an old recording?) and no "
            "meta/sampling was given to recompute it from.")
    return recompute_visible_cells(meta, sampling)


def recompute_visible_cells(meta, sampling):
    """Recompute the texture-visible cell set from geometry alone.

    Should agree with the protocol's own `one_d_visible_cells` bit-for-bit;
    kept as a fallback/cross-check, not the primary path -- see
    `epoch_visible_cells`.
    """
    gen, _ = _generator()
    stim = gen.OneDTernaryStimulus.from_metadata(meta)
    return np.asarray(stim.visible_cells(
        int(sampling['n_rows']), int(sampling['n_cols']),
        float(sampling['width_deg']), float(sampling['height_deg']),
        float(sampling.get('center_az_deg', 0.0)),
        float(sampling.get('center_el_deg', 0.0))), dtype=np.int64)


# ─────────────────────────────────────────────────────────────────────────────
# reconstruction
# ─────────────────────────────────────────────────────────────────────────────

def load_one_d_cells(meta, n_frames, frame_indices=None):
    """Exact cell values at raw display-frame resolution.

    Returns (n_frames, n_cells) float32, mean-subtracted. Rows are RAW display
    frames, matching the convention octave_stim.load_octave_cells and
    noise_stim.load_bar_stim_1d both use, which is what lets the caller ask
    for whichever indices actually rendered (dropped frames included) rather
    than repairing an array afterwards.
    """
    gen, _ = _generator()
    idx = np.arange(int(n_frames)) if frame_indices is None else np.asarray(frame_indices)
    cells = gen.OneDTernaryStimulus.reproduce_cells(meta, idx)   # (n_cells, T) int8
    out = np.ascontiguousarray(cells.T, dtype=np.float32)        # (T, n_cells)
    # Zero-mean by construction, but subtract the realised mean so the STA is
    # not biased by a finite-sample offset in a short trial.
    out -= out.mean(axis=0, keepdims=True)
    return out


# ─────────────────────────────────────────────────────────────────────────────
# geometry: one coordinate per cell, not a direction
# ─────────────────────────────────────────────────────────────────────────────

def cell_positions_deg(meta):
    """(n_cells,) the coordinate at each cell's CENTRE, along the stimulus's
    own binned axis.

    beta (0 at the pole, 180 at the antipole) for bars; alpha (azimuth about
    the pole, wraps at 360, 0 at the recorded `wedge_reference`) for wedges.
    This is deliberately NOT an (az, el) pair: a bar or wedge is a whole ring
    or sector of directions, not a point, so unlike octave_stim.cell_azel_deg
    there is no single sphere location to give it. Use `one_d_strf.centroid_1d`
    to average positions -- wedges wrap and a plain mean of angles is wrong at
    the seam.
    """
    _, tess = _generator()
    b = tess.cell_bounds(meta['mode'], meta['cell_width_deg'], meta['phase_deg'])
    return 0.5 * (b[:, 0] + b[:, 1])


def cell_bounds_deg(meta):
    """(n_cells, 2) [start, end) of each cell along the binned axis. See
    one_d_tessellation.cell_bounds for the exact convention at each end."""
    _, tess = _generator()
    return tess.cell_bounds(meta['mode'], meta['cell_width_deg'], meta['phase_deg'])


def cell_widths_deg(meta):
    """Each cell's true extent along the binned axis (narrower at the ends --
    see `runt_cell_ids`)."""
    return np.asarray(meta['cell_widths_deg'], dtype=float)


def runt_cell_ids(meta):
    """Ids of the leftover, narrower-than-nominal cells (recorded at write
    time; see OneDTernaryStimulus.runt_cells). Being a runt is about binning,
    not about the aperture -- cross-reference with `epoch_visible_cells`."""
    return np.asarray(meta['runt_cells'], dtype=np.int64)


def cell_solid_angles(meta):
    """Full solid angle of each cell over the WHOLE sphere, steradians.

    The area a cell WOULD have if the whole sphere were lit. Bars are far from
    equal-area -- the cell at the pole is a small cap, the one at the equator
    a wide band -- which is inherent to "iso-elevation" and not a defect, but
    it means raw per-cell coefficients are not comparable between cells without
    accounting for it. Compare against the area actually displayed (the bowl
    lights a circular cap) to find cells the aperture cuts hardest.
    """
    _, tess = _generator()
    return tess.cell_solid_angles(meta['mode'], meta['cell_width_deg'], meta['phase_deg'])


def cell_variance(meta):
    """Per-cell variance of the regressor, 2p.

    UNITS. reproduce_cells returns the RAW ternary values in {-1, 0, +1}, the
    same values (not display bytes) the octave stimulus's cells return, so
    there is no separate clip_k here the way octave_stim.display_scale has --
    a recovered filter is already in units of full stimulus contrast.

    C_ss = 2p * I exactly, the same reasoning as octave_stim.cell_variance:
    every cell has the same variance and cells are mutually independent, so
    the plain STA is the ML estimate up to this one known constant.
    """
    return 2.0 * float(meta['p'])


def family_constraint(meta, coord_deg):
    """(a, c) with `d . a = c` for a direction measured into one 1D family.

    Both parcellations are the same kind of statement about a direction, which
    is what lets any two of them be combined:

      bars   bin the angle FROM the pole, so a cell is the cone d . p = cos(beta)
      wedges bin the azimuth ABOUT it, so a cell is the half plane through the
             pole at azimuth alpha, whose great circle has normal
             n = v cos(alpha) - u sin(alpha), giving d . n = 0

    `coord_deg` may be an array; `a` is then (..., 3) and `c` (...,).
    """
    pole, u, v, _ref = pole_frame(meta)
    x = np.radians(np.asarray(coord_deg, dtype=float))
    if str(meta['mode']) == 'bars':
        a = np.broadcast_to(pole, x.shape + (3,))
        return np.array(a, dtype=float), np.cos(x)
    a = (v[None, :] * np.cos(x)[..., None] - u[None, :] * np.sin(x)[..., None])
    return a, np.zeros(x.shape)


def direction_from_two_families(meta_a, coord_a, meta_b, coord_b, prefer=None):
    """Direction(s) satisfying two 1D measurements at once.

    Each family pins the direction to a curve -- a ring for bars, a great
    circle for wedges -- and two curves meet in (generally) two antipodal
    points. Writing both as `d . a = c` makes the solve one small linear system
    whatever the two families are, so bars+wedges about one pole, wedges+wedges
    about different poles, and bars+bars all go through the same code:

        d = x a1 + y a2 + z (a1 x a2),  x, y from the two dot products,
        z from |d| = 1, which is where the two solutions come from.

    `prefer` picks between them -- pass the patch center, since the lit region
    is a cap about it and the other root is generally behind the fly. Returns
    NaN where the two constraints cannot be met (inconsistent measurements) or
    where the families are parallel and pin nothing down.
    """
    a1, c1 = family_constraint(meta_a, coord_a)
    a2, c2 = family_constraint(meta_b, coord_b)
    a1 = np.atleast_2d(a1)
    a2 = np.atleast_2d(a2)
    c1 = np.atleast_1d(c1)
    c2 = np.atleast_1d(c2)
    n = max(a1.shape[0], a2.shape[0])
    a1 = np.broadcast_to(a1, (n, 3))
    a2 = np.broadcast_to(a2, (n, 3))
    c1 = np.broadcast_to(c1, (n,))
    c2 = np.broadcast_to(c2, (n,))

    g11 = np.einsum('ij,ij->i', a1, a1)
    g22 = np.einsum('ij,ij->i', a2, a2)
    g12 = np.einsum('ij,ij->i', a1, a2)
    det = g11 * g22 - g12 ** 2
    out = np.full((n, 3), np.nan)
    ok = np.abs(det) > 1e-12                      # parallel families pin nothing
    if not np.any(ok):
        return out
    x = np.where(ok, (c1 * g22 - c2 * g12) / np.where(ok, det, 1.0), np.nan)
    y = np.where(ok, (c2 * g11 - c1 * g12) / np.where(ok, det, 1.0), np.nan)
    p = x[:, None] * a1 + y[:, None] * a2
    w = np.cross(a1, a2)
    w2 = np.einsum('ij,ij->i', w, w)
    rem = 1.0 - np.einsum('ij,ij->i', p, p)
    good = ok & (rem >= -1e-9) & (w2 > 1e-12)     # negative => no intersection
    z = np.sqrt(np.clip(rem, 0.0, None) / np.where(w2 > 0, w2, 1.0))
    d1 = p + z[:, None] * w
    d2 = p - z[:, None] * w
    if prefer is None:
        pick = d1
    else:
        pref = np.asarray(prefer, dtype=float)
        pick = np.where((d1 @ pref)[:, None] >= (d2 @ pref)[:, None], d1, d2)
        # WHEN THE CHOICE IS SAFE, AND WHEN IT IS NOT. Two great circles (a
        # wedges/wedges pair) meet at ANTIPODAL points, so one root is always
        # behind the fly and `prefer` decides correctly. Two rings (bars/bars)
        # meet at two MIRROR-IMAGE points that can both sit in the lit region,
        # and no amount of preference distinguishes them -- a pair of 1D ring
        # measurements genuinely does not name one direction. Refuse those
        # rather than return the wrong one with a straight face; measured, the
        # naive pick was up to 118 deg out.
        n1 = np.linalg.norm(d1, axis=1)
        n2 = np.linalg.norm(d2, axis=1)
        safe = np.einsum('ij,ij->i', d1, d2) <= -0.999 * n1 * n2
        # Otherwise the pick is only trustworthy when exactly ONE root is in
        # front of the fly at all. A cone threshold around the patch center is
        # not enough: the patch reaches ~80 deg into its corners, so a true
        # direction can sit outside any cone tight enough to be useful.
        f1 = (d1 @ pref) > 0
        f2 = (d2 @ pref) > 0
        good = good & (safe | (f1 ^ f2))
    nrm = np.linalg.norm(pick, axis=1, keepdims=True)
    out[good] = (pick / np.where(nrm > 0, nrm, 1.0))[good]
    return out


def azel_of(directions):
    """(azimuth, elevation) in degrees for (..., 3) FLY-FRAME unit vectors."""
    d = np.asarray(directions, dtype=float)
    az = np.degrees(np.arctan2(d[..., 1], d[..., 0]))
    el = np.degrees(np.arcsin(np.clip(d[..., 2], -1.0, 1.0)))
    return az, el


def screen_azel_of(directions, center_az_deg, center_el_deg):
    """(azimuth, elevation) in degrees RELATIVE TO THE PATCH CENTRE: zero is the
    middle of the display on both axes.

    By ROTATION into the frame whose forward direction is the centre, not by
    subtracting the centre's az/el. Subtraction is only a rotation when the
    centre sits on the horizon; on a mount pitched by 38 degrees it stretches
    azimuth by 1/cos(38), bends elevation, and reports a direction AT the
    screen centre as 38 degrees high.
    """
    _gen, tess = _generator()
    R = tess.frame_rotation(float(center_az_deg), float(center_el_deg))
    return azel_of(np.asarray(directions, dtype=float) @ R)


def azel_from_pole_coords(meta, beta_deg, alpha_deg):
    """Fly-frame (azimuth, elevation) from a BARS coordinate and a WEDGES
    coordinate measured about the SAME pole.

    Bars bin the angle FROM the pole and wedges the azimuth ABOUT it, so the
    pair is a complete spherical coordinate system on that axis -- whatever the
    pole is. Together they name a point; either alone names a curve (a ring, or
    a meridian), which is why a single 1D measurement cannot give a position.

        d = cos(beta) * p + sin(beta) * (cos(alpha) * u + sin(alpha) * v)

    with (u, v, p) the frame wedge_reference_frame fixes, so alpha is measured
    from the same reference the stimulus binned it from.

    For a zenith pole this reduces to the identity -- elevation = 90 - beta and
    azimuth = alpha -- which is why that pair reads directly off the display.
    Any other pole measures the same point in a rotated frame, and converting
    is this function rather than a reason to exclude it.
    """
    pole, u, v, _ref = pole_frame(meta)
    b = np.radians(np.asarray(beta_deg, dtype=float))
    a = np.radians(np.asarray(alpha_deg, dtype=float))
    d = (np.cos(b)[..., None] * pole
         + np.sin(b)[..., None] * (np.cos(a)[..., None] * u
                                   + np.sin(a)[..., None] * v))
    az = np.degrees(np.arctan2(d[..., 1], d[..., 0]))
    el = np.degrees(np.arcsin(np.clip(d[..., 2], -1.0, 1.0)))
    return az, el


def screen_center_coord(meta, sampling):
    """The binned coordinate of the PATCH CENTRE. Returns (coord_deg, defined).

    This is the anchor that turns the stimulus's own coordinate into a
    screen-referenced one: beta of the patch center for bars, alpha of it for
    wedges. Subtracting it puts zero where the fly is looking.

    `defined` is False in one case, and it is a real configuration rather than
    an error: a wedge family whose POLE IS the patch center -- the dartboard,
    pole_el_deg = -5.5 on this bowl. Every wedge meets the center there, so the
    center has no azimuth about the pole and there is nothing to subtract. The
    caller falls back to the stimulus's own reference direction.
    """
    _gen, tess = _generator()
    pole, u, v, _ref = pole_frame(meta)
    c = tess.pole_vector(float(sampling.get('center_az_deg', 0.0)),
                         float(sampling.get('center_el_deg', 0.0)))
    if str(meta['mode']) == 'bars':
        return float(np.degrees(np.arccos(np.clip(float(np.dot(c, pole)), -1.0, 1.0)))), True
    if abs(float(np.dot(c, pole))) > np.cos(np.radians(1.0)):
        return 0.0, False
    return float(np.degrees(np.arctan2(float(np.dot(c, v)), float(np.dot(c, u))))), True


def cell_offsets_deg(meta, sampling):
    """(n_cells,) signed offset of each cell from the PATCH CENTRE, and whether
    that center was defined.

    WHY THIS RATHER THAN cell_positions_deg. The raw binned coordinate is
    anchored to the POLE and runs `phase_deg .. phase_deg + 360` for wedges, so
    two cells that are neighbours on the sphere can sit at opposite ends of the
    axis -- at phase 0 the first and last wedge are 345 deg apart on the axis
    and 15 deg apart in the world. Anchoring to the patch center instead puts
    the branch cut antipodal to the screen, which is never lit, so the visible
    cells are contiguous by construction and the wrap stops being a problem to
    work around. It also makes the number mean something physical: 0 is the
    middle of the display.

    SIGN. For bars the offset is `center - beta`, so POSITIVE is toward the
    pole: with the zenith pole that is above the patch center, i.e. exactly
    elevation relative to it. For wedges it increases anticlockwise about the
    pole, the same sense `wedge_reference_frame` fixes.

    In the degenerate wedge case (pole at the patch center) the offsets are the
    raw alpha wrapped to [-180, 180), so zero is the stimulus's own reference
    direction rather than the screen center; `defined` is False to say so.
    """
    pos = cell_positions_deg(meta)
    c0, ok = screen_center_coord(meta, sampling)
    if str(meta['mode']) == 'bars':
        return c0 - pos, ok
    if not ok:
        return ((pos + 180.0) % 360.0) - 180.0, False
    return ((pos - c0 + 180.0) % 360.0) - 180.0, True


def screen_axis_label(meta, sampling):
    """Axis label for `cell_offsets_deg`, naming what the coordinate IS.

    The same offset means different things depending on the pole, and a reader
    cannot tell which from the numbers, so the label carries it.
    """
    mode = str(meta['mode'])
    el = float(meta['pole_el_deg'])
    az = float(meta['pole_az_deg'])
    _c0, ok = screen_center_coord(meta, sampling)
    at_zenith = abs(el - 90.0) < 1e-6
    at_forward = abs(el) < 1e-6 and abs(az) < 1e-6
    if mode == 'bars':
        if at_zenith:
            return 'Elevation re. patch center (deg, + = up)'
        if at_forward:
            return 'Eccentricity re. patch center (deg)'
        return 'Ring offset from patch center (deg), pole az {:g} el {:g}'.format(az, el)
    if not ok:
        return 'Polar angle about patch center (deg, 0 = stimulus reference)'
    if at_zenith:
        return 'Azimuth re. patch center (deg)'
    return 'Angle about pole re. patch center (deg), pole az {:g} el {:g}'.format(az, el)


def pole_label(meta):
    """Human axis label for the binned coordinate, analogous in spirit to
    strf_plot._pole_label but for THIS pole convention -- a fly-frame
    (az, el) pole, not the noizone rtz orientation about the optical axis.
    """
    mode = str(meta['mode'])
    az, el = float(meta['pole_az_deg']), float(meta['pole_el_deg'])
    at_zenith = abs(el - 90.0) < 1e-6
    at_forward = abs(el) < 1e-6 and abs(az) < 1e-6
    if mode == 'bars':
        if at_zenith:
            return 'Elevation from zenith (deg from pole, 0=up)'
        if at_forward:
            return 'Eccentricity from straight ahead (deg from pole)'
        return 'Angle from pole at (az={:g}, el={:g}) (deg)'.format(az, el)
    if at_zenith:
        return 'Azimuth about the zenith (deg, wraps at 360)'
    if at_forward:
        return 'Radial angle about straight ahead (deg, wraps at 360)'
    return 'Azimuth about pole at (az={:g}, el={:g}) (deg, wraps at 360)'.format(az, el)
