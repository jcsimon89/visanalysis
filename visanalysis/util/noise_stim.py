"""
Helpers for locating and parsing flymax noizone bar-noise stimulus files.

https://github.com/ClandininLab/visanalysis
"""
import json
import pathlib

import numpy as np


def _decode(x):
    return x.decode('utf-8') if isinstance(x, bytes) else str(x)


def is_flymax_movie(run_parameters):
    """
    True if a series' run_parameters indicate ANY flymax movie-playback
    protocol -- i.e. stimpack is just playing back a pre-rendered video
    (has a movie_bin_name candidate list at all), rather than rendering the
    stimulus live. Content-agnostic: matches noizone, the Edge/search
    stimulus, or anything else built the same way.
    """
    movie_bin_names = run_parameters.get('movie_bin_name', [])
    if isinstance(movie_bin_names, (str, bytes)):
        movie_bin_names = [movie_bin_names]
    return len(movie_bin_names) > 0


def is_noizone(run_parameters):
    """
    True if a series' run_parameters indicate a flymax noizone protocol,
    based on whether any of the series' candidate movie filenames
    (run_parameters['movie_bin_name']) contain "noizone". protocol_ID is
    not useful for this -- it's a generic protocol class ("FlymaxMovie")
    shared by every movie-playback stimulus, noise or not.
    """
    movie_bin_names = run_parameters.get('movie_bin_name', [])
    if isinstance(movie_bin_names, (str, bytes)):
        movie_bin_names = [movie_bin_names]
    return any('noizone' in _decode(name).lower() for name in movie_bin_names)


def _values_equal(a, b):
    """Numpy-aware equality: hdf5-read values come back as numpy scalars/
    arrays, json-derived ones as plain Python types."""
    try:
        return bool(np.array_equal(np.asarray(a), np.asarray(b)))
    except Exception:
        return a == b


def resolve_attrs_to_write(existing_params, new_attrs, label):
    """
    Compare new_attrs (e.g. resolved from a movie's params.json) against
    existing_params (an epoch's or series' current hdf5 attrs) before
    writing, to guard against silently overwriting a native stimpack
    parameter -- or anything else already present -- with a same-named but
    different-meaning cached movie param.

    - Key not present yet: included in the returned dict (safe to write).
    - Key present with the SAME value: excluded from the returned dict
      (nothing to write -- covers idempotent re-runs of process_data.py).
    - Key present with a DIFFERENT value: raises ValueError immediately,
      naming the exact key/values, rather than overwriting silently.

    Returns
    -------
    dict of just the keys that actually need to be written.
    """
    to_write = {}
    conflicts = []
    for key, new_val in new_attrs.items():
        if key in existing_params:
            old_val = existing_params[key]
            if _values_equal(old_val, new_val):
                continue  # already there with the same value -- nothing to do
            conflicts.append((key, old_val, new_val))
        else:
            to_write[key] = new_val
    if conflicts:
        lines = '\n'.join(f'  {k}: existing={o!r}  new={n!r}' for k, o, n in conflicts)
        raise ValueError(
            f'{label}: cached movie params would overwrite existing attribute(s) '
            f'with a DIFFERENT value -- naming collision between the flymax '
            f'params.json and native/existing fields:\n{lines}\n'
            f'Rename the conflicting field(s) in the flymax stimulus generator '
            f'params.json to avoid this.'
        )
    return to_write


def unique_values(values):
    """
    De-duplicate values (numpy-aware, via _values_equal, so unhashable
    numpy arrays don't break a plain set()), sorted if the values happen to
    be sortable, in original order otherwise.
    """
    seen = []
    for v in values:
        if not any(_values_equal(v, s) for s in seen):
            seen.append(v)
    try:
        return sorted(seen)
    except TypeError:
        return seen


def get_epoch_bin_filename(epoch_params):
    """
    Extract the noise .bin filename from one epoch's parameter dict
    (one entry of ID.getEpochParameters()).

    Stimpack saves 'movie_bin_name' (the bin file's stem, no extension) and
    'movie_bin_path' (the full path on the stimulus computer, not valid
    here) on every epoch. The .bin file on disk is that stem + '.bin'.
    """
    return epoch_params['movie_bin_name'] + '.bin'


NOISE_PARAM_KEYS = ('bin_path', 'rtz', 'noisetau')


def resolve_noizone_epoch_noise_params(epoch_params):
    """
    Return this epoch's noise-stimulus params, cached onto the epoch by
    process_data.py.

    Only three are needed now: bin_path locates the movie's master, and
    rtz/noisetau group epochs into filters. Everything else about the stimulus
    -- bar width, phase, pole, pose, phimx, dm, timebase -- is read from the
    master itself, which is authoritative. In particular the legacy 'width'
    epoch attr is deliberately NOT required or used: it is a square-wave
    PERIOD, while the master's bar_geometry.width_deg is the width of one BAR,
    and conflating them is what produced a silent 2x error before.

    Raises KeyError if any are missing -- process_data.py (with
    --flymax_movies_root) must be run for this fly first. No fallback
    resolution strategy here on purpose.
    """
    missing = [k for k in NOISE_PARAM_KEYS if k not in epoch_params]
    if missing:
        raise KeyError(
            f'Epoch is missing cached noise params {missing}. '
            f'Run process_data.py (with --flymax_movies_root) for this fly first.'
        )
    return {k: epoch_params[k] for k in NOISE_PARAM_KEYS}


def movie_params_from_master(bin_path):
    """
    Epoch attrs for a movie that ships a noiz1d master instead of a params.json.

    noiz1d_make.m writes <base>_master.json / _values.bin / _barindex.bin and
    NO _params.json, so process_data.py cannot resolve such a movie the legacy
    way. This maps the master onto the same attr names downstream code expects.

    Deliberately does NOT emit a 'width' key. The legacy params.json 'width' is
    a square-wave PERIOD; the master's is a single BAR. Reusing the name would
    hand the old 2x trap to anyone reading epoch attrs, so the bar width is
    exported under its own unambiguous name and the analysis takes it from the
    master regardless.

    Returns None if there is no master next to the .bin.
    """
    try:
        mp = find_master(bin_path)
    except FileNotFoundError:
        return None
    m = json.loads(mp.read_text())
    dmp, dmt = (int(v) for v in m['render']['dm'])
    return {
        'stimtype': 'NoizOne',
        'master_file': mp.name,
        'dmp': dmp,
        'dmt': dmt,
        'rtz': float(m['bar_geometry']['orientation_deg']),
        'noisetau': float(m['timebase']['tau_actual_ms']),
        'bar_width_deg': float(m['bar_geometry']['width_deg']),
        'bar_phase_deg': float(m['bar_geometry']['phase_deg']),
        'num_bars': int(m['bar_geometry']['num_bars']),
        'noisedur': float(m['timebase']['duration_actual_s']) / 60.0,
    }


def resolve_movie_params(bin_path):
    """
    Stimulus params for any flymax movie, from whichever sidecar it ships.

    noiz1d movies carry a master; legacy noizone and the Edge/Flash movies
    carry a _params.json. Master wins when both exist, since it is the
    generator's own complete specification rather than a summary.
    """
    from_master = movie_params_from_master(bin_path)
    if from_master is not None:
        return from_master
    params_json = load_stim_params_json(bin_path)
    if params_json is None:
        raise FileNotFoundError(
            'No sidecar found next to {}: expected either a noiz1d '
            '<base>_master.json or a legacy <base>_params.json.'.format(
                pathlib.Path(bin_path).name))
    return params_json


def mesh_shape_from_path(mesh_path):
    """
    (n_phi, n_theta) from a display mesh filename, e.g.
    'mesh_7_31_stimpack.txt' -> (7, 31).

    Parsed exactly the way GlMesh does it at display time (clandinin_labpack
    shapes.py, via stimuli.py:1224-1225), so the analysis reconstructs the same
    grid the stimulus was actually drawn through.
    """
    name = pathlib.Path(_decode(mesh_path)).name
    parts = name.split('_')
    try:
        return int(parts[1]), int(parts[2])
    except (IndexError, ValueError):
        raise ValueError(
            f'Could not parse mesh dimensions from {name!r}; expected '
            f'mesh_<nphi>_<ntheta>_*.txt'
        )


def resolve_noizone_series_mesh_shape(run_parameters):
    """
    The display mesh a series was shown through, as (n_phi, n_theta).

    stimpack writes protocol_parameters as attributes on the SERIES group
    (experiment/data.py:107-108) and 'mesh_path' is one of them, so it can be
    read once per series rather than per epoch. (It is also written onto every
    epoch -- twice, since FlymaxMovie puts it in both epoch_stim_parameters and
    epoch_protocol_parameters, which land as the same attr -- but the series
    copy is the natural place to read it.)

    The stored path points at the stimulus computer and will not exist here;
    only the basename is used.

    Raises KeyError if absent and ValueError if the series was configured with
    several different meshes, deliberately: silently defaulting would apply the
    wrong position correction with no indication anything was guessed.
    """
    if 'mesh_path' not in run_parameters:
        raise KeyError(
            "Series has no 'mesh_path' attribute, so the display mesh cannot be "
            "determined and the bar-position correction would have to guess. "
            "Pass mesh_shape explicitly if you know which mesh was used."
        )
    value = run_parameters['mesh_path']
    if isinstance(value, (str, bytes)):
        candidates = [value]
    else:
        candidates = list(np.atleast_1d(value))

    shapes = unique_values([mesh_shape_from_path(v) for v in candidates])
    if len(shapes) > 1:
        raise ValueError(
            f'Series lists more than one display mesh {shapes}; the mesh would '
            f'differ between epochs, so it cannot be resolved once at the series '
            f'level. Resolve per epoch from each epoch\'s own mesh_path attr.'
        )
    return tuple(shapes[0])


def find_bin_files(root_dir, filenames):
    """
    Recursively search root_dir for each bin filename in filenames.

    Parameters
    ----------
    root_dir  : str or Path -- root directory to search (e.g. flymax_movies)
    filenames : iterable of str -- exact .bin file basenames to locate

    Returns
    -------
    dict {filename: full_path (str)}

    Raises
    ------
    FileNotFoundError if any filename is not found anywhere under root_dir.
    ValueError if any filename is found in more than one location under
    root_dir (ambiguous -- requires manual resolution).
    """
    wanted = set(filenames)
    found = {}  # filename -> list of matching paths
    for p in pathlib.Path(root_dir).rglob('*.bin'):
        if p.name in wanted:
            found.setdefault(p.name, []).append(p)

    missing = wanted - found.keys()
    if missing:
        raise FileNotFoundError(
            f'Could not find {len(missing)} bin file(s) under {root_dir}: {sorted(missing)}'
        )

    duplicates = {name: paths for name, paths in found.items() if len(paths) > 1}
    if duplicates:
        lines = '\n'.join(
            f'  {name}: {[str(p) for p in paths]}' for name, paths in duplicates.items()
        )
        raise ValueError(
            f'Found multiple files under {root_dir} matching the same bin filename '
            f'(ambiguous, needs manual resolution):\n{lines}'
        )

    return {name: str(paths[0]) for name, paths in found.items()}


def load_stim_params_json(bin_path):
    """
    Load the _params.json sidecar file that accompanies a noizone .bin file.

    The .bin filename's stem already ends in '_' (e.g. '..._001_.bin'), so the
    sidecar is just that stem + 'params.json' (e.g. '..._001_params.json').
    Returns the full parsed dict verbatim (all fields, e.g. dmp, dmt, rtz,
    noisetau, width, noisedur, bf, stimtype), or None if the file does not
    exist.
    """
    bin_path = pathlib.Path(bin_path)
    params_path = bin_path.parent / (bin_path.stem + 'params.json')
    if not params_path.exists():
        return None
    with open(params_path, 'r') as f:
        return json.load(f)


# ── projector cap geometry ──────────────────────────────────────────────────
# dmp indexes phi, the POLAR ANGLE from the projector's optical axis, spanning
# [0, phimx]. dmt indexes theta, the AZIMUTHAL angle around that axis, spanning
# the full [-pi, +pi]. phi is NOT elevation and theta is NOT azimuth: the axis
# is the projector's, not the fly's, and theta covers 360 deg regardless of how
# much of the sphere is actually lit.
#
# phimx is where the projector's rays run tangent to the sphere. It follows from
# the flymax mesh geometry (pmeshdf.m defaults -> sphere2plane.m -> circleToLine.m:11):
#     half_angle = c / r / 2
#     d          = 2*r*t*a*sin(half_angle) + r*cos(half_angle)
#     phimx      = pi/2 - asin(r / d)
# Cross-checked against the movies themselves: fitting phimx purely from the
# requirement that recovered bar boundaries land on exact multiples of width/2
# gives 76.335 deg, against 76.308 deg from this formula (agreement < 2 arcmin).
SCREEN_GEOMETRY = {'s': dict(r=71.5, c=110.0),     # new small screen
                   'l': dict(r=77.5, c=105.0)}     # old large screen
_THROW = 1.57523511    # projector throw ratio, measured (pmeshdf.m:23)
_ASPECT = 1.6          # projector aspect ratio, width/height
_HALF_H = 1.0          # half-height of the projection plane (pmeshdf.m:25)


def compute_phimx_deg(screen='s'):
    """Max phi (deg) of the projector cap, for a given flymax spherical screen."""
    g = SCREEN_GEOMETRY[screen]
    r, c = g['r'], g['c']
    half_angle = c / r / 2.0
    d = 2 * r * _THROW * _ASPECT * np.sin(half_angle) + r * np.cos(half_angle)
    return float(np.degrees(np.pi / 2 - np.arcsin(r / d)))


PHIMX_DEG = compute_phimx_deg('s')     # 76.308 deg -- previously hardcoded 77.0


# ────────────────────────────────────────────────────────────────────────────
# bar-identity stimulus regressor
#
# A noizone movie is not really a 36000-frame movie. cmake_noizone.m renders a
# 1D bar template rotated once onto the projector cap, then looks up one
# ternary noise value per bar per update. The entire stimulus is therefore
# (vals, bar_idx) -- see cmake_noizone.m:52 -- and the exact 1D regressor is
# just vals indexed by bar. No projection, binning, or interpolation is
# involved, so nothing below approximates anything.
#
# Bars are bands of constant psi about the pole axis a_hat = (cos rtz, sin rtz,
# 0), i.e. cos(psi) = sin(phi) * cos(theta - rtz), so iso-bar sets are circles
# on the viewing sphere. Do not treat them as straight parallel stripes in any
# flattened (x, y) map: at phi = phimx, theta - rtz = 45 deg the two differ by
# 10.6 deg, over two bar widths.
# ────────────────────────────────────────────────────────────────────────────

# ────────────────────────────────────────────────────────────────────────────
# display mesh warp
#
# stimpack does not warp the movie analytically. It texture-maps it onto GlMesh
# (clandinin_labpack shapes.py:142), a coarse vertex grid read from
# mesh_<nphi>_<ntheta>_stimpack.txt, and each quad is split into two triangles
# (stimpack shapes.py:83, diagonal (ii,jj)->(ii+1,jj+1)) that are AFFINE in
# texture coordinates. The true sphere->image-plane warp
#     e(phi) = b*r*sin(phi) / (d - r*cos(phi))
# is nonlinear and concave near the field edge, so linear interpolation across
# a quad falls below it and content is drawn at too small a radius -- i.e.
# displaced toward the middle of the field.
#
# With the 7 x 31 mesh used for the 2026-07-29 recordings this shifts the
# outermost bars inward by up to 4.07 deg, 81% of a 5 deg bar. Both ends shift
# toward the centre so they cancel in a global mean while being strongly
# systematic per bar -- which is exactly what a position axis must not do.
#
# The displacement is almost entirely radial, but do NOT read that as "only the
# phi sampling matters". The coarse theta sampling walks a polygon instead of a
# circle, and that chordal contraction is itself radial: at 12 deg steps the
# mid-edge sits at e*cos(6 deg), 0.55% inside the true radius. Near the cap edge
# that tiny radial error explodes, because phimx is the tangency angle where
# de/dphi -> 0: at phi = 76.0 deg a radial error of 0.001 is 8 deg of phi, and
# at 76.3 deg it is 330 deg. So refining phi alone saturates (13 x 31 -> 39 x 31
# barely moves the max), and at equal triangle count refining BOTH wins:
#     7 x 31   (360 tri)   median 0.717   max 4.072 deg
#    26 x 31  (1500 tri)   median 0.312   max 3.289
#    13 x 61  (1440 tri)   median 0.202   max 2.231   <- better, same cost
#    26 x 121 (6000 tri)   median 0.052   max 1.056
# The residual max always sits on the outermost bar, which is ill-conditioned
# by geometry rather than by mesh resolution and stays uncertain regardless.
# ────────────────────────────────────────────────────────────────────────────

_MESH_WARP_CACHE = {}


def mesh_warp_directions(dmp, dmt, mesh_shape=(7, 31), screen='s'):
    """
    Where each texel of a (dmp, dmt) movie actually lands on the sphere, given
    the display mesh it is drawn through.

    Returns (phi_i, th_i, phi_r, th_r, lit), each (dmp, dmt); angles in radians.
    phi_i/th_i are where a texel should land, phi_r/th_r where it does, and lit
    is the projector frustum mask (verified pixel-exact against the movies).
    """
    key = (dmp, dmt, tuple(mesh_shape), screen)
    if key in _MESH_WARP_CACHE:
        return _MESH_WARP_CACHE[key]

    g = SCREEN_GEOMETRY[screen]
    r, c = g['r'], g['c']
    half = c / r / 2.0
    d = 2 * r * _THROW * _ASPECT * np.sin(half) + r * np.cos(half)
    b = 2 * _HALF_H * _THROW * _ASPECT
    phimx = np.pi / 2 - np.arcsin(r / d)

    def e_of_phi(phi):
        return b * r * np.sin(phi) / (d - r * np.cos(phi))

    _pg = np.linspace(0, phimx, 200_001)
    _eg = e_of_phi(_pg)

    n_phi, n_theta = mesh_shape
    phi_v = np.linspace(0, phimx, n_phi)
    th_v = np.linspace(-np.pi, np.pi, n_theta)
    PX = e_of_phi(phi_v)[:, None] * np.cos(th_v)[None, :]
    PY = e_of_phi(phi_v)[:, None] * np.sin(th_v)[None, :]

    phi_i = np.linspace(0, phimx, dmp)[:, None] * np.ones((1, dmt))
    th_i = np.ones((dmp, 1)) * np.linspace(-np.pi, np.pi, dmt)[None, :]

    u = phi_i / phimx * (n_phi - 1)          # texture coord in mesh-vertex units
    v = (th_i + np.pi) / (2 * np.pi) * (n_theta - 1)
    ii = np.clip(np.floor(u).astype(int), 0, n_phi - 2)
    jj = np.clip(np.floor(v).astype(int), 0, n_theta - 2)
    aa, bb = u - ii, v - jj

    lower = bb <= aa    # GlQuad triangle (0,0),(1,0),(1,1); upper is (0,0),(1,1),(0,1)
    rx = np.where(lower,
                  (1 - aa) * PX[ii, jj] + (aa - bb) * PX[ii + 1, jj] + bb * PX[ii + 1, jj + 1],
                  (1 - bb) * PX[ii, jj] + aa * PX[ii + 1, jj + 1] + (bb - aa) * PX[ii, jj + 1])
    ry = np.where(lower,
                  (1 - aa) * PY[ii, jj] + (aa - bb) * PY[ii + 1, jj] + bb * PY[ii + 1, jj + 1],
                  (1 - bb) * PY[ii, jj] + aa * PY[ii + 1, jj + 1] + (bb - aa) * PY[ii, jj + 1])

    phi_r = np.interp(np.hypot(rx, ry), _eg, _pg)
    th_r = np.arctan2(ry, rx)

    e_i = e_of_phi(phi_i)
    lit = (np.abs(e_i * np.cos(th_i)) < _ASPECT * _HALF_H) & (np.abs(e_i * np.sin(th_i)) < _HALF_H)

    _MESH_WARP_CACHE[key] = (phi_i, th_i, phi_r, th_r, lit)
    return _MESH_WARP_CACHE[key]


def bar_psi_from_direction(phi, theta, pole_out):
    """
    psi (deg): the angle between each movie direction and the bar pole axis.

    pole_out is the pole as a unit vector in the OUTPUT (movie) frame -- i.e.
    P @ pole_fly, using the master's own pose matrix. Taking the pole as a
    vector rather than as an rtz angle means any screen pose works, not just
    the default one. The angle itself is frame-independent, so this equals the
    master's `beta`.
    """
    p = np.asarray(pole_out, dtype=float).ravel()
    sp = np.sin(phi)
    cosb = sp * np.cos(theta) * p[0] + sp * np.sin(theta) * p[1] + np.cos(phi) * p[2]
    return np.degrees(np.arccos(np.clip(cosb, -1.0, 1.0)))


def _weighted_bar_centroid(bar_idx, psi, weight):
    flat_bar = bar_idx.astype(int).ravel()
    n = flat_bar.max() + 1
    tot_w = np.bincount(flat_bar, weights=weight.ravel(), minlength=n)
    tot_wp = np.bincount(flat_bar, weights=(weight * psi).ravel(), minlength=n)
    return {int(bi): float(tot_wp[bi] / tot_w[bi])
            for bi in range(1, n) if tot_w[bi] > 0}


def effective_bar_psi(bar_idx, pole_out, mesh_shape=(7, 31), screen='s'):
    """
    Solid-angle-weighted psi centroid of each bar, in degrees -- both as
    actually displayed and as it would be with a perfect display mesh.

    Uses the movie's own bar_idx (which already encodes how cmake_noizone
    assigned pixels to bars, and is already zeroed outside the frustum) rather
    than nominal 5 deg boundaries, and weights each texel by the solid angle it
    covers on the sphere, so the result is the centroid of the region the eye
    integrated over rather than of the texture grid.

    pole_out is the bar pole in the output frame (P @ pole_fly from the master).

    Two separate effects move a bar away from its nominal phase + (bar-0.5)*width:

      * CLIPPING -- bars at the edge of the projector cap are only partly
        displayed, so the centroid of the visible part is pulled inward. This
        is geometry, not a display artefact, and no mesh can fix it.
      * MESH WARP -- affine interpolation across coarse quads draws content
        inward of where it belongs (see the section header above).

    Returns
    -------
    (psi_displayed, psi_geometric) : two dicts {bar_index: psi_deg}
        psi_displayed includes both effects; psi_geometric includes only
        clipping. Their difference isolates the mesh contribution.
    """
    dmp, dmt = bar_idx.shape
    phi_i, th_i, phi_r, th_r, lit = mesh_warp_directions(dmp, dmt, mesh_shape, screen)

    # solid angle per texel after warping. theta is unwrapped via (th_r - th_i),
    # which is tiny and continuous, to avoid the atan2 branch cut at +-pi.
    dth = np.arctan2(np.sin(th_r - th_i), np.cos(th_r - th_i))
    dphi_di, dphi_dj = np.gradient(phi_r, axis=0), np.gradient(phi_r, axis=1)
    dth_di, dth_dj = np.gradient(dth, axis=0), np.gradient(dth, axis=1) + 2 * np.pi / (dmt - 1)
    det_j = np.abs(dphi_di * dth_dj - dphi_dj * dth_di)

    w_disp = np.where(lit, np.sin(phi_r) * det_j, 0.0)
    w_geom = np.where(lit, np.sin(phi_i), 0.0)      # identity warp, det_j = 1

    return (_weighted_bar_centroid(bar_idx, bar_psi_from_direction(phi_r, th_r, pole_out), w_disp),
            _weighted_bar_centroid(bar_idx, bar_psi_from_direction(phi_i, th_i, pole_out), w_geom))


def find_master(bin_path, bardata_dir=None):
    """
    Locate a movie's noiz1d master JSON.

    noiz1d_make.m names it fn_base + '_master.json', where fn_base lacks the
    trailing underscore the .bin filename carries (namefile.m:41 strips it), so
    '..._001_.bin' -> '..._001_master.json'.

    Looked for next to the .bin first, so a natively generated master is picked
    up with no configuration. bardata_dir is an optional fallback.
    """
    bin_path = pathlib.Path(bin_path)
    name = bin_path.stem.rstrip('_') + '_master.json'
    candidates = [bin_path.parent / name]
    if bardata_dir is not None:
        candidates.append(pathlib.Path(bardata_dir) / name)
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        'No noiz1d master found for {}; looked in {}. Movies generated before '
        'the stimulus code wrote a master can be converted from the rendered '
        '.bin with legacy_to_noiz1d.py.'.format(
            bin_path.name, ', '.join(str(c.parent) for c in candidates))
    )


def load_noiz1d_master(bin_path, bardata_dir=None):
    """
    Read a noiz1d master and its two binary sidecars.

    Returns (master_dict, values, bar_index):
        values     uint8  (num_bars, num_updates)  -- one column per noise update
        bar_index  uint16 (dmp, dmt)               -- 0 means "not displayed"
    """
    mp = find_master(bin_path, bardata_dir)
    m = json.loads(mp.read_text())
    nb, nu = (int(v) for v in m['values']['shape'])
    values = np.fromfile(mp.parent / m['values']['file'], dtype=np.uint8).reshape(nb, nu)
    dmp, dmt = (int(v) for v in m['render']['dm'])
    bar_index = np.fromfile(mp.parent / m['render']['barindex_file'],
                            dtype=np.dtype(m['render']['barindex_dtype'])).reshape(dmp, dmt)
    return m, values, bar_index


def load_bar_stim_1d(bin_path, bardata_dir=None,
                     mesh_correct=True, mesh_shape=None, screen='s'):
    """
    Exact 1D noise regressor for a noizone / noiz1d movie, from its master.

    Everything comes from the master -- bar width, phase, pole, pose, phimx and
    the timebase -- so nothing is re-derived from filenames or params.json, and
    the legacy "width is a PERIOD, halve it" trap cannot recur: the master's
    bar_geometry.width_deg is already the width of ONE BAR.

    values is stored per noise UPDATE; it is expanded here by frames_per_update
    so stim_1d stays at raw display-frame resolution, which is what per-trial
    dropped-frame correction indexes against.

    Only bars that actually land on the projector cap are returned -- 32 of 36
    at orientation 0, 22 of 36 at orientation 90 on the standard flymax screen.

    Position convention
    -------------------
    psi is the angle from this movie's pole axis, so psi = 0 sits at the pole,
    outside the lit field. The returned axis is centred on the field:

        bar_pos_deg = 90 - psi

    putting 0 at the projector's optical axis and making the axis signed and
    symmetric. Positive is anterior at orientation 0 and dorsal at orientation
    90; at 90 it is exactly elevation.

    Mesh correction
    ---------------
    With mesh_correct=True each bar is labelled with where it was ACTUALLY
    displayed -- its solid-angle-weighted psi centroid after the display mesh
    warp -- rather than its nominal phase + (bar - 0.5) * width. mesh_shape is
    REQUIRED when it is on; read it from the series metadata with
    resolve_noizone_series_mesh_shape().

    Returns
    -------
    stim_1d     : (n_frames, n_bars_on_screen) float32, NOT mean-subtracted
    bar_pos_deg : (n_bars_on_screen,) ascending degrees, signed as above
    meta        : dict -- bar_index, psi_centre_deg, psi_nominal_deg,
                  psi_geometric_deg, clip_shift_deg, mesh_shift_deg,
                  bar_width_deg, orientation_deg, pole_fly, mesh_correct,
                  mesh_shape, master_path
    """
    m, values, bar_index = load_noiz1d_master(bin_path, bardata_dir)

    fpu = int(m['timebase']['frames_per_update'])
    stim_full = np.repeat(values, fpu, axis=1)              # -> raw display frames
    n_frames = int(m['timebase']['total_frames'])
    if stim_full.shape[1] != n_frames:
        raise ValueError(
            'master timebase inconsistent for {}: values {} x {} expanded by {} '
            'gives {} frames, but total_frames says {}'.format(
                pathlib.Path(bin_path).name, values.shape[0], values.shape[1],
                fpu, stim_full.shape[1], n_frames))

    # Which bars the fly actually saw. The master's bars_used is authoritative:
    # noiz1d_make renders the FULL projector rectangle, so its bar_index is
    # 1..num_bars everywhere and "index > 0" would wrongly admit bars that fall
    # outside the lit frustum -- they would enter the design matrix as
    # regressors with no corresponding visual input. Legacy converted movies are
    # masked and carry 0 off-screen, so the fallback is correct for those.
    bars_used = m['bar_geometry'].get('bars_used')
    if bars_used:
        on_screen = np.unique(np.asarray(bars_used, dtype=int))
    else:
        on_screen = np.unique(bar_index[bar_index > 0]).astype(int)
    stim_1d = stim_full[on_screen - 1, :].T.astype(np.float32)

    geom = m['bar_geometry']
    bar_width_deg = float(geom['width_deg'])                # ONE bar, not a period
    phase_deg = float(geom['phase_deg'])
    psi_nominal = phase_deg + (on_screen - 0.5) * bar_width_deg

    # the pole, rotated into the movie frame using the master's own pose, so any
    # screen pose works rather than only the default one
    P = np.asarray(m['pose']['matrix_P'], dtype=float).reshape(3, 3)
    pole_out = P @ np.asarray(geom['pole_fly'], dtype=float).ravel()

    if mesh_correct:
        if mesh_shape is None:
            raise ValueError(
                'mesh_correct=True requires mesh_shape -- the correction depends on '
                'which display mesh the epoch was shown through. Read it from the '
                'series metadata with resolve_noizone_series_mesh_shape(), or pass '
                'mesh_correct=False to skip the correction entirely.')
        # the master records the phimx the stimulus was built with; if it
        # disagrees with this screen's geometry the correction would be applied
        # in the wrong frame, so fail rather than silently proceed
        phimx_master = float(m['render']['phimx_deg'])
        phimx_here = compute_phimx_deg(screen)
        if abs(phimx_master - phimx_here) > 0.05:
            raise ValueError(
                "master phimx {:.3f} deg does not match screen '{}' ({:.3f} deg) -- "
                'wrong screen, or the master came from different rig geometry.'
                .format(phimx_master, screen, phimx_here))
        disp, geo = effective_bar_psi(bar_index, pole_out, mesh_shape, screen)
        psi_centre = np.array([disp[int(b)] for b in on_screen])
        psi_geom = np.array([geo[int(b)] for b in on_screen])
    else:
        psi_centre = psi_nominal.astype(float)
        psi_geom = psi_nominal.astype(float)

    bar_pos_deg = 90.0 - psi_centre
    order = np.argsort(bar_pos_deg)          # keep the position axis ascending
    return (stim_1d[:, order],
            bar_pos_deg[order],
            {'bar_index': on_screen[order],
             'psi_centre_deg': psi_centre[order],
             'psi_nominal_deg': psi_nominal[order],
             'psi_geometric_deg': psi_geom[order],
             'clip_shift_deg': (psi_geom - psi_nominal)[order],
             'mesh_shift_deg': (psi_centre - psi_geom)[order],
             'bar_width_deg': bar_width_deg,
             'orientation_deg': float(geom['orientation_deg']),
             'pole_fly': list(geom['pole_fly']),
             'mesh_correct': mesh_correct,
             'mesh_shape': tuple(mesh_shape) if mesh_shape else None,
             'master_path': str(find_master(bin_path, bardata_dir))})
