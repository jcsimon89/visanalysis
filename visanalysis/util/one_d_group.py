"""
one_d_group.py

Registering 1D ternary STRFs to each ROI's own receptive-field center and
averaging them, to answer "what does this cell type's filter look like".

WHY REGISTRATION IS CHEAP HERE. The per-fly analysis already reports every cell
in a SCREEN-REFERENCED coordinate (one_d_stim.cell_offsets_deg), so every ROI
arrives on a common axis before anything is aligned. The only shift left is each
ROI's own centroid. average_strf_noizone and average_strf_octave both have to do
considerably more -- rotate on the sphere, or reconcile tessellations that differ
between flies -- and none of that applies.

WHY BINNING AND NOT INTERPOLATION. Cell values ARE the measurement; a cell is
15 degrees wide and nothing was measured inside it. Interpolating onto a finer
grid would invent structure at a scale the stimulus cannot resolve. Registered
samples are therefore binned, the same choice average_strf_noizone makes. Because
the centroids are now sub-cell (see one_d_strf.centroid_1d), the shifted samples
do NOT land on a common grid, and with enough ROIs the bins can be finer than one
cell -- the population fills in what one ROI cannot.

WHY z UNITS. Filters are averaged after z-scoring against each ROI's own noise
(the analysis saves strf_z). Averaging raw filters lets the loudest ROI dominate
the mean; normalising each ROI by its own PEAK biases the shape, because the peak
is where that ROI's noise is largest and dividing by it shrinks exactly the
feature being measured. Dividing by the noise level does neither.

WHY NOT SIGN-ALIGN. Never flip an ROI to match the others before averaging:
aligning each to its own extremum makes noise add constructively there, so the
mean grows a peak even when nothing carries signal. util.nullcal says the same at
more length. Mixed populations are SPLIT by polarity, not folded together.
"""

import numpy as np


def auto_grid_step(cell_width_deg, rel_deg, target_per_bin=8, max_subdiv=4):
    """Bin width for the registered grid, in degrees, chosen from the DATA.

    One bin per cell is always safe. Sub-dividing is licensed by how many ROIs
    actually land in each finer bin, which depends on how the centroids scatter
    and not just on how many ROIs there are: a hundred ROIs whose fields sit on
    top of each other fill one bin, not four. So the candidate steps are tried
    from finest to coarsest and the finest is taken whose typical CENTRAL bin
    still holds `target_per_bin` samples.

    An earlier version keyed the subdivision on the ROI count alone, which is
    the wrong quantity: it cannot know whether those ROIs' fields overlap.

    THE TRADE-OFF IS REAL AND IS NOT SETTLED BY CORRELATION WITH A KNOWN
    ANSWER. Finer bins resolve more and average fewer samples each, so a
    recovered profile scored by correlation against ground truth gets BETTER as
    the bins get coarser -- measured, 0.968 at one bin per cell against 0.935 at
    a quarter of one. That metric rewards smoothing, and a 15 degree bin cannot
    represent a 14 degree center at all, so following it would throw away the
    structure the average exists to show. The criterion here is therefore
    occupancy, not agreement: subdivide while the bins are still filled, and let
    the resolution be whatever the population supports.

    `rel_deg` is the registered relative positions, any shape; NaNs ignored.
    """
    rel = np.asarray(rel_deg, dtype=float).ravel()
    rel = rel[np.isfinite(rel)]
    if rel.size == 0:
        return float(cell_width_deg)
    # the central half of the registered range, where the population actually is
    lo, hi = np.percentile(rel, [25, 75])
    for subdiv in range(int(max_subdiv), 0, -1):
        step = float(cell_width_deg) / subdiv
        edges = np.arange(lo, hi + step, step)
        if edges.size < 2:
            continue
        counts = np.histogram(rel, bins=edges)[0]
        if counts.size and np.median(counts) >= target_per_bin:
            return step
    return float(cell_width_deg)


def register_average(filters, positions_deg, centers_deg, cell_width_deg,
                     mode='bars', grid_step=None, min_roi_per_bin=3,
                     half_width_deg=None):
    """Shift each ROI's filter onto its own RF center and average.

    filters      : (n_roi, n_cells, n_lags), already z-scored per ROI
    positions_deg: (n_cells,) screen-referenced offset of each cell
    centers_deg  : (n_roi,) each ROI's RF center in that same frame
    cell_width_deg : nominal cell width, which sets the default bin width
    mode         : 'wedges' wraps the relative coordinate to +-180, 'bars' does
                   not -- beta runs pole to antipole and a difference of two
                   positions is already signed and bounded
    grid_step    : bin width; None asks auto_grid_step
    min_roi_per_bin : bins with fewer contributing ROIs are returned as NaN.
                   A mean over two ROIs is not a population mean, and leaving
                   it in makes the tails of the figure look like measurements.
    half_width_deg : half extent of the grid; None uses the widest registered
                   sample actually present

    Returns dict with 'grid_deg', 'mean', 'sem', 'n_roi_per_bin', 'n_roi',
    'grid_step_deg'.
    """
    filters = np.asarray(filters, dtype=float)
    n_roi, n_cells, n_lags = filters.shape
    positions_deg = np.asarray(positions_deg, dtype=float)
    centers_deg = np.asarray(centers_deg, dtype=float)

    rel = positions_deg[None, :] - centers_deg[:, None]      # (n_roi, n_cells)
    if str(mode) == 'wedges':
        rel = ((rel + 180.0) % 360.0) - 180.0

    ok_roi = np.isfinite(centers_deg)
    if not np.any(ok_roi):
        raise ValueError('no roi has a finite center to register to')

    step = (auto_grid_step(cell_width_deg, rel[ok_roi])
            if grid_step is None else float(grid_step))
    if half_width_deg is None:
        half_width_deg = float(np.nanmax(np.abs(rel[ok_roi]))) + 0.5 * step
    n_half = int(np.ceil(half_width_deg / step))
    grid = np.arange(-n_half, n_half + 1) * step

    total = np.zeros((grid.size, n_lags))
    total_sq = np.zeros((grid.size, n_lags))
    count = np.zeros((grid.size, n_lags))

    for r in np.flatnonzero(ok_roi):
        idx = np.round(rel[r] / step).astype(int) + n_half
        keep = (idx >= 0) & (idx < grid.size)
        f = filters[r][keep]
        good = np.isfinite(f)
        np.add.at(total, idx[keep], np.where(good, f, 0.0))
        np.add.at(total_sq, idx[keep], np.where(good, f * f, 0.0))
        np.add.at(count, idx[keep], good.astype(float))

    with np.errstate(invalid='ignore', divide='ignore'):
        mean = total / count
        var = total_sq / count - mean ** 2
        sem = np.sqrt(np.clip(var, 0.0, None) / count)
    thin = count < float(min_roi_per_bin)
    mean[thin] = np.nan
    sem[thin] = np.nan
    return dict(grid_deg=grid, mean=mean, sem=sem,
                n_roi_per_bin=count, n_roi=int(ok_roi.sum()),
                grid_step_deg=step)


def polarity_groups(signs, minority_tol=0.1):
    """Split ROI indices by filter polarity, keeping the pooled group too.

    Returns (groups, minority_fraction) with groups a dict of name -> indices.
    'all' is ALWAYS present. Above `minority_tol` the two signs are additionally
    reported on their own.

    Why keep 'all' when the population is mixed. A mean over both signs does
    describe neither half, and that is the reason the split exists -- but a fly
    whose minority happens to land at 11% would otherwise have no figure that
    any other fly could be compared against, and whether a fly splits is partly
    an accident of where it falls relative to one threshold. Measured across one
    dataset, two flies of seven split, and on one of them the two orientations
    of the SAME channel came out at 9.8% and 17.1% -- split in one and not the
    other. Carrying 'all' everywhere keeps the cross-fly comparison on one
    footing; the split groups are there to show what it is hiding.
    """
    signs = np.asarray(signs, dtype=int)
    if signs.size == 0:
        return {'all': np.array([], dtype=int)}, 0.0
    pos = np.flatnonzero(signs > 0)
    neg = np.flatnonzero(signs < 0)
    minor = min(pos.size, neg.size) / float(signs.size)
    out = {'all': np.arange(signs.size)}
    if minor > minority_tol:
        out['positive'] = pos
        out['negative'] = neg
    return out, minor


def coverage_bounds(count, frac_of_max=0.5):
    """Grid indices whose ROI count is at least `frac_of_max` of the best bin.

    The registered grid always tapers at its edges -- only the ROIs whose
    receptive field sat near the edge of the lit region contribute out there --
    so the tails are averages over a handful of unusual ROIs. Cropping to where
    the population actually is keeps the figure honest about what it shows.
    """
    per_bin = np.nanmax(np.asarray(count, dtype=float), axis=1)
    if not np.any(per_bin > 0):
        return np.array([], dtype=int)
    return np.flatnonzero(per_bin >= frac_of_max * per_bin.max())


# ─────────────────────────────────────────────────────────────────────────────
# radial symmetry: recover h(r) from the registered projection
# ─────────────────────────────────────────────────────────────────────────────
#
# WHAT THE MEASUREMENT IS. A bar or wedge cell integrates the receptive field
# over a whole strip, so a registered 1D profile is a PROJECTION of the 2D
# field, not a slice through it. Recovering the field from projections is
# tomography; with two families there are two angles, which is far too few for
# a general inversion, so one assumption is needed. Radial symmetry is the
# strongest available and reduces the field to a single radial profile h(r).
#
# WHY FORWARD-MODEL RATHER THAN INVERT. The textbook route is the inverse Abel
# transform, which differentiates the measured profile -- the same noise
# amplifying move that made the temporal deconvolution useless on real data.
# Instead the forward operator is built explicitly and a regularised least
# squares is solved. It never differentiates data, it is the same machinery
# octave_strf.deconvolve_octaves uses, and it is far better conditioned.
#
# WHY IT BEATS FITTING A SHAPE. A difference-of-Gaussians fit would also avoid
# differentiation, but it assumes the answer's form. This does not: h comes out
# of the data. It also gets the CELL WIDTH right -- strips are 15 degrees wide,
# not lines, and the overlap areas below account for that exactly, where a fit
# to cell-center values quietly ignores it.
#
# WHY THE TWO FAMILIES STACK. Under radial symmetry both measure the same h, so
# they are not two reconstructions to compare but twice the constraints on one.
# Stack them and the disagreement becomes a residual to report rather than an
# assumption to hope about.
#
# WHY EACH FAMILY NEEDS ITS OWN GAIN. They measure the same h only up to SCALE.
# Each family's registered mean comes from its own kept-roi set with its own
# SNR, and z units normalise per roi, not across populations -- so the two mean
# amplitudes are independently arbitrary. Forcing one h to explain both at face
# value turns that arbitrary difference into residual: measured, on a field
# that is exactly radially symmetric, a 1.5x amplitude difference alone gives
# 0.196, past the 0.15 at which the residual is read as evidence of
# elongation. One scalar per family, shared across lags, removes it -- shared,
# not per lag, because a per-lag gain would flatten each family's temporal
# profile, which is part of the answer rather than a nuisance.
#
# AND WHY THAT GAIN MUST NOT GO NEGATIVE. Only the MAGNITUDE is arbitrary. A
# negative gain would let the solve invert one family to match the other, which
# is not a rescaling -- it is the two families disagreeing about the sign of the
# receptive field, the single most important thing the residual exists to
# report. Unconstrained, two families of opposite polarity fit at a residual of
# 0.001 instead of 1.0. Clamped at zero, such a family simply stops
# constraining h and its own residual goes to 1, which is the honest answer.
# The usual cause is upstream: the dominant polarity GROUP differs between the
# two families, so they describe different roi populations.


def _disc_strip_area(radius, a, b):
    """Area of the infinite vertical strip x in [a, b] inside a disc of
    `radius` centerd at the origin. Closed form, so no quadrature."""
    radius = float(radius)
    if radius <= 0:
        return np.zeros(np.shape(a))

    def g(x):
        xc = np.clip(x, -radius, radius)
        out = (xc * np.sqrt(np.maximum(radius ** 2 - xc ** 2, 0.0))
               + radius ** 2 * np.arcsin(np.clip(xc / radius, -1.0, 1.0)))
        out = np.where(np.asarray(x) > radius, radius ** 2 * np.pi / 2.0, out)
        out = np.where(np.asarray(x) < -radius, -radius ** 2 * np.pi / 2.0, out)
        return out

    return g(np.asarray(b, dtype=float)) - g(np.asarray(a, dtype=float))


def projection_matrix(strip_centers_deg, strip_width_deg, r_edges_deg):
    """(n_strips, n_annuli) area of overlap between each strip and each annulus.

    This IS the forward model: a radially symmetric field that is h_j on
    annulus j projects to A @ h. Rows are the registered grid bins, columns the
    radial bins.
    """
    c = np.asarray(strip_centers_deg, dtype=float)
    a, b = c - float(strip_width_deg) / 2.0, c + float(strip_width_deg) / 2.0
    edges = np.asarray(r_edges_deg, dtype=float)
    A = np.zeros((c.size, edges.size - 1))
    prev = _disc_strip_area(edges[0], a, b)
    for j in range(1, edges.size):
        cur = _disc_strip_area(edges[j], a, b)
        A[:, j - 1] = cur - prev
        prev = cur
    return A


def radial_from_projection(profile, strip_centers_deg, strip_width_deg,
                           r_edges_deg=None, ridge=1e-2):
    """Radial profile h(r) behind a registered projection.

    profile : (n_bins,) or (n_bins, n_lags); NaN bins are dropped, which is how
              the thin edges of a registered grid are excluded
    ridge   : second-difference penalty, RELATIVE to trace(AtA)/trace(DtD), so
              it means the same thing whatever the units and binning. Small and
              fixed on purpose -- it is conditioning, not a fitted parameter.
              Report `residual` rather than tuning this to taste.

    Returns dict with 'r_centers_deg', 'h', 'residual', 'A'. `residual` is the
    relative misfit ‖A h − p‖ / ‖p‖ per lag.

    THAT RESIDUAL IS NOT A TEST OF RADIAL SYMMETRY, and it is worth being blunt
    about why. ONE projection cannot falsify the assumption: some radial h
    reproduces any 1D profile you hand it. Measured, projections of isotropic,
    2:1 and 4:1 elongated fields all fit to a residual of 0.0007 here. The test
    requires BOTH families constraining a single h -- see
    radial_from_projections, whose residual does discriminate.
    """
    p = np.asarray(profile, dtype=float)
    single = p.ndim == 1
    if single:
        p = p[:, None]
    c = np.asarray(strip_centers_deg, dtype=float)

    if r_edges_deg is None:
        step = float(np.median(np.diff(np.sort(c)))) if c.size > 1 else 1.0
        r_max = float(np.nanmax(np.abs(c))) + step / 2.0
        r_edges_deg = np.arange(0.0, r_max + step, step)
    r_edges_deg = np.asarray(r_edges_deg, dtype=float)

    A_full = projection_matrix(c, strip_width_deg, r_edges_deg)
    n_r = A_full.shape[1]
    D2 = np.diff(np.eye(n_r), 2, axis=0) if n_r > 2 else np.eye(n_r)

    h = np.full((n_r, p.shape[1]), np.nan)
    resid = np.full(p.shape[1], np.nan)
    for k in range(p.shape[1]):
        ok = np.isfinite(p[:, k])
        if ok.sum() < 3:
            continue
        A = A_full[ok]
        y = p[ok, k]
        AtA, DtD = A.T @ A, D2.T @ D2
        lam = float(ridge) * np.trace(AtA) / max(np.trace(DtD), 1e-12)
        h[:, k] = np.linalg.solve(AtA + lam * DtD, A.T @ y)
        pred = A @ h[:, k]
        denom = float(np.linalg.norm(y))
        resid[k] = float(np.linalg.norm(pred - y) / denom) if denom > 0 else np.nan

    r_centers = 0.5 * (r_edges_deg[:-1] + r_edges_deg[1:])
    if single:
        h, resid = h[:, 0], float(resid[0])
    return dict(r_centers_deg=r_centers, h=h, residual=resid, A=A_full)


def radial_to_image(r_centers_deg, h, extent_deg=None, n_pix=101):
    """Render h(r) as a 2D map, for display.

    Under radial symmetry the image carries no information the radial profile
    does not, so this is a rendering rather than a reconstruction -- worth
    saying plainly, because a 2D picture invites the reader to believe 2D
    structure was measured.
    """
    r_centers_deg = np.asarray(r_centers_deg, dtype=float)
    h = np.asarray(h, dtype=float)
    if extent_deg is None:
        extent_deg = float(np.nanmax(r_centers_deg))
    ax = np.linspace(-extent_deg, extent_deg, int(n_pix))
    rr = np.hypot(*np.meshgrid(ax, ax, indexing='xy'))
    if h.ndim == 1:
        return ax, np.interp(rr, r_centers_deg, h, left=h[0], right=np.nan)
    out = np.empty(rr.shape + (h.shape[1],))
    for k in range(h.shape[1]):
        out[..., k] = np.interp(rr, r_centers_deg, h[:, k],
                                left=h[0, k], right=np.nan)
    return ax, out


def radial_from_projections(profiles, strip_centers, strip_width_deg,
                            r_edges_deg=None, ridge=1e-2, fit_gain=True,
                            n_gain_iter=6):
    """One radial profile h(r) constrained by SEVERAL families at once.

    Under radial symmetry every family measures the same h, so they are not
    separate reconstructions to compare afterwards -- they are twice the
    constraints on one answer. Stacking them is also what makes the assumption
    testable: an elongated field projects differently along different axes, and
    no single h can satisfy both, so the residual rises. With one family the
    residual says nothing (see radial_from_projection).

    profiles      : list of (n_bins,) or (n_bins, n_lags) registered profiles
    strip_centers : list of (n_bins,) registered grids, one per profile
    fit_gain      : fit one amplitude per family, shared across lags, so that an
                    arbitrary difference in mean amplitude between two roi
                    populations is not read as a failure of radial symmetry.
                    Reported as 'gain', relative to the first family.
    Returns dict with 'r_centers_deg', 'h', 'residual' (over everything),
    'residual_per_family' and 'gain'.
    """
    profs = [np.atleast_2d(np.asarray(p, dtype=float).T).T for p in profiles]
    n_lags = profs[0].shape[1]
    assert all(p.shape[1] == n_lags for p in profs), 'lag axes must match'
    grids = [np.asarray(c, dtype=float) for c in strip_centers]

    if r_edges_deg is None:
        allc = np.concatenate(grids)
        step = float(np.median(np.diff(np.sort(np.unique(allc)))))
        r_max = float(np.nanmax(np.abs(allc))) + step / 2.0
        r_edges_deg = np.arange(0.0, r_max + step, step)
    r_edges_deg = np.asarray(r_edges_deg, dtype=float)

    mats = [projection_matrix(g, strip_width_deg, r_edges_deg) for g in grids]
    n_r = mats[0].shape[1]
    D2 = np.diff(np.eye(n_r), 2, axis=0) if n_r > 2 else np.eye(n_r)
    DtD = D2.T @ D2

    # p_i ~= gain_i * (A_i h), bilinear in (h, gain), so alternate: each half
    # is a least squares with the other held fixed. Initialised at gain = 1,
    # which is the old behaviour, and converged in a handful of passes.
    gain = np.ones(len(profs))
    h = np.full((n_r, n_lags), np.nan)
    for _ in range(int(n_gain_iter) if fit_gain else 1):
        h = np.full((n_r, n_lags), np.nan)
        for k in range(n_lags):
            rows, ys = [], []
            for A, p, g in zip(mats, profs, gain):
                ok = np.isfinite(p[:, k])
                if ok.sum() < 3:
                    continue
                rows.append(g * A[ok])
                ys.append(p[ok, k])
            if not rows:
                continue
            As = np.vstack(rows)
            y = np.concatenate(ys)
            AtA = As.T @ As
            lam = float(ridge) * np.trace(AtA) / max(np.trace(DtD), 1e-12)
            h[:, k] = np.linalg.solve(AtA + lam * DtD, As.T @ y)
        if not fit_gain:
            break
        # one scalar per family, pooled over every lag it has data for
        new = np.ones(len(profs))
        for i, (A, p) in enumerate(zip(mats, profs)):
            num = den = 0.0
            for k in range(n_lags):
                if not np.all(np.isfinite(h[:, k])):
                    continue
                ok = np.isfinite(p[:, k])
                if ok.sum() < 3:
                    continue
                pred = A[ok] @ h[:, k]
                num += float(pred @ p[ok, k])
                den += float(pred @ pred)
            # clamped at zero: a gain rescales a family, it never inverts it
            new[i] = max(num / den, 0.0) if den > 0 else 1.0
        if np.allclose(new, gain, rtol=1e-6, atol=1e-12):
            gain = new
            break
        gain = new

    # h and gain trade scale freely; report gains relative to the first family
    # so they read as 'family i is this much louder than family 0'. If family 0
    # is the one that was clamped out, normalise by the loudest instead, so the
    # numbers stay readable rather than dividing by zero.
    _ref = gain[0] if (np.isfinite(gain[0]) and gain[0] > 0) else np.nanmax(gain)
    if np.isfinite(_ref) and _ref > 0:
        h = h * _ref
        gain = gain / _ref

    resid = np.full(n_lags, np.nan)
    per_fam = np.full((len(profs), n_lags), np.nan)
    for k in range(n_lags):
        if not np.all(np.isfinite(h[:, k])):
            continue
        num_sq = den_sq = 0.0
        for i, (A, p, g) in enumerate(zip(mats, profs, gain)):
            ok = np.isfinite(p[:, k])
            if ok.sum() < 3:
                continue
            yi = p[ok, k]
            pi = g * (A[ok] @ h[:, k])
            di = float(np.linalg.norm(yi))
            per_fam[i, k] = (float(np.linalg.norm(pi - yi) / di)
                             if di > 0 else np.nan)
            num_sq += float(np.sum((pi - yi) ** 2))
            den_sq += float(np.sum(yi ** 2))
        resid[k] = np.sqrt(num_sq / den_sq) if den_sq > 0 else np.nan

    r_centers = 0.5 * (r_edges_deg[:-1] + r_edges_deg[1:])
    single = np.asarray(profiles[0]).ndim == 1
    if single:
        h, resid, per_fam = h[:, 0], float(resid[0]), per_fam[:, 0]
    return dict(r_centers_deg=r_centers, h=h, residual=resid,
                residual_per_family=per_fam, gain=gain)
