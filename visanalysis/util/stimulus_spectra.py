"""Spatial and temporal frequency analysis of visual stimuli.

Written for Zebra noise but deliberately stimulus-agnostic: everything here
takes a movie array plus a geometry description, so it applies equally to
NoizOne/NoiseTwo movies, drifting gratings (useful as a known-answer test),
or anything else.

Reproduces the characterisations in Figure 1f-1g of

    Skriabine, S., Shinn, M., Picard, S., Harris, K. D., & Carandini, M.
    (2026). Mapping the visual cortex with Zebra noise and wavelets.
    Journal of Vision, 26(1):1, 1-16.

which computes the power spectral density of stimulus frames by FFT,
converts to polar coordinates to separate radial (spatial frequency) from
angular (orientation) content, and reports mean and standard deviation of
power at each. Figure 1e's spatial autocorrelation is also implemented.


Generalising to a sphere
------------------------
The intended future case is noise generated as 4D Perlin (x, y, z, t) and
sampled on the sphere's surface, which avoids the seams and pole distortion
of warping a 2D field.

Only ONE thing changes in that case: how spatial structure is decomposed.

    planar     2D FFT           -> coefficients indexed by (fx, fy)
    spherical  harmonic transform -> coefficients indexed by (l, m)

Everything downstream is shared, because both produce "power at a scalar
spatial frequency". For the plane that is sqrt(fx^2 + fy^2) in cycles/degree;
for the sphere, a harmonic of degree l oscillates l times around a great
circle, i.e. l cycles per 360 degrees, so the comparable quantity is

    cycles_per_degree = l / 360

That makes planar and spherical spatial spectra directly comparable on the
same axis, which is the whole point of the abstraction.

The seam is `FrequencyCoords`: a spatial transform supplies the power array
plus per-coefficient radial and orientation coordinates, and the binning and
temporal machinery below never needs to know which geometry produced them.

Two honest caveats for the spherical case, recorded now so they are not
discovered late:

1. Partial coverage. The rig's screen is a cap (max phi ~76 deg), not the
   full sphere. Spherical harmonics are orthogonal only over the FULL
   sphere, so decomposing partial data leaks power between degrees -- the
   same problem as the partial-sky/mode-coupling issue in CMB analysis.
   Options are to apodise and accept broadened resolution in l, to use a
   local tangent-plane projection over the cap, or to compute the true
   coupling matrix and deconvolve. Do not silently run a full-sphere SHT on
   cap data and report C_l.

2. Orientation on a sphere is not a single global quantity the way it is on
   a plane; m indexes orientation relative to the chosen pole, so a naive
   "power vs m" is pole-dependent, not an orientation tuning curve. Local
   analysis in a tangent plane is the meaningful analogue.
"""

from dataclasses import dataclass

import warnings

import numpy as np


# --------------------------------------------------------------------------
# geometry / the pluggable seam
# --------------------------------------------------------------------------

@dataclass
class PlanarGrid:
    """Sampling geometry of a flat-screen stimulus.

    deg_per_px_x, deg_per_px_y : visual angle subtended by one pixel.
    These are usually NOT equal on a rig with a portrait-mounted projector,
    and getting them wrong rescales every frequency axis, so they are
    required rather than defaulted.
    """
    deg_per_px_x: float
    deg_per_px_y: float


@dataclass
class FrequencyCoords:
    """Per-coefficient coordinates for a spatial decomposition.

    This is what a spatial transform must provide. A spherical transform
    would populate the same fields (radial = l / 360, orientation = None).

    radial : ndarray
        Spatial frequency of each coefficient, cycles/degree.
    orientation : ndarray or None
        Orientation of each coefficient in [0, 180) degrees, or None if the
        geometry does not define a global orientation.
    shape : tuple
        Shape of the coefficient array these coordinates describe.
    nyquist : float
        Highest frequency at which the coefficient set covers ALL orientations.
        For a 2D FFT this is the smaller of the two axis Nyquists, NOT
        max(radial): the frequency array is a square, so radial reaches
        Nyquist*sqrt(2) at the corners, and every bin above Nyquist is sampled
        only near 45/135 deg. Averaging over such a bin mixes an orientation
        subset with a radial offset, which on any non-flat spectrum reads as
        anisotropy that is not in the stimulus. A spherical transform would set
        this to its own max fully-sampled degree / 360.
    """
    radial: np.ndarray
    orientation: np.ndarray
    shape: tuple
    nyquist: float


def planar_frequency_coords(shape, grid):
    """Frequency coordinates for a 2D FFT on a regular pixel grid.

    Uses fftfreq with the pixel size as the sample spacing, so the returned
    frequencies are already in cycles/degree rather than cycles/pixel.
    """
    nx, ny = shape
    fx = np.fft.fftfreq(nx, d=grid.deg_per_px_x)
    fy = np.fft.fftfreq(ny, d=grid.deg_per_px_y)
    FX, FY = np.meshgrid(fx, fy, indexing="ij")

    radial = np.sqrt(FX ** 2 + FY ** 2)
    # Orientation of the sinusoid's modulation direction, folded to [0, 180)
    # since a grating at theta and theta+180 are the same grating.
    orientation = np.degrees(np.arctan2(FY, FX)) % 180.0
    nyquist = min(1.0 / (2.0 * grid.deg_per_px_x), 1.0 / (2.0 * grid.deg_per_px_y))
    return FrequencyCoords(radial=radial, orientation=orientation,
                           shape=(nx, ny), nyquist=nyquist)


# --------------------------------------------------------------------------
# windowing and periodograms
# --------------------------------------------------------------------------

def _hann2d(nx, ny):
    """Separable 2D Hann window.

    A stimulus frame is not periodic across its edges, so an unwindowed FFT
    smears power from the discontinuity across all frequencies (spectral
    leakage) -- which for a 1/f-like spectrum masquerades as excess
    high-frequency power. The paper does not mention windowing; it matters.
    """
    wx = np.hanning(nx).astype(np.float32)
    wy = np.hanning(ny).astype(np.float32)
    return np.outer(wx, wy)


def _prepare(frames):
    """(X, Y, T) of any dtype -> float32, guarding against integer input."""
    a = np.asarray(frames)
    if a.ndim == 2:
        a = a[:, :, None]
    return a.astype(np.float32, copy=False)


# --------------------------------------------------------------------------
# invalid (unlit / padded) regions
#
# A rectangular array containing a non-rectangular stimulus -- a spherical
# rig's projector cap, a circular aperture, a letterboxed movie -- has a hard
# mask edge running through it. An FFT convolves the true spectrum with the
# mask's transform, and the mask usually wins. Measured on an ISOTROPIC 1/f
# field zero-padded outside a barrel-shaped lit region: radial slope went from
# -2.04 to -2.50 and apparent anisotropy from 0.18 to 3.25, out of a field with
# none. Nothing downstream can detect this from the numbers, so detect it here.
# --------------------------------------------------------------------------

def valid_mask(frames, n_probe=32):
    """Pixels that carry real stimulus, as a 2D bool array.

    A padded / unlit pixel holds the same value in every frame. A live one
    almost surely does not -- but "almost surely" is weaker for a BINARY
    stimulus (zebra noise is 0/255) than for a continuous one, and weaker
    still if the probes are not independent. Both failure modes were measured
    on live zebra noise where every pixel was valid:

        probes  probe window   tau    falsely flagged invalid
            16        10.0 s  0.25 s              0.0000%
            16         1.0 s  0.25 s              0.90%
            16         0.5 s  0.25 s             11.28%
            16        0.25 s  0.25 s             25.42%
             8        10.0 s  0.25 s              0.67%
            32        10.0 s  0.25 s              0.0000%

    With independent probes a binary pixel is constant with probability
    2^-(n-1): 0.78% at n=8, 0.003% at n=16, 5e-10 at n=32. The n=8 row above
    matches that prediction, so the default is 32 -- the extra probes cost
    almost nothing and buy five orders of magnitude.

    The second hazard is not fixed by more probes. When the probe window is
    comparable to the stimulus's correlation time the probes are not
    independent, and the rate explodes regardless of n. Analysing a 0.3-0.5 s
    epoch of a stimulus with tau = 0.25 s is squarely in that regime.
    `crop_to_valid` therefore checks the GEOMETRY of what it found rather
    than trusting the count.
    """
    a = _prepare(frames)
    nt = a.shape[2]
    idx = np.unique(np.linspace(0, nt - 1, min(n_probe, nt)).round().astype(int))
    sample = a[:, :, idx]
    return ~np.all(sample == sample[:, :, :1], axis=2)


def split_half_agreement(frames, n_probe=32):
    """Jaccard agreement between invalid masks from two disjoint frame halves.

    The discriminator between real padding and chance-constant pixels.
    An unlit pixel is constant in EVERY subset of frames, so the two masks
    coincide and this returns ~1. A pixel that merely happened to hold still
    is unlikely to do so in both halves, so the masks disagree and this drops.

    Note what does NOT work here. Geometric contiguity was the obvious idea --
    padding is a solid region, chance ought to be scattered pixels -- but
    measured on live zebra noise it fails outright:

        case                                    contiguity
        real barrel-shaped unlit region              0.94
        false positives, 0.25 s window               0.79
        false positives, 0.50 s window               0.84

    Zebra noise is spatially smooth, so over a short window it is whole
    STRIPES that hold still, and a stripe is just as contiguous as a padded
    region. The split-half test does not care about shape, only about whether
    the same pixels hold still in independent samples, which is exactly the
    property that separates the two.

    Returns (agreement, invalid_both) where invalid_both is the intersection,
    i.e. the conservative padding estimate.
    """
    a = _prepare(frames)
    nt = a.shape[2]
    if nt < 4:
        return float("nan"), np.zeros(a.shape[:2], dtype=bool)
    mid = nt // 2
    inv1 = ~valid_mask(a[:, :, :mid], n_probe=n_probe)
    inv2 = ~valid_mask(a[:, :, mid:], n_probe=n_probe)
    both = inv1 & inv2
    either = inv1 | inv2
    n_either = int(either.sum())
    if n_either == 0:
        return 1.0, both
    return float(both.sum()) / float(n_either), both


def largest_valid_rect(mask):
    """Largest all-valid axis-aligned rectangle. Returns (r0, r1, c0, c1) slices.

    Classic maximal-rectangle-in-a-binary-matrix, O(rows*cols) via the
    largest-rectangle-in-histogram stack. Cropping to this is the cheap exact
    fix: no mask edge survives, so no leakage to correct for, at the cost of
    the data outside the rectangle.
    """
    mask = np.asarray(mask, dtype=bool)
    nrows, ncols = mask.shape
    heights = np.zeros(ncols, dtype=np.int64)
    best = (0, 0, 0, 0, 0)                       # area, r0, r1, c0, c1
    for r in range(nrows):
        heights = np.where(mask[r], heights + 1, 0)
        stack = []                               # (start_col, height)
        for c in range(ncols + 1):
            h = heights[c] if c < ncols else 0
            start = c
            while stack and stack[-1][1] > h:
                s, hh = stack.pop()
                area = hh * (c - s)
                if area > best[0]:
                    best = (area, r - hh + 1, r, s, c - 1)
                start = s
            stack.append((start, h))
    _, r0, r1, c0, c1 = best
    return r0, r1 + 1, c0, c1 + 1


def crop_to_valid(frames, min_fraction=0.999, min_agreement=0.9,
                  verbose=True):
    """Detect a padded/unlit region and crop to the largest clean rectangle.

    Returns (frames_cropped, info). If the array is already fully valid the
    frames come back untouched and info['cropped'] is False, so this is safe to
    call unconditionally.

    For a SPHERICAL stimulus, prefer extracting your own tangent-plane patch:
    this crops in array coordinates, and for a (phi, theta) movie those are
    polar, not a uniform angular grid -- a rectangle there is not a rectangle
    on the sphere, and the sampling density still varies across it. Use this
    for flat-screen movies and as a safety net.
    """
    a = _prepare(frames)
    m = valid_mask(a)
    frac = float(m.mean())
    info = {"valid_fraction": frac, "cropped": False, "slices": None,
            "retained_fraction": 1.0, "agreement": None}
    if frac >= min_fraction:
        return a, info

    # Is what we found actually a padded region, or pixels that merely held
    # still? Cropping on the latter is destructive and silent -- in the worst
    # measured case it retained 2.7% of the array and moved the spectral slope
    # from -1.815 to -1.744 while printing a routine-looking message.
    agreement, both = split_half_agreement(a)
    info["agreement"] = agreement
    if not np.isfinite(agreement) or agreement < min_agreement:
        if verbose:
            print(f"  [mask] {100*(1-frac):.1f}% of pixels are constant across "
                  f"frames, but two disjoint halves of the clip disagree about "
                  f"WHICH pixels (agreement {agreement:.2f} < {min_agreement}). "
                  f"Real padding is constant in every subset; this is the "
                  f"signature of a binary stimulus sampled over too short a "
                  f"window relative to its correlation time -- NOT cropping. "
                  f"Pass a longer clip, or an explicit mask if you know the "
                  f"lit region.")
        return a, info

    # Use the intersection of the two halves: a pixel must have held still in
    # both to count as unlit. Strictly more conservative than the full-clip
    # mask, and it costs nothing since it is already computed.
    m = ~both
    r0, r1, c0, c1 = largest_valid_rect(m)
    out = a[r0:r1, c0:c1, :]
    info.update(cropped=True, slices=(r0, r1, c0, c1),
                retained_fraction=out[:, :, 0].size / a[:, :, 0].size)
    if verbose:
        print(f"  [mask] {100*(1-frac):.1f}% of the array is constant across frames "
              f"(padded/unlit, split-half agreement {agreement:.2f}). Cropped to the largest "
              f"clean rectangle "
              f"{out.shape[0]}x{out.shape[1]} = {100*info['retained_fraction']:.0f}% "
              f"of the array, {100*out[:,:,0].size/max(m.sum(),1):.0f}% of the valid pixels.")
    return out, info


def frame_power_spectrum(frames, window=True, remove_mean=True, grid=None,
                         frame_stride=1):
    """Mean 2D periodogram over frames.

    Averaging per-frame 2D spectra (rather than taking one 3D FFT) keeps
    memory flat in the number of frames, which is what makes it practical
    to average over the thousands of frames needed for a stable estimate.
    The paper averages over 1500 frames.

    grid : PlanarGrid or None
        If given, the result is scaled by the pixel area to a true power
        spectral DENSITY, in signal^2 per (cycles/degree)^2, which is what
        makes levels comparable across different pixel sizes. Left as None the
        result is power per coefficient -- a constant factor away, so spectral
        slopes and anisotropy are identical either way.

    frame_stride : int
        Use every Nth frame. Consecutive frames of a slowly-evolving stimulus
        are highly correlated, so averaging all nt of them buys far fewer than
        nt independent samples; the variance of the estimate is set by the
        independent count, not the raw one. Set this to at least
        tau_correlation * fps to make the averaged frames near-independent.
        See effective_frame_count().

    Returns power in fft (unshifted) layout, matching planar_frequency_coords.
    """
    a = _prepare(frames)[:, :, ::max(1, int(frame_stride))]
    nx, ny, nt = a.shape
    w = _hann2d(nx, ny) if window else None
    # Correct for the power the window removes, so absolute levels stay
    # comparable between windowed and unwindowed estimates.
    wcorr = np.mean(w ** 2) if w is not None else 1.0

    acc = np.zeros((nx, ny), dtype=np.float64)
    for t in range(nt):
        f = a[:, :, t]
        if remove_mean:
            f = f - f.mean()
        if w is not None:
            f = f * w
        F = np.fft.fft2(f)
        acc += (F.real ** 2 + F.imag ** 2)
    acc /= (nt * nx * ny * wcorr)
    if grid is not None:
        acc *= grid.deg_per_px_x * grid.deg_per_px_y   # -> PSD per (c/deg)^2
    return acc


def effective_frame_count(frames, fps):
    """(n_independent, tau_seconds) for a frame stack.

    Estimates the temporal correlation time as the 1/e point of the pixel-
    averaged temporal autocorrelation, then reports how many effectively
    independent frames the stack contains. Statistics whose variance you care
    about -- notably the std/mean anisotropy index -- are set by THIS number,
    not by the raw frame count. A single realisation's angular spectrum is
    chi-squared speckle with a couple of degrees of freedom per bin, which is
    easily mistaken for real structure.
    """
    lags, ac = temporal_autocorrelation(frames, fps)
    below = np.flatnonzero(ac < np.exp(-1.0))
    tau = float(lags[below[0]]) if below.size else float(lags[-1])
    nt = np.asarray(frames).shape[2]
    n_eff = nt / max(tau * fps, 1.0)
    return n_eff, tau


# --------------------------------------------------------------------------
# generic binning -- geometry agnostic, works for planar or spherical coords
# --------------------------------------------------------------------------

def bin_average(power, coord, bins, mask=None):
    """Average `power` into bins of `coord`, returning mean, std and count.

    Both the radial (spatial frequency) and angular (orientation) reductions
    are this same operation on different coordinates, and a spherical
    transform reuses it unchanged.
    """
    p = np.asarray(power).ravel()
    c = np.asarray(coord).ravel()
    if mask is not None:
        m = np.asarray(mask).ravel()
        p, c = p[m], c[m]

    idx = np.digitize(c, bins) - 1
    n = len(bins) - 1
    valid = (idx >= 0) & (idx < n)
    idx, p = idx[valid], p[valid]

    count = np.bincount(idx, minlength=n).astype(float)
    total = np.bincount(idx, weights=p, minlength=n)
    total_sq = np.bincount(idx, weights=p ** 2, minlength=n)

    with np.errstate(invalid="ignore", divide="ignore"):
        mean = total / count
        var = total_sq / count - mean ** 2
        std = np.sqrt(np.maximum(var, 0.0))
    mean[count == 0] = np.nan
    std[count == 0] = np.nan

    centres = 0.5 * (bins[:-1] + bins[1:])
    return centres, mean, std, count


def radial_spectrum(power, coords, n_bins=48, log_bins=True, f_min=None,
                    f_max=None):
    """Power vs spatial frequency (cycles/degree). Figure 1f.

    Log-spaced bins by default: a 1/f-like spectrum spans decades, and
    linear bins would put almost every sample in the highest-frequency bins.
    The DC term is excluded -- it carries the mean luminance, not structure.

    f_max defaults to coords.nyquist, NOT to max(radial). The frequency array
    is a square, so radial runs out to Nyquist*sqrt(2) in the corners; bins
    above Nyquist are populated only by those corners, i.e. only near 45/135
    deg. Including them mixes an orientation subset into the radial average and
    leaves the top bins with a fraction of the samples (measured: 10066 -> 1440
    coefficients per bin across the Nyquist boundary on a 256^2 grid). Pass
    f_max explicitly if you really want the corner region.

    This matters more here than for orientation_spectrum: it is the radial
    axis whose top bins become both undersampled and orientation-restricted.
    """
    r = coords.radial
    positive = r > 0
    if f_min is None:
        f_min = np.min(r[positive])
    if f_max is None:
        f_max = coords.nyquist

    if log_bins:
        bins = np.logspace(np.log10(f_min), np.log10(f_max), n_bins + 1)
    else:
        bins = np.linspace(f_min, f_max, n_bins + 1)

    return bin_average(power, r, bins, mask=positive)


def orientation_spectrum(power, coords, n_bins=36, f_range=None):
    """Power vs orientation in [0, 180). Figure 1g.

    f_range restricts the analysis to an annulus of spatial frequency. When
    None it defaults to (4 * grid frequency step, coords.nyquist), because BOTH
    ends of the full square are hostile -- but very unequally so. Measured on an
    ISOTROPIC 1/f field, where the true anisotropy is 0:

        no limits at all                    0.687
        corners excluded only (r <= Nyq)    0.687   <- unchanged
        low frequencies excluded only       0.171
        both excluded                       0.171

    So essentially ALL of the ~4x inflation comes from the handful of
    coefficients nearest DC, not from the corners. On a 1/f-like spectrum those
    carry overwhelmingly the most power while being sparse in angle, so they
    dominate the per-bin mean and scatter it. The corner coefficients sit where
    power is already tiny and barely move the average -- excluding them matters
    for the RADIAL spectrum (where their bins are undersampled and reach only
    45/135 deg; see radial_spectrum) far more than for this one.

    The Nyquist cap is applied even if f_range asks for more, since a bin that
    samples only two orientations cannot inform an orientation average.
    """
    if coords.orientation is None:
        raise ValueError(
            "This geometry does not define a global orientation; use a local "
            "tangent-plane analysis instead (see module docstring)."
        )
    r = coords.radial
    if f_range is None:
        f_range = (4.0 * np.min(r[r > 0]), coords.nyquist)
    mask = (r > 0) & (r <= coords.nyquist)
    mask &= (r >= f_range[0]) & (r <= f_range[1])

    bins = np.linspace(0.0, 180.0, n_bins + 1)
    return bin_average(power, coords.orientation, bins, mask=mask)


# --------------------------------------------------------------------------
# temporal and joint spatiotemporal
# --------------------------------------------------------------------------

def temporal_spectrum(frames, fps, nperseg=256, noverlap=None,
                      remove_mean=True):
    """Power vs temporal frequency (Hz), averaged over pixels.

    Welch's method: split into overlapping segments, window each, average the
    periodograms. Trades frequency resolution for a lower-variance estimate,
    which is the right trade for noise.
    """
    a = _prepare(frames)
    nx, ny, nt = a.shape
    if noverlap is None:
        noverlap = nperseg // 2
    step = nperseg - noverlap
    if nt < nperseg:
        raise ValueError(f"need at least {nperseg} frames, got {nt}")

    w = np.hanning(nperseg).astype(np.float32)
    wcorr = np.mean(w ** 2)

    starts = range(0, nt - nperseg + 1, step)
    acc = np.zeros(nperseg // 2 + 1, dtype=np.float64)
    n_seg = 0
    for s in starts:
        seg = a[:, :, s:s + nperseg]
        if remove_mean:
            seg = seg - seg.mean(axis=2, keepdims=True)
        seg = seg * w
        F = np.fft.rfft(seg, axis=2)
        acc += np.mean(F.real ** 2 + F.imag ** 2, axis=(0, 1))
        n_seg += 1
    acc /= (n_seg * nperseg * wcorr)

    freqs = np.fft.rfftfreq(nperseg, d=1.0 / fps)
    return freqs, acc


def spatiotemporal_spectrum(frames, grid, fps, nperseg=64, spatial_stride=1,
                            remove_mean=True):
    """Joint power over (spatial frequency, temporal frequency).

    For a motion-sensitive system this is more informative than the two
    marginal spectra: a rigidly drifting pattern puts power on the line
    ft = speed * fs, so the joint distribution is what encodes the range of
    speeds present. Relevant for fly T4/T5, which are tuned to speed rather
    than to spatial or temporal frequency alone.

    spatial_stride subsamples pixels before the 3D transform; a full
    640x1280x64 complex FFT is several GB, and the joint spectrum rarely
    needs full spatial resolution. Striding lowers the spatial Nyquist
    accordingly.
    """
    a = _prepare(frames)[::spatial_stride, ::spatial_stride, :]
    nx, ny, nt = a.shape
    if nt < nperseg:
        raise ValueError(f"need at least {nperseg} frames, got {nt}")

    eff = PlanarGrid(grid.deg_per_px_x * spatial_stride,
                     grid.deg_per_px_y * spatial_stride)
    coords = planar_frequency_coords((nx, ny), eff)

    wt = np.hanning(nperseg).astype(np.float32)
    ws = _hann2d(nx, ny)
    wcorr = np.mean(ws ** 2) * np.mean(wt ** 2)

    n_ft = nperseg // 2 + 1
    acc = np.zeros((nx, ny, n_ft), dtype=np.float64)
    n_seg = 0
    for s in range(0, nt - nperseg + 1, nperseg // 2):
        seg = a[:, :, s:s + nperseg]
        if remove_mean:
            seg = seg - seg.mean()
        seg = seg * ws[:, :, None] * wt[None, None, :]
        F = np.fft.rfftn(seg, axes=(0, 1, 2))
        acc += (F.real ** 2 + F.imag ** 2)
        n_seg += 1
    acc /= (n_seg * nx * ny * nperseg * wcorr)

    ft = np.fft.rfftfreq(nperseg, d=1.0 / fps)
    return coords, ft, acc


def speed_spectrum(coords, ft, power, n_bins=32, speed_range=(1.0, 1000.0),
                   f_min=None):
    """Power vs speed (degrees/second), from the joint spectrum.

    Each coefficient at spatial frequency fs and temporal frequency ft
    corresponds to a drifting component of speed ft / fs. Coefficients near
    fs = 0 are excluded: their speed is undefined (a spatially uniform
    flicker has no direction), and dividing by a near-zero fs would scatter
    enormous spurious speeds through the distribution.
    """
    r = coords.radial[:, :, None]
    ftg = ft[None, None, :]

    if f_min is None:
        f_min = np.min(r[r > 0])
    valid = (r > f_min) & (ftg > 0)

    with np.errstate(divide="ignore", invalid="ignore"):
        speed = ftg / r
    speed = np.broadcast_to(speed, power.shape)

    bins = np.logspace(np.log10(speed_range[0]), np.log10(speed_range[1]),
                       n_bins + 1)
    return bin_average(power, speed, bins,
                       mask=np.broadcast_to(valid, power.shape))


# --------------------------------------------------------------------------
# autocorrelation
# --------------------------------------------------------------------------

def spatial_autocorrelation(frames, grid, max_lag_deg=None, remove_mean=True):
    """2D spatial autocorrelation averaged over frames. Figure 1e.

    Computed via the Wiener-Khinchin theorem (inverse FFT of the power
    spectrum) rather than by direct pairwise correlation, which is far
    cheaper and equivalent. Normalised so zero lag is 1.

    Returns (lag_x_deg, lag_y_deg, correlation) with zero lag at the centre.
    """
    a = _prepare(frames)
    nx, ny, nt = a.shape

    acc = np.zeros((nx, ny), dtype=np.float64)
    for t in range(nt):
        f = a[:, :, t]
        if remove_mean:
            f = f - f.mean()
        F = np.fft.fft2(f)
        acc += np.fft.ifft2(F * np.conj(F)).real
    acc /= nt
    acc = np.fft.fftshift(acc)
    if acc.max() > 0:
        acc /= acc.max()

    lag_x = (np.arange(nx) - nx // 2) * grid.deg_per_px_x
    lag_y = (np.arange(ny) - ny // 2) * grid.deg_per_px_y

    if max_lag_deg is not None:
        kx = np.abs(lag_x) <= max_lag_deg
        ky = np.abs(lag_y) <= max_lag_deg
        acc = acc[np.ix_(kx, ky)]
        lag_x, lag_y = lag_x[kx], lag_y[ky]

    return lag_x, lag_y, acc


def temporal_autocorrelation(frames, fps, max_lag_s=None, remove_mean=True,
                             max_pixels=4096, seed=0):
    """Pixelwise temporal autocorrelation, averaged over pixels.

    Gives the stimulus's correlation time, which sets how quickly the
    pattern refreshes -- directly relevant to choosing tscale.

    max_pixels : int
        Estimate from at most this many randomly chosen pixels. Averaging over
        every pixel is not just wasteful but a memory hazard: np.fft upcasts to
        complex128 whatever it is given, so a full 640x1280x512 stack costs
        ~1.7 GB for the reshape and ~6.7 GB EACH for the rfft and irfft
        buffers. That is ~15 GB to estimate one scalar. A few thousand pixels
        give the same tau to well within its own uncertainty, since neighbouring
        pixels of a smooth stimulus are highly redundant anyway. Pass None to
        use all pixels.
    """
    a = _prepare(frames)
    nx, ny, nt = a.shape
    x = a.reshape(-1, nt)
    if max_pixels is not None and x.shape[0] > max_pixels:
        rng = np.random.default_rng(seed)
        x = x[rng.choice(x.shape[0], max_pixels, replace=False)]
    x = x.astype(np.float32, copy=False)
    if remove_mean:
        x = x - x.mean(axis=1, keepdims=True)

    n_pad = 1 << int(np.ceil(np.log2(2 * nt)))
    F = np.fft.rfft(x, n=n_pad, axis=1)
    ac = np.fft.irfft(F * np.conj(F), n=n_pad, axis=1)[:, :nt]
    # Unbiased: each lag k is averaged over only nt-k overlapping products.
    ac /= np.arange(nt, 0, -1)[None, :]

    with np.errstate(invalid="ignore", divide="ignore"):
        ac = ac / ac[:, :1]
    # A pixel that never changes has zero variance, so its normalised autocorrelation is 0/0 and
    # the whole row is NaN. That is legitimate -- a coarse, slow stimulus really does hold some
    # pixels constant over a short record -- and dropping those rows is the right answer, which
    # nanmean already gives. The warning is silenced here rather than left to fire on every run
    # of a perfectly valid measurement.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Mean of empty slice",
                                category=RuntimeWarning)
        ac = np.nanmean(ac, axis=0)

    lags = np.arange(nt) / fps
    if max_lag_s is not None:
        keep = lags <= max_lag_s
        lags, ac = lags[keep], ac[keep]
    return lags, ac


# --------------------------------------------------------------------------
# convenience
# --------------------------------------------------------------------------

def summarize(frames, grid, fps, f_range=None, verbose=True, frame_stride=1,
              handle_mask=False):
    """Run the standard battery and return a dict of results.

    handle_mask : bool, default False
        Detect a padded / unlit region and crop to the largest clean
        rectangle. OPT-IN, because for stimuli generated by this project
        there is no padding to find and the detector can only do harm:

            tangent-plane patch   none by construction -- the reason that is
                                  the recommended method for a sphere
            4D polar texture      none; every texel within max_phi is a real
                                  viewing direction
            flat-screen movie     none; the whole frame is lit
            projector-space frame padded, but nobody analyses in projector
                                  coordinates

        Turn it on for FOREIGN data whose sampling you did not control --
        legacy flymax movies assign zero to offscreen texels, letterboxed
        recordings, circular apertures. There the mask's hard edge really
        does dominate the transform: measured on an isotropic 1/f field
        padded outside a barrel-shaped lit region, the mask alone steepened
        the radial slope from -2.04 to -2.50 and inflated apparent anisotropy
        from 0.18 to 3.25, out of a field with none.

        The cost of leaving it on by default is not zero. On a BINARY
        stimulus the detector can misfire (see valid_mask); in the worst
        measured case it cropped to 2.7% of the array and shifted the slope
        from -1.815 to -1.744 while printing an ordinary-looking message.
        crop_to_valid now refuses when its split-half check fails, but the
        cleanest guard is not to invoke it on data that cannot be padded.

    For a SPHERICAL stimulus the crop is never the right answer even when
    enabled: it crops in array coordinates, and a (phi, theta) movie's array
    coordinates are polar, so a rectangle there is not a rectangle on the
    sphere. Extract a tangent-plane patch inside the lit region and pass that.

    f_range defaults to an annulus from 4 grid frequencies up to Nyquist,
    rather than the whole square -- see orientation_spectrum, which now applies
    that default itself.
    """
    if handle_mask:
        frames, mask_info = crop_to_valid(frames, verbose=verbose)
    else:
        mask_info = {"valid_fraction": None, "cropped": False,
                     "checked": False}
    power = frame_power_spectrum(frames, grid=grid, frame_stride=frame_stride)
    coords = planar_frequency_coords(power.shape, grid)

    if f_range is None:
        f_range = (4.0 * np.min(coords.radial[coords.radial > 0]), coords.nyquist)

    f, p_f, s_f, _ = radial_spectrum(power, coords)
    o, p_o, s_o, _ = orientation_spectrum(power, coords, f_range=f_range)

    n_eff, tau = effective_frame_count(frames, fps)

    out = {
        "power2d": power,
        "coords": coords,
        "spatial_freq": f, "spatial_power": p_f, "spatial_std": s_f,
        "orientation": o, "orientation_power": p_o,
        "orientation_std": s_o,
        "f_range": f_range,
        "n_effective_frames": n_eff, "correlation_time_s": tau,
        "mask": mask_info,
    }

    if frames.shape[2] >= 64:
        nper = min(256, 1 << int(np.floor(np.log2(frames.shape[2]))))
        t_f, t_p = temporal_spectrum(frames, fps, nperseg=nper)
        out["temporal_freq"] = t_f
        out["temporal_power"] = t_p

    if verbose:
        valid = np.isfinite(p_f)
        # Slope of log power vs log frequency: -2 would be classic 1/f^2
        # (pink in 2D), 0 would be white.
        slope = np.polyfit(np.log10(f[valid]), np.log10(p_f[valid]), 1)[0]
        aniso = np.nanstd(p_o) / np.nanmean(p_o)
        # chi-squared speckle floor: with n_eff independent samples an
        # ISOTROPIC field still measures an anisotropy of roughly 1/sqrt(n_eff).
        floor = 1.0 / np.sqrt(max(n_eff, 1.0))
        print(f"  spatial slope (log-log)   : {slope:+.2f}")
        print(f"  analysis band (c/deg)     : {f_range[0]:.3f} - {f_range[1]:.3f} "
              f"(Nyquist {coords.nyquist:.3f})")
        print(f"  effective independent frames: {n_eff:.0f} of "
              f"{np.asarray(frames).shape[2]} (tau = {tau*1e3:.0f} ms)")
        print(f"  orientation anisotropy    : {aniso:.3f} "
              f"(0 = isotropic; noise floor ~{floor:.3f})"
              + ("   <- AT THE NOISE FLOOR, not measurable"
                 if aniso < 2 * floor else ""))
        if "temporal_power" in out:
            tv = np.isfinite(t_p) & (t_f > 0)
            tslope = np.polyfit(np.log10(t_f[tv]),
                                np.log10(t_p[tv]), 1)[0]
            print(f"  temporal slope (log-log)  : {tslope:+.2f}")
        out["spatial_slope"] = slope
        out["anisotropy"] = aniso

    return out
