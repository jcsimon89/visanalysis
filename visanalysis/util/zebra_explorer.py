"""Offline parameter exploration for the Zebra noise stimulus.

Set the parameters you would set in the protocol, render the stimulus the display would render,
and look at what its spatial and temporal frequency content actually is -- without a rig, a fly,
or a trigger.

    from visanalysis.util import zebra_explorer as ze

    r = ze.analyze(scale=1.09, n_teeth=4, w_per_second=0.056)
    ze.plot(r, 'zebra.png')

    ze.compare([dict(scale=0.8), dict(scale=1.09), dict(scale=1.5)], 'scale_sweep.png')

or from a shell::

    python -m visanalysis.util.zebra_explorer --scale 1.09 --n-teeth 4
    python -m visanalysis.util.zebra_explorer --sweep scale 0.8,1.09,1.5

WHAT IS MEASURED, AND ON WHAT
-----------------------------
The stimulus is rendered exactly as the protocol would render it -- same generator, same
resolution, same upsampling, same comb -- and the spectra are then computed on a CENTRAL WINDOW
of that texture, not the whole thing.

The window is there because the spectral tools are planar: they assume one pixel subtends one
constant solid angle. A GlSphericalTexturedRect's texel grid is linear in (theta, phi), which is
locally equal-angle at the patch centre and increasingly compressed away from it -- a step of
d(theta) subtends d(theta)*sin(phi). Over the default 40 degree window the worst corner is off by
1 - sin(70 deg) = 6%, which is small against the things being read off these plots. Over the full
136 x 102 degree patch it would be 36%, which is not, and it would show up as spurious anisotropy.

So: the stimulus is the real one, the analysis is honest for the middle of it, and anything read
here is a statement about the stimulus near the centre of gaze.
"""
import argparse
import sys
import time
import warnings

import numpy as np

from visanalysis.util import stimulus_spectra as ss

# The reference categorical order from the dataviz palette, used unchanged and in slot order.
# The ordering is itself the colorblind-safety mechanism, so slots are taken from the front and
# never reordered or cycled; past 5 curves this tool facets instead of inventing a 6th hue.
SERIES_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']

# Fraction of an autocorrelation's lag range that is worth showing. Past this the unbiased
# estimator has too few overlapping products left to mean anything -- see full_patch_measures.
LAG_TRUST = 0.25

# Radial spectrum binning. Fewer bins than the 40 first used, because a log spacing that fine
# out-resolves the discrete Fourier grid at low frequency and leaves bins empty; and a floor on
# how many coefficients a bin needs before its average is worth plotting.
RADIAL_BINS = 28
MIN_BIN_COUNT = 4
INK = '#1a1a1a'
INK_MUTED = '#6b6b6b'
GRID = '#d8d8d8'

DEFAULTS = dict(
    # Post-comb: these describe the displayed stimulus, matching JCS_protocol's defaults.
    scale=1.09,
    n_teeth=4,
    w_per_second=0.056,
    seed=0,
    octaves=6,
    persistence=0.2,
    width=136.0,
    height=102.0,
    n_rows=384,
    n_cols=512,
    gen_n_rows=96,
    gen_n_cols=128,
    fps=120.0,
)


def _zebra_noise():
    try:
        from labpack.visual_stim.clandinin import zebra_noise
    except ImportError as exc:                                     # pragma: no cover
        raise ImportError(
            'zebra_explorer needs the labpack generator. Put clandinin_labpack on the path, '
            'e.g. PYTHONPATH=/path/to/clandinin_labpack.') from exc
    return zebra_noise


def render(n_frames=1024, frame_stride=1, verbose=True, **params):
    """Render the stimulus as the display would. Returns (frames, resolved parameters).

    frames is (rows, cols, T) uint8 of 0/255, in the texture's own (elevation, azimuth) layout.
    """
    zn = _zebra_noise()
    p = dict(DEFAULTS)
    p.update(params)

    stim = zn.ZebraSphereStimulus(
        scale=p['scale'], w_per_second=p['w_per_second'], n_teeth=p['n_teeth'],
        octaves=p['octaves'], persistence=p['persistence'], seed=p['seed'])

    indices = np.arange(0, n_frames * frame_stride, frame_stride)
    t0 = time.perf_counter()
    frames = stim.render(indices, p['n_rows'], p['n_cols'], p['width'], p['height'],
                         fps=p['fps'], gen_n_rows=p['gen_n_rows'],
                         gen_n_cols=p['gen_n_cols'])
    if verbose:
        print(f'  rendered {frames.shape[2]} frames in {time.perf_counter() - t0:.1f} s '
              f'({stim.scale:.4f} scale, {stim.w_per_second:.4f} w/s)')
    p['scale'] = float(stim.scale)
    p['w_per_second'] = float(stim.w_per_second)
    p['measured_feature_deg'] = (None if stim.measured_feature_deg is None
                                 else float(stim.measured_feature_deg))
    return frames, p


def analysis_windows(frames, p, window_deg=24.0, n_az=4, n_el=3, max_gradient=0.10):
    """Tile the patch with small windows, each with its own correctly scaled grid.

    Yields (cropped frames, PlanarGrid, (azimuth, elevation)) per window.

    WHY MANY SMALL WINDOWS RATHER THAN ONE BIG ONE

    The spectral tools are planar: they assume a pixel subtends a constant solid angle. On a
    (theta, phi) texel grid a texel spans a constant d(phi) in elevation but d(theta)*sin(phi)
    in azimuth, so the azimuthal scale depends on WHERE you look. One 40 degree window at the
    centre handled that by staying where sin(phi) ~ 1; small windows handle it better, because
    each gets the sin(phi) of its own centre and the leftover error is only the variation ACROSS
    the window.

    That leftover goes as 2*h*tan(elevation) for a window of half-height h. It is exactly zero at
    the equator and grows away from it, which is why the elevation rows are capped by
    `max_gradient` rather than by a fixed extent -- shrink the windows and more rows qualify.

    Azimuth is free: sin(phi) does not depend on theta at all, so tiling across azimuth adds
    samples at no geometric cost whatever.

    More windows also buy statistics. The orientation spread has a sampling floor set by the
    number of independent samples, and -- more usefully -- real anisotropy is CONSISTENT across
    windows while sampling noise is not, which is a far stronger test than splitting one window.
    """
    deg_x_equator = p['width'] / p['n_cols']
    deg_y = p['height'] / p['n_rows']
    n_r = max(int(round(window_deg / deg_y)), 16)
    n_c = max(int(round(window_deg / deg_x_equator)), 16)
    half_h = 0.5 * n_r * deg_y

    rows, cols = frames.shape[0], frames.shape[1]

    # How far from the equator a window of this height may sit and still be locally planar to
    # within max_gradient. Solving 2*h*tan(e) <= max_gradient for e.
    #
    # Rows are placed INSIDE this band rather than spread over the whole patch and then filtered.
    # Spreading them first put the outer rows at +-34 degrees, where the gradient is 0.28, so
    # every row but the equator was discarded and the elevation sampling silently collapsed to
    # one. Note the trade: smaller windows reach higher elevations, so window_deg buys either
    # frequency resolution or elevation coverage, not both.
    max_elev = float(np.degrees(np.arctan(max_gradient / (2.0 * np.radians(half_h)))))
    max_elev = min(max_elev, 0.5 * p['height'] - half_h)

    out = []
    for ei in range(n_el):
        frac = 0.5 if n_el == 1 else ei / (n_el - 1)
        elev_target = (2.0 * frac - 1.0) * max_elev
        r0 = int(round((elev_target / p['height'] + 0.5) * rows - n_r / 2))
        r0 = min(max(r0, 0), rows - n_r)
        elev = ((r0 + n_r / 2) / rows - 0.5) * p['height']
        # A texel's azimuthal extent at this elevation. phi is the polar angle, so
        # sin(phi) = cos(elevation).
        deg_x = deg_x_equator * float(np.cos(np.radians(elev)))
        for ai in range(n_az):
            c0 = int(round(((ai + 0.5) / n_az) * cols - n_c / 2))
            c0 = min(max(c0, 0), cols - n_c)
            azim = ((c0 + n_c / 2) / cols - 0.5) * p['width']
            crop = np.ascontiguousarray(frames[r0:r0 + n_r, c0:c0 + n_c, :])
            out.append((crop, ss.PlanarGrid(deg_per_px_x=deg_x, deg_per_px_y=deg_y),
                        (azim, elev)))
    if not out:
        raise ValueError('no analysis window met the flatness limit; raise max_gradient '
                         'or lower window_deg')
    return out


def _median_of(x, weight):
    """Frequency below which half the (weighted) power lies. Robust to NaN bins."""
    ok = np.isfinite(x) & np.isfinite(weight) & (x > 0) & (weight > 0)
    if ok.sum() < 3:
        return float('nan')
    cum = np.cumsum(weight[ok])
    cum = cum / cum[-1]
    return float(x[ok][np.searchsorted(cum, 0.5)])


def _crossing(x, y, level=0.5):
    """Where a falling curve first crosses `level`, linearly interpolated."""
    below = np.flatnonzero(np.asarray(y) < level)
    if below.size == 0:
        return float('nan')
    i = int(below[0])
    if i == 0:
        return float(x[0])
    y0, y1 = y[i - 1], y[i]
    if y1 == y0:
        return float(x[i])
    return float(x[i - 1] + (level - y0) * (x[i] - x[i - 1]) / (y1 - y0))


def _mean_curve(curves):
    """Average binned curves that share x but were measured on different windows.

    A bin that no window sampled is all-NaN across the stack, and NaN is the right answer for it
    -- the speed and radial binnings both span wider ranges than any one window populates, so
    empty bins at the extremes are normal rather than a symptom. numpy warns about the all-NaN
    slice anyway, so the warning is silenced HERE, narrowly, rather than left to print on every
    run until people stop reading warnings.
    """
    x = curves[0][0]
    stack = np.vstack([c[1] for c in curves])
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Mean of empty slice',
                                category=RuntimeWarning)
        y = np.nanmean(stack, axis=0)
    return x, y


def analyze(n_frames=1024, window_deg=24.0, n_az=4, n_el=3, verbose=True,
            min_snr=3.0, max_frames=8192, n_seeds=8, **params):
    """Characterise one parameter set, adding samples until the measurement can resolve.

    The orientation statistic has a sampling floor that depends on how much independent data went
    into it, and that in turn depends on the STIMULUS: a coarse, slow stimulus puts far fewer
    independent samples in a fixed number of frames than a fine, fast one. Measured across a
    feature-size sweep, the floor ran 0.090 at 2 degree clusters up to 0.350 at 20 -- at which
    point the spread itself was only 0.49, so "spread exceeds floor" meant almost nothing.

    So the frame count is not a fixed number but a target: keep doubling until the spread is
    `min_snr` times the floor, or until `max_frames`. Doubling frames rather than windows because
    the windows are already tiled as densely as the flatness limit allows -- more of them would
    overlap, which adds pixels but not independent samples.

    A run that stops at max_frames without reaching min_snr is reported rather than silently
    accepted: at that point the honest reading is "too few independent samples to tell", and the
    difference matters when the number being read off is whether the stimulus is isotropic.

    `n_seeds` is the other lever, and usually the better one. These are properties of the noise
    PROCESS, not of one realisation, so independent seeds are independent samples of the same
    thing and averaging their spectra cuts the sampling floor by sqrt(n_seeds). Frames bought
    within one seed are correlated over tau and stop helping once the record is many tau long;
    seeds never are.
    """
    n = int(n_frames)
    while True:
        result = _analyze_over_seeds(n, window_deg, n_az, n_el, verbose,
                                     int(n_seeds), **params)
        s = result['summary']
        floor = s['anisotropy_floor']
        snr = (s['anisotropy'] / floor) if (floor and np.isfinite(floor) and floor > 0) else 0.0
        result['summary']['anisotropy_snr'] = float(snr)
        result['summary']['n_frames_used'] = n
        if snr >= min_snr or n * 2 > max_frames:
            if snr < min_snr and verbose:
                print(f'  NOTE: orientation spread is only {snr:.1f}x its sampling floor at '
                      f'{n} frames (wanted {min_snr:g}x). Too few independent samples to call '
                      f'this stimulus isotropic or not.')
            result['summary']['anisotropy_resolved'] = bool(snr >= min_snr)
            return result
        n *= 2
        if verbose:
            print(f'  spread {s["anisotropy"]:.3f} is only {snr:.1f}x the floor '
                  f'{floor:.3f}; retrying at {n} frames')


def full_patch_measures(frames, p, frame_stride=8):
    """Stripe width, correlation lengths and speed, measured over the WHOLE patch.

    These four do not need the windows, and are wrong when taken from them. The windows exist so
    the planar spectral tools see a locally equal-angle grid; but counting stripe transitions
    along a meridian needs no such thing, because the ELEVATION axis has a constant degrees per
    texel everywhere (a texel spans d(phi) in elevation and only d(theta)*sin(phi) in azimuth).
    Taking them from a 24 degree window instead undersamples badly once the stimulus is coarse:
    at the current defaults that read 38 degree stripes against a true 16, and 2.8 deg/s against
    a true 90, simply because fewer than one stripe fits inside a window.

    So: spectra from the windows, these from the whole patch. Measured along elevation only, and
    the azimuth axis is left alone precisely because its scale is not constant.
    """
    deg_y = p['height'] / p['n_rows']
    idx = range(0, frames.shape[2], max(int(frame_stride), 1))

    # Stripe width from edge density down the elevation axis.
    edge = np.mean([np.mean(np.abs(np.diff(frames[:, :, k].astype(np.int16), axis=0)) > 0)
                    for k in idx])
    stripe = (deg_y / edge) if edge > 0 else float('inf')

    # Spatial correlation along the same axis, as a 1-D autocorrelation per column.
    x = frames[:, ::4, ::max(int(frame_stride), 1)].astype(np.float32)
    x = x.reshape(x.shape[0], -1)
    x = x - x.mean(axis=0, keepdims=True)
    n = x.shape[0]
    pad = 1 << int(np.ceil(np.log2(2 * n)))
    F = np.fft.rfft(x, n=pad, axis=0)
    ac = np.fft.irfft(F * np.conj(F), n=pad, axis=0)[:n]
    ac /= np.arange(n, 0, -1)[:, None]
    with np.errstate(invalid='ignore', divide='ignore'):
        ac = ac / ac[:1]
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Mean of empty slice',
                                category=RuntimeWarning)
        ac_space = np.nanmean(ac, axis=1)
    lags_space = np.arange(n) * deg_y
    keep_s = lags_space <= LAG_TRUST * lags_space[-1]
    lags_space, ac_space = lags_space[keep_s], ac_space[keep_s]
    corr_deg = _crossing(lags_space, ac_space, 0.5)

    # Temporal correlation over the whole patch, on a subsampled pixel set.
    #
    # Truncated to a quarter of the record, because past that the estimate is not usable. The
    # unbiased normalisation divides lag k by the (n - k) overlapping products available, so its
    # variance explodes as k approaches the record length; and removing a mean estimated from the
    # same finite record forces the autocovariance to sum to about zero, which drags mid-lags
    # negative and lifts the tail back up. That upturn is entirely an artifact -- its position
    # tracks the record length (95% of a 2.1 s record, 51% of a 17 s one) rather than sitting at
    # any fixed lag -- so plotting it invites reading an estimator's failure as structure.
    sub = np.ascontiguousarray(frames[::5, ::5, :])
    lags_t, ac_t = ss.temporal_autocorrelation(sub, p['fps'])
    keep_t = lags_t <= LAG_TRUST * lags_t[-1]
    lags_t, ac_t = lags_t[keep_t], ac_t[keep_t]
    corr_s = _crossing(lags_t, ac_t, 0.5)

    # Speed as one stripe crossed per correlation time. Robust because both terms are.
    speed = (stripe / corr_s) if (np.isfinite(corr_s) and corr_s > 0) else float('nan')
    return dict(stripe_deg=float(stripe), corr_deg=float(corr_deg),
                corr_s=float(corr_s), speed_dps=float(speed),
                lags_space=lags_space, ac_space=ac_space,
                lags_t=lags_t, ac_t=ac_t)


def _analyze_over_seeds(n_frames, window_deg, n_az, n_el, verbose, n_seeds, **params):
    """Average the curves over independent seeds, which are independent samples of the process.

    Only the curves and the scalars derived from them are averaged; the sample frame comes from
    the first seed, since an average of binary images is not a stimulus anyone will display.
    """
    base_seed = int(params.pop('seed', DEFAULTS['seed']))
    runs = []
    for k in range(max(int(n_seeds), 1)):
        runs.append(_analyze_once(n_frames, window_deg, n_az, n_el,
                                  verbose and k == 0, seed=base_seed + 1000 * k, **params))
    if len(runs) == 1:
        return runs[0]

    out = dict(runs[0])
    # Bins no seed populated stay NaN, which is the right answer; see _mean_curve for why the
    # warning is silenced rather than the case avoided.
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='Mean of empty slice',
                                category=RuntimeWarning)
        for key in ('radial_p', 'orient_p', 'temporal_p', 'speed_p', 'acorr', 'acorr_t'):
            out[key] = np.nanmean(np.vstack([r[key] for r in runs]), axis=0)

    # The orientation floor has to be recomputed from the seed-to-seed scatter, not averaged:
    # averaging the per-run floors would report the floor of a single run while the spread is
    # now that of the mean, which would overstate the anisotropy by sqrt(n_seeds).
    per_seed = np.vstack([r['orient_p'] for r in runs])
    fin = np.all(np.isfinite(per_seed), axis=0)
    summary = dict(runs[0]['summary'])
    if fin.sum() > 4:
        norm = per_seed[:, fin] / np.nanmean(per_seed[:, fin], axis=1, keepdims=True)
        summary['anisotropy'] = float(np.nanstd(out['orient_p'][fin])
                                      / np.nanmean(out['orient_p'][fin]))
        summary['anisotropy_floor'] = float(np.nanmean(np.nanstd(norm, axis=0))
                                            / np.sqrt(norm.shape[0]))
        cc = np.corrcoef(norm)
        iu = np.triu_indices(cc.shape[0], 1)
        summary['anisotropy_consistency'] = float(np.nanmean(cc[iu]))
    for key in ('stripe_deg', 'stimulus_corr_deg', 'stimulus_corr_s', 'speed_dps',
                'median_spatial_cpd', 'median_temporal_hz', 'white_fraction'):
        summary[key] = float(np.nanmean([r['summary'][key] for r in runs]))
    summary['n_seeds'] = len(runs)
    out['summary'] = summary
    out['params'] = dict(runs[0]['params'], seed=base_seed)
    return out


def _analyze_once(n_frames, window_deg, n_az, n_el, verbose, **params):
    """One pass at a fixed frame count. See analyze() for the adaptive wrapper.

    Spectra are averaged over a tile of small windows across azimuth and elevation rather than
    computed on one big central one -- see analysis_windows for the geometry, and for why the
    binned 1D curves are averaged rather than the raw 2D power (windows at different elevations
    have different frequency axes, so their power arrays are not commensurable).
    """
    frames, p = render(n_frames=n_frames, verbose=verbose, **params)
    fps = p['fps']
    windows = analysis_windows(frames, p, window_deg, n_az, n_el)

    if verbose:
        crop0, grid0, _ = windows[0]
        els = sorted({round(w[2][1]) for w in windows})
        print(f'  {len(windows)} windows of {crop0.shape[1]}x{crop0.shape[0]} texels '
              f'({crop0.shape[1] * grid0.deg_per_px_x:.0f} x '
              f'{crop0.shape[0] * grid0.deg_per_px_y:.0f} deg), elevations {els}')

    radial, orient, temporal, speed, sp_ac, t_ac, stripes = [], [], [], [], [], [], []
    for crop, grid, _pos in windows:
        power = ss.frame_power_spectrum(crop, grid=grid)
        coords = ss.planar_frequency_coords(power.shape, grid)

        # Blank out radial bins with too few Fourier coefficients to mean anything. At the
        # low-frequency end a log-spaced binning is finer than the discrete Fourier grid, so
        # bins land BETWEEN adjacent modes: measured on a 90x90 window, 7 of 40 bins held zero
        # coefficients and the first held two. The empty ones plot as gaps -- the visible breaks
        # in the spatial spectrum -- and the nearly-empty ones plot as points made of noise,
        # which is worse because they look like data.
        f_r, p_r, _, n_r = ss.radial_spectrum(power, coords, n_bins=RADIAL_BINS)
        p_r = np.where(n_r >= MIN_BIN_COUNT, p_r, np.nan)
        radial.append((f_r, p_r))
        ang, p_ang, _, _ = ss.orientation_spectrum(power, coords, n_bins=36)
        orient.append((ang, p_ang))
        temporal.append(ss.temporal_spectrum(crop, fps, nperseg=min(512, crop.shape[2] // 2)))

        st_c, st_ft, st_p = ss.spatiotemporal_spectrum(
            crop, grid, fps, nperseg=min(256, crop.shape[2] // 4), spatial_stride=1)
        sp, p_sp, _, _ = ss.speed_spectrum(st_c, st_ft, st_p, n_bins=40,
                                           speed_range=(0.2, 2000.0))
        speed.append((sp, p_sp))

        lag_x, lag_y, acorr = ss.spatial_autocorrelation(crop, grid)
        LX, LY = np.meshgrid(lag_x, lag_y, indexing='xy')
        rad = np.hypot(LX, LY)
        edges = np.linspace(0, float(np.nanmax(rad)) / np.sqrt(2), 40)
        r_c, r_ac, _, _ = ss.bin_average(acorr, rad, edges)
        sp_ac.append((r_c, r_ac))
        t_ac.append(ss.temporal_autocorrelation(crop, fps))

        e = float(np.mean(np.abs(np.diff(crop[:, :, 0].astype(np.int16), axis=1)) > 0))
        stripes.append(grid.deg_per_px_x / e if e > 0 else np.nan)

    f_r, p_r = _mean_curve(radial)
    ang, p_ang = _mean_curve(orient)
    f_t, p_t = _mean_curve(temporal)
    spd, p_spd = _mean_curve(speed)
    r_c, r_ac = _mean_curve(sp_ac)
    lags_t, ac_t = _mean_curve(t_ac)

    # Orientation, judged against the two things that matter.
    #
    # SPREAD is the coefficient of variation of the averaged angular spectrum. On its own it is
    # not interpretable: a bin's mean is dominated by its few lowest-frequency members on a
    # steeply falling spectrum, so even an isotropic field gives a non-zero spread. Two earlier
    # attempts at a floor were both wrong -- a derived 1/sqrt(N) gave 0.005, forty times too
    # small; the pre-comb field, isotropic by construction, gave 0.71, because its spectrum is
    # steeper and the floor depends on spectral slope.
    #
    # FLOOR is the scatter of the per-window means, which needs no model: independent windows of
    # the same isotropic stimulus differ only by sampling.
    #
    # CONSISTENCY is the stronger test and the reason for tiling. Real anisotropy prefers the
    # same orientations in every window and correlates across them; sampling noise does not. A
    # spread above the floor with consistency near zero is noise that happened to be lumpy.
    per_window = np.vstack([c[1] for c in orient])
    fin = np.all(np.isfinite(per_window), axis=0)
    if fin.sum() > 4 and per_window.shape[0] > 1:
        norm = per_window[:, fin] / np.nanmean(per_window[:, fin], axis=1, keepdims=True)
        anisotropy = float(np.nanstd(p_ang[fin]) / np.nanmean(p_ang[fin]))
        anisotropy_floor = float(np.nanmean(np.nanstd(norm, axis=0)) / np.sqrt(norm.shape[0]))
        cc = np.corrcoef(norm)
        iu = np.triu_indices(cc.shape[0], 1)
        anisotropy_consistency = float(np.nanmean(cc[iu]))
    else:
        anisotropy = anisotropy_floor = anisotropy_consistency = float('nan')

    full = full_patch_measures(frames, p)
    stripe_deg = full['stripe_deg']
    r_c, r_ac = full['lags_space'], full['ac_space']
    lags_t, ac_t = full['lags_t'], full['ac_t']

    # The field's own correlation lengths, for reference only. feature_deg and tau_s are defined
    # on the DISPLAYED stimulus now, so these are no longer the calibration targets -- they are
    # here because the ratio between them and the displayed values is what n_teeth controls.
    zn = _zebra_noise()
    field_deg, _, _ = zn.correlation_length_deg(p['scale'], p['octaves'], p['persistence'],
                                                p['seed'])
    field_tau, _, _ = zn.temporal_correlation_s(p['scale'], p['w_per_second'], p['octaves'],
                                                p['persistence'], p['seed'],
                                                max_s=max(2.0, 60.0 * zn.cluster_seconds(p['w_per_second'])))

    grid = windows[0][1]
    crop_shape = windows[0][0].shape
    coords = ss.planar_frequency_coords(
        ss.frame_power_spectrum(windows[0][0], grid=grid).shape, grid)

    return dict(
        params=p, frame=frames[:, :, 0], crop_shape=crop_shape, grid=grid,
        n_windows=len(windows),
        radial_f=f_r, radial_p=p_r,
        orient_deg=ang, orient_p=p_ang,
        temporal_f=f_t, temporal_p=p_t,
        speed=spd, speed_p=p_spd,
        acorr_lag_deg=r_c, acorr=r_ac,
        acorr_lag_s=lags_t, acorr_t=ac_t,
        summary=dict(
            median_spatial_cpd=_median_of(f_r, p_r * f_r),
            median_temporal_hz=_median_of(f_t, p_t),
            median_speed_dps=_median_of(spd, p_spd),
            field_feature_deg=float(field_deg),
            field_tau_s=float(field_tau),
            stimulus_corr_deg=full['corr_deg'],
            stimulus_corr_s=full['corr_s'],
            speed_dps=full['speed_dps'],
            anisotropy_consistency=anisotropy_consistency,
            stripe_deg=stripe_deg,
            white_fraction=float(np.mean(frames > 0)),
            anisotropy=anisotropy,
            anisotropy_floor=anisotropy_floor,
            nyquist_cpd=float(coords.nyquist),
        ),
    )


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------

def _populated(x, y):
    """Drop bins with no data, so a curve plots as one line rather than several.

    The empty bins have to survive the averaging -- every window shares one bin grid, and
    dropping them per-window would misalign the stack -- so they are only removed here, on the
    way to the axes. A gap in a log-binned radial spectrum is a statement about the Fourier grid,
    not about the stimulus, and drawing it as a break in the curve reads as the latter.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(x) & np.isfinite(y)
    return x[ok], y[ok]


def _style(ax):
    """Recessive axes and grid: the data should be the only assertive thing on the panel."""
    ax.grid(True, color=GRID, linewidth=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_MUTED, labelsize=8, length=3)
    ax.xaxis.label.set_color(INK_MUTED)
    ax.yaxis.label.set_color(INK_MUTED)
    ax.title.set_color(INK)


def _marker_line(ax, x, label, color=INK_MUTED):
    """A vertical reference at a summary value, labelled in ink rather than in a series color."""
    if not np.isfinite(x):
        return
    ax.axvline(x, color=color, linewidth=1.0, linestyle='--', alpha=0.8)
    ax.annotate(label, xy=(x, 1.0), xycoords=('data', 'axes fraction'),
                xytext=(3, -10), textcoords='offset points',
                fontsize=7.5, color=INK_MUTED, ha='left', va='top')


def plot(result, path=None, title=None):
    """One parameter set: the stimulus, its spectra, and the numbers read off them."""
    import matplotlib
    if path is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    p, s = result['params'], result['summary']
    fig = plt.figure(figsize=(13.5, 8.6), facecolor='white')
    gs = fig.add_gridspec(2, 4, hspace=0.42, wspace=0.30,
                          left=0.055, right=0.985, top=0.86, bottom=0.08)

    head = title or (f"Zebra noise  |  scale {p['scale']:g}, {p['n_teeth']:g} teeth, "
                     f"w {p['w_per_second']:g}/s")
    fig.text(0.055, 0.955, head, fontsize=14, color=INK, weight='bold')
    fig.text(0.055, 0.918,
             f"{p['width']:g} x {p['height']:g} deg patch at {p['n_rows']}x{p['n_cols']} texels "
             f"from {p['gen_n_rows']}x{p['gen_n_cols']}  |  {p['fps']:g} Hz  |  seed {p['seed']}"
             f"  |  spectra averaged over {result['n_windows']} windows of "
             f"{result['crop_shape'][1] * result['grid'].deg_per_px_x:.0f} deg",
             fontsize=9, color=INK_MUTED)

    # 1. the stimulus itself. Grayscale because the data IS luminance -- a color map here would
    # invent a dimension the stimulus does not have.
    ax = fig.add_subplot(gs[0, 0])
    ax.imshow(result['frame'], cmap='gray', vmin=0, vmax=255, aspect='equal',
              extent=[-p['width'] / 2, p['width'] / 2, -p['height'] / 2, p['height'] / 2])
    ax.set_title('stimulus, frame 0', fontsize=10, loc='left')
    ax.set_xlabel('azimuth (deg)')
    ax.set_ylabel('elevation (deg)')
    ax.tick_params(colors=INK_MUTED, labelsize=8)
    ax.title.set_color(INK)

    # 2. spatial
    ax = fig.add_subplot(gs[0, 1])
    ax.loglog(*_populated(result['radial_f'], result['radial_p']),
              color=SERIES_COLORS[0], linewidth=1.8)
    _marker_line(ax, s['median_spatial_cpd'], f"median {s['median_spatial_cpd']:.3f}")
    _marker_line(ax, s['nyquist_cpd'], 'Nyquist')
    ax.set_title('spatial power spectrum', fontsize=10, loc='left')
    ax.set_xlabel('spatial frequency (cycles/deg)')
    ax.set_ylabel('power')
    _style(ax)

    # 3. temporal
    ax = fig.add_subplot(gs[0, 2])
    m = result['temporal_f'] > 0
    ax.loglog(*_populated(result['temporal_f'][m], result['temporal_p'][m]),
              color=SERIES_COLORS[0], linewidth=1.8)
    _marker_line(ax, s['median_temporal_hz'], f"median {s['median_temporal_hz']:.2f} Hz")
    ax.set_title('temporal power spectrum', fontsize=10, loc='left')
    ax.set_xlabel('temporal frequency (Hz)')
    ax.set_ylabel('power')
    _style(ax)

    # 4. speed
    ax = fig.add_subplot(gs[0, 3])
    ax.semilogx(*_populated(result['speed'], result['speed_p']),
                color=SERIES_COLORS[0], linewidth=1.8)
    _marker_line(ax, s['median_speed_dps'], f"median {s['median_speed_dps']:.0f} deg/s")
    ax.set_title('speed distribution', fontsize=10, loc='left')
    ax.set_xlabel('speed (deg/s)')
    ax.set_ylabel('power')
    _style(ax)

    # 5. spatial autocorrelation
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(result['acorr_lag_deg'], result['acorr'], color=SERIES_COLORS[0], linewidth=1.8)
    ax.axhline(0.5, color=GRID, linewidth=1.0)
    _marker_line(ax, s['stimulus_corr_deg'], f"{s['stimulus_corr_deg']:.2f} deg")
    ax.set_title('spatial autocorrelation, full patch', fontsize=10, loc='left')
    ax.set_xlabel('lag (deg)')
    ax.set_ylabel('correlation')
    _style(ax)

    # 6. temporal autocorrelation
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(result['acorr_lag_s'], result['acorr_t'], color=SERIES_COLORS[0], linewidth=1.8)
    ax.axhline(0.5, color=GRID, linewidth=1.0)
    _marker_line(ax, s['stimulus_corr_s'], f"{s['stimulus_corr_s']:.3f} s")
    ax.set_title('temporal autocorrelation, full patch', fontsize=10, loc='left')
    ax.set_xlabel('lag (s)')
    ax.set_ylabel('correlation')
    _style(ax)

    # 7. orientation
    ax = fig.add_subplot(gs[1, 2])
    ax.plot(result['orient_deg'], result['orient_p'], color=SERIES_COLORS[0], linewidth=1.8)
    ax.set_ylim(bottom=0)
    ax.axhline(np.nanmean(result['orient_p']), color=GRID, linewidth=1.0)
    ax.set_title(f"orientation  (spread {s['anisotropy']:.3f} vs floor "
                 f"{s['anisotropy_floor']:.3f}, consistency {s['anisotropy_consistency']:+.2f})",
                 fontsize=9.5, loc='left')
    ax.set_xlabel('orientation (deg)')
    ax.set_ylabel('power')
    _style(ax)

    # 8. the numbers, as text rather than as a chart of eight unrelated scalars
    ax = fig.add_subplot(gs[1, 3])
    ax.axis('off')
    rows = [
        ('asked for', ''),
        ('  scale', f"{p['scale']:g}"),
        ('  w_per_second', f"{p['w_per_second']:g}"),
        ('  n_teeth', f"{p['n_teeth']:g}"),
        ('measured on the frames', ''),
        # The same quantities feature_deg and tau_s name, but read off the rendered frames rather
        # than from the calibration's own direction-based estimator. Agreement is a check that the
        # two independent measurements of the same thing land together; they differ by ~10% because
        # this one is limited by the window size and the radial binning.
        ('  feature_deg', f"{s['stimulus_corr_deg']:.2f}"),
        ('  tau_s', f"{s['stimulus_corr_s']:.3f}"),
        ('  stripe width (deg)', f"{s['stripe_deg']:.2f}"),
        ('  median spatial (c/deg)', f"{s['median_spatial_cpd']:.3f}"),
        ('  median temporal (Hz)', f"{s['median_temporal_hz']:.2f}"),
        ('  speed (deg/s)', f"{s['speed_dps']:.0f}"),
        ('  median speed, windows', f"{s['median_speed_dps']:.0f}"),
        ('  white fraction', f"{s['white_fraction']:.3f}"),
        ('  orientation spread', f"{s['anisotropy']:.3f}"),
        ('  sampling floor', f"{s['anisotropy_floor']:.3f}"),
        ('  window consistency', f"{s['anisotropy_consistency']:+.3f}"),
        ('  spread / floor', f"{s.get('anisotropy_snr', float('nan')):.1f}x"
                             f"{'' if s.get('anisotropy_resolved', True) else '  UNRESOLVED'}"),
        ('  frames used', f"{s.get('n_frames_used', 0)}"),
        ('  seeds averaged', f"{s.get('n_seeds', 1)}"),
        ('field, for reference', ''),
        ('  correlation (deg)', f"{s['field_feature_deg']:.1f}"),
        ('  correlation (s)', f"{s['field_tau_s']:.3f}"),
    ]
    y = 0.99
    for name, val in rows:
        weight = 'bold' if not name.startswith(' ') else 'normal'
        ax.text(0.0, y, name, fontsize=8.5, color=INK if weight == 'bold' else INK_MUTED,
                weight=weight, va='top', transform=ax.transAxes)
        ax.text(1.0, y, val, fontsize=8.5, color=INK, va='top', ha='right',
                family='monospace', transform=ax.transAxes)
        y -= 0.058

    if path:
        fig.savefig(path, dpi=130, facecolor='white')
        plt.close(fig)
        print(f'wrote {path}')
    return fig


def compare(param_sets, path=None, labels=None, n_frames=512, window_deg=24.0,
            title=None, verbose=True, min_snr=0.0, max_frames=2048, n_seeds=2):
    """Overlay several parameter sets on the four spectral panels.

    Capped at five sets: the categorical order is only validated for colorblind separation over
    its published slots, and a sixth curve would mean inventing a hue. Run two figures instead.

    Cheaper defaults than analyze() on purpose. A comparison is about the SHAPE of the curves --
    where the spectrum sits, which way it moved -- not about resolving whether one stimulus is
    isotropic, so min_snr defaults to 0 and the frame escalation is off. With the analyze()
    defaults a three-set sweep took over ten minutes, since each set independently climbed to
    8192 frames chasing an anisotropy figure nobody reads off a comparison. Raise min_snr if you
    do want the orientation panel to be trustworthy here.
    """
    if len(param_sets) > len(SERIES_COLORS):
        raise ValueError(
            f'compare() takes at most {len(SERIES_COLORS)} parameter sets; got '
            f'{len(param_sets)}. Split them across two figures rather than adding a hue.')

    import matplotlib
    if path is not None:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    results = []
    for i, ps in enumerate(param_sets):
        if verbose:
            print(f'[{i + 1}/{len(param_sets)}] {ps}')
        results.append(analyze(n_frames=n_frames, window_deg=window_deg,
                               verbose=verbose, min_snr=min_snr,
                               max_frames=max_frames, n_seeds=n_seeds, **ps))

    if labels is None:
        varying = sorted({k for ps in param_sets for k in ps})
        labels = [', '.join(f'{k}={ps.get(k, DEFAULTS.get(k))!s}' for k in varying)
                  for ps in param_sets]

    fig = plt.figure(figsize=(13.5, 7.4), facecolor='white')
    gs = fig.add_gridspec(2, 3, hspace=0.40, wspace=0.28,
                          left=0.06, right=0.985, top=0.84, bottom=0.09)
    fig.text(0.06, 0.94, title or 'Zebra noise parameter comparison',
             fontsize=14, color=INK, weight='bold')

    panels = [
        (gs[0, 0], 'radial_f', 'radial_p', 'spatial power spectrum',
         'spatial frequency (cycles/deg)', 'loglog'),
        (gs[0, 1], 'temporal_f', 'temporal_p', 'temporal power spectrum',
         'temporal frequency (Hz)', 'loglog'),
        (gs[0, 2], 'speed', 'speed_p', 'speed distribution', 'speed (deg/s)', 'semilogx'),
        (gs[1, 0], 'acorr_lag_deg', 'acorr', 'spatial autocorrelation', 'lag (deg)', 'plot'),
        (gs[1, 1], 'acorr_lag_s', 'acorr_t', 'temporal autocorrelation', 'lag (s)', 'plot'),
        (gs[1, 2], 'orient_deg', 'orient_p', 'orientation', 'orientation (deg)', 'plot'),
    ]
    for cell, xk, yk, name, xlabel, scale in panels:
        ax = fig.add_subplot(cell)
        for r, colour, label in zip(results, SERIES_COLORS, labels):
            x, y = np.asarray(r[xk], dtype=float), np.asarray(r[yk], dtype=float)
            if scale == 'loglog':
                m = x > 0
                ax.loglog(*_populated(x[m], y[m]), color=colour, linewidth=1.8, label=label)
            elif scale == 'semilogx':
                ax.semilogx(*_populated(x, y), color=colour, linewidth=1.8, label=label)
            else:
                ax.plot(*_populated(x, y), color=colour, linewidth=1.8, label=label)
        if 'autocorrelation' in name:
            ax.axhline(0.5, color=GRID, linewidth=1.0)
        if name == 'orientation':
            ax.set_ylim(bottom=0)
        ax.set_title(name, fontsize=10, loc='left')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('correlation' if 'autocorrelation' in name else 'power')
        _style(ax)

    # A legend is always present for two or more series, so identity is never carried by color
    # alone for someone who cannot separate the hues.
    handles, lbls = fig.axes[0].get_legend_handles_labels()
    leg = fig.legend(handles, lbls, loc='upper right', bbox_to_anchor=(0.985, 0.965),
                     frameon=False, fontsize=9, ncol=min(len(labels), 3))
    for text in leg.get_texts():
        text.set_color(INK)

    print()
    print('feature/tau measured on the rendered frames; field = the pre-comb reference')
    print(f"{'set':<26} {'feature':>8} {'tau':>8} {'stripe':>8} {'c/deg':>8} "
          f"{'Hz':>7} {'deg/s':>7} {'field deg':>10}")
    print('-' * 92)
    for label, r in zip(labels, results):
        s = r['summary']
        print(f"{label[:25]:<26} {s['stimulus_corr_deg']:>8.2f} "
              f"{s['stimulus_corr_s']:>8.3f} {s['stripe_deg']:>8.2f} "
              f"{s['median_spatial_cpd']:>8.3f} {s['median_temporal_hz']:>7.2f} "
              f"{s['median_speed_dps']:>7.0f} {s['field_feature_deg']:>10.1f}")

    if path:
        fig.savefig(path, dpi=130, facecolor='white')
        plt.close(fig)
        print(f'\nwrote {path}')
    return fig, results


def main(argv=None):
    ap = argparse.ArgumentParser(
        description='Render Zebra noise offline and plot its spatial and temporal spectra.')
    for name, default in DEFAULTS.items():
        ap.add_argument(f"--{name.replace('_', '-')}", type=type(default), default=default)
    ap.add_argument('--n-frames', type=int, default=1024,
                    help='frames to render; more gives finer temporal resolution')
    ap.add_argument('--window-deg', type=float, default=24.0,
                    help='size of each analysis window, in degrees')
    ap.add_argument('--n-az', type=int, default=4, help='analysis windows across azimuth')
    ap.add_argument('--n-el', type=int, default=3, help='analysis window rows in elevation')
    ap.add_argument('--min-snr', type=float, default=3.0,
                    help='keep doubling frames until orientation spread is this '
                         'many times its sampling floor')
    ap.add_argument('--max-frames', type=int, default=8192,
                    help='give up adding frames past this')
    ap.add_argument('--n-seeds', type=int, default=8,
                    help='independent seeds to average; cuts the sampling floor '
                         'by sqrt of this')
    ap.add_argument('-o', '--out', default='zebra_spectra.png')
    ap.add_argument('--sweep', nargs=2, metavar=('PARAM', 'VALUES'),
                    help='compare one parameter across comma-separated values, '
                         'e.g. --sweep scale 0.8,1.09,1.5')
    args = ap.parse_args(argv)

    params = {k: getattr(args, k) for k in DEFAULTS}

    if args.sweep:
        name, raw = args.sweep
        # Accept the dashed spelling too: every other flag is --n-teeth, so requiring
        # --sweep n_teeth here is a trap the error message alone does not undo.
        name = name.replace('-', '_')
        if name not in DEFAULTS:
            ap.error(f'unknown parameter {name!r}; choose from {sorted(DEFAULTS)}')
        cast = type(DEFAULTS[name])
        values = [cast(v) for v in raw.split(',')]
        sets = [dict(params, **{name: v}) for v in values]
        compare(sets, args.out, labels=[f'{name}={v:g}' for v in values],
                n_frames=min(args.n_frames, 512), window_deg=args.window_deg,
                title=f'Zebra noise: {name} sweep',
                min_snr=args.min_snr if args.min_snr != 3.0 else 0.0,
                max_frames=args.max_frames, n_seeds=args.n_seeds)
        return 0

    print(f'rendering {args.n_frames} frames...')
    result = analyze(n_frames=args.n_frames, window_deg=args.window_deg,
                     n_az=args.n_az, n_el=args.n_el, min_snr=args.min_snr,
                     max_frames=args.max_frames, n_seeds=args.n_seeds, **params)
    plot(result, args.out)
    print()
    for k, v in result['summary'].items():
        print(f'  {k:<24} {v:.4f}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
