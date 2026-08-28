# -*- coding: utf-8 -*-
"""Bleaching correction for long-epoch noise recordings.

WHY MULTIPLICATIVE. Bleaching scales the whole signal, so raw fluorescence is
F(t) = B(t) * (1 + r(t)) with B the bleach envelope and r the fractional
response. Against a FIXED pre-stimulus baseline F0 = B(0),

    dF/F = [b(t) - 1] + b(t) * r(t),        b = B / B(0)

which injects two distinct artifacts: an additive trend b-1, AND a declining
gain b(t) on the response itself. Subtracting a trend -- polynomial, or any
high-pass -- fixes only the first. Dividing by a slowly varying envelope fixes
both, and matches the physics rather than approximating it.

WHY THE OWN TRACE. Two other envelope sources were considered and rejected as
the standard method:

  * across-ROI mean: relies on ROIs being decorrelated. True for 1D bar noise,
    FALSE for the full-field stimulus, where every ROI sees identical input and
    the common mode IS the signal. Disqualified as a uniform method.
  * background ROI: better in principle -- external information, so it costs no
    signal bandwidth -- but it must then exist in every dataset, drawn
    consistently, forever, and a method that silently degrades when one is
    missing cannot be applied consistently. Use it to VALIDATE this one
    instead: its slow envelope should match what this removes.

Choosing a source per dataset was rejected outright: it would let the
correction be picked for the filter it produces.

WHAT IT COSTS, AND WHY THAT IS RECOVERABLE. For small modulations
F / lowpass(F) ~ 1 + r - lowpass(r), so the correction high-passes the response
by a KNOWN transfer function,

    HP(f) = 1 - exp(-2 pi^2 f^2 sigma^2)

which `high_pass_response` returns so it can be divided back out. At sigma =
25 s that is 1.000 at 0.1 Hz, 0.993 at 0.02 Hz, 0.709 at 0.01 Hz. So the
low-frequency floor is ~0.01 Hz, not 1/T. Note 1/T is NOT a resolution limit
for a filter that decays to zero inside the window: such an h has a compactly
supported transform, so H(f) is exactly determined at every f and zero-padding
interpolates it exactly. Check `returns_to_zero` before trusting low
frequencies -- that is the condition, not the window length.

SIGMA is fixed by timescale separation, not tuned: bleaching is monotonic over
a 300 s epoch so its energy sits below ~0.01 Hz, while the filter band of
interest starts near 0.1 Hz. sigma = 25 s sits between them with >10x margin on
both sides. Gaussian rather than boxcar because a running mean has sinc nulls
-- blind spots at f = k/W where it removes nothing at all.

NOT FOR SHORT EPOCHS. Flashes, steps and moving patches have no room for the
window, and there the sustained level is the measurement rather than an
artifact. `correct` refuses rather than leaving that to the caller.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d

SIGMA_S = 25.0          # fixed; see module docstring
MIN_SPAN_FACTOR = 4.0   # epoch must be at least this many sigma long


def _interp_nans(y):
    """Fill NaNs by interpolation so the smoother does not spread them."""
    y = np.asarray(y, dtype=float).copy()
    bad = ~np.isfinite(y)
    if bad.all():
        return y, bad
    if bad.any():
        idx = np.arange(len(y))
        y[bad] = np.interp(idx[bad], idx[~bad], y[~bad])
    return y, bad


def envelope(F, dt, sigma_s=SIGMA_S):
    """Slow multiplicative envelope of a raw fluorescence trace.

    Edges use mode='nearest'; 'reflect' would fold the epoch's own trend back
    on itself and flatten the envelope exactly where bleaching is steepest.
    """
    y, bad = _interp_nans(F)
    env = gaussian_filter1d(y, sigma=float(sigma_s) / float(dt), mode='nearest')
    return env, bad


def correct(F, dt, sigma_s=SIGMA_S, min_span_factor=MIN_SPAN_FACTOR,
            allow_short=False):
    """Divide out the slow bleach envelope of a RAW fluorescence trace.

    Parameters
    ----------
    F : (n_t,) raw fluorescence. NOT dF/F -- the correction is multiplicative
        and must act before any baseline normalisation.
    dt : sample period (s)
    sigma_s : Gaussian sigma (s). Leave at the default; it is fixed by
        timescale separation and varying it per dataset defeats the point.
    allow_short : bypass the epoch-length guard. For tests only.

    Returns
    -------
    F_corrected : F / (env / env[0]), so the trace keeps its original scale
    env : the removed envelope
    info : dict, recorded so an analysis can be traced to its correction
    """
    F = np.asarray(F, dtype=float)
    span = len(F) * float(dt)
    if not allow_short and span < min_span_factor * sigma_s:
        raise ValueError(
            'epoch is {:.1f} s but the correction needs at least {:.0f} s '
            '({:g} x sigma={:g}s). Short-epoch protocols (flashes, steps, '
            'moving patches) should NOT be bleach-corrected: there the '
            'sustained level is the measurement, not an artifact.'
            .format(span, min_span_factor * sigma_s, min_span_factor, sigma_s))
    env, bad = envelope(F, dt, sigma_s)
    ref = env[np.isfinite(env)][0] if np.isfinite(env).any() else 1.0
    if not np.isfinite(ref) or ref == 0:
        raise ValueError('envelope reference is {!r}; trace is unusable'
                         .format(ref))
    out = F / (env / ref)
    out[bad] = np.nan
    return out, env, dict(method='gaussian_own_trace', sigma_s=float(sigma_s),
                          dt=float(dt), span_s=float(span),
                          n_nan=int(bad.sum()), env_ref=float(ref),
                          frac_removed=float(1.0 - env[-1] / ref))


def high_pass_response(f, sigma_s=SIGMA_S):
    """Transfer function this correction imposes on the response.

    Divide an estimated filter's spectrum by this to undo the attenuation. It
    falls to 0.709 at 0.01 Hz and 0.265 at 0.005 Hz for the default sigma, so
    inverting below ~0.01 Hz amplifies noise -- report a ceiling with it, the
    same way the sample-and-hold sinc^2 is handled.
    """
    f = np.asarray(f, dtype=float)
    return 1.0 - np.exp(-2.0 * np.pi ** 2 * f ** 2 * float(sigma_s) ** 2)


def roi_responses(ID, name, roi_prefix='aligned', dff='pre', sigma_s=SIGMA_S,
                  correct_bleach='auto', quiet=False):
    """Drop-in for ID.getRoiResponses with the bleach envelope divided out.

    The correction is multiplicative and must act on raw F BEFORE any baseline
    normalisation -- correcting dF/F would fix the additive trend and leave the
    declining gain -- so this asks for dff='none' and reproduces
    imaging_data.get_dff's own 'pre' and 'mean' conventions afterwards.

    correct_bleach : True | False | 'auto'
        'auto' (default) applies the correction when the epoch is long enough
        and skips it otherwise, so a short-epoch protocol degrades to the
        uncorrected result instead of raising. What happened is recorded in
        the returned dict under 'bleach_info'; check it rather than assuming.
    """
    rd = ID.getRoiResponses(name, roi_prefix=roi_prefix, dff='none')
    er = np.asarray(rd['epoch_response'], dtype=float)
    # sample_period, not diff(time_vector_by_epoch): those vectors carry
    # per-epoch acquisition jitter (they start at 0.026 / 0.072 / 0.032 s on
    # this series), and sample_period is what imaging_data itself uses.
    dt = float(ID.getResponseTiming()['sample_period'])
    span = er.shape[2] * dt

    do = bool(correct_bleach) and correct_bleach != 'auto'
    if correct_bleach == 'auto':
        do = span >= MIN_SPAN_FACTOR * sigma_s
    info = dict(applied=bool(do), sigma_s=float(sigma_s), span_s=float(span),
                reason='' if do else
                ('disabled by caller' if not correct_bleach else
                 'epoch {:.0f}s < {:.0f}s, too short to separate bleaching '
                 'from signal'.format(span, MIN_SPAN_FACTOR * sigma_s)))
    if do:
        frac = []
        for r in range(er.shape[0]):
            for e in range(er.shape[1]):
                try:
                    er[r, e], _, inf = correct(er[r, e], dt, sigma_s,
                                               allow_short=True)
                    frac.append(inf['frac_removed'])
                except ValueError:
                    pass
        info['frac_removed_mean'] = float(np.mean(frac)) if frac else float('nan')
        # ONLY epoch_response is corrected. 'roi_response' is the continuous
        # series trace, so correcting it would apply a DIFFERENT envelope --
        # whole-session decay rather than within-epoch -- and the two would not
        # be comparable: on Mi1 the per-ROI medians ended up a factor 1.22
        # apart. A caller wanting a baseline must take it from
        # epoch_response, not from roi_response, whenever this correction ran.
        rd['bleach_roi_response_corrected'] = False
    elif not quiet:
        print('    [bleach] skipped: {}'.format(info['reason']))

    if dff != 'none':
        # imaging_data builds time_vector as arange(n) * sample_period, i.e.
        # starting at 0 rather than at -pre_time, so the pre period CANNOT be
        # found as tv < 0. Match get_dff exactly: pre_frames is derived per
        # epoch from pre_time / sample_period.
        sp = float(ID.getResponseTiming()['sample_period'])
        pre_frames = np.array(
            [int(x) for x in np.asarray(ID.getEpochParameters('pre_time'),
                                        dtype=float) / sp])
        for e in range(er.shape[1]):
            blk = er[:, e, :]
            npre = int(pre_frames[e]) if e < len(pre_frames) else 0
            if dff == 'pre' and npre > 0:
                base = np.mean(blk[:, 0:npre], axis=1, keepdims=True)
            else:                                   # 'mean', or no pre period
                base = np.mean(blk, axis=1, keepdims=True)
            with np.errstate(invalid='ignore', divide='ignore'):
                er[:, e, :] = (blk - base) / base
    rd['epoch_response'] = er
    rd['bleach_info'] = info
    return rd


def returns_to_zero(h, dt, tail_frac=0.2, tol_frac=0.02):
    """Does an estimated filter decay to zero inside its window?

    This is the condition under which low frequencies are trustworthy -- not
    the window length. Returns (ok, tail_mean_as_fraction_of_peak).
    """
    h = np.asarray(h, dtype=float)
    n = max(1, int(round(tail_frac * len(h))))
    peak = np.max(np.abs(h))
    if peak <= 0:
        return False, np.nan
    tail = float(np.nanmean(h[-n:])) / peak
    return bool(abs(tail) < tol_frac), tail
