"""Space-time separability of a (cells x lags) filter, judged across trial folds.

A receptive field is separable when every position carries the same time
course up to a scale factor, i.e. the (cells x lags) block has rank one.
The singular values of a noisy estimate cannot say so directly: noise fills
every singular value, so an exactly separable filter at realistic SNR puts
only about half its energy in the first one, and the fraction measures SNR
rather than shape (see average_strf_octave.separability for that lesson).

Everything here therefore uses TWO estimates of the same filter with
independent noise -- the two trial folds the noise pipelines already compute --
and asks the question across them, so noise can never manufacture structure:

  cross-validated fraction (descriptive)
      Rank-1 of fold A, scored against fold B: <A_1, B> / <A, B>. Fold B's
      noise is independent of A, so the ratio estimates how much of the
      filter's energy one spatial profile times one time course accounts
      for. The rank-2 version is returned as well; its increment over rank 1
      is the energy a reproducible SECOND component carries. Both are
      averaged over the two directions. Biased LOW at low SNR (rank-1 of a
      noisy block captures less of the truth than rank-1 of the truth
      would): on synthetic separable filters the median is 0.98 at peak z
      5-7 and 1.00 above. Read it against z, or against ROIs of similar z.

  second-component test (the yes/no)
      Take fold A's SECOND spatial profile, project out the full block's
      rank-1 profile (an imperfect fold profile leaves some of the true
      profile in the second singular vector, which would then correlate
      with the other fold's main component), and project fold B onto it:
      one time course, which for a separable filter is pure noise and for
      anything else holds the second component. Do the same with the folds
      swapped, sum the two traces sign-aligned (independent noise, shared
      signal), and take the peak |z| of the sum (nullcal.peak_z: late-lag
      noise scale, peak in the lag window). The null is the SAME statistic
      on the shift-null filters the ROI selection already computes, via
      nullcal.calibrate's extra_stats hook, and ROIs are called inseparable
      by BH at q = 0.05 over the kept set.

      Measured on synthetic 1D-noise data through the real estimator path
      (C:/Users/Public/flymax_analysis/separability_test, 2026-11-18):
      false positives on planted separable filters 0-1.2% of kept ROIs;
      detection of a centre plus a delayed, wider surround 20-50% overall,
      0 below peak z 7, 20-75% at z 10-15 and 75-90% above z 15, the lower
      numbers for a surround lagging the centre by 0.1 s (smeared time
      courses correlated 0.81), the higher for 0.2 s (0.40). Alternatives
      measured and rejected there: a per-cell residual test (fold B with
      fold A's rank-1 profile removed, both directions required) had zero
      false positives and under 10% power; the un-orthogonalised sum had the
      same power as this one but 13-14% false positives above z 15; the
      rank-2 increment against a no-signal null never fired.

Convolution along the lag axis -- the stimulus autocorrelation, the indicator
kernel -- acts identically at every cell, so it preserves rank. These
statistics are computed on the RAW STA, never the deconvolved one, and the
hold tau cannot create or hide inseparability.

Geometry-agnostic: blocks are (n_cells, n_lags) with any cell axis, so the 1D
and octave pipelines share it.
"""
import numpy as np

from . import nullcal


def window(block, n_win):
    """The leading n_win lags of a (..., n_lags) block; the whole block when
    n_win is None."""
    return block if n_win is None else block[..., :int(n_win)]


def rank_k(block_w, k):
    """Rank-k truncation of a (n_cells, n_win) block via SVD.

    Returns (approx, profiles, time_courses, singular_values) with profiles
    (k, n_cells) unit vectors, each oriented so its dominant cell is positive,
    and time_courses (k, n_win) carrying the singular value and the polarity.
    """
    u, s, vt = np.linalg.svd(np.asarray(block_w, dtype=float), full_matrices=False)
    k = int(min(k, s.size))
    prof = u[:, :k].T.copy()
    tc = (s[:k, None] * vt[:k]).copy()
    for j in range(k):
        sgn = np.sign(prof[j][int(np.argmax(np.abs(prof[j])))]) or 1.0
        prof[j] *= sgn
        tc[j] *= sgn
    approx = prof.T @ tc
    return approx, prof, tc, s


def _corr(a, b):
    a = a - a.mean()
    b = b - b.mean()
    d = np.linalg.norm(a) * np.linalg.norm(b)
    return float(a @ b / d) if d > 0 else np.nan


def cv_fractions(a_w, b_w, max_rank=2):
    """<A_k, B> / <A, B> for k = 1..max_rank, averaged over both directions.

    Returns (fractions (max_rank,), cv_corr). NaN when the folds do not agree
    on a filter at all (<A, B> <= 0), where the ratio has no meaning.
    """
    num = np.zeros(max_rank)
    den = 0.0
    for x, y in ((a_w, b_w), (b_w, a_w)):
        _, prof, tc, _ = rank_k(x, max_rank)
        den += float(np.sum(x * y))
        acc = np.zeros_like(x)
        for k in range(prof.shape[0]):
            acc = acc + np.outer(prof[k], tc[k])
            num[k] += float(np.sum(acc * y))
        for k in range(prof.shape[0], max_rank):
            num[k] = num[prof.shape[0] - 1]
    cv_corr = _corr(a_w.ravel(), b_w.ravel())
    if not den > 0:
        return np.full(max_rank, np.nan), cv_corr
    return num / den, cv_corr


def profiles(block, n_win, k=1):
    """Spatial profile k (0-based) of each ROI's filter, from the SVD of its
    lag window. block (n_roi, n_cells, n_lags) -> (n_roi, n_cells), zero
    where the block is empty or has fewer than k + 1 components."""
    block = np.asarray(block, dtype=float)
    n_roi, n_cells = block.shape[:2]
    out = np.zeros((n_roi, n_cells))
    for r in range(n_roi):
        w = np.nan_to_num(window(block[r], n_win))
        if not np.any(w):
            continue
        prof = rank_k(w, k + 1)[1]
        if prof.shape[0] > k:
            out[r] = prof[k]
    return out


def fold_profiles(fold_block, n_win):
    """Rank-1 spatial profile per ROI; kept for callers of the first version."""
    return profiles(fold_block, n_win, 0)


def orthogonalise(vec, against):
    """Rows of vec with the matching rows of against projected out, renormalised."""
    v = np.asarray(vec, dtype=float) - np.asarray(against, dtype=float) * \
        np.einsum('rc,rc->r', vec, against)[:, None]
    n = np.linalg.norm(v, axis=1)
    n[n == 0] = 1.0
    return v / n[:, None]


def test_vectors(fold_a, fold_b, full, n_win):
    """The two unit spatial vectors the second-component test projects onto:
    each fold's second profile with the full block's rank-1 profile removed.
    Returns (v_a, v_b), each (n_roi, n_cells)."""
    p1 = profiles(full, n_win, 0)
    return (orthogonalise(profiles(fold_a, n_win, 1), p1),
            orthogonalise(profiles(fold_b, n_win, 1), p1))


def second_trace(filt_for_a, filt_for_b, v_a, v_b):
    """Sign-aligned sum of v_a . filt_for_a and v_b . filt_for_b, per ROI:
    (n_roi, n_lags). For the real statistic filt_for_a is fold B (A's vector
    tests B) and filt_for_b is fold A; for a null draw both are the same
    shifted filter, which is equivalent because peak_z is self-normalising
    and the lag correlation of the noise is the same either way."""
    ta = np.einsum('rc,rcl->rl', v_a, filt_for_a)
    tb = np.einsum('rc,rcl->rl', v_b, filt_for_b)
    sgn = np.sign(np.einsum('rc,rc->r', v_a, v_b))
    sgn[sgn == 0] = 1.0
    return ta + sgn[:, None] * tb


def second_peak_z(filt_for_a, filt_for_b, v_a, v_b, noise_frac=0.3, peak_lags=None):
    """The statistic: peak |z| of the summed second-component trace, plus the
    two directional traces' peak |z| as (2, n_roi) diagnostics."""
    z_sum = nullcal.peak_z(second_trace(filt_for_a, filt_for_b, v_a, v_b),
                           noise_frac=noise_frac, peak_lags=peak_lags)
    z_dir = np.stack([
        nullcal.peak_z(np.einsum('rc,rcl->rl', v_a, filt_for_a), noise_frac, peak_lags),
        nullcal.peak_z(np.einsum('rc,rcl->rl', v_b, filt_for_b), noise_frac, peak_lags)])
    return z_sum, z_dir


def extra_stats(v_a, v_b, noise_frac=0.3, peak_lags=None):
    """The second-component statistics as callables for nullcal.calibrate's
    extra_stats: each takes one shifted (n_roi, n_cells, n_lags) filter."""
    def _sum(filt):
        return second_peak_z(filt, filt, v_a, v_b, noise_frac, peak_lags)[0]

    def _ab(filt):
        return nullcal.peak_z(np.einsum('rc,rcl->rl', v_a, filt), noise_frac, peak_lags)

    def _ba(filt):
        return nullcal.peak_z(np.einsum('rc,rcl->rl', v_b, filt), noise_frac, peak_lags)
    return {'sep_second': _sum, 'sep_second_ab': _ab, 'sep_second_ba': _ba}


def analyze(fold_a, fold_b, full, keep, n_win, v_a=None, v_b=None, null=None,
            noise_frac=0.3, q=0.05):
    """Separability of every kept ROI from its two fold filters.

    fold_a, fold_b, full : (n_roi, n_cells, n_lags), raw STAs over the SAME
                           cells (the visible ones)
    keep                 : indices of the ROIs to test
    n_win                : lags in the window (nullcal.peak_lags_for) or None
    v_a, v_b             : test_vectors(...) if already computed
    null                 : calibrate(...)['null_extra'] holding 'sep_second'
                           (and optionally the directional draws); None
                           skips the test and leaves p NaN

    Returns a dict of per-ROI arrays, NaN outside keep:
      frac         (n_roi, 2)  cross-validated fractions for rank 1 and 2
      cv_corr      (n_roi,)    correlation of the two folds in the window
      second_z     (n_roi,)    peak |z| of the summed second-component trace
      second_p     (n_roi,)    its p against the null
      second_z_dir (2, n_roi)  the two directional traces' peak |z|
      second_p_dir (2, n_roi)  their p (NaN when the null lacks them)
      inseparable  (n_roi,)    1 where BH (q) passes, 0 tested and not,
                               -1 untested
      second_trace (n_roi, n_win)  the summed trace that was tested
      test_vector  (2, n_roi, n_cells)  v_a, v_b
      energy_frac  (n_roi, 2)  s_k^2 / sum s^2 of the full block, descriptive
      profile      (2, n_roi, n_cells)  first two spatial profiles of the full block
      time_course  (2, n_roi, n_win)    their time courses (window lags)
    plus 'bh_cutoff' and 'n_tested'.
    """
    fold_a = np.asarray(fold_a, dtype=float)
    fold_b = np.asarray(fold_b, dtype=float)
    full = np.asarray(full, dtype=float)
    n_roi, n_cells, n_lags = full.shape
    n_w = n_lags if n_win is None else int(n_win)
    keep = np.asarray(keep, dtype=int) if keep is not None else np.array([], int)
    if v_a is None or v_b is None:
        v_a, v_b = test_vectors(fold_a, fold_b, full, n_win)

    out = dict(frac=np.full((n_roi, 2), np.nan), cv_corr=np.full(n_roi, np.nan),
               second_z=np.full(n_roi, np.nan), second_p=np.full(n_roi, np.nan),
               second_z_dir=np.full((2, n_roi), np.nan),
               second_p_dir=np.full((2, n_roi), np.nan),
               inseparable=np.full(n_roi, -1, dtype=int),
               second_trace=np.full((n_roi, n_w), np.nan),
               test_vector=np.stack([v_a, v_b]),
               energy_frac=np.full((n_roi, 2), np.nan),
               profile=np.full((2, n_roi, n_cells), np.nan),
               time_course=np.full((2, n_roi, n_w), np.nan),
               bh_cutoff=np.nan, n_tested=int(keep.size))
    if keep.size == 0:
        return out

    for r in keep:
        a_w = np.nan_to_num(window(fold_a[r], n_win))
        b_w = np.nan_to_num(window(fold_b[r], n_win))
        out['frac'][r], out['cv_corr'][r] = cv_fractions(a_w, b_w, 2)
        _, prof, tc, s = rank_k(np.nan_to_num(window(full[r], n_win)), 2)
        tot = float(np.sum(s ** 2))
        if tot > 0:
            out['energy_frac'][r, :s[:2].size] = (s[:2] ** 2) / tot
        out['profile'][:prof.shape[0], r] = prof
        out['time_course'][:tc.shape[0], r] = tc

    # A's vector tests fold B and B's tests fold A
    fa_k, fb_k = np.nan_to_num(fold_a[keep]), np.nan_to_num(fold_b[keep])
    z_sum, z_dir = second_peak_z(fb_k, fa_k, v_a[keep], v_b[keep], noise_frac, n_win)
    out['second_z'][keep] = z_sum
    out['second_z_dir'][:, keep] = z_dir
    out['second_trace'][keep] = window(second_trace(fb_k, fa_k, v_a[keep], v_b[keep]),
                                       n_win)

    if null is not None and np.asarray(null.get('sep_second', [])).size:
        ok = np.isfinite(z_sum)
        p = np.full(z_sum.size, np.nan)
        p[ok] = nullcal.pvalues(z_sum[ok], null['sep_second'])[0]
        out['second_p'][keep] = p
        idx, cut = nullcal.benjamini_hochberg(np.where(ok, p, 1.0), q)
        out['bh_cutoff'] = float(cut)
        flag = np.zeros(keep.size, dtype=int)
        flag[idx] = 1
        out['inseparable'][keep] = flag
        for d, name in enumerate(('sep_second_ab', 'sep_second_ba')):
            draws = np.asarray(null.get(name, []))
            if draws.size:
                okd = np.isfinite(z_dir[d])
                pd = np.full(z_dir.shape[1], np.nan)
                pd[okd] = nullcal.pvalues(z_dir[d][okd], draws)[0]
                out['second_p_dir'][d, keep] = pd
    return out
