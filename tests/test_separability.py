"""Unit tests for visanalysis.util.separability on planted (cells x lags) blocks.

Fast, block-level checks. The end-to-end calibration through the STA
estimator, the shift null and realistic SNR lives in
C:/Users/Public/flymax_analysis/separability_test/ (validate_separability.py
for the module's statistic, validate_variants.py for the alternatives it was
chosen over).
"""
import numpy as np

from visanalysis.util import separability as S
from visanalysis.util import nullcal

N_CELLS, N_LAGS, N_WIN = 15, 240, 120
T = np.arange(N_LAGS) / 120.0


def _kernel():
    return np.exp(-T / 0.18) - np.exp(-T / 0.05)


def _centre():
    return np.exp(-0.5 * ((np.arange(N_CELLS) - 7) / 1.5) ** 2)


def separable_block():
    return np.outer(_centre(), _kernel())


def rank2_block(surround=0.5):
    u2 = np.exp(-0.5 * ((np.arange(N_CELLS) - 7) / 4.0) ** 2) - 0.6 * _centre()
    g2 = np.exp(-(T - 0.2) / 0.3) * (T > 0.2)
    return separable_block() - surround * np.outer(u2, g2)


def noisy_pair(block, sd, rng):
    return (block + sd * rng.standard_normal(block.shape),
            block + sd * rng.standard_normal(block.shape))


def test_rank_k_orientation_and_reconstruction():
    blk = separable_block()[:, :N_WIN]
    approx, prof, tc, s = S.rank_k(blk, 1)
    assert np.allclose(approx, blk, atol=1e-10)
    assert prof[0][int(np.argmax(np.abs(prof[0])))] > 0
    assert np.isclose(np.linalg.norm(prof[0]), 1.0)


def test_separable_fraction_is_one_at_high_snr():
    rng = np.random.default_rng(1)
    a, b = noisy_pair(separable_block(), 0.02, rng)
    frac, r = S.cv_fractions(a[:, :N_WIN], b[:, :N_WIN])
    assert abs(frac[0] - 1.0) < 0.03
    assert abs(frac[1] - frac[0]) < 0.02
    assert r > 0.9


def test_rank2_fraction_below_one_and_increment_positive():
    rng = np.random.default_rng(2)
    a, b = noisy_pair(rank2_block(), 0.02, rng)
    frac, _ = S.cv_fractions(a[:, :N_WIN], b[:, :N_WIN])
    assert frac[0] < 0.85
    assert frac[1] - frac[0] > 0.1


def test_fraction_shape_without_signal():
    rng = np.random.default_rng(3)
    a = rng.standard_normal((N_CELLS, N_WIN))
    b = rng.standard_normal((N_CELLS, N_WIN))
    frac, _ = S.cv_fractions(a, b)
    assert frac.shape == (2,)


def test_profiles_orthogonal_and_test_vectors_unit():
    rng = np.random.default_rng(4)
    a, b = noisy_pair(rank2_block(), 0.03, rng)
    full = 0.5 * (a + b)
    p1 = S.profiles(full[None], N_WIN, 0)
    p2 = S.profiles(full[None], N_WIN, 1)
    assert abs(p1[0] @ p2[0]) < 1e-9
    va, vb = S.test_vectors(a[None], b[None], full[None], N_WIN)
    assert np.isclose(np.linalg.norm(va[0]), 1.0) and np.isclose(np.linalg.norm(vb[0]), 1.0)
    assert abs(va[0] @ p1[0]) < 1e-9 and abs(vb[0] @ p1[0]) < 1e-9


def test_second_z_separates_rank1_from_rank2():
    rng = np.random.default_rng(5)
    sd = 0.03
    a1, b1 = noisy_pair(separable_block(), sd, rng)
    a2, b2 = noisy_pair(rank2_block(), sd, rng)
    for (a, b), expect_big in (((a1, b1), False), ((a2, b2), True)):
        full = 0.5 * (a + b)
        va, vb = S.test_vectors(a[None], b[None], full[None], N_WIN)
        z, zdir = S.second_peak_z(b[None], a[None], va, vb, peak_lags=N_WIN)
        if expect_big:
            assert z[0] > 8.0
        else:
            assert z[0] < 6.0
        assert zdir.shape == (2, 1)


def test_calibrate_extra_stats_and_analyze_flags():
    rng = np.random.default_rng(6)
    sd = 0.03
    blocks = [separable_block()] * 4 + [rank2_block()] * 4
    a = np.stack([blk + sd * rng.standard_normal(blk.shape) for blk in blocks])
    b = np.stack([blk + sd * rng.standard_normal(blk.shape) for blk in blocks])
    full = 0.5 * (a + b)
    va, vb = S.test_vectors(a, b, full, N_WIN)

    def shifted(_s):
        return sd * rng.standard_normal(a.shape)
    cal = nullcal.calibrate(shifted, list(range(1, 41)), peak_lags=N_WIN,
                            extra_stats=S.extra_stats(va, vb, peak_lags=N_WIN))
    assert set(cal['null_extra']) == {'sep_second', 'sep_second_ab', 'sep_second_ba'}
    assert cal['null_extra']['sep_second'].size == 40 * 8
    out = S.analyze(a, b, full, np.arange(8), N_WIN, va, vb, null=cal['null_extra'])
    assert out['inseparable'].tolist() == [0, 0, 0, 0, 1, 1, 1, 1]
    assert np.all(out['frac'][:4, 0] > 0.9)
    assert np.all(out['frac'][4:, 0] < 0.9)
    assert np.all(out['energy_frac'][:, 0] > out['energy_frac'][:, 1])
    assert np.isfinite(out['second_p_dir']).all()
    assert out['second_trace'].shape == (8, N_WIN)


def test_analyze_untested_rois_are_nan():
    rng = np.random.default_rng(7)
    a = rng.standard_normal((3, N_CELLS, N_LAGS))
    b = rng.standard_normal((3, N_CELLS, N_LAGS))
    out = S.analyze(a, b, a, [1], N_WIN)
    assert np.isnan(out['frac'][0]).all() and np.isnan(out['frac'][2]).all()
    assert out['inseparable'].tolist() == [-1, -1, -1]   # no null given: untested
    assert np.isnan(out['second_p']).all()
    assert np.isfinite(out['second_z'][1])


def test_analyze_empty_keep():
    a = np.zeros((2, N_CELLS, N_LAGS))
    out = S.analyze(a, a, a, [], N_WIN)
    assert out['n_tested'] == 0 and np.isnan(out['frac']).all()


def test_rank_preserved_by_lag_convolution():
    # the smear from the hold / indicator acts identically at every cell
    blk = separable_block()
    R = np.exp(-np.abs(np.arange(N_LAGS)[:, None] - np.arange(N_LAGS)[None, :]) / 12.0)
    sm = blk @ R
    s = np.linalg.svd(sm[:, :N_WIN], compute_uv=False)
    assert s[1] / s[0] < 1e-8
