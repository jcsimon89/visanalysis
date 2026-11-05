"""Tests for visanalysis.model. Run from the repo root:
    python -m pytest tests/test_model.py -q
Needs clandinin_labpack beside the visanalysis checkout for the ContrastSwitch stimulus tests."""
import os

import numpy as np
import pytest

from visanalysis.model import spec, stimulus, measurement, readouts
from visanalysis.model.filters import TabulatedFilter, GammaDiffFilter, ExpFilter, bm_scale, filter_from_spec
from visanalysis.model.nonlinearities import LinearNL, TabulatedNL, ReluNL, nl_from_spec
from visanalysis.model.models import LN, SwitchedLN, model_from_spec, bm_normalize
from visanalysis.util import contrast_switch as cs


def _cs_stim(duration=60.0, seed=1):
    return stimulus.from_protocol("ContrastSwitch", dict(low_contrast=0.2, high_contrast=0.92, modulator_dwell=10,
                                                         probe_tau=0.05), duration=duration, seed=seed,
                                  pre_time=0.5, tail_time=0.5)


# ── filters ─────────────────────────────────────────────────────────────────

def test_tabulated_filter_matches_pipeline_taps():
    """A filter built from pipeline taps gives the pipeline's own A @ h at imaging times."""
    stim = stimulus.white_noise(20.0, 0.3, fps=120.0, seed=0)
    rng = np.random.default_rng(0)
    lag_dt, n_lag = 0.02, 50
    h = rng.standard_normal(n_lag) * np.exp(-np.arange(n_lag) / 10.0)
    f = TabulatedFilter.from_taps(h, lag_dt, output_units="dF/F")
    y = f.apply(stim.c, stim.dt)
    t = np.arange(2.0, 19.0, 0.1153) + 0.0001          # off the .5 ms rounding boundary
    st = dict(t0=stim.t[0], cs=np.r_[0.0, np.cumsum(stim.c)])
    A = cs.taps(st, t, n_lag, lag_dt)
    i = stim.index(t)
    assert np.array_equal(i, np.round((t - stim.t[0]) * cs.FS).astype(int))
    assert np.allclose(y[i], A @ h, atol=1e-9)


def test_parametric_filters_scale_and_spec():
    f = GammaDiffFilter(amp1=-1.0, n1=4, tau1=0.03, amp2=-0.3, n2=6, tau2=0.1, duration=3.0, output_units="contrast")
    k = f.kernel(0.001)
    assert k.size == 3001 and k[0] == 0.0
    assert abs(k.sum() * 0.001 - (-1.0 - -0.3)) < 0.02        # amplitudes are lobe integrals
    f2 = f.scaled(2.5)
    assert np.allclose(f2.kernel(0.001), 2.5 * k)
    f3 = filter_from_spec(f2.to_spec())
    assert np.allclose(f3.kernel(0.001), f2.kernel(0.001))
    e = ExpFilter(amp=1.0, tau=0.1)
    assert abs(e.kernel(0.0005).sum() * 0.0005 - 1.0) < 0.01


def test_bm_scale_makes_unit_variance_generator():
    stim = _cs_stim(40.0)
    f = GammaDiffFilter(amp1=-1.0, n1=4, tau1=0.03, duration=0.5, output_units="contrast")
    lab = stim.labels["contrast_state"]
    for v in (-1, 1):
        k = bm_scale(f, stim, lab == v)
        g = f.scaled(k).apply(stim.c, stim.dt)
        # rows whose whole history is in state v
        inside = np.convolve((lab == v).astype(float), np.ones(501), "full")[:lab.size] >= 501
        assert abs(g[inside].std() / stim.c[lab == v].std() - 1.0) < 1e-6


# ── nonlinearities ──────────────────────────────────────────────────────────

def test_tabulated_nl_extrapolation_and_range():
    nl = TabulatedNL([-1, 0, 1], [0, 1, 3], input_units="contrast", offset=0.5)
    assert nl.range == (-1.0, 1.0)
    assert np.allclose(nl([-1, 0, 1]), [0.5, 1.5, 3.5])
    assert np.isclose(nl(2.0), 3.5 + 2.0)                       # end slope continues
    assert nl.out_of_range_fraction([-2, 0, 0, 2]) == 0.5
    nl2 = TabulatedNL([-1, 0, 1], [0, 1, 3], extrapolation="clamp")
    assert nl2(5.0) == 3.0
    nl3 = nl_from_spec(nl.to_spec())
    assert np.allclose(nl3([-0.5, 0.7]), nl([-0.5, 0.7]))


def test_pair_units_must_match():
    f = ExpFilter(1.0, 0.1, output_units="contrast")
    with pytest.raises(ValueError):
        LN(f, LinearNL(gain=1.0, input_units="dF/F"))


# ── models ──────────────────────────────────────────────────────────────────

def test_switched_ln_selects_state_and_nan_elsewhere():
    stim = _cs_stim(30.0)
    f = ExpFilter(1.0, 0.05, output_units="contrast")
    m = SwitchedLN("contrast_state", {1: LN(f, LinearNL(1.0, offset=10.0), name="high"),
                                      -1: LN(f, LinearNL(1.0, offset=-10.0), name="low")}, name="sw")
    r = m.respond(stim)
    lab = stim.labels["contrast_state"]
    assert np.all(np.isnan(r.y[lab == 0]))
    assert np.all(r.y[lab == 1] > 5) and np.all(r.y[lab == -1] < -5)
    # the generator is the same filter on the same history, whichever state picks it
    g_high = f.apply(stim.c, stim.dt)
    assert np.allclose(r.g[lab != 0], g_high[lab != 0])


def test_recover_ln_returns_the_model_without_noise():
    stim = _cs_stim(240.0, seed=5)
    f = GammaDiffFilter(amp1=-1.0, n1=4, tau1=0.025, amp2=-0.3, n2=6, tau2=0.08, duration=1.0, output_units="contrast")
    base = SwitchedLN("contrast_state", {1: LN(f, LinearNL(0.1, offset=0.15), name="high"),
                                         -1: LN(f, LinearNL(0.25, offset=0.05), name="low")}, name="m")
    m = bm_normalize(base, stim)
    r = m.respond(stim)
    planes = measurement.sample(r, volume_period=0.11524, plane_offsets=np.arange(10) * 0.00549, noise_sd=0.0)
    rec = readouts.recover_ln(planes, stim, label="contrast_state", filter_length=1.0, lag_dt=0.02, ridge=1e-4)
    for v in (-1, 1):
        lag, d_true = m.states[v].filter.tabulate(0.02)
        _, d_rec = rec[v]["filter"].tabulate(0.02)
        assert np.corrcoef(d_true, d_rec)[0, 1] > 0.995
        assert rec[v]["r"] > 0.99
        slope_true = m.states[v].nonlinearity.gain
        assert abs(rec[v]["nonlinearity"].slope(-0.3, 0.3) / slope_true - 1.0) < 0.05
        assert abs(rec[v]["nonlinearity"].offset - m.states[v].nonlinearity.offset) < 0.01


def test_switch_triggered_average_sees_the_step():
    stim = _cs_stim(60.0)
    f = ExpFilter(1.0, 0.05, output_units="contrast")
    m = SwitchedLN("contrast_state", {1: LN(f, LinearNL(0.0, offset=1.0), name="high"),
                                      -1: LN(f, LinearNL(0.0, offset=0.0), name="low")})
    r = m.respond(stim)
    sta = readouts.switch_triggered_average(r, "contrast_state", pre=1, post=2, bin_s=0.5)
    up = sta[(-1, 1)]
    assert np.allclose(up["mean"][up["lag"] < 0], 0.0) and np.allclose(up["mean"][up["lag"] > 0], 1.0)
    assert (1, -1) in sta


def test_indicator_refuses_imaged_and_keeps_steps():
    stim = stimulus.step(10.0, [(2.0, 1.0)])
    f = ExpFilter(1.0, 0.02, output_units="dF/F")
    neural = LN(f, LinearNL(1.0, input_units="dF/F"), output=dict(signal="neural")).respond(stim)
    imaged = LN(f, LinearNL(1.0, input_units="dF/F")).respond(stim)
    ind = measurement.Indicator.preset("jRGECO1a")
    with pytest.raises(ValueError):
        ind.apply(imaged)
    out = ind.apply(neural)
    assert out.signal == "imaged"
    assert abs(out.y[-1] - neural.y[-1]) < 1e-3                 # unit-area kernel keeps the plateau


# ── spec / library ──────────────────────────────────────────────────────────

def test_spec_roundtrip_and_hash(tmp_path, monkeypatch):
    monkeypatch.setenv(spec.ENV_VAR, str(tmp_path))
    f = TabulatedFilter.from_taps(np.linspace(1, 0, 25), 0.02, output_units="contrast", normalization="BM2002")
    nl = TabulatedNL(np.linspace(-1, 1, 12), np.linspace(-1, 1, 12) ** 2, input_units="contrast", offset=0.1)
    m = SwitchedLN("contrast_state", {1: LN(f, nl, name="high"), -1: LN(f.scaled(2.0), nl, name="low")},
                   name="rt", output=dict(channel="calcium"), provenance=dict(source="hypothetical"))
    s = m.to_spec()
    p = spec.save_model(s, "grp/rt")
    assert p.exists() and p.with_suffix(".npz").exists()
    s2 = spec.load_model("grp/rt")
    m2 = model_from_spec(s2)
    stim = _cs_stim(20.0)
    assert np.allclose(np.nan_to_num(m.respond(stim).y), np.nan_to_num(m2.respond(stim).y))
    assert spec.model_hash(s) == spec.model_hash(m2.to_spec())
    assert spec.list_models() == ["grp/rt"]
    with pytest.raises(spec.LibraryError):
        spec.save_model(s, "grp/rt")                             # no silent overwrite


def test_library_resolution_order(tmp_path, monkeypatch):
    monkeypatch.delenv(spec.ENV_VAR, raising=False)
    monkeypatch.setattr(spec, "CONFIG_PATH", tmp_path / "nowhere.yaml")
    with pytest.raises(spec.LibraryError):
        spec.resolve_library()
    cfg = tmp_path / "config.yaml"
    lib_cfg = tmp_path / "from_cfg"
    lib_cfg.mkdir()
    cfg.write_text("model_library: {}\n".format(lib_cfg.as_posix()))
    monkeypatch.setattr(spec, "CONFIG_PATH", cfg)
    assert spec.resolve_library() == lib_cfg
    lib_env = tmp_path / "from_env"
    lib_env.mkdir()
    monkeypatch.setenv(spec.ENV_VAR, str(lib_env))
    assert spec.resolve_library() == lib_env
    lib_arg = tmp_path / "from_arg"
    assert spec.resolve_library(lib_arg, create=True) == lib_arg


def test_contrast_switch_export_builds_a_runnable_model(tmp_path, monkeypatch):
    """The exporter on a minimal stand-in for the pipeline's outputs gives a SwitchedLN that runs."""
    from visanalysis.model import export
    monkeypatch.setenv(spec.ENV_VAR, str(tmp_path))
    n_lag, lag_dt = 150, 0.02
    h = -np.exp(-np.arange(n_lag) / 5.0) * 0.01
    nl = np.c_[np.linspace(-1, 1, 12), np.linspace(-0.05, 0.05, 12), np.full(12, 0.002)]
    ch = dict(n_roi=10, n_rows=100, n_planes=2, ridge=1.0, filter_length=3.0, lag_dt=lag_dt, steady_start=3.0,
              first_dwell=2, cv_r_low=0.1, cv_r_high=0.2, null_max_low=0.0, null_max_high=0.0,
              has_filter_both_states=True, gain=dict(low=1, high=1, ratio=1.0, ci=[0.9, 1.1]))
    res = dict(channels=dict(calcium=ch), gate=dict(method="high"), settings=dict(tag="raw"),
               switch_triggered=dict(calcium=dict(to_low=dict(level_late=0.06), to_high=dict(level_late=0.15))))
    arrs = dict(calcium=dict(h_low=h, h_high=h, nl_low=nl, nl_high=nl, bm_scale=np.array([5.0, 1.2]),
                             sd_state=np.array([0.117, 0.55])))
    ctx = dict(metas=[dict(seed=1, fps=120.0)], seeds=[1], frozen=True, n_ep=2, stim_dur=[300.0], dwell_s=10.0,
               realized=(0.2, 0.92), contrast_ratio=4.6)
    p = export.contrast_switch_model(ctx, res, arrs, "calcium", "D:/data/JS256_x_JS261/fly_011", "series_002")
    assert p.name == "JS256_x_JS261_fly_011_series_002_calcium.yaml"
    m = model_from_spec(spec.load_model("contrast_switch/JS256_x_JS261_fly_011_series_002_calcium"))
    assert m.states[1].nonlinearity.offset == 0.15 and m.states[-1].nonlinearity.offset == 0.06
    assert np.isclose(m.states[-1].filter.density[0], h[0] * 5.0 / lag_dt)
    r = m.respond(_cs_stim(30.0))
    assert np.isfinite(r.y[r.labels["contrast_state"] != 0]).all()


def test_example_model_loads_and_runs():
    here = os.path.dirname(os.path.abspath(__file__))
    p = os.path.join(here, "..", "visanalysis", "model", "examples", "tm9_like_off_switched.yaml")
    m = model_from_spec(spec.load_model(p))
    stim = _cs_stim(30.0)
    r = bm_normalize(m, stim).respond(stim)
    assert np.isfinite(r.y[stim.labels["contrast_state"] != 0]).all()
