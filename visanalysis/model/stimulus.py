"""Stimuli as a model sees them: contrast on a fine grid, with labels.

A Stimulus is Weber contrast c(t) on a grid of step DT (1 ms, the grid the
noise pipelines hold the displayed stimulus on), plus LABELS: named arrays
aligned with t that say what else was true at each sample -- which contrast
state the display was in, which dwell, how long since it changed. Models key
on labels by name, so a stimulus source is free to invent labels and a model
written against one stimulus works on any other that carries the same label.

Sources
-------
from_protocol(protocol, params, duration, ...)
    The stimulus a protocol WOULD show, before any fly is run, built through
    the labpack generator from protocol-style parameters. Only the parameters
    given are set; the rest are the generator's own defaults.
from_recording(ID)
    One Stimulus per epoch of a recorded series: the frames the display
    actually showed, timed from the photodiode (first edge = stimulus time 0),
    with dropped frames held, exactly as the analysis pipelines rebuild it.
step, sinusoid, white_noise, ternary_renewal
    Synthetic signals for sanity checks and hypothetical experiments.

Both protocol sources go through a registry keyed on the protocol name, so a
new protocol is one entry here and nothing in the models.
"""
import numpy as np

from visanalysis.util import contrast_switch as _cs

DT = 1.0 / _cs.FS

__all__ = ["DT", "Stimulus", "from_protocol", "from_recording", "step", "sinusoid", "white_noise",
           "ternary_renewal", "PROTOCOLS"]


class Stimulus(object):
    def __init__(self, t, c, labels=None, name="", meta=None, mean_intensity=None):
        self.t = np.asarray(t, dtype=float)
        self.c = np.asarray(c, dtype=float)
        if self.t.shape != self.c.shape or self.t.ndim != 1:
            raise ValueError("t and c must be the same 1-d shape")
        # rounded so a 1 ms grid is exactly 0.001 and index() rounds as the pipelines do (t * FS)
        self.dt = float(np.round(np.median(np.diff(self.t)), 12)) if self.t.size > 1 else DT
        self.labels = {}
        for k, v in (labels or {}).items():
            v = np.asarray(v)
            if v.shape != self.t.shape:
                raise ValueError("label {!r} has shape {}, stimulus {}".format(k, v.shape, self.t.shape))
            self.labels[k] = v
        self.name = name
        self.meta = dict(meta or {})
        self.mean_intensity = mean_intensity

    @property
    def n(self):
        return self.t.size

    @property
    def duration(self):
        return float(self.t[-1] - self.t[0] + self.dt)

    def index(self, times):
        """Nearest sample index for each time (clipped)."""
        return np.clip(np.round((np.asarray(times, float) - self.t[0]) / self.dt).astype(np.int64), 0, self.n - 1)

    def label_at(self, name, times):
        return self.labels[name][self.index(times)]

    def transitions(self, label):
        """[(index, from_value, to_value)] where the label changes value."""
        v = np.asarray(self.labels[label])
        if v.dtype.kind == "f":
            same = np.isclose(v[1:], v[:-1]) | (np.isnan(v[1:]) & np.isnan(v[:-1]))
        else:
            same = v[1:] == v[:-1]
        idx = np.flatnonzero(~same) + 1
        return [(int(i), v[i - 1].item(), v[i].item()) for i in idx]

    def window(self, t0, t1):
        m = (self.t >= t0) & (self.t < t1)
        return Stimulus(self.t[m], self.c[m], {k: v[m] for k, v in self.labels.items()}, self.name, self.meta,
                        self.mean_intensity)

    def rms_contrast(self, label=None):
        """RMS contrast overall, or per value of a label."""
        if label is None:
            return float(np.sqrt(np.mean(self.c ** 2)))
        v = self.labels[label]
        return {u.item(): float(np.sqrt(np.mean(self.c[v == u] ** 2))) for u in np.unique(v)}

    def __repr__(self):
        return "Stimulus({!r}, {:.1f} s at {:.0f} Hz grid, labels {})".format(
            self.name, self.duration, 1 / self.dt, sorted(self.labels))


# ── synthetic ───────────────────────────────────────────────────────────────

def _grid(duration, t0=0.0):
    n = int(round(duration / DT))
    return t0 + np.arange(n) * DT


def step(duration, steps, t0=0.0, name="step"):
    """steps: [(time, contrast)] -- the contrast from each time on."""
    t = _grid(duration, t0)
    c = np.zeros_like(t)
    for at, v in sorted(steps):
        c[t >= at] = v
    return Stimulus(t, c, dict(stim_time=t - t0), name)


def sinusoid(duration, freq_hz, amplitude, t0=0.0, phase=0.0, name="sinusoid"):
    t = _grid(duration, t0)
    return Stimulus(t, amplitude * np.sin(2 * np.pi * freq_hz * (t - t0) + phase), dict(stim_time=t - t0), name)


def white_noise(duration, rms, fps=120.0, t0=0.0, seed=0, name="white_noise"):
    """Gaussian contrast, a new value every display frame, held between frames."""
    rng = np.random.default_rng(seed)
    t = _grid(duration, t0)
    k = np.floor((t - t0) * fps + 1e-9).astype(int)
    vals = rng.standard_normal(k.max() + 1) * rms
    return Stimulus(t, vals[k], dict(stim_time=t - t0), name)


def ternary_renewal(duration, contrast, tau, fps=120.0, p=1.0 / 6, t0=0.0, seed=0, name="ternary"):
    """Ternary noise (-c, 0, +c) redrawn at each frame with probability 1 - exp(-1/(tau fps))."""
    rng = np.random.default_rng(seed)
    t = _grid(duration, t0)
    k = np.floor((t - t0) * fps + 1e-9).astype(int)
    n = k.max() + 1
    q = 1.0 - np.exp(-1.0 / (tau * fps)) if tau > 0 else 1.0
    redraw = rng.random(n) < q
    redraw[0] = True
    draws = rng.choice([-1.0, 0.0, 1.0], n, p=[p, 1 - 2 * p, p])
    idx = np.maximum.accumulate(np.where(redraw, np.arange(n), 0))
    return Stimulus(t, contrast * draws[idx][k], dict(stim_time=t - t0), name)


# ── protocols ───────────────────────────────────────────────────────────────

_CS_KEYS = {"centre_code", "modulator_schedule", "modulator_dwell", "modulator_start", "modulator_tau",
            "probe_tau", "probe_p", "fps"}


def _cs_generator(params, seed):
    from visanalysis.util.octave_stim import _find_labpack
    _find_labpack()
    from labpack.visual_stim.clandinin.contrast_switch_ternary import ContrastSwitchTernaryStimulus
    p = dict(params or {})
    kw = {k: p.pop(k) for k in list(p) if k in _CS_KEYS}
    if "contrasts" in p:
        kw["contrasts"] = tuple(float(x) for x in p.pop("contrasts"))
    elif "low_contrast" in p or "high_contrast" in p:
        g0 = ContrastSwitchTernaryStimulus(seed=0)
        lo, hi = g0.contrasts if hasattr(g0, "contrasts") else g0.realized_contrasts
        kw["contrasts"] = (float(p.pop("low_contrast", lo)), float(p.pop("high_contrast", hi)))
    if p:
        raise KeyError("unknown ContrastSwitch parameters {}".format(sorted(p)))
    return ContrastSwitchTernaryStimulus(seed=int(seed), **kw)


def _from_epoch_stimulus(st, name, meta):
    """Wrap contrast_switch.epoch_stimulus's dict as a Stimulus with the standard labels."""
    t = st["t0"] + np.arange(st["c"].size) * DT
    labels = dict(contrast_state=st["state"].astype(int), dwell=st["dwell"].astype(int),
                  since_switch=st["since"], stim_time=np.where(np.isfinite(st["pos"]), t, np.nan))
    return Stimulus(t, st["c"], labels, name, meta)


def _cs_from_protocol(params, duration, seed=0, pre_time=0.0, tail_time=0.0, name=None):
    g = _cs_generator(params, seed)
    meta = g.as_metadata()
    st = _cs.epoch_stimulus(meta, float(duration), float(pre_time), float(tail_time))
    name = name or "ContrastSwitch c={:.2f}/{:.2f} dwell {:g} s tau {:g} s seed {}".format(
        *g.realized_contrasts, g.modulator_dwell if g.fixed_modulator else g.modulator_tau, g.probe_tau, seed)
    return _from_epoch_stimulus(st, name, dict(protocol="ContrastSwitchTernaryNoise", generator=meta, seed=int(seed),
                                               realized_contrasts=list(g.realized_contrasts),
                                               rms_contrasts=list(g.rms_contrasts), idle_intensity=g.idle_intensity,
                                               pre_time=pre_time, tail_time=tail_time))


def _cs_from_recording(ID):
    _, _, sample_rate = ID.getVoltageData()
    timing = ID.getStimulusTiming()
    starts = np.asarray(timing["stimulus_start_times"], dtype=float)
    ends = np.asarray(timing["stimulus_end_times"], dtype=float)
    drops = np.asarray(timing["dropped_frame_times"], dtype=float) / sample_rate
    eps = ID.getEpochParameters()
    n_ep = min(len(starts), len(eps))
    pre = np.asarray(ID.getEpochParameters("pre_time"), dtype=float)[:n_ep]
    tail = np.asarray(ID.getEpochParameters("tail_time"), dtype=float)[:n_ep]
    out, cache = [], {}
    for i in range(n_ep):
        meta = _cs.resolve_meta(eps[i])
        d = float(ends[i] - starts[i])
        dd = drops[(drops >= starts[i]) & (drops <= ends[i])]
        di = np.round((dd - starts[i]) * float(meta["fps"])).astype(int)
        st = _cs.epoch_stimulus(meta, d, pre[i], tail[i], di, cache)
        out.append(_from_epoch_stimulus(st, "epoch {} seed {}".format(i, meta["seed"]),
                                        dict(protocol="ContrastSwitchTernaryNoise", epoch=i, generator=meta,
                                             seed=int(meta["seed"]), stim_dur=d, n_dropped=int(dd.size),
                                             realized_contrasts=list(st["realized"]), pre_time=float(pre[i]),
                                             tail_time=float(tail[i]), start_time=float(starts[i]))))
    return out


PROTOCOLS = {
    "ContrastSwitch": dict(from_protocol=_cs_from_protocol, from_recording=_cs_from_recording,
                           match=lambda pid: "ContrastSwitch" in str(pid)),
}


def _lookup(protocol):
    for k, v in PROTOCOLS.items():
        if k == protocol or v["match"](protocol):
            return v
    raise KeyError("no stimulus source for protocol {!r}; known: {}".format(protocol, sorted(PROTOCOLS)))


def from_protocol(protocol, params=None, duration=60.0, seed=0, pre_time=0.0, tail_time=0.0, name=None):
    """The stimulus a protocol would show for `duration` s with these protocol parameters."""
    return _lookup(protocol)["from_protocol"](params or {}, duration, seed=seed, pre_time=pre_time,
                                              tail_time=tail_time, name=name)


def from_recording(ID):
    """One Stimulus per epoch of the series an ImagingDataObject points at."""
    pid = ID.getRunParameters().get("protocol_ID", "")
    return _lookup(pid)["from_recording"](ID)
