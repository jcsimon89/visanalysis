"""From a model's output to what the microscope would record.

Indicator   a calcium (or cAMP) indicator's impulse response, unit area so a
            sustained step keeps its size. Applied only to a NEURAL response;
            a model exported from imaging data is already AS IMAGED and
            applying one would count the indicator twice, so apply() refuses.
sample      the fine-grid trace at imaging times -- a volume rate and the z
            offset of each plane, or explicit times from a recording -- with
            additive Gaussian noise. Each plane is one Response (its own
            times), as the pipelines treat them.

Indicator presets hold half-rise / half-decay times and their source; they are
converted to exponential time constants (t_half / ln 2) for a difference of
exponentials. They are starting points for a hypothetical model, not
measurements on this rig.
"""
import numpy as np
from scipy.signal import fftconvolve

from visanalysis.model.models import Response

__all__ = ["Indicator", "INDICATORS", "sample", "plane_times"]

LN2 = np.log(2.0)

# name: (half-rise s, half-decay s, source)
INDICATORS = {
    "jRGECO1a": (0.025, 0.250, "2 AP, jcsimon's numbers (2026-10-28)"),
    "GCaMP8m": (0.0074, 0.134, "1 AP, cultured neurons at 37 C (Zhang et al. 2023); no fly measurement"),
    "GCaMP6m": (0.080, 0.300, "approximate, 1 AP in vivo (Chen et al. 2013)"),
}


class Indicator(object):
    def __init__(self, tau_rise, tau_decay, name=None, duration=None):
        self.tau_rise, self.tau_decay = float(tau_rise), float(tau_decay)
        self.name = name
        self.duration = float(duration) if duration else self.tau_decay * 8.0 + self.tau_rise * 8.0

    @classmethod
    def preset(cls, name):
        if name not in INDICATORS:
            raise KeyError("unknown indicator {!r}; known: {}".format(name, sorted(INDICATORS)))
        hr, hd, _ = INDICATORS[name]
        return cls(hr / LN2, hd / LN2, name=name)

    def kernel(self, dt):
        t = np.arange(int(np.ceil(self.duration / dt)) + 1) * dt
        k = np.exp(-t / self.tau_decay) - np.exp(-t / self.tau_rise) if self.tau_rise > 0 else np.exp(-t / self.tau_decay)
        s = k.sum() * dt
        return k / s if s > 0 else k

    def apply(self, resp, force=False):
        """The response as the indicator would report it. Refuses an 'imaged' response unless forced."""
        if resp.signal == "imaged" and not force:
            raise ValueError("{!r} is already as imaged; applying {} would count the indicator twice "
                             "(force=True if you mean it)".format(resp.name, self.name or "an indicator"))
        if resp.dt is None:
            raise ValueError("apply the indicator before sampling: the response has no fine grid")
        y = np.nan_to_num(resp.y)
        out = fftconvolve(y, self.kernel(resp.dt), mode="full")[:y.size] * resp.dt
        out[np.isnan(resp.y)] = np.nan
        return resp.copy(y=out, signal="imaged", valid_from=resp.valid_from + self.duration,
                         notes=resp.notes + ["indicator {} (rise {:.0f} ms, decay {:.0f} ms)".format(
                             self.name or "", self.tau_rise * 1e3, self.tau_decay * 1e3)])

    def __repr__(self):
        return "Indicator({}, tau rise {:.3f} s, decay {:.3f} s)".format(self.name or "", self.tau_rise, self.tau_decay)


def plane_times(t_start, t_end, volume_period, plane_offsets=(0.0,), scan_time=0.0):
    """Sample times per plane: volume starts every volume_period, plane p at offset[p] + scan/2."""
    starts = np.arange(t_start, t_end, volume_period)
    return [starts + float(o) + scan_time / 2.0 for o in plane_offsets]


def sample(resp, times=None, volume_period=None, plane_offsets=(0.0,), scan_time=0.0, noise_sd=0.0,
           seed=0, t_start=None, t_end=None):
    """Sample a fine-grid Response at imaging times. Returns a list of Responses, one per plane.

    Give either `times` (a list of arrays, one per plane) or `volume_period` with plane_offsets.
    The sample is the mean of the fine-grid response over the fine step containing the time."""
    if resp.dt is None:
        raise ValueError("sample() needs a fine-grid response")
    rng = np.random.default_rng(seed)
    if times is None:
        if volume_period is None:
            raise ValueError("give times or volume_period")
        t0 = resp.t[0] if t_start is None else t_start
        t1 = resp.t[-1] if t_end is None else t_end
        times = plane_times(t0, t1, volume_period, plane_offsets, scan_time)
    out = []
    for p, tp in enumerate(times):
        tp = np.asarray(tp, dtype=float)
        tp = tp[(tp >= resp.t[0]) & (tp <= resp.t[-1])]
        i = np.clip(np.round((tp - resp.t[0]) / resp.dt).astype(np.int64), 0, resp.t.size - 1)
        y = resp.y[i] + (noise_sd * rng.standard_normal(tp.size) if noise_sd > 0 else 0.0)
        labels = {k: np.asarray(v)[i] for k, v in resp.labels.items()}
        out.append(Response(tp, y, dt=None, g=None if resp.g is None else resp.g[i], labels=labels,
                            name=resp.name, signal=resp.signal, units=resp.units, valid_from=resp.valid_from,
                            out_of_range=resp.out_of_range, plane=p,
                            notes=resp.notes + ["sampled, noise sd {:g}".format(noise_sd)]))
    return out
