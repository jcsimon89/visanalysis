"""Linear filters, stored as densities on lag intervals.

A filter here is response per unit contrast PER SECOND of stimulus. Applying it
to a contrast signal c on a grid of step dt is

    y[n] = dt * sum_m k[m] * c[n - m],        k[m] = density at lag m * dt,

so the same filter gives the same response on any grid fine enough to resolve
it. The pipelines fit filters "per tap of lag_dt" (each tap the mean contrast
over one lag interval); that is density * lag_dt, and the export divides it
out. A tabulated filter keeps the intervals it was fitted on: interval j is
(edge_j, edge_{j+1}], the same half-open convention as contrast_switch.taps, so
a model built from a pipeline filter and the pipeline itself agree to the
sample.

Kinds
-----
tabulated   lag_edges (n + 1,), density (n,)
gamma_diff  a difference of two gamma-function lobes, each a unit-area pdf
            scaled by an amplitude (the lobe's integral, response per unit
            contrast): amp1, n1, tau1, amp2, n2, tau2, delay, duration.
            Monophasic when amp2 = 0.
exp         a single exponential, amp * exp(-(t - delay) / tau) / tau.

Every kind has a `scale` multiplying its density. output_units says what the
filter's output is: 'contrast' for a Baccus & Meister filter (the generator is
in units of stimulus contrast and the nonlinearity holds the gain), or the
response's own units for a raw filter followed by the identity.
bm_scale() gives the scale that puts a filter in the Baccus & Meister
convention on a given stimulus: var(F * s) = var(s). A pipeline export applies
it once; a hypothetical filter can be given it so its nonlinearity is in the
same units as a measured one.
"""
import numpy as np
from scipy.signal import fftconvolve
from scipy.special import gammaln

__all__ = ["Filter", "TabulatedFilter", "GammaDiffFilter", "ExpFilter", "filter_from_spec", "bm_scale"]


class Filter(object):
    kind = None

    def __init__(self, output_units="contrast", normalization=None, name=None, scale=1.0):
        self.output_units = str(output_units)
        self.normalization = normalization
        self.name = name
        self.scale = float(scale)

    # -- what subclasses provide ------------------------------------------------
    @property
    def length(self):
        """Support in seconds."""
        raise NotImplementedError

    def _kernel(self, dt):
        """Unscaled density at lags m * dt for m = 0 .. ceil(length / dt). k[0] is the lag-0 point."""
        raise NotImplementedError

    def _spec_body(self):
        raise NotImplementedError

    # -- shared -----------------------------------------------------------------
    def kernel(self, dt):
        return self.scale * self._kernel(dt)

    def apply(self, c, dt):
        """Causal convolution of contrast c (fine grid, step dt) with the filter."""
        c = np.asarray(c, dtype=float)
        k = self.kernel(dt)
        if k.size == 0:
            return np.zeros_like(c)
        return fftconvolve(c, k, mode="full")[:c.size] * dt

    def scaled(self, factor, output_units=None, normalization=None):
        """A copy with the density multiplied by factor (and, optionally, new units)."""
        d = self.to_spec()
        d["scale"] = self.scale * float(factor)
        if output_units is not None:
            d["output_units"] = output_units
        if normalization is not None:
            d["normalization"] = normalization
        return filter_from_spec(d)

    def tabulate(self, lag_dt):
        """(lag_mid, density) on intervals of lag_dt covering the support: for plots and comparison."""
        n = int(np.ceil(self.length / lag_dt - 1e-9))
        edges = np.arange(n + 1) * lag_dt
        fine = max(lag_dt / 20.0, 1e-4)
        k = self.kernel(fine)
        lag = np.arange(k.size) * fine
        j = np.clip(np.searchsorted(edges, lag - 1e-9, side="left") - 1, -1, n - 1)
        dens = np.array([k[j == i].mean() if np.any(j == i) else 0.0 for i in range(n)])
        return 0.5 * (edges[:-1] + edges[1:]), dens

    def peak(self, dt=0.001):
        """(lag, value) of the largest |density|."""
        k = self.kernel(dt)
        i = int(np.argmax(np.abs(k)))
        return i * dt, float(k[i])

    def to_spec(self):
        d = dict(kind=self.kind, output_units=self.output_units)
        if self.normalization:
            d["normalization"] = self.normalization
        if self.name:
            d["name"] = self.name
        if self.scale != 1.0:
            d["scale"] = self.scale
        d.update(self._spec_body())
        return d

    def __repr__(self):
        lag, v = self.peak()
        return "{}({} s, peak {:+.3g} at {:.0f} ms, -> {})".format(
            type(self).__name__, self.length, v, lag * 1e3, self.output_units)


class TabulatedFilter(Filter):
    kind = "tabulated"

    def __init__(self, lag_edges, density, **kw):
        Filter.__init__(self, **kw)
        self.lag_edges = np.asarray(lag_edges, dtype=float)
        self.density = np.asarray(density, dtype=float)
        if self.lag_edges.ndim != 1 or self.density.shape != (self.lag_edges.size - 1,):
            raise ValueError("lag_edges must be (n+1,) and density (n,); got {} and {}".format(
                self.lag_edges.shape, self.density.shape))
        if np.any(np.diff(self.lag_edges) <= 0) or self.lag_edges[0] < 0:
            raise ValueError("lag_edges must be increasing and non-negative")

    @classmethod
    def from_taps(cls, h, lag_dt, scale=1.0, **kw):
        """From a pipeline filter h (response per unit contrast per tap of lag_dt), times scale."""
        h = np.asarray(h, dtype=float)
        return cls(np.arange(h.size + 1) * lag_dt, h * scale / lag_dt, **kw)

    @property
    def length(self):
        return float(self.lag_edges[-1])

    def _kernel(self, dt):
        m = int(np.ceil(self.length / dt - 1e-9)) + 1
        lag = np.arange(m) * dt
        # interval j is (edge_j, edge_{j+1}]; the epsilon keeps a lag that is an edge up to
        # rounding (m * dt against j * lag_dt) in the lower interval, as exact arithmetic would
        j = np.searchsorted(self.lag_edges, lag - 1e-9, side="left") - 1
        return np.where((j >= 0) & (j < self.density.size), self.density[np.clip(j, 0, self.density.size - 1)], 0.0)

    def _spec_body(self):
        return dict(lag_edges=self.lag_edges, density=self.density)


def _gamma_pdf(t, n, tau):
    t = np.asarray(t, dtype=float)
    out = np.zeros_like(t)
    ok = t > 0
    out[ok] = np.exp((n - 1) * np.log(t[ok] / tau) - t[ok] / tau - gammaln(n)) / tau
    return out


class GammaDiffFilter(Filter):
    kind = "gamma_diff"

    def __init__(self, amp1, n1, tau1, amp2=0.0, n2=2.0, tau2=0.1, delay=0.0, duration=None, **kw):
        Filter.__init__(self, **kw)
        self.amp1, self.n1, self.tau1 = float(amp1), float(n1), float(tau1)
        self.amp2, self.n2, self.tau2 = float(amp2), float(n2), float(tau2)
        self.delay = float(delay)
        if duration is None:
            duration = self.delay + 8.0 * max(self.n1 * self.tau1, self.n2 * self.tau2 if self.amp2 else 0.0)
        self.duration = float(duration)

    @property
    def length(self):
        return self.duration

    def _kernel(self, dt):
        m = int(np.ceil(self.duration / dt - 1e-9)) + 1
        t = np.arange(m) * dt - self.delay
        return self.amp1 * _gamma_pdf(t, self.n1, self.tau1) - self.amp2 * _gamma_pdf(t, self.n2, self.tau2)

    def _spec_body(self):
        return dict(amp1=self.amp1, n1=self.n1, tau1=self.tau1, amp2=self.amp2, n2=self.n2, tau2=self.tau2,
                    delay=self.delay, duration=self.duration)


class ExpFilter(Filter):
    kind = "exp"

    def __init__(self, amp, tau, delay=0.0, duration=None, **kw):
        Filter.__init__(self, **kw)
        self.amp, self.tau, self.delay = float(amp), float(tau), float(delay)
        self.duration = float(self.delay + 8.0 * self.tau if duration is None else duration)

    @property
    def length(self):
        return self.duration

    def _kernel(self, dt):
        m = int(np.ceil(self.duration / dt - 1e-9)) + 1
        t = np.arange(m) * dt - self.delay
        return np.where(t > 0, self.amp * np.exp(-np.clip(t, 0, None) / self.tau) / self.tau, 0.0)

    def _spec_body(self):
        return dict(amp=self.amp, tau=self.tau, delay=self.delay, duration=self.duration)


_KINDS = {c.kind: c for c in (TabulatedFilter, GammaDiffFilter, ExpFilter)}


def filter_from_spec(d):
    d = dict(d)
    kind = d.pop("kind")
    if kind not in _KINDS:
        raise ValueError("unknown filter kind {!r}; known: {}".format(kind, sorted(_KINDS)))
    return _KINDS[kind](**d)


def bm_scale(filt, stim, mask=None):
    """The factor that puts filt in the Baccus & Meister convention on this stimulus:
    var(F * s) = var(s) over the samples in mask (all valid samples by default).

    The filtered stimulus is compared where the whole history window is inside the stimulus
    (and inside mask), so a mask selecting one contrast state measures that state alone."""
    g = filt.apply(stim.c, stim.dt)
    ok = stim.t >= stim.t[0] + filt.length
    if mask is not None:
        m = np.asarray(mask, bool)
        # the row and its whole history must be in mask
        n_hist = int(np.ceil(filt.length / stim.dt))
        inside = np.convolve(m.astype(float), np.ones(n_hist + 1), mode="full")[:m.size] >= n_hist + 1
        ok &= inside
        sd_s = stim.c[m].std()
    else:
        sd_s = stim.c[ok].std()
    sd_g = g[ok].std()
    return float(sd_s / sd_g) if sd_g > 0 else np.nan
