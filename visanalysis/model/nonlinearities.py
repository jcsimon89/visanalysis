"""Static nonlinearities: generator -> response.

A nonlinearity is the second half of a matched pair; its input_units must be
the filter's output_units. Every kind has an additive `offset`, the response
level the nonlinearity sits on. A pipeline's tabulated nonlinearity is
measured on a demeaned response, so the export puts the state's measured level
in `offset` rather than folding it into y: what was measured and what was
added stay separate.

Kinds
-----
tabulated  x (n,), y (n,), optional y_se (n,): binned means, interpolated
           linearly between points. Beyond the measured range the response is
           extrapolated by `extrapolation`: 'linear' continues the end slopes,
           'clamp' holds the end values, 'nan' refuses. The fraction of inputs
           outside the range is reported by out_of_range_fraction and by every
           Response, because a prediction that lives mostly outside the range
           the nonlinearity was measured on is a guess about the extrapolation,
           not a prediction of the model.
linear     gain * g
relu       gain * max(g - threshold, 0)
power      gain * max(g - threshold, 0) ** exponent
sigmoid    amplitude / (1 + exp(-(g - x0) / width))
"""
import numpy as np

__all__ = ["Nonlinearity", "TabulatedNL", "LinearNL", "ReluNL", "PowerNL", "SigmoidNL", "nl_from_spec"]


class Nonlinearity(object):
    kind = None

    def __init__(self, input_units="contrast", offset=0.0, name=None):
        self.input_units = str(input_units)
        self.offset = float(offset)
        self.name = name

    @property
    def range(self):
        """(lo, hi) of generator the nonlinearity is trusted on, or None."""
        return None

    def _f(self, g):
        raise NotImplementedError

    def _spec_body(self):
        raise NotImplementedError

    def __call__(self, g):
        g = np.asarray(g, dtype=float)
        return self._f(g) + self.offset

    def out_of_range_fraction(self, g):
        r = self.range
        if r is None:
            return 0.0
        g = np.asarray(g, dtype=float)
        g = g[np.isfinite(g)]
        return float(np.mean((g < r[0]) | (g > r[1]))) if g.size else 0.0

    def slope(self, lo, hi, n=200):
        """Mean slope over [lo, hi]: the Baccus & Meister 'sensitivity' on that generator range."""
        x = np.linspace(lo, hi, n)
        y = self(x)
        return float(np.polyfit(x, y, 1)[0])

    def to_spec(self):
        d = dict(kind=self.kind, input_units=self.input_units, offset=self.offset)
        if self.name:
            d["name"] = self.name
        d.update(self._spec_body())
        return d

    def __repr__(self):
        r = self.range
        rs = " on [{:.3g}, {:.3g}]".format(*r) if r else ""
        return "{}(offset {:+.3g}{}, <- {})".format(type(self).__name__, self.offset, rs, self.input_units)


class TabulatedNL(Nonlinearity):
    kind = "tabulated"

    def __init__(self, x, y, y_se=None, extrapolation="linear", **kw):
        Nonlinearity.__init__(self, **kw)
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        ok = np.isfinite(x) & np.isfinite(y)
        o = np.argsort(x[ok])
        self.x, self.y = x[ok][o], y[ok][o]
        self.y_se = None if y_se is None else np.asarray(y_se, dtype=float)[ok][o]
        if self.x.size < 2:
            raise ValueError("a tabulated nonlinearity needs at least two points")
        if extrapolation not in ("linear", "clamp", "nan"):
            raise ValueError("extrapolation must be linear, clamp or nan")
        self.extrapolation = extrapolation

    @property
    def range(self):
        return float(self.x[0]), float(self.x[-1])

    def _f(self, g):
        out = np.interp(g, self.x, self.y)
        lo, hi = g < self.x[0], g > self.x[-1]
        if self.extrapolation == "linear":
            s0 = (self.y[1] - self.y[0]) / (self.x[1] - self.x[0])
            s1 = (self.y[-1] - self.y[-2]) / (self.x[-1] - self.x[-2])
            out = np.where(lo, self.y[0] + s0 * (g - self.x[0]), out)
            out = np.where(hi, self.y[-1] + s1 * (g - self.x[-1]), out)
        elif self.extrapolation == "nan":
            out = np.where(lo | hi, np.nan, out)
        return out

    def _spec_body(self):
        d = dict(x=self.x, y=self.y, extrapolation=self.extrapolation)
        if self.y_se is not None:
            d["y_se"] = self.y_se
        return d


class LinearNL(Nonlinearity):
    kind = "linear"

    def __init__(self, gain=1.0, **kw):
        Nonlinearity.__init__(self, **kw)
        self.gain = float(gain)

    def _f(self, g):
        return self.gain * g

    def _spec_body(self):
        return dict(gain=self.gain)


class ReluNL(Nonlinearity):
    kind = "relu"

    def __init__(self, gain=1.0, threshold=0.0, **kw):
        Nonlinearity.__init__(self, **kw)
        self.gain, self.threshold = float(gain), float(threshold)

    def _f(self, g):
        return self.gain * np.clip(g - self.threshold, 0, None)

    def _spec_body(self):
        return dict(gain=self.gain, threshold=self.threshold)


class PowerNL(Nonlinearity):
    kind = "power"

    def __init__(self, gain=1.0, exponent=2.0, threshold=0.0, **kw):
        Nonlinearity.__init__(self, **kw)
        self.gain, self.exponent, self.threshold = float(gain), float(exponent), float(threshold)

    def _f(self, g):
        return self.gain * np.clip(g - self.threshold, 0, None) ** self.exponent

    def _spec_body(self):
        return dict(gain=self.gain, exponent=self.exponent, threshold=self.threshold)


class SigmoidNL(Nonlinearity):
    kind = "sigmoid"

    def __init__(self, amplitude=1.0, x0=0.0, width=1.0, **kw):
        Nonlinearity.__init__(self, **kw)
        self.amplitude, self.x0, self.width = float(amplitude), float(x0), float(width)

    def _f(self, g):
        return self.amplitude / (1.0 + np.exp(-(g - self.x0) / self.width))

    def _spec_body(self):
        return dict(amplitude=self.amplitude, x0=self.x0, width=self.width)


_KINDS = {c.kind: c for c in (TabulatedNL, LinearNL, ReluNL, PowerNL, SigmoidNL)}


def nl_from_spec(d):
    d = dict(d)
    kind = d.pop("kind")
    if kind not in _KINDS:
        raise ValueError("unknown nonlinearity kind {!r}; known: {}".format(kind, sorted(_KINDS)))
    return _KINDS[kind](**d)
