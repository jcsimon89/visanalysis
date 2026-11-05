"""Models and their responses.

Model.respond(stimulus) -> Response. A Response carries the fine-grid trace,
the generator(s) that produced it, the stimulus labels (so readouts can key on
them without the stimulus), and a record of how much of the prediction ran
outside the range any nonlinearity was measured on.

LN          one filter, one nonlinearity.
SwitchedLN  a label name and an LN per label value. At every sample the LN of
            the CURRENT label value is applied to the WHOLE stimulus history,
            including the part before the label last changed. That is the
            simplest model of adaptation -- every parameter changes the instant
            the label does -- and it is a hypothesis, not a convenience: the
            switch-triggered average is where it will be wrong if it is wrong,
            and the models that would improve on it (state parameters that
            relax after a switch; a feedback stage) belong beside it here.
            Samples whose label value has no LN (the grey, say) are NaN unless a
            default LN is given.

Specs: to_spec() / model_from_spec() round-trip through spec.py's YAML + npz.
"""
import numpy as np

from visanalysis.model.filters import filter_from_spec
from visanalysis.model.nonlinearities import nl_from_spec

__all__ = ["Response", "Model", "LN", "SwitchedLN", "model_from_spec"]


class Response(object):
    """A predicted (or measured) trace with its context.

    t          sample times (s); evenly spaced with step dt on the fine grid, or irregular
               after measurement.sample (dt None)
    y          the response
    g          the generator (fine grid) or None
    labels     dict name -> array aligned with t
    signal     'imaged' or 'neural'
    valid_from time before which the filter's history window reaches into unknown stimulus
    out_of_range  dict state name -> fraction of samples whose generator was outside the
               nonlinearity's measured range (0 when it has none)
    """

    def __init__(self, t, y, dt=None, g=None, labels=None, name="", signal="imaged", units="dF/F",
                 valid_from=None, out_of_range=None, components=None, plane=None, notes=None):
        self.t = np.asarray(t, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.dt = dt
        self.g = None if g is None else np.asarray(g, dtype=float)
        self.labels = {} if labels is None else dict(labels)
        self.name = name
        self.signal = signal
        self.units = units
        self.valid_from = self.t[0] if valid_from is None else float(valid_from)
        self.out_of_range = {} if out_of_range is None else dict(out_of_range)
        self.components = {} if components is None else dict(components)   # name -> (g, y) per state
        self.plane = plane
        self.notes = [] if notes is None else list(notes)

    @property
    def valid(self):
        return self.t >= self.valid_from

    def copy(self, **kw):
        d = dict(t=self.t, y=self.y, dt=self.dt, g=self.g, labels=self.labels, name=self.name,
                 signal=self.signal, units=self.units, valid_from=self.valid_from,
                 out_of_range=self.out_of_range, components=self.components, plane=self.plane,
                 notes=self.notes)
        d.update(kw)
        return Response(**d)

    def __repr__(self):
        oor = ", ".join("{} {:.0%} outside NL range".format(k, v) for k, v in self.out_of_range.items() if v > 0)
        return "Response({!r}, {} samples, {} .. {} s, {}{})".format(
            self.name, self.t.size, self.t[0], self.t[-1], self.signal, "; " + oor if oor else "")


class Model(object):
    kind = None

    def __init__(self, name="", output=None, provenance=None):
        self.name = name
        self.output = dict(channel=None, signal="imaged", units="dF/F", indicator=None)
        if output:
            self.output.update(output)
        self.provenance = dict(provenance or {})

    def respond(self, stim):
        raise NotImplementedError

    def _spec_body(self):
        raise NotImplementedError

    def to_spec(self):
        d = dict(spec_version=1, name=self.name, kind=self.kind, output=dict(self.output),
                 provenance=dict(self.provenance))
        d.update(self._spec_body())
        return d

    def describe(self):
        return repr(self)


def _check_pair(f, nl, where):
    if f.output_units != nl.input_units:
        raise ValueError("{}: filter outputs {!r} but the nonlinearity takes {!r}; they are not a pair".format(
            where, f.output_units, nl.input_units))


class LN(Model):
    kind = "LN"

    def __init__(self, filter, nonlinearity, name="", **kw):
        Model.__init__(self, name=name, **kw)
        _check_pair(filter, nonlinearity, name or "LN")
        self.filter = filter
        self.nonlinearity = nonlinearity

    def generator(self, stim):
        return self.filter.apply(stim.c, stim.dt)

    def respond(self, stim):
        g = self.generator(stim)
        y = self.nonlinearity(g)
        return Response(stim.t, y, dt=stim.dt, g=g, labels=stim.labels, name=self.name,
                        signal=self.output["signal"], units=self.output["units"],
                        valid_from=stim.t[0] + self.filter.length,
                        out_of_range={self.name or "LN": self.nonlinearity.out_of_range_fraction(g)})

    def _spec_body(self):
        return dict(filter=self.filter.to_spec(), nonlinearity=self.nonlinearity.to_spec())

    def __repr__(self):
        return "LN({!r}: {!r} -> {!r})".format(self.name, self.filter, self.nonlinearity)


class SwitchedLN(Model):
    kind = "SwitchedLN"

    def __init__(self, label, states, default=None, name="", **kw):
        """states: dict label value -> LN (or list of (value, LN)). default: LN for other values."""
        Model.__init__(self, name=name, **kw)
        self.label = str(label)
        items = states.items() if isinstance(states, dict) else list(states)
        self.states = {}
        for v, ln in items:
            if not isinstance(ln, LN):
                raise TypeError("each state of a SwitchedLN is an LN")
            self.states[_key(v)] = ln
        self.default = default

    def respond(self, stim):
        if self.label not in stim.labels:
            raise KeyError("stimulus {!r} has no label {!r} (it has {})".format(
                stim.name, self.label, sorted(stim.labels)))
        lab = np.asarray(stim.labels[self.label])
        y = np.full(stim.n, np.nan)
        g = np.full(stim.n, np.nan)
        comps, oor, length = {}, {}, 0.0
        used = np.zeros(stim.n, bool)
        for v, ln in self.states.items():
            m = _match(lab, v)
            gi = ln.generator(stim)
            yi = ln.nonlinearity(gi)
            comps[ln.name or str(v)] = (gi, yi)
            y[m], g[m] = yi[m], gi[m]
            used |= m
            oor[ln.name or str(v)] = ln.nonlinearity.out_of_range_fraction(gi[m]) if m.any() else 0.0
            length = max(length, ln.filter.length)
        if self.default is not None and not used.all():
            m = ~used
            gi = self.default.generator(stim)
            y[m], g[m] = self.default.nonlinearity(gi)[m], gi[m]
            length = max(length, self.default.filter.length)
        return Response(stim.t, y, dt=stim.dt, g=g, labels=stim.labels, name=self.name,
                        signal=self.output["signal"], units=self.output["units"],
                        valid_from=stim.t[0] + length, out_of_range=oor, components=comps)

    def _spec_body(self):
        return dict(label=self.label,
                    states=[dict(value=_plain(v), name=ln.name, filter=ln.filter.to_spec(),
                                 nonlinearity=ln.nonlinearity.to_spec()) for v, ln in self.states.items()],
                    default=None if self.default is None else dict(
                        name=self.default.name, filter=self.default.filter.to_spec(),
                        nonlinearity=self.default.nonlinearity.to_spec()))

    def __repr__(self):
        return "SwitchedLN({!r} on {!r}: {})".format(
            self.name, self.label, "; ".join("{}={!r}".format(v, ln) for v, ln in self.states.items()))


def bm_normalize(model, stim, label=None):
    """A copy of the model whose filters are Baccus & Meister scaled on this stimulus.

    For a SwitchedLN each state's filter is scaled on the samples in its own state (and only
    where the whole history window is inside that state), so var(F * s) = var(s) there; for an LN
    the scale is taken over the whole stimulus (or the samples where `label` is non-zero). Use it
    on a hypothetical model whose nonlinearity is meant in the units of a measured one; a model
    exported from data is already scaled and must not be scaled again."""
    from visanalysis.model.filters import bm_scale
    if isinstance(model, SwitchedLN):
        lab = np.asarray(stim.labels[model.label])
        states = {}
        for v, ln in model.states.items():
            k = bm_scale(ln.filter, stim, _match(lab, v))
            states[v] = LN(ln.filter.scaled(k, normalization="BM2002"), ln.nonlinearity, name=ln.name)
        return SwitchedLN(model.label, states, default=model.default, name=model.name, output=model.output,
                          provenance=dict(model.provenance, bm_normalized_on=stim.name))
    if isinstance(model, LN):
        mask = None if label is None else np.asarray(stim.labels[label]) != 0
        k = bm_scale(model.filter, stim, mask)
        return LN(model.filter.scaled(k, normalization="BM2002"), model.nonlinearity, name=model.name,
                  output=model.output, provenance=dict(model.provenance, bm_normalized_on=stim.name))
    raise TypeError("bm_normalize knows LN and SwitchedLN")


def _key(v):
    """Label values as dict keys: integers stay integers, integral floats become integers."""
    if isinstance(v, str):
        try:
            v = float(v)
        except ValueError:
            return v
    if isinstance(v, (float, np.floating)) and float(v).is_integer():
        return int(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    return v


def _plain(v):
    return v.item() if isinstance(v, np.generic) else v


def _match(lab, v):
    if isinstance(v, str):
        return lab.astype(str) == v
    return np.isclose(lab.astype(float), float(v))


def _ln_from(d, name_default=""):
    return LN(filter_from_spec(d["filter"]), nl_from_spec(d["nonlinearity"]), name=d.get("name") or name_default)


def model_from_spec(spec):
    kind = spec.get("kind")
    common = dict(name=spec.get("name", ""), output=spec.get("output"), provenance=spec.get("provenance"))
    if kind == "LN":
        m = LN(filter_from_spec(spec["filter"]), nl_from_spec(spec["nonlinearity"]), **common)
    elif kind == "SwitchedLN":
        states = [(s["value"], _ln_from(s, str(s["value"]))) for s in spec["states"]]
        default = spec.get("default")
        m = SwitchedLN(spec["label"], states, default=None if not default else _ln_from(default, "default"),
                       **common)
    else:
        raise ValueError("unknown model kind {!r}".format(kind))
    return m
