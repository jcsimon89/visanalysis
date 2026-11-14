"""Idealise a measured LN model into a few parameters: a toy model.

fit_gamma_diff     a difference of gamma lobes fitted to a tabulated filter's
                   density from `min_lag` on. The first taps of a measured
                   filter carry the display light leak, which the pipelines
                   leave in on purpose; fitting from 40 ms replaces them with
                   the lobe's own rise.
fit_nonlinearity   every parametric family (linear, relu, power, sigmoid)
                   fitted by weighted least squares to the binned
                   nonlinearity; ranked by AICc so the extra parameters have
                   to earn their place on 12 points. Returns all candidates.
idealize_switched  both, for every state of a SwitchedLN, with the fit
                   residuals and the rejected candidates kept in provenance.
"""
import numpy as np
from scipy.optimize import least_squares

from visanalysis.model.filters import GammaDiffFilter, _gamma_pdf
from visanalysis.model.nonlinearities import LinearNL, ReluNL, PowerNL, SigmoidNL
from visanalysis.model.models import LN, SwitchedLN

__all__ = ["fit_gamma_diff", "fit_nonlinearity", "idealize_switched", "NL_FAMILIES"]


# ── filters ─────────────────────────────────────────────────────────────────

def _gd(p, t):
    amp1, n1, tau1, amp2, n2, tau2 = p
    return amp1 * _gamma_pdf(t, n1, tau1) - amp2 * _gamma_pdf(t, n2, tau2)


def fit_gamma_diff(filt, min_lag=0.04, max_lag=None, lag_dt=None, biphasic=True):
    """Fit a GammaDiffFilter to a tabulated filter's density on lags >= min_lag.

    Returns (GammaDiffFilter, info) with info = dict(r, rel_rms, lag, data, fit, params)."""
    lag_dt = lag_dt or float(np.min(np.diff(filt.lag_edges)))
    lag, d = filt.tabulate(lag_dt)
    max_lag = filt.length if max_lag is None else max_lag
    m = (lag >= min_lag) & (lag <= max_lag)
    t, y = lag[m], d[m]
    i = int(np.argmax(np.abs(y)))
    s = np.sign(y[i]) or 1.0
    n1 = 4.0
    tau1 = max(lag[m][i] / (n1 - 1), 0.005)
    amp1 = s * np.abs(y[i]) * tau1 * 2.5
    p0 = [amp1, n1, tau1, 0.3 * amp1 if biphasic else 0.0, 6.0, 2.5 * tau1]
    # n >= 2 on both lobes: a gamma lobe with n < 2 is finite (n = 1) or infinite (n < 1) at
    # lag 0, and a filter that does not start from zero would put a tap on the sample itself
    lo = [-np.inf, 2.0, 0.003, -np.inf, 2.0, 0.003]
    hi = [np.inf, 30.0, 2.0, np.inf, 30.0, 3.0]
    if not biphasic:
        lo[3], hi[3] = -1e-12, 1e-12

    def resid(p):
        return _gd(p, t) - y
    best = None
    for scale2 in ((1.0, 2.5), (0.5, 4.0), (0.1, 1.5)):
        p_init = list(p0)
        if biphasic:
            p_init[3] = scale2[0] * 0.3 * amp1
            p_init[5] = scale2[1] * tau1
        try:
            r = least_squares(resid, p_init, bounds=(lo, hi), max_nfev=4000)
        except ValueError:
            continue
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    f = GammaDiffFilter(*p, delay=0.0, duration=filt.length, output_units=filt.output_units,
                        normalization=filt.normalization, name=filt.name)
    fit = _gd(p, t)
    rr = float(np.corrcoef(fit, y)[0, 1]) if y.std() > 0 else np.nan
    info = dict(r=rr, rel_rms=float(np.sqrt(np.mean((fit - y) ** 2)) / (np.abs(y).max() or 1.0)),
                lag=lag, data=d, fit=_gd(p, lag), min_lag=min_lag,
                params=dict(zip(("amp1", "n1", "tau1", "amp2", "n2", "tau2"), [float(v) for v in p])),
                peak_lag=float(f.peak()[0]))
    return f, info


# ── nonlinearities ──────────────────────────────────────────────────────────

def _families(x, y):
    """name -> (constructor(params) -> NL, p0, lower, upper). Parameters end with the offset."""
    sx, sy = (x.max() - x.min()) or 1.0, (y.max() - y.min()) or 1.0
    slope = np.polyfit(x, y, 1)[0]
    return {
        "linear": (lambda p: LinearNL(gain=p[0], offset=p[1]), [slope, y.mean()],
                   [-np.inf, -np.inf], [np.inf, np.inf]),
        "relu": (lambda p: ReluNL(gain=p[0], threshold=p[1], offset=p[2]), [slope, x.min(), y.min()],
                 [-np.inf, x.min() - sx, -np.inf], [np.inf, x.max(), np.inf]),
        "power": (lambda p: PowerNL(gain=p[0], exponent=p[1], threshold=p[2], offset=p[3]),
                  [slope, 1.5, x.min() - 0.1 * sx, y.min()],
                  [-np.inf, 0.3, x.min() - 2 * sx, -np.inf], [np.inf, 5.0, x.max(), np.inf]),
        "sigmoid": (lambda p: SigmoidNL(amplitude=p[0], x0=p[1], width=p[2], offset=p[3]),
                    [sy * 1.5, x.mean(), sx / 4, y.min() - 0.25 * sy],
                    [0.0, x.min() - sx, sx / 100, -np.inf], [np.inf, x.max() + sx, 10 * sx, np.inf]),
    }


NL_FAMILIES = ("linear", "relu", "power", "sigmoid")


def fit_nonlinearity(x, y, y_se=None, families=NL_FAMILIES, input_units="contrast", name=None, se_floor=None):
    """Fit each family to binned (x, y +- y_se). Returns a list of dicts sorted by AICc (best first):
    dict(family, nl, params, rss, wrss, aicc, k, r). The chosen NL keeps the measured x range."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if y_se is None or not np.any(np.isfinite(y_se)):
        w = np.ones_like(y)
    else:
        se = np.asarray(y_se, float)[ok]
        floor = np.nanmedian(se) * 0.5 if se_floor is None else se_floor
        w = 1.0 / np.maximum(np.nan_to_num(se, nan=np.nanmedian(se)), floor)
    n = x.size
    out = []
    for fam in families:
        make, p0, lo, hi = _families(x, y)[fam]
        k = len(p0)

        def resid(p):
            return (make(p)(x) - y) * w
        try:
            r = least_squares(resid, p0, bounds=(lo, hi), max_nfev=4000)
        except ValueError:
            continue
        nl = make(r.x)
        nl.input_units, nl.name = input_units, name
        nl._range = (float(x.min()), float(x.max()))
        pred = nl(x)
        rss = float(np.sum((pred - y) ** 2))
        wrss = float(np.sum(((pred - y) * w) ** 2))
        aicc = n * np.log(max(wrss, 1e-30) / n) + 2 * k + (2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else np.inf)
        out.append(dict(family=fam, nl=nl, params={kk: float(v) for kk, v in nl.to_spec().items()
                                                   if isinstance(v, (int, float))},
                        rss=rss, wrss=wrss, aicc=float(aicc), k=k,
                        r=float(np.corrcoef(pred, y)[0, 1]) if y.std() > 0 else np.nan))
    out.sort(key=lambda d: d["aicc"])
    return out


# ── whole models ────────────────────────────────────────────────────────────

def idealize_switched(model, min_lag=0.04, families=NL_FAMILIES, nl_family=None, biphasic=True, name=None):
    """A SwitchedLN whose every state is a GammaDiffFilter and a parametric nonlinearity.

    nl_family: force one family (for all states), else the AICc winner per state. Returns
    (toy SwitchedLN, report) with report[state] = dict(filter=info, nl_candidates=[...], chosen)."""
    states, report = {}, {}
    for v, ln in model.states.items():
        f, finfo = fit_gamma_diff(ln.filter, min_lag=min_lag, biphasic=biphasic)
        x, y, se = ln.nonlinearity.total()
        cands = fit_nonlinearity(x, y, se, families=families, input_units=ln.nonlinearity.input_units, name=ln.name)
        chosen = next((c for c in cands if c["family"] == nl_family), cands[0]) if nl_family else cands[0]
        states[v] = LN(f, chosen["nl"], name=ln.name)
        report[v] = dict(filter=finfo, nl_candidates=cands, chosen=chosen["family"],
                         measured=dict(x=x, y=y, se=se))
    toy = SwitchedLN(model.label, states, name=name or model.name + "_toy", output=dict(model.output),
                     provenance=dict(source="toy", idealized_from=model.name))
    return toy, report
