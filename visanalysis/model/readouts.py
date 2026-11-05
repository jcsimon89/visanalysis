"""What to look at in a predicted response.

switch_triggered_average   the response around every change of a label,
                           grouped by (from, to): the readout where a
                           SwitchedLN's instant switch shows itself.
by_label                   mean, SD and count of the response per label value.
recover_ln                 fit a filter and nonlinearity back from a (sampled)
                           response per label value, Baccus & Meister scaled,
                           with ridge and a smoothness penalty as the
                           pipelines use. Run on a model's own prediction it
                           must return the model; run on a prediction with
                           noise at the real imaging rate it says what an
                           experiment of that length could resolve.
"""
import numpy as np

from visanalysis.model.filters import TabulatedFilter
from visanalysis.model.nonlinearities import TabulatedNL

__all__ = ["switch_triggered_average", "by_label", "recover_ln", "concat"]


def concat(responses):
    """Flatten a list of sampled Responses (planes) into (t, y, labels, plane) sorted by time."""
    t = np.concatenate([r.t for r in responses])
    y = np.concatenate([r.y for r in responses])
    pl = np.concatenate([np.full(r.t.size, -1 if r.plane is None else r.plane) for r in responses])
    labels = {k: np.concatenate([np.asarray(r.labels[k]) for r in responses]) for k in responses[0].labels}
    o = np.argsort(t, kind="stable")
    return t[o], y[o], {k: v[o] for k, v in labels.items()}, pl[o]


def switch_triggered_average(resp, label, pre=2.0, post=10.0, bin_s=0.1, ignore=(0,), subtract_pre=False):
    """Mean response around each change of `label`, by (from, to).

    resp may be a fine-grid Response or a list of sampled ones. Transitions from or to a value in
    `ignore` (the grey) are skipped. Returns {(from, to): dict(lag, mean, se, n)}."""
    rs = resp if isinstance(resp, (list, tuple)) else [resp]
    t, y, labels, _ = concat(rs)
    lab = labels[label]
    change = np.flatnonzero(lab[1:] != lab[:-1]) + 1
    edges = np.arange(-pre, post + bin_s / 2, bin_s)
    lag = 0.5 * (edges[:-1] + edges[1:])
    groups = {}
    for i in change:
        a, b = lab[i - 1].item(), lab[i].item()
        if a in ignore or b in ignore:
            continue
        # the switch time is between samples i-1 and i; on the fine grid that is t[i] - dt/2
        ts = t[i] - 0.5 * (t[i] - t[i - 1]) if rs[0].dt is None else t[i] - 0.5 * rs[0].dt
        m = (t >= ts - pre) & (t < ts + post)
        if not m.any():
            continue
        yy = y[m]
        if subtract_pre:
            yy = yy - np.nanmean(yy[t[m] < ts]) if np.any(t[m] < ts) else yy
        b_ = np.clip(np.searchsorted(edges, t[m] - ts, "right") - 1, 0, lag.size - 1)
        row = np.full(lag.size, np.nan)
        for k in np.unique(b_):
            vals = yy[b_ == k]
            if np.isfinite(vals).any():
                row[k] = np.nanmean(vals)
        groups.setdefault((a, b), []).append(row)
    out = {}
    for k, rows in groups.items():
        R = np.array(rows)
        n = np.sum(np.isfinite(R), axis=0)
        with np.errstate(invalid="ignore"):
            mean = np.nanmean(R, axis=0)
            se = np.nanstd(R, axis=0) / np.sqrt(np.maximum(n, 1))
        out[k] = dict(lag=lag, mean=mean, se=se, n=int(R.shape[0]), rows=R)
    return out


def by_label(resp, label, ignore=(0,), valid_only=True):
    rs = resp if isinstance(resp, (list, tuple)) else [resp]
    t, y, labels, _ = concat(rs)
    lab = labels[label]
    ok = np.isfinite(y)
    if valid_only:
        ok &= t >= rs[0].valid_from
    out = {}
    for v in np.unique(lab):
        if v.item() in ignore:
            continue
        m = ok & (lab == v)
        if m.sum() > 1:
            out[v.item()] = dict(mean=float(y[m].mean()), sd=float(y[m].std()), n=int(m.sum()))
    return out


def _penalty(n):
    D = np.diff(np.eye(n), n=2, axis=0)
    return D.T @ D


def _taps(stim, t, n_lag, lag_dt):
    """(len(t), n_lag): tap j = mean contrast over (t - (j+1) lag_dt, t - j lag_dt]."""
    cs = np.r_[0.0, np.cumsum(stim.c)]
    edges = np.round((t[:, None] - np.arange(n_lag + 1)[None, :] * lag_dt - stim.t[0]) / stim.dt) + 1
    edges = np.clip(edges.astype(np.int64), 0, cs.size - 1)
    hi, lo = edges[:, :-1], edges[:, 1:]
    return (cs[hi] - cs[lo]) / np.maximum(hi - lo, 1)


def recover_ln(resp, stim, label=None, filter_length=3.0, lag_dt=0.02, ridge=1.0, steady_label="since_switch",
               steady_start=None, n_bins=12, ignore=(0,), demean=True):
    """Fit an LN per label value from a sampled response, Baccus & Meister scaled.

    Returns {value: dict(filter, nonlinearity, h_tap, lags, scale, sd_c, r, n_rows)} where `filter`
    is a TabulatedFilter in contrast units and `nonlinearity` a TabulatedNL on the generator
    (offset = the mean response in that state, so demeaning is undone). label None fits one LN
    over everything (value 'all'). steady_start (default filter_length) drops rows whose history
    reaches across a change of state, using the `steady_label` if the stimulus has it."""
    pairs = _pairs(resp, stim)
    n_lag = int(round(filter_length / lag_dt))
    steady_start = filter_length if steady_start is None else steady_start
    A_, y_, lab_, ok_, c_, cl_ = [], [], [], [], [], []
    for rs, st in pairs:
        t, y, labels, _ = concat(rs)
        ok = np.isfinite(y) & (t >= st.t[0] + filter_length)
        if steady_label in labels:
            with np.errstate(invalid="ignore"):
                ok &= np.nan_to_num(labels[steady_label], nan=np.inf) >= steady_start
        A_.append(_taps(st, t, n_lag, lag_dt))
        y_.append(y)
        ok_.append(ok)
        lab_.append(None if label is None else labels[label])
        c_.append(st.c)
        cl_.append(None if label is None else np.asarray(st.labels[label]))
    A, y, ok = np.concatenate(A_), np.concatenate(y_), np.concatenate(ok_)
    lab = None if label is None else np.concatenate(lab_)
    c_all = np.concatenate(c_)
    cl_all = None if label is None else np.concatenate(cl_)
    P = _penalty(n_lag)
    lags = (np.arange(n_lag) + 0.5) * lag_dt
    groups = {"all": ok} if label is None else {
        v.item(): ok & (lab == v) for v in np.unique(lab) if v.item() not in ignore}
    # the stimulus SD per state, for the B&M scale, from the stimulus itself (not the rows)
    out = {}
    for v, m in groups.items():
        if m.sum() < 2 * n_lag:
            continue
        Am, ym = A[m], y[m]
        y0 = ym.mean() if demean else 0.0
        ax = Am.mean(0)
        X = Am - ax
        G = X.T @ X
        h = np.linalg.solve(G + ridge * np.trace(G) / n_lag * P, X.T @ (ym - y0))
        sd_c = float(c_all[np.isfinite(c_all)].std()) if label is None else float(c_all[cl_all == v].std())
        raw = Am @ h
        scale = sd_c / raw.std() if raw.std() > 0 else np.nan
        g = raw * scale
        pred = (Am - ax) @ h + y0
        r = float(np.corrcoef(pred, ym)[0, 1]) if ym.std() > 0 else np.nan
        e = np.quantile(g, np.linspace(0, 1, n_bins + 1))
        b = np.clip(np.searchsorted(e, g, "right") - 1, 0, n_bins - 1)
        nl = np.array([[g[b == i].mean(), ym[b == i].mean() - y0, ym[b == i].std() / np.sqrt(max((b == i).sum(), 1))]
                       for i in range(n_bins) if (b == i).any()])
        out[v] = dict(filter=TabulatedFilter.from_taps(h, lag_dt, scale, output_units="contrast",
                                                       normalization="BM2002", name=str(v)),
                      nonlinearity=TabulatedNL(nl[:, 0], nl[:, 1], nl[:, 2], input_units="contrast", offset=y0,
                                               name=str(v)),
                      h_tap=h, lags=lags, scale=float(scale), sd_c=sd_c, r=r, n_rows=int(m.sum()), ridge=ridge)
    return out


def _pairs(resp, stim):
    """Normalise (resp, stim) to [(list of Responses, Stimulus)]: one stimulus with one Response or a
    list of planes, or a list of stimuli with a matching list of Response(s) each."""
    if isinstance(stim, (list, tuple)):
        if not isinstance(resp, (list, tuple)) or len(resp) != len(stim):
            raise ValueError("with a list of stimuli, resp must be a list of the same length")
        return [(r if isinstance(r, (list, tuple)) else [r], s) for r, s in zip(resp, stim)]
    return [(resp if isinstance(resp, (list, tuple)) else [resp], stim)]
