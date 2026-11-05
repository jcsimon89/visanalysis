"""Figures for a prediction: what the model is, what it predicts, what comes back."""
import numpy as np
import matplotlib.pyplot as plt

from visanalysis.model.models import LN, SwitchedLN

__all__ = ["summary_figure", "state_colors"]

STATE_COLORS = {-1: "#1f77b4", 1: "#d62728", "all": "#333333"}
_CYCLE = ["#2ca02c", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]


def state_colors(values):
    out, k = {}, 0
    for v in values:
        if v in STATE_COLORS:
            out[v] = STATE_COLORS[v]
        else:
            out[v] = _CYCLE[k % len(_CYCLE)]
            k += 1
    return out


def _lns(model):
    if isinstance(model, SwitchedLN):
        return dict(model.states)
    return {"all": model}


def summary_figure(model, stim, resp, sampled=None, sta=None, recovered=None, title=None, window=None,
                   label=None):
    """2 x 3 panels: stimulus and prediction in a window; filters and nonlinearities (model solid,
    recovered dashed/points); switch-triggered averages; a text summary."""
    lns = _lns(model)
    colors = state_colors(list(lns))
    label = label or (model.label if isinstance(model, SwitchedLN) else None)
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle(title or "{} on {}".format(model.name, stim.name), fontsize=11)

    # window: around the second transition if there is one, else the first 20 s after valid_from
    if window is None:
        tr = stim.transitions(label) if label and label in stim.labels else []
        tr = [x for x in tr if 0 not in (x[1], x[2])]
        if len(tr) >= 2:
            tc = stim.t[tr[1][0]]
            window = (tc - 4.0, tc + 8.0)
        else:
            window = (resp.valid_from, min(resp.valid_from + 20.0, stim.t[-1]))
    w = (stim.t >= window[0]) & (stim.t < window[1])

    a = ax[0, 0]
    a.plot(stim.t[w], stim.c[w], color="0.3", lw=0.6)
    _shade_states(a, stim, label, w, colors)
    a.set_xlim(window)
    a.set_ylabel("Weber contrast")
    a.set_title("stimulus", fontsize=10)

    a = ax[1, 0]
    a.plot(resp.t[w], resp.y[w], color="k", lw=0.8, label="model (fine grid)")
    if sampled:
        t, y = np.concatenate([s.t for s in sampled]), np.concatenate([s.y for s in sampled])
        m = (t >= window[0]) & (t < window[1])
        a.plot(t[m], y[m], ".", color="#ff7f0e", ms=3, label="sampled{}".format(
            " + noise" if any("noise sd 0" not in n for s in sampled for n in s.notes if "noise" in n) else ""))
    _shade_states(a, stim, label, w, colors)
    a.set_xlim(window)
    a.set_xlabel("time (s)")
    a.set_ylabel("response ({})".format(resp.units))
    a.set_title("prediction", fontsize=10)
    a.legend(fontsize=8, loc="upper right")

    a = ax[0, 1]
    for v, ln in lns.items():
        lag, dens = ln.filter.tabulate(0.01)
        a.plot(lag, dens, color=colors[v], lw=1.5, label="{} (model)".format(ln.name or v))
        if recovered and v in recovered:
            lag_r, dens_r = recovered[v]["filter"].tabulate(0.01)
            a.plot(lag_r, dens_r, "--", color=colors[v], lw=1, label="{} recovered, r {:.2f}".format(
                ln.name or v, recovered[v]["r"]))
    a.axhline(0, color="0.7", lw=0.5)
    a.set_xlabel("lag (s)")
    a.set_ylabel("density ({} / contrast / s)".format(next(iter(lns.values())).filter.output_units))
    a.set_title("filters", fontsize=10)
    a.legend(fontsize=8)
    a.set_xlim(0, min(max(ln.filter.length for ln in lns.values()), 1.5))

    a = ax[1, 1]
    for v, ln in lns.items():
        nl = ln.nonlinearity
        r = nl.range
        if r is None:
            g = resp.components[ln.name or str(v)][0] if resp.components else resp.g
            g = g[np.isfinite(g)]
            r = (np.percentile(g, 0.5), np.percentile(g, 99.5)) if g.size else (-1, 1)
        x = np.linspace(r[0], r[1], 200)
        a.plot(x, nl(x), color=colors[v], lw=1.5, label="{} (model)".format(ln.name or v))
        if recovered and v in recovered:
            rn = recovered[v]["nonlinearity"]
            a.errorbar(rn.x, rn.y + rn.offset, yerr=rn.y_se, fmt="o", color=colors[v], ms=3, lw=0.8,
                       label="{} recovered".format(ln.name or v))
    a.set_xlabel("generator ({})".format(next(iter(lns.values())).nonlinearity.input_units))
    a.set_ylabel("response ({})".format(resp.units))
    a.set_title("nonlinearities (over measured / visited range)", fontsize=10)
    a.legend(fontsize=8)

    a = ax[0, 2]
    if sta:
        for (f, to), d in sorted(sta.items()):
            c = colors.get(to, "k")
            a.plot(d["lag"], d["mean"], color=c, lw=1.5, label="{} -> {} (n={})".format(f, to, d["n"]))
            a.fill_between(d["lag"], d["mean"] - d["se"], d["mean"] + d["se"], color=c, alpha=0.2, lw=0)
        a.axvline(0, color="0.6", lw=0.6)
        a.set_xlabel("time from switch (s)")
        a.set_ylabel("response ({})".format(resp.units))
        a.legend(fontsize=8)
    else:
        a.text(0.5, 0.5, "no label transitions", ha="center", va="center", transform=a.transAxes)
    a.set_title("switch-triggered average", fontsize=10)

    a = ax[1, 2]
    a.axis("off")
    lines = [repr(stim), ""]
    for v, ln in lns.items():
        lines.append("{}: {!r}".format(ln.name or v, ln.filter))
        lines.append("    {!r}".format(ln.nonlinearity))
    lines.append("")
    for k, f in resp.out_of_range.items():
        lines.append("{}: {:.1%} of samples outside the NL's measured range".format(k, f))
    if recovered:
        lines.append("")
        for v, d in recovered.items():
            lines.append("recovered {}: r {:.3f}, B&M scale {:.3g}, {} rows, peak {:.0f} ms".format(
                v, d["r"], d["scale"], d["n_rows"], d["filter"].peak()[0] * 1e3))
    for n in resp.notes:
        lines.append(n)
    a.text(0.0, 1.0, "\n".join(_wrap(l, 70) for l in lines), va="top", ha="left", fontsize=7.5, family="monospace",
           transform=a.transAxes)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def _shade_states(a, stim, label, w, colors):
    if not label or label not in stim.labels:
        return
    lab = stim.labels[label][w]
    t = stim.t[w]
    if t.size == 0:
        return
    change = np.r_[0, np.flatnonzero(lab[1:] != lab[:-1]) + 1, t.size]
    for i0, i1 in zip(change[:-1], change[1:]):
        v = lab[i0].item()
        if v in colors:
            a.axvspan(t[i0], t[min(i1, t.size - 1)], color=colors[v], alpha=0.08, lw=0)


def _wrap(s, n):
    if len(s) <= n:
        return s
    out, cur = [], ""
    for word in s.split(" "):
        if len(cur) + len(word) + 1 > n and cur:
            out.append(cur)
            cur = "    " + word
        else:
            cur = (cur + " " + word) if cur else word
    out.append(cur)
    return "\n".join(out)
