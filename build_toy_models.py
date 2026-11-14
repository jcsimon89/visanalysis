"""
build_toy_models.py
Idealise measured LN models into toy models (a gamma-difference filter and a
parametric nonlinearity per contrast state) and compare them on one stimulus.

    python build_toy_models.py --sources "toy/tm9_control=contrast_switch/JS256_x_JS261_fly_011_series_002_calcium,\
toy/tm9_pde=contrast_switch/JS256_x_JS261_fly_012_series_002_calcium,toy/mi1_control=group/mi1_control" \
        --output_directory DIR [--model_library DIR] [--nl_family auto|linear|relu|power|sigmoid]

For every source: the filter of each state is fitted from --min_lag on (40 ms,
past the display leak tap), every nonlinearity family is fitted and ranked by
AICc, and the winner (or the forced family) is used. The toy is saved to the
library with source 'toy', the source model's name and hash, every candidate's
parameters and AICc, and the filter fit residual. Figures: one per toy showing
measured vs idealised, and one comparing all toys on the same stimulus.
"""
import argparse
import json
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from visanalysis.model import spec as mspec, stimulus as mstim, measurement, readouts, idealize
from visanalysis.model.models import model_from_spec
from visanalysis.model.plots import state_colors

try:
    from visanalysis.util.strf_plot import apply_figure_style
    apply_figure_style()
except Exception:
    pass


def _fit_figure(name, src, toy, report, colors):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    fig.suptitle("{}  (from {})".format(name, src.name), fontsize=11)
    for v, rep in report.items():
        c = colors[v]
        fi = rep["filter"]
        ax[0].plot(fi["lag"], fi["data"], color=c, lw=1, alpha=0.6, label="{} measured".format(toy.states[v].name))
        ax[0].plot(fi["lag"], fi["fit"], "--", color=c, lw=1.6, label="{} gamma-diff, r {:.3f} (lags >= {:.0f} ms)".format(
            toy.states[v].name, fi["r"], fi["min_lag"] * 1e3))
        m = rep["measured"]
        ax[1].errorbar(m["x"], m["y"], yerr=m["se"], fmt="o", color=c, ms=3.5, lw=0.8, label="{} measured".format(toy.states[v].name))
        xx = np.linspace(m["x"].min(), m["x"].max(), 200)
        for k, cand in enumerate(rep["nl_candidates"]):
            style = "-" if cand["family"] == rep["chosen"] else ":"
            ax[1].plot(xx, cand["nl"](xx), style, color=c, lw=2 if style == "-" else 0.9, alpha=1 if style == "-" else 0.6,
                       label="{} {}{} AICc {:.1f}".format(toy.states[v].name, cand["family"],
                                                        " *" if cand["family"] == rep["chosen"] else "", cand["aicc"]))
    ax[0].axhline(0, color="0.7", lw=0.5)
    ax[0].set_xlim(0, 1.0)
    ax[0].set_xlabel("lag (s)")
    ax[0].set_ylabel("density (contrast / contrast / s)")
    ax[0].set_title("filters", fontsize=10)
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("generator (contrast)")
    ax[1].set_ylabel("response (dF/F)")
    ax[1].set_title("nonlinearities: measured, candidates (* chosen)", fontsize=10)
    ax[1].legend(fontsize=6.5)
    ax[2].axis("off")
    lines = []
    for v, rep in report.items():
        p = rep["filter"]["params"]
        lines.append("{}: filter peak {:.0f} ms; amp1 {:+.3g} n1 {:.2f} tau1 {:.3f}; amp2 {:+.3g} n2 {:.2f} tau2 {:.3f}".format(
            toy.states[v].name, rep["filter"]["peak_lag"] * 1e3, p["amp1"], p["n1"], p["tau1"], p["amp2"], p["n2"], p["tau2"]))
        for cand in rep["nl_candidates"]:
            lines.append("   {:8s} AICc {:7.1f}  r {:.3f}  {}".format(
                cand["family"], cand["aicc"], cand["r"],
                ", ".join("{} {:.3g}".format(k, x) for k, x in cand["params"].items())))
    ax[2].text(0, 1, "\n".join(lines), va="top", ha="left", fontsize=7, family="monospace", transform=ax[2].transAxes)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def _compare_figure(toys, stim, args):
    n = len(toys)
    fig, ax = plt.subplots(3, n, figsize=(5 * n, 10), squeeze=False)
    tr = [x for x in stim.transitions("contrast_state") if 0 not in (x[1], x[2])]
    tc = stim.t[tr[1][0]] if len(tr) > 1 else stim.t[0] + 5
    win = (tc - 3.0, tc + 6.0)
    w = (stim.t >= win[0]) & (stim.t < win[1])
    for j, (name, toy) in enumerate(toys.items()):
        colors = state_colors(list(toy.states))
        r = toy.respond(stim)
        planes = measurement.sample(r, volume_period=args.volume_period, plane_offsets=np.arange(args.n_planes) * args.plane_dt,
                                    scan_time=args.plane_dt, noise_sd=0.0)
        sta = readouts.switch_triggered_average(planes, "contrast_state", pre=2, post=min(10, args.sta_post), bin_s=0.2)
        a = ax[0, j]
        for v, ln in toy.states.items():
            lag, d = ln.filter.tabulate(0.005)
            a.plot(lag, d, color=colors[v], lw=1.6, label=ln.name)
        a.axhline(0, color="0.7", lw=0.5)
        a.set_xlim(0, 0.6)
        a.set_title(name, fontsize=10)
        a.set_xlabel("lag (s)")
        if j == 0:
            a.set_ylabel("filter density (contrast / contrast / s)")
        a.legend(fontsize=8)
        a = ax[1, j]
        for v, ln in toy.states.items():
            rg = ln.nonlinearity.range or (-1, 1)
            xx = np.linspace(rg[0], rg[1], 200)
            a.plot(xx, ln.nonlinearity(xx), color=colors[v], lw=1.6, label="{} ({})".format(ln.name, ln.nonlinearity.kind))
        a.set_xlabel("generator (contrast)")
        if j == 0:
            a.set_ylabel("response (dF/F)")
        a.legend(fontsize=8)
        a = ax[2, j]
        for (f, to), d in sorted(sta.items()):
            a.plot(d["lag"], d["mean"], color=colors[to], lw=1.6, label="{} -> {}".format(f, to))
        a.axvline(0, color="0.6", lw=0.6)
        a.set_xlabel("time from switch (s)")
        if j == 0:
            a.set_ylabel("switch-triggered average (dF/F)")
        a.legend(fontsize=8)
    fig.suptitle("toy models on {}".format(stim.name), fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sources", required=True, help="toy_name=source_name, comma separated")
    ap.add_argument("--model_library", default=None)
    ap.add_argument("--min_lag", type=float, default=0.04)
    ap.add_argument("--nl_family", default="auto")
    ap.add_argument("--monophasic", default="", help="comma-separated toy names whose filters get a single lobe")
    ap.add_argument("--stim_params", default="low_contrast=0.2,high_contrast=0.92,modulator_dwell=10,probe_tau=0.05")
    ap.add_argument("--duration", type=float, default=200.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--volume_period", type=float, default=0.11524)
    ap.add_argument("--n_planes", type=int, default=10)
    ap.add_argument("--plane_dt", type=float, default=0.00549)
    ap.add_argument("--sta_post", type=float, default=10.0)
    ap.add_argument("--output_directory", required=True)
    ap.add_argument("--overwrite", default="True")
    ap.add_argument("--show_figs", default="False")
    args = ap.parse_args()
    if args.show_figs != "True":
        matplotlib.use("Agg")
    os.makedirs(args.output_directory, exist_ok=True)
    mono = {x.strip() for x in args.monophasic.split(",") if x.strip()}
    toys, summary = {}, {}
    for item in [x for x in args.sources.split(",") if x.strip()]:
        toy_name, src_name = [x.strip() for x in item.split("=", 1)]
        s = mspec.load_model(src_name, library=args.model_library)
        src = model_from_spec(s)
        toy, report = idealize.idealize_switched(src, min_lag=args.min_lag, nl_family=None if args.nl_family == "auto" else args.nl_family,
                                                 biphasic=toy_name not in mono, name=toy_name.rsplit("/", 1)[-1])
        spec_out = toy.to_spec()
        spec_out["provenance"] = dict(
            source="toy", created=mspec.new_spec("x", "LN", "toy", "")["provenance"]["created"],
            description="Toy model: gamma-difference filters fitted from {:.0f} ms on and the AICc-best parametric "
                        "nonlinearity per state, idealised from {}".format(args.min_lag * 1e3, src.name),
            idealized_from=dict(name=s.get("name"), path=s.get("_path"), hash=mspec.model_hash(s),
                                source=s.get("provenance", {}).get("source"), description=s.get("provenance", {}).get("description")),
            min_lag=args.min_lag, nl_family=args.nl_family,
            per_state={str(v): dict(filter_fit=dict(r=rep["filter"]["r"], rel_rms=rep["filter"]["rel_rms"],
                                                    peak_lag=rep["filter"]["peak_lag"], params=rep["filter"]["params"]),
                                    nl_chosen=rep["chosen"],
                                    nl_candidates=[dict(family=c["family"], aicc=c["aicc"], r=c["r"], rss=c["rss"], params=c["params"])
                                                   for c in rep["nl_candidates"]])
                       for v, rep in report.items()})
        spec_out["output"] = dict(src.output)
        p = mspec.save_model(spec_out, toy_name, library=args.model_library, overwrite=args.overwrite == "True")
        toys[toy_name] = toy
        summary[toy_name] = spec_out["provenance"]["per_state"]
        print("{} <- {}".format(toy_name, src.name))
        for v, rep in report.items():
            print("  {}: filter r {:.3f} peak {:.0f} ms; NL {} (AICc {})".format(
                toy.states[v].name, rep["filter"]["r"], rep["filter"]["peak_lag"] * 1e3, rep["chosen"],
                ", ".join("{} {:.1f}".format(c["family"], c["aicc"]) for c in rep["nl_candidates"])))
        print("  saved {}".format(p))
        fig = _fit_figure(toy_name, src, toy, report, state_colors(list(toy.states)))
        fig.savefig(os.path.join(args.output_directory, toy_name.replace("/", "_") + "_fit.png"), dpi=140)
        plt.close(fig)
    p = {}
    for kv in args.stim_params.split(","):
        k, v = kv.split("=")
        p[k.strip()] = float(v)
    stim = mstim.from_protocol("ContrastSwitch", p, duration=args.duration, seed=args.seed, pre_time=1.0, tail_time=1.0)
    fig = _compare_figure(toys, stim, args)
    fig.savefig(os.path.join(args.output_directory, "toy_models_compare.png"), dpi=140)
    with open(os.path.join(args.output_directory, "toy_models.json"), "w") as f:
        json.dump(summary, f, indent=1, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else float(v))
    if args.show_figs == "True":
        plt.show()


if __name__ == "__main__":
    main()
