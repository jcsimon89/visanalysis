"""
summarize_model_fits.py
One table and one figure from the .json files compare_model_to_data.py writes:
how well each model predicts the recording it was compared with, per contrast state.

    python summarize_model_fits.py --inputs DIR_OR_JSON [DIR_OR_JSON ...] --output_directory DIR

Columns: r on single-trial rows (all rows, and steady rows >= steady_start after a switch);
slope of data on prediction (1 = gain right); for frozen-noise recordings the r against the
trial-averaged data, the noise ceiling (sqrt of the Spearman-Brown reliability of that average,
from odd/even trial halves) and the fraction of it reached; response level data vs model; the
switch-triggered step data vs model.
"""
import argparse
import glob
import json
import os

import numpy as np
import matplotlib
import matplotlib.pyplot as plt


def _load(inputs):
    files = []
    for p in inputs:
        if os.path.isdir(p):
            files += sorted(glob.glob(os.path.join(p, "*.json")))
        else:
            files.append(p)
    return [json.load(open(f)) for f in files]


def _step(sta, key_hi="-1->1", key_lo="1->-1"):
    """late level after the switch minus the level before, for each direction."""
    out = {}
    for key, nm in ((key_hi, "to_high"), (key_lo, "to_low")):
        d = sta.get(key)
        if not d:
            continue
        lag, mean = np.asarray(d["lag"]), np.asarray(d["mean"], float)
        pre = np.nanmean(mean[lag < 0])
        late = np.nanmean(mean[lag > lag.max() - 2.0])
        out[nm] = float(late - pre)
    return out


def rows_from(results):
    rows = []
    for r in results:
        fly = os.path.basename(os.path.normpath(r["fly"]))
        geno = os.path.basename(os.path.dirname(os.path.normpath(r["fly"])))
        sd, sm = _step(r["sta_data"]), _step(r["sta_model"])
        for v, nm in (("-1", "low"), ("1", "high")):
            d = r["per_state"].get(v, {})
            fz = (r.get("frozen_average") or {}).get(v, {})
            rows.append(dict(model=r["model"].split("/")[-1], fly="{} {} s{}".format(geno, fly, r["series"]),
                             protocol="switch" if "Switch" in r["protocol"] else "blocks", n_roi=r["n_keep"], state=nm,
                             r_all=d.get("r_all", np.nan), r_steady=d.get("r_steady", np.nan),
                             slope_all=d.get("slope_all", np.nan), slope_steady=d.get("slope_steady", np.nan),
                             r_avg=fz.get("r", np.nan), ceiling=fz.get("ceiling", np.nan),
                             frac=fz.get("fraction_of_ceiling", np.nan), n_trials=fz.get("n_trials", 1),
                             level_data=d.get("mean_data_steady", np.nan), level_model=d.get("mean_model_steady", np.nan),
                             sd_data=d.get("sd_data_steady", np.nan), sd_model=d.get("sd_model_steady", np.nan),
                             step_data=sd.get("to_high" if nm == "high" else "to_low", np.nan),
                             step_model=sm.get("to_high" if nm == "high" else "to_low", np.nan)))
    return rows


def table_figure(rows, title):
    cols = ["model", "fly", "ROIs", "state", "r single-trial\nall / steady", "slope\nall / steady",
            "r vs trial avg\n(n trials)", "ceiling", "% of\nceiling", "level dF/F\ndata / model", "step dF/F\ndata / model"]
    cells = []
    for r in rows:
        cells.append([r["model"], r["fly"], str(r["n_roi"]), r["state"],
                      "{:.3f} / {:.3f}".format(r["r_all"], r["r_steady"]),
                      "{:.2f} / {:.2f}".format(r["slope_all"], r["slope_steady"]),
                      "{:.3f} ({})".format(r["r_avg"], r["n_trials"]) if np.isfinite(r["r_avg"]) else "-",
                      "{:.3f}".format(r["ceiling"]) if np.isfinite(r["ceiling"]) else "-",
                      "{:.0f}%".format(100 * r["frac"]) if np.isfinite(r["frac"]) else "-",
                      "{:+.3f} / {:+.3f}".format(r["level_data"], r["level_model"]),
                      "{:+.3f} / {:+.3f}".format(r["step_data"], r["step_model"])])
    fig, ax = plt.subplots(figsize=(17, 0.5 * len(rows) + 1.8))
    ax.axis("off")
    tb = ax.table(cellText=cells, colLabels=cols, loc="center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 1.6)
    for (i, j), c in tb.get_celld().items():
        if i == 0:
            c.set_text_props(weight="bold")
            c.set_height(c.get_height() * 1.6)
    ax.set_title(title, fontsize=11, pad=12)
    fig.tight_layout()
    return fig


def bar_figure(rows, title):
    labels = ["{}\n{} {}".format(r["model"], r["fly"].split()[-2], r["fly"].split()[-1]) for r in rows]
    uniq = []
    for l in labels:
        if l not in uniq:
            uniq.append(l)
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.4))
    x = np.arange(len(uniq))
    w = 0.36
    colors = {"low": "#1f77b4", "high": "#d62728"}
    for k, st in enumerate(("low", "high")):
        col = colors[st]
        rr = [next(r for r, l in zip(rows, labels) if l == u and r["state"] == st) for u in uniq]
        xx = x + (k - 0.5) * w
        # bar: r between the prediction and the measured single-trial responses
        ax[0].bar(xx, [r["r_steady"] for r in rr], w, color=col, alpha=0.85)
        # dot: r between the prediction and the data averaged over repeated (frozen-noise) trials
        ax[0].scatter(xx, [r["r_avg"] for r in rr], color=col, edgecolor="k", zorder=3, s=40)
        # tick: the r a perfect model could reach against that trial average (noise ceiling)
        ax[0].scatter(xx, [r["ceiling"] for r in rr], marker="_", color="k", s=320, zorder=4, linewidths=2)
        ax[1].bar(xx, [r["slope_steady"] for r in rr], w, color=col, alpha=0.85)
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(uniq, fontsize=8)
    ax[0].set_ylabel("r, prediction vs data (rows >= 2 s after a switch)")
    ax[0].set_title("how much of the response the model predicts", fontsize=10)
    ax[0].legend(handles=[Patch(color=colors["low"], label="low contrast"),
                          Patch(color=colors["high"], label="high contrast"),
                          Patch(color="0.6", label="single trials"),
                          Line2D([], [], marker="o", color="w", markerfacecolor="0.6", markeredgecolor="k", markersize=7,
                                 label="repeat average"),
                          Line2D([], [], marker="_", color="k", markersize=14, markeredgewidth=2, linestyle="none",
                                 label="noise ceiling")],
                 loc="upper center", bbox_to_anchor=(0.5, -0.16), fontsize=8, ncol=5, frameon=False)
    ax[0].set_ylim(0, 1.0)
    reps = sorted({(r["model"], r["fly"].split()[-1], r["n_trials"]) for r in rows if r["n_trials"] > 1})
    fig.text(0.01, 0.01, "repeat average and ceiling only where the series repeated one seed: " +
             ", ".join("{} {}: {} repeats".format(m, s, n) for m, s, n in reps), fontsize=8, color="0.35")
    ax[1].axhline(1, color="0.5", lw=0.8)
    ax[1].set_xticks(x)
    ax[1].set_xticklabels(uniq, fontsize=8)
    ax[1].set_ylabel("slope of data on prediction")
    ax[1].set_title("gain calibration (1 = right)", fontsize=10)
    ax[1].legend(handles=[Patch(color=colors["low"], label="low contrast"), Patch(color=colors["high"], label="high contrast")],
                 fontsize=8, frameon=False)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    return fig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--inputs", nargs="+", required=True)
    ap.add_argument("--output_directory", required=True)
    ap.add_argument("--title", default="toy LN models vs data")
    ap.add_argument("--show_figs", default="False")
    args = ap.parse_args()
    if args.show_figs != "True":
        matplotlib.use("Agg")
    os.makedirs(args.output_directory, exist_ok=True)
    rows = rows_from(_load(args.inputs))
    fig = table_figure(rows, args.title)
    fig.savefig(os.path.join(args.output_directory, "model_fit_table.png"), dpi=150)
    fig2 = bar_figure(rows, args.title)
    fig2.savefig(os.path.join(args.output_directory, "model_fit_bars.png"), dpi=150)
    with open(os.path.join(args.output_directory, "model_fit_table.csv"), "w") as f:
        keys = list(rows[0].keys())
        f.write(",".join(keys) + "\n")
        for r in rows:
            f.write(",".join("{:.4f}".format(r[k]) if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    for r in rows:
        print("{:14s} {:28s} {:4s} r {:.3f}/{:.3f} slope {:.2f}/{:.2f} r_avg {:.3f} ceil {:.3f} frac {:.2f}".format(
            r["model"], r["fly"], r["state"], r["r_all"], r["r_steady"], r["slope_all"], r["slope_steady"],
            r["r_avg"], r["ceiling"], r["frac"]))
    if args.show_figs == "True":
        plt.show()


if __name__ == "__main__":
    main()
