"""
compare_model_to_data.py
Run a library model on the stimulus a recording actually showed and compare the
prediction with the measured responses of the ROIs the pipeline kept.

    python compare_model_to_data.py --model toy/tm9_control \
        --experiment_file_directory ".../JS256_x_JS261/fly_011" --series 2 \
        --keep ".../contrast_switch_export/JS256_x_JS261/fly_011/series_002/keep_calcium.npy" \
        --output_directory DIR

    python compare_model_to_data.py --model toy/mi1_control \
        --experiment_file_directory ".../JS140_x_JS252/fly_013" --series 3 --tag final \
        --keep ".../contrast_blocks/JS140_x_JS252/fly_013/contrast_blocks_final.hdf5" --channel 1 \
        --output_directory DIR

The protocol is read from the series; contrast-switch and contrast-blocks recordings are
supported (visanalysis.model.recording). Writes <model>__<fly>_<series>.png and .json.
"""
import argparse
import json
import os
import re

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from visanalysis.analysis import imaging_data
from visanalysis.model import spec as mspec, recording, compare
from visanalysis.model.models import model_from_spec

try:
    from visanalysis.util.strf_plot import apply_figure_style
    apply_figure_style()
except Exception:
    pass


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--model_library", default=None)
    ap.add_argument("--experiment_file_directory", required=True)
    ap.add_argument("--series", type=int, required=True)
    ap.add_argument("--tag", default="raw", help="reads fly_<tag>.hdf5 if present, else fly.hdf5")
    ap.add_argument("--keep", required=True, help="kept-ROI .npy, or a contrast-blocks result .hdf5")
    ap.add_argument("--keep_series", type=int, default=None,
                    help="series whose kept ROIs to read from a result .hdf5 (default: --series); a frozen-noise "
                         "series has no gate of its own and uses the random-seed series' ROIs")
    ap.add_argument("--channel", default="calcium", help="contrast switch: channel name; contrast blocks: channel number")
    ap.add_argument("--roi_set", default="mask_ch1")
    ap.add_argument("--steady_start", type=float, default=2.0)
    ap.add_argument("--window_s", type=float, default=16.0)
    ap.add_argument("--trace_plane", type=int, default=None,
                    help="example trace from this plane only; default pools all planes by stimulus position")
    ap.add_argument("--trace_bin", type=float, default=None, help="s; bin for the pooled example trace")
    ap.add_argument("--output_directory", required=True)
    ap.add_argument("--show_figs", default="False")
    args = ap.parse_args()
    if args.show_figs != "True":
        matplotlib.use("Agg")
    os.makedirs(args.output_directory, exist_ok=True)

    efd = args.experiment_file_directory
    tagged = "fly.hdf5" if args.tag == "raw" else "fly_" + args.tag + ".hdf5"
    efn = tagged if os.path.exists(os.path.join(efd, tagged)) else "fly.hdf5"
    ID = imaging_data.ImagingDataObject(os.path.join(efd, efn), args.series, quiet=True)
    pid = str(ID.getRunParameters().get("protocol_ID", ""))
    s = mspec.load_model(args.model, library=args.model_library)
    model = model_from_spec(s)
    print("model {!r}\nrecording {} series {} ({})".format(model, efd, args.series, pid))
    if "ContrastBlock" in pid:
        ch = int(args.channel) if str(args.channel).isdigit() else 1
        keep = recording.load_keep(args.keep, args.keep_series or args.series, ch, tag=args.tag)
        segs = recording.contrast_blocks_segments(ID, keep, channel=ch)
    else:
        keep = recording.load_keep(args.keep)
        segs = recording.contrast_switch_segments(ID, keep, channel=args.channel, roi_set=args.roi_set)
    print("{} ROIs kept, {} segment(s), {} planes".format(keep.size, len(segs), len(segs[0]["planes"])))
    res = compare.compare(model, segs, steady_start=args.steady_start)
    print("overall r {:.3f}".format(res["r_all"]))
    for v, d in res["per_state"].items():
        print("  state {:+d}: r all {:.3f} steady {:.3f}; sd data {:.4f} model {:.4f}; mean data {:+.4f} model {:+.4f}".format(
            v, d.get("r_all", np.nan), d.get("r_steady", np.nan), d.get("sd_data_all", np.nan), d.get("sd_model_all", np.nan),
            d.get("mean_data_all", np.nan), d.get("mean_model_all", np.nan)))
    for v, d in res.get("frozen_average", {}).items():
        print("  frozen-noise average over {} trials, state {:+d}: r {:.3f} slope {:.2f} sd data {:.4f} model {:.4f}".format(
            d["n_trials"], v, d["r"], d["slope"], d["sd_data"], d["sd_model"]))
    fly = os.path.basename(os.path.normpath(efd))
    base = os.path.join(args.output_directory, "{}__{}_series_{:03d}".format(
        re.sub(r"[^A-Za-z0-9._-]+", "_", args.model), fly, args.series))
    title = "{} vs {} series {} ({} ROIs)".format(args.model, fly, args.series, keep.size)
    fig = compare.summary_figure(model, segs, res, title=title, window_s=args.window_s, plane=args.trace_plane,
                                 trace_bin_s=args.trace_bin)
    fig.savefig(base + ".png", dpi=140)
    if args.show_figs != "True":
        plt.close(fig)
    out = dict(model=args.model, model_hash=mspec.model_hash(s), fly=efd, series=args.series, protocol=pid,
               n_keep=int(keep.size), r_all=res["r_all"], per_state=res["per_state"], steady_start=args.steady_start,
               frozen_average=res.get("frozen_average", {}),
               sta_data={"{}->{}".format(a, b): dict(lag=d["lag"].tolist(), mean=d["mean"].tolist(), se=d["se"].tolist(), n=d["n"])
                         for (a, b), d in res["sta_data"].items()},
               sta_model={"{}->{}".format(a, b): dict(lag=d["lag"].tolist(), mean=d["mean"].tolist(), n=d["n"])
                          for (a, b), d in res["sta_model"].items()})
    with open(base + ".json", "w") as f:
        json.dump(out, f, indent=1, default=lambda v: v.tolist() if isinstance(v, np.ndarray) else float(v))
    print("saved {}.png/.json".format(base))
    if args.show_figs == "True":
        plt.show()


if __name__ == "__main__":
    main()
