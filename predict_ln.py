"""
predict_ln.py
Run LN models from the model library forward on a stimulus and show what they
predict: the trace, the switch-triggered average, the response per state, and
what an experiment of this length and noise would recover of the model.

    python predict_ln.py --models contrast_switch/JS256_x_JS261_fly_011_series_002_calcium \
        --stimulus ContrastSwitch --stim_params "low_contrast=0.2,high_contrast=0.92,modulator_dwell=10,probe_tau=0.05" \
        --duration 300 --seed 0 --output_directory "path/to/out"

    python predict_ln.py --models <name>[,<name>] --stimulus recording \
        --experiment_file_directory "path/to/fly_folder" --series 2 --output_directory ...

--models        library names (subdirectory/name) or paths to .yaml files, comma separated.
--stimulus      a protocol known to visanalysis.model.stimulus (ContrastSwitch), 'recording'
                (the frames a recorded series showed, from its photodiode timing), or a
                synthetic signal: step, sinusoid, white_noise, ternary.
--stim_params   "key=value,key=value" for the protocol or synthetic source. A protocol takes its
                own parameter names (low_contrast, high_contrast, modulator_dwell, probe_tau,
                fps, ...); unspecified ones are the generator's defaults.
--bm_normalize  True: scale each model filter so var(F*s) = var(s) on THIS stimulus, per state.
                For hypothetical models whose nonlinearities are meant in measured units. Never
                for a model exported from data: that one is already scaled.
--indicator     name from measurement.INDICATORS, applied to models whose output is 'neural'.
Imaging: --volume_period, --n_planes, --plane_dt, --noise_sd give the sampling the recovery
sees; a recording supplies its own volume period.

Outputs, per model, in <output_directory>: <model>__<stimulus>.png / .json / .npz.
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from visanalysis.model import spec as mspec
from visanalysis.model import stimulus as mstim
from visanalysis.model import measurement, readouts, plots
from visanalysis.model.models import model_from_spec, bm_normalize, SwitchedLN

try:
    from visanalysis.util.strf_plot import apply_figure_style
    apply_figure_style()
except Exception:
    pass


def _parse_params(s):
    out = {}
    for item in [x for x in (s or "").split(",") if x.strip()]:
        k, v = item.split("=", 1)
        v = v.strip()
        for cast in (int, float):
            try:
                v = cast(v)
                break
            except ValueError:
                continue
        out[k.strip()] = v
    return out


def _slug(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s).strip("_")[:80]


def _recording_stimuli(args):
    from visanalysis.analysis import imaging_data
    efd = args.experiment_file_directory
    tagged = "fly.hdf5" if args.tag == "raw" else "fly_" + args.tag + ".hdf5"
    efn = tagged if os.path.exists(os.path.join(efd, tagged)) else "fly.hdf5"
    efp = os.path.join(efd, efn)
    if not os.path.exists(efp):
        sys.exit("No such file: " + efp)
    if args.rig == "Bruker":
        from visanalysis.plugin import bruker
        plug = bruker.BrukerPlugin()
    else:
        from visanalysis.plugin import base as base_plugin
        plug = base_plugin.BasePlugin()
    series = [int(s) for s in plug.getSeriesNumbers(efp)]
    if args.series:
        series = [int(args.series)]
    for s in series:
        ID = imaging_data.ImagingDataObject(efp, s, quiet=True)
        pid = ID.getRunParameters().get("protocol_ID", "")
        try:
            mstim._lookup(pid)
        except KeyError:
            if args.series:
                sys.exit("series {} is {} and has no stimulus source".format(s, pid))
            continue
        stims = mstim.from_recording(ID)
        # one timeline across epochs: shift each epoch to its start in the series
        out = []
        for st in stims:
            t0 = st.meta.get("start_time", 0.0)
            out.append(mstim.Stimulus(st.t + t0, st.c, st.labels, st.name, st.meta))
        vp = None
        try:
            vp = float(ID.getResponseTiming()["sample_period"])
        except Exception:
            pass
        try:
            offs = np.asarray(ID.getVolumeFrameOffsets(), float)
        except Exception:
            offs = None
        name = "{}_series_{:03d}".format(os.path.basename(os.path.normpath(efd)), s)
        return out, name, vp, offs
    sys.exit("no series in {} has a stimulus source".format(efp))


def _build_stimuli(args):
    """-> (list of Stimulus, name, volume_period or None, plane offsets or None)"""
    p = _parse_params(args.stim_params)
    if args.stimulus == "recording":
        return _recording_stimuli(args)
    if args.stimulus == "step":
        steps = p.pop("steps", None)
        if steps is None:
            steps = [(args.pre_time, p.pop("contrast", 1.0))]
        st = mstim.step(args.duration + args.pre_time + args.tail_time, steps, t0=-args.pre_time)
    elif args.stimulus == "sinusoid":
        st = mstim.sinusoid(args.duration, p.get("freq_hz", 1.0), p.get("amplitude", 0.5))
    elif args.stimulus == "white_noise":
        st = mstim.white_noise(args.duration, p.get("rms", 0.3), fps=p.get("fps", 120.0), seed=args.seed)
    elif args.stimulus == "ternary":
        st = mstim.ternary_renewal(args.duration, p.get("contrast", 0.5), p.get("tau", 0.05), fps=p.get("fps", 120.0),
                                   p=p.get("p", 1.0 / 6), seed=args.seed)
    else:
        st = mstim.from_protocol(args.stimulus, p, duration=args.duration, seed=args.seed, pre_time=args.pre_time,
                                 tail_time=args.tail_time)
    return [st], st.name, None, None


def _predict(model, stims, args, vp, offs):
    """Fine-grid responses and sampled planes for each stimulus."""
    ind = None if args.indicator in (None, "", "none", "None") else measurement.Indicator.preset(args.indicator)
    vp = vp or args.volume_period
    offs = np.arange(args.n_planes) * args.plane_dt if offs is None else offs
    fine, sampled = [], []
    for k, st in enumerate(stims):
        r = model.respond(st)
        if r.signal == "neural" and ind is not None:
            r = ind.apply(r)
        fine.append(r)
        sampled.append(measurement.sample(r, volume_period=vp, plane_offsets=offs, scan_time=args.plane_dt,
                                          noise_sd=args.noise_sd, seed=args.seed + k))
    return fine, sampled


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--models", required=True)
    ap.add_argument("--model_library", default=None)
    ap.add_argument("--stimulus", default="ContrastSwitch")
    ap.add_argument("--stim_params", default="")
    ap.add_argument("--duration", type=float, default=300.0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--pre_time", type=float, default=1.0)
    ap.add_argument("--tail_time", type=float, default=1.0)
    ap.add_argument("--experiment_file_directory", default=None)
    ap.add_argument("--series", default=None)
    ap.add_argument("--tag", default="raw")
    ap.add_argument("--rig", default="Bruker")
    ap.add_argument("--bm_normalize", default="False")
    ap.add_argument("--indicator", default="none")
    ap.add_argument("--volume_period", type=float, default=0.11524, help="s; a recording supplies its own")
    ap.add_argument("--n_planes", type=int, default=10)
    ap.add_argument("--plane_dt", type=float, default=0.00549, help="s between planes within a volume")
    ap.add_argument("--noise_sd", type=float, default=0.0, help="additive noise on each sample, response units")
    ap.add_argument("--recover", default="True")
    ap.add_argument("--filter_length", type=float, default=3.0)
    ap.add_argument("--lag_dt", type=float, default=0.02)
    ap.add_argument("--ridge", type=float, default=1.0)
    ap.add_argument("--sta_pre", type=float, default=2.0)
    ap.add_argument("--sta_post", type=float, default=None,
                    help="s after a switch; default = the typical time to the next switch, at most 10 s")
    ap.add_argument("--sta_bin", type=float, default=0.25)
    ap.add_argument("--output_directory", required=True)
    ap.add_argument("--save_figs", default="True")
    ap.add_argument("--show_figs", default="False")
    args = ap.parse_args()
    if args.show_figs != "True":
        matplotlib.use("Agg")
    os.makedirs(args.output_directory, exist_ok=True)

    stims, stim_name, vp, offs = _build_stimuli(args)
    if args.sta_post is None:
        args.sta_post = _typical_dwell(stims)
    print("stimulus: {} ({} segment{}, {:.0f} s)".format(stim_name, len(stims), "s" if len(stims) > 1 else "",
                                                       sum(s.duration for s in stims)))
    for name in [m.strip() for m in args.models.split(",") if m.strip()]:
        s = mspec.load_model(name, library=args.model_library)
        model = model_from_spec(s)
        if args.bm_normalize == "True":
            model = bm_normalize(model, stims[0])
        print("model: {!r}".format(model))
        fine, sampled = _predict(model, stims, args, vp, offs)
        label = model.label if isinstance(model, SwitchedLN) else None
        flat = [p for planes in sampled for p in planes]
        sta = readouts.switch_triggered_average(flat, label, args.sta_pre, args.sta_post, args.sta_bin) if label else {}
        per_state = readouts.by_label(flat, label) if label else {}
        rec = {}
        if args.recover == "True":
            rec = readouts.recover_ln(sampled, stims, label=label, filter_length=args.filter_length,
                                      lag_dt=args.lag_dt, ridge=args.ridge)
        # report
        summary = dict(model=name, model_hash=mspec.model_hash(s), stimulus=stim_name,
                       stimulus_meta=[_jsonable(st.meta) for st in stims],
                       settings=vars(args), out_of_range=fine[0].out_of_range, per_state=per_state,
                       switch_triggered={"{}->{}".format(a, b): dict(lag=d["lag"].tolist(), mean=d["mean"].tolist(),
                                                                    se=d["se"].tolist(), n=d["n"])
                                         for (a, b), d in sta.items()},
                       recovered={str(v): dict(r=d["r"], scale=d["scale"], n_rows=d["n_rows"],
                                               peak_lag=d["filter"].peak()[0], peak=d["filter"].peak()[1],
                                               nl_offset=d["nonlinearity"].offset)
                                  for v, d in rec.items()})
        for v, d in per_state.items():
            print("  state {}: mean {:+.4f}, sd {:.4f} ({} samples)".format(v, d["mean"], d["sd"], d["n"]))
        for (a, b), d in sorted(sta.items()):
            early = d["mean"][(d["lag"] > 0) & (d["lag"] < 1.0)]
            late = d["mean"][d["lag"] > args.sta_post - 2.0]
            print("  switch {} -> {}: n {}, pre {:+.4f}, 0-1 s {:+.4f}, late {:+.4f}".format(
                a, b, d["n"], np.nanmean(d["mean"][d["lag"] < 0]), np.nanmean(early), np.nanmean(late)))
        for v, d in rec.items():
            print("  recovered {}: r {:.3f}, peak {:+.3g} at {:.0f} ms, B&M scale {:.3g}".format(
                v, d["r"], d["filter"].peak()[1], d["filter"].peak()[0] * 1e3, d["scale"]))
        for k, f in fine[0].out_of_range.items():
            if f > 0.05:
                print("  NOTE {}: {:.0%} of samples lie outside the nonlinearity's measured range".format(k, f))
        base = os.path.join(args.output_directory, "{}__{}".format(_slug(name.replace("/", "_")), _slug(stim_name)))
        with open(base + ".json", "w") as f:
            json.dump(summary, f, indent=1, default=_jsonable)
        np.savez(base + ".npz", **{"t_{}".format(k): r.t for k, r in enumerate(fine)},
                 **{"y_{}".format(k): r.y for k, r in enumerate(fine)},
                 **{"g_{}".format(k): r.g for k, r in enumerate(fine)},
                 **{"sampled_t_{}".format(k): np.concatenate([p.t for p in pl]) for k, pl in enumerate(sampled)},
                 **{"sampled_y_{}".format(k): np.concatenate([p.y for p in pl]) for k, pl in enumerate(sampled)})
        if args.save_figs == "True" or args.show_figs == "True":
            fig = plots.summary_figure(model, stims[0], fine[0], sampled[0], sta, rec,
                                       title="{} on {}".format(name, stim_name))
            if args.save_figs == "True":
                fig.savefig(base + ".png", dpi=140)
            if args.show_figs != "True":
                plt.close(fig)
        print("  saved {}.*".format(base))
    if args.show_figs == "True":
        plt.show()


def _typical_dwell(stims, cap=10.0):
    """Median time between consecutive non-grey label changes over all labelled stimuli, capped."""
    gaps = []
    for st in stims:
        for lab in st.labels:
            if lab in ("stim_time", "since_switch", "dwell"):
                continue
            tr = [i for i, a, b in st.transitions(lab) if 0 not in (a, b)]
            if len(tr) > 1:
                gaps.extend(np.diff(tr) * st.dt)
    return float(min(cap, np.median(gaps))) if gaps else cap


def _jsonable(v):
    if isinstance(v, dict):
        return {str(k): _jsonable(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_jsonable(x) for x in v]
    if isinstance(v, np.ndarray):
        return v.tolist()
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, (str, int, float, bool)) or v is None:
        return v
    return str(v)


if __name__ == "__main__":
    main()
