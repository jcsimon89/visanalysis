"""Exporters: from a pipeline's fitted results to a model in the library.

Each pipeline that fits an LN model gets one function here that turns what it
fitted into a spec and saves it. The pipelines keep writing their own outputs;
this is in addition, and it is the only way a measured model enters the
library. The exporter records enough provenance to trace the model back: fly,
series, channel, ROI gate, the fit settings, the stimulus the fit was made on,
and the headline numbers of the fit.

contrast_switch_model(ctx, res, arrs, channel, ...)
    SwitchedLN on 'contrast_state' with the steady-state filter and
    nonlinearity of each contrast. The filter is the pipeline's h (per tap of
    lag_dt) times its Baccus & Meister scale (arrs['bm_scale']), stored as a
    density in contrast units. The nonlinearity is the pipeline's binned
    (generator, response) with the response's mean level in that state as
    `offset` -- the pipeline fits on a response demeaned per state, so the
    switch-triggered average's late level is the only record of where each
    state sits.
"""
import os

import numpy as np

from visanalysis.model import spec as mspec

__all__ = ["contrast_switch_model", "contrast_switch_model_name"]


def _fly_series(fly_dir, series_name):
    fly_dir = os.path.normpath(str(fly_dir))
    fly = os.path.basename(fly_dir)
    geno = os.path.basename(os.path.dirname(fly_dir))
    return geno, fly, str(series_name)


def contrast_switch_model_name(fly_dir, series_name, channel):
    geno, fly, sn = _fly_series(fly_dir, series_name)
    return "contrast_switch/{}_{}_{}_{}".format(geno, fly, sn, channel)


def contrast_switch_model(ctx, res, arrs, channel, fly_dir, series_name, library=None, overwrite=True,
                          indicator=None, save=True):
    """Build (and save) the SwitchedLN spec for one channel of a contrast-switch fit.

    ctx, res, arrs are build_context's context and analyse_channel's two outputs as the
    analyze_data_trf_contrast_switch script holds them (res = the whole results dict, arrs[channel]
    = the channel's arrays, arrs['sta_' + channel] its switch-triggered arrays)."""
    ch = res["channels"][channel]
    a = arrs[channel]
    if "bm_scale" not in a:
        raise KeyError("analyse_channel did not return bm_scale; the pipeline is older than the exporter")
    sta = res.get("switch_triggered", {}).get(channel, {})
    lag_dt = float(ch["lag_dt"])
    geno, fly, sn = _fly_series(fly_dir, series_name)
    name = contrast_switch_model_name(fly_dir, series_name, channel)
    meta0 = ctx["metas"][0]
    levels = {-1: sta.get("to_low", {}).get("level_late", 0.0), 1: sta.get("to_high", {}).get("level_late", 0.0)}
    states = []
    for q, v, nm in ((0, -1, "low"), (1, 1, "high")):
        h = np.asarray(a["h_" + nm], float)
        nl = np.asarray(a["nl_" + nm], float)
        states.append(dict(
            value=v, name=nm,
            filter=dict(kind="tabulated", lag_edges=np.arange(h.size + 1) * lag_dt,
                        density=h * float(a["bm_scale"][q]) / lag_dt, output_units="contrast",
                        normalization="BM2002", name=nm),
            nonlinearity=dict(kind="tabulated", x=nl[:, 0], y=nl[:, 1], y_se=nl[:, 2], input_units="contrast",
                              offset=float(levels[v]), extrapolation="linear", name=nm)))
    desc = ("Steady-state LN per contrast state, {} {} {} {}: pooled filter over the gated ROIs, Baccus & "
            "Meister scaled; nonlinearity = binned response on the out-of-sample generator, offset = late "
            "level after the switch into that state (switch-triggered average)").format(geno, fly, sn, channel)
    s = mspec.new_spec(
        name=name.rsplit("/", 1)[-1], kind="SwitchedLN", source="pipeline", description=desc, channel=channel,
        signal="imaged", units="dF/F vs epoch median", indicator=indicator,
        pipeline="analyze_data_trf_contrast_switch", fly=str(fly_dir), genotype=geno, fly_id=fly, series=sn,
        gate=dict(res.get("gate", {})),
        settings=dict(res.get("settings", {})),
        fit=dict(n_roi=ch["n_roi"], n_rows=ch["n_rows"], n_planes=ch["n_planes"], ridge=ch["ridge"],
                 filter_length=ch["filter_length"], lag_dt=lag_dt, steady_start=ch["steady_start"],
                 first_dwell=ch["first_dwell"], cv_r_low=ch["cv_r_low"], cv_r_high=ch["cv_r_high"],
                 null_max_low=ch["null_max_low"], null_max_high=ch["null_max_high"],
                 has_filter_both_states=ch.get("has_filter_both_states"),
                 bm_scale=[float(x) for x in a["bm_scale"]], stimulus_sd=[float(x) for x in a["sd_state"]],
                 gain=_sub(ch.get("gain"), ("low", "high", "ratio", "ci", "ratio_vs_static")),
                 shared_range=_sub(ch.get("shared_range"), ("low", "high", "ratio", "ci", "ratio_vs_static")),
                 stretch_low_vs_high=_sub(ch.get("stretch_low_vs_high"), ("value", "ci")),
                 levels=dict(low=float(levels[-1]), high=float(levels[1]))),
        fit_stimulus=dict(protocol="ContrastSwitchTernaryNoise", generator=_plain(meta0), seeds=list(ctx["seeds"]),
                          frozen=bool(ctx["frozen"]), n_trials=int(ctx["n_ep"]),
                          trial_duration=float(np.median(ctx["stim_dur"])), dwell_s=float(ctx["dwell_s"]),
                          realized_contrasts=[float(x) for x in ctx["realized"]],
                          contrast_ratio=float(ctx["contrast_ratio"])))
    s["label"] = "contrast_state"
    s["states"] = states
    s["default"] = None
    if save:
        return mspec.save_model(s, name, library=library, overwrite=overwrite)
    return s


def _sub(d, keys):
    if not isinstance(d, dict):
        return None
    return {k: _plain(d[k]) for k in keys if k in d}


def _plain(v):
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    if isinstance(v, np.generic):
        return v.item()
    if isinstance(v, np.ndarray):
        return v.tolist()
    return v
