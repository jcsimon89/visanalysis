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


def contrast_blocks_model_name(fly_dir, series_name, channel):
    geno, fly, sn = _fly_series(fly_dir, series_name)
    return "contrast_blocks/{}_{}_{}_{}".format(geno, fly, sn, channel)


def contrast_blocks_model(h5_path, fly_dir, series_name, channel="ch1", tag="final", library=None,
                          overwrite=True, indicator=None, save=True):
    """SwitchedLN on 'contrast_state' from a contrast-blocks result file (analyze_data_trf_contrast_blocks).

    The file holds, per kept ROI, each condition's Baccus & Meister filter (ln_filters, per tap of
    lag_dt) and binned nonlinearity (nl_x, nl_y, response NOT demeaned, so the level is in y).
    The model is the ROI average of both: the mean normalized filter (whose generator variance is
    then slightly below the stimulus's when ROI filters differ in shape; recorded, not corrected)
    and the mean nonlinearity, with the SEM over ROIs as y_se. Condition 0 is the low contrast
    (value -1), condition 1 the high (value +1)."""
    import h5py
    geno, fly, sn = _fly_series(fly_dir, series_name)
    name = contrast_blocks_model_name(fly_dir, series_name, channel)
    with h5py.File(h5_path, "r") as f:
        top = "TRF_CONTRAST_BLOCKS_" + tag
        g = f["{}/{}/{}".format(top, sn, channel)]
        A = dict(g.attrs)
        R = {k: g[k][...] for k in ("ln_filters", "nl_x", "nl_y", "lags", "cv_r", "ln_gain", "ln_shared_gain",
                                    "intercepts", "stim_sd", "rms_contrast", "kept_rois", "peak_lag_s",
                                    "shared_gain_ratio", "gain_ratio")}
        lag_dt = float(f[top].attrs["lag_dt"])
        flen = float(f[top].attrs["filter_length_s"])
    n_roi = int(R["ln_filters"].shape[1])
    states = []
    for q, v, nm in ((0, -1, "low"), (1, 1, "high")):
        h = np.nanmean(R["ln_filters"][q], axis=0)                 # (n_lag,), per tap of lag_dt
        x = np.nanmean(R["nl_x"][q], axis=1)
        y = np.nanmean(R["nl_y"][q], axis=1)
        n_ok = np.sum(np.isfinite(R["nl_y"][q]), axis=1)
        se = np.nanstd(R["nl_y"][q], axis=1) / np.sqrt(np.maximum(n_ok, 1))
        states.append(dict(
            value=v, name=nm,
            filter=dict(kind="tabulated", lag_edges=np.arange(h.size + 1) * lag_dt, density=h / lag_dt,
                        output_units="contrast", normalization="BM2002", name=nm),
            nonlinearity=dict(kind="tabulated", x=x, y=y, y_se=se, input_units="contrast", offset=0.0,
                              extrapolation="linear", name=nm)))
    cond = str(A.get("condition_keys", ""))
    desc = ("Steady-state LN per contrast block, {} {} {} {}: ROI-averaged Baccus & Meister filter and binned "
            "nonlinearity (response level included in y) over the {} ROIs the pipeline kept. Measured with "
            "ContrastBlockUniformNoise: uniform noise at 30 Hz, 5 s blocks, conditions {}").format(
                geno, fly, sn, channel, n_roi, cond)
    s = mspec.new_spec(
        name=name.rsplit("/", 1)[-1], kind="SwitchedLN", source="pipeline", description=desc, channel=channel,
        signal="imaged", units="dF/F", indicator=indicator,
        pipeline="analyze_data_trf_contrast_blocks", fly=str(fly_dir), genotype=geno, fly_id=fly, series=sn,
        result_file=str(h5_path), tag=tag,
        fit=dict(n_roi=n_roi, kept_rois=[int(r) for r in R["kept_rois"]], ridge=_plain(A.get("ridge")),
                 filter_length=flen, lag_dt=lag_dt, steady_start=_plain(A.get("steady_start_s")),
                 has_filter=_plain(A.get("has_filter")),
                 cv_r=dict(low=float(np.nanmean(R["cv_r"][0])), high=float(np.nanmean(R["cv_r"][1]))),
                 gain=dict(low=float(np.nanmean(R["ln_gain"][0])), high=float(np.nanmean(R["ln_gain"][1])),
                           ratio=float(np.nanmean(R["ln_gain"][0]) / np.nanmean(R["ln_gain"][1]))),
                 shared_range=dict(low=float(np.nanmean(R["ln_shared_gain"][0])),
                                   high=float(np.nanmean(R["ln_shared_gain"][1])),
                                   ratio=float(np.nanmedian(R["shared_gain_ratio"]))),
                 levels=dict(low=float(np.nanmean(R["intercepts"][0])), high=float(np.nanmean(R["intercepts"][1]))),
                 peak_lag=dict(low=float(np.nanmedian(R["peak_lag_s"][0])), high=float(np.nanmedian(R["peak_lag_s"][1]))),
                 stimulus_sd=[float(x) for x in R["stim_sd"]]),
        fit_stimulus=dict(protocol="ContrastBlockUniformNoise", conditions=cond, block_s=5.0, update_hz=30.0,
                          rms_contrasts=[float(x) for x in R["rms_contrast"]],
                          contrast_ratio=_plain(A.get("contrast_ratio"))))
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
