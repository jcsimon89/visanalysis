"""Measured responses alongside the stimulus that drove them, for comparing a model with data.

Each adapter returns a list of SEGMENTS, one per continuous stretch of
stimulus the model can be run on:

    dict(stim=Stimulus (absolute series time), planes=[dict(t, y, plane, n_roi)], ...)

`y` is the mean dF/F over the kept ROIs of one z plane at that plane's own
sample times, which is the quantity the pipelines fitted their pooled models
to; the model is run on `stim` and sampled at `t`. Two protocols so far:

contrast_switch   one segment per epoch (trial), timing and dF/F exactly as
                  contrast_switch.build_context has them (per-ROI dF/F against
                  the epoch median; dark and faded ROI-epochs dropped).
contrast_blocks   one segment for the whole series, through
                  contrast_blocks.build_design (continuous bleach correction,
                  dF/F against the series median), labels from stimulus_trace.

The kept ROIs come from the pipeline that fitted the model: a keep .npy from
the contrast-switch output directory, or the kept_rois of a contrast-blocks
result file.
"""
import numpy as np

from visanalysis.model import stimulus as mstim

__all__ = ["load_segments", "load_keep", "contrast_switch_segments", "contrast_blocks_segments"]


def load_keep(path, series=None, channel=None, tag="final"):
    """Kept ROI indices from a .npy (contrast switch) or a contrast-blocks result .hdf5."""
    path = str(path)
    if path.endswith(".npy"):
        return np.sort(np.load(path).astype(int))
    import h5py
    with h5py.File(path, "r") as f:
        g = f["TRF_CONTRAST_BLOCKS_{}/sn{}/ch{}".format(tag, int(series), int(channel))]
        return np.sort(np.asarray(g["kept_rois"][...]).astype(int))


def contrast_switch_segments(ID, keep, channel="calcium", roi_set="mask_ch1", correct_bleach="auto"):
    from visanalysis.util import contrast_switch as cs
    ctx = cs.build_context(ID, {channel: roi_set}, correct_bleach=correct_bleach)
    R = ctx["resp"][channel]
    stims = mstim.from_recording(ID)
    keep = np.asarray(keep, int)
    out = []
    for e in range(ctx["n_ep"]):
        st = stims[e]
        t0 = st.meta.get("start_time", 0.0)
        st_abs = mstim.Stimulus(st.t + t0, st.c, st.labels, st.name, st.meta)
        planes = []
        for p in np.unique(R["plane"][keep]):
            rois = keep[R["plane"][keep] == p]
            off = float(np.mean(R["off"][rois])) + R["scan"] / 2.0
            t = R["tvecs"][e] - ctx["pre"][e] + off + t0
            with np.errstate(invalid="ignore"):
                y = np.nanmean(R["resp"][e][rois], axis=0)
            planes.append(dict(t=t, y=y, plane=int(p), n_roi=int(rois.size)))
        out.append(dict(stim=st_abs, planes=planes, epoch=e, steady_label="since_switch",
                        dwell_s=float(ctx["dwell_s"]), label="contrast_state", t0=float(t0),
                        seed=int(st.meta.get("seed", -1))))
    return out


def contrast_blocks_segments(ID, keep, channel=1, filter_length_s=2.0, roi_prefix="aligned", correct_bleach=True,
                             split_cycles=True):
    """channel: the number (1 -> ROI set 'mask_ch1', as the pipeline names it) or a ROI set name.
    A frozen-noise series (every epoch the same seed) is returned as one segment per cycle of the
    conditions when split_cycles is True, otherwise as one whole-series segment."""
    from visanalysis.util import contrast_blocks as cb
    roi_set = "mask_ch{}".format(channel) if str(channel).isdigit() else str(channel)
    d = cb.build_design(ID, roi_set, filter_length_s, roi_prefix=roi_prefix, correct_bleach=correct_bleach)
    eps = ID.getEpochParameters()
    tr = cb.stimulus_trace(eps, d["onsets"], d["tvol"][-1] + 1.0, d["idle"], fs=d["fs"])
    cond = np.asarray(tr["cond"])
    state = np.where(cond < 0, 0, np.where(cond == 0, -1, 1)).astype(int)     # 0 low -> -1, 1 high -> +1
    labels = dict(contrast_state=state, since_switch=np.where(cond < 0, np.nan, np.asarray(tr["since_switch"], float)),
                  epoch=np.asarray(tr["epoch"]), stim_time=np.asarray(tr["t"]))
    stim = mstim.Stimulus(np.asarray(tr["t"]), np.asarray(tr["c"]), labels,
                          "ContrastBlocks series {} (conditions {})".format(ID.series_number, tr.get("keys")),
                          dict(protocol="ContrastBlockUniformNoise", conditions=[list(k) for k in tr.get("keys", [])],
                               rms_contrast=[float(x) for x in np.asarray(tr.get("rms_contrast", []))]))
    keep = np.asarray(keep, int)
    planes = []
    for g in np.unique(d["roi_group"][keep]):
        rois = keep[d["roi_group"][keep] == g]
        t = d["tvol"] + float(d["z_offsets"][g]) + d["dt"] / (2.0 * max(len(d["z_offsets"]), 1))
        with np.errstate(invalid="ignore"):
            y = np.nanmean(d["Y_full"][rois], axis=0)
        planes.append(dict(t=t, y=y, plane=int(g), n_roi=int(rois.size)))
    dwell = float(np.median(np.diff(d["onsets"])))
    whole = dict(stim=stim, planes=planes, epoch=0, steady_label="since_switch", label="contrast_state",
                 dwell_s=dwell, t0=float(d["onsets"][0]), seed=-1)
    if not (split_cycles and d.get("repeated")):
        return [whole]
    # frozen noise: every epoch of a condition showed the same noise, so each cycle of the conditions
    # (low epoch + high epoch) is a repeat of the same stimulus. One segment per cycle, so that
    # compare() can average the data over cycles at the same stimulus position.
    onsets = np.asarray(d["onsets"], float)
    ec = np.asarray(d["epoch_cond"])
    period = len(d["cond_keys"])
    first = int(np.argmax(ec == ec.min())) if ec.size else 0
    seed = int(eps[0].get("seed", 0)) if isinstance(eps[0], dict) and "seed" in eps[0] else 0
    out = []
    for k in range(first, len(onsets) - period + 1, period):
        t_a = onsets[k]
        t_b = onsets[k + period] if k + period < len(onsets) else onsets[k] + period * dwell
        # the segment's stimulus starts a filter length early so the first samples have history
        st_k = stim.window(t_a - filter_length_s, t_b)
        pk = []
        for p in planes:
            m = (p["t"] >= t_a) & (p["t"] < t_b)
            pk.append(dict(t=p["t"][m], y=p["y"][m], plane=p["plane"], n_roi=p["n_roi"]))
        out.append(dict(stim=mstim.Stimulus(st_k.t, st_k.c, st_k.labels, "cycle {} seed {}".format(k // period, seed),
                                            dict(stim.meta, seed=seed, cycle=k // period)),
                        planes=pk, epoch=k, steady_label="since_switch", label="contrast_state", dwell_s=dwell,
                        t0=float(t_a), seed=seed))
    return out


def load_segments(ID, keep, **kw):
    """Dispatch on the series' protocol."""
    pid = str(ID.getRunParameters().get("protocol_ID", ""))
    if "ContrastSwitch" in pid:
        return contrast_switch_segments(ID, keep, **kw)
    if "ContrastBlock" in pid:
        return contrast_blocks_segments(ID, keep, **kw)
    raise KeyError("no recording adapter for protocol {!r}".format(pid))
