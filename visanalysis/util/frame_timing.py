"""
Generic display-frame drop reconstruction.

https://github.com/ClandininLab/visanalysis
"""
import numpy as np


def reconstruct_displayed_frame_indices(n_raw_frames, dropped_frame_inds):
    """
    Reconstruct which raw display-frame index was actually on screen at each
    raw display-frame slot, given a set of dropped (skipped) frame indices.

    The display recomputes its intended frame index fresh from wall-clock
    time on every successful draw, rather than incrementing a counter from
    wherever the last successful draw left off. So a dropped frame is never
    shown later -- it's permanently skipped, and the display just continues
    showing whatever was drawn last until the next successful draw, which
    lands exactly on its own correct nominal index (not delayed). This is
    content-agnostic: it works the same whether every raw frame is distinct
    or many consecutive raw frames happen to show identical content.

    Parameters
    ----------
    n_raw_frames       : int -- number of raw display-frame slots to
                         reconstruct (e.g. epoch_duration * frame_rate)
    dropped_frame_inds : array-like of int -- raw frame indices (0-based)
                         that were dropped (not successfully drawn).
                         Indices outside [0, n_raw_frames) are ignored.

    Returns
    -------
    actual_idx : (n_raw_frames,) int array -- actual_idx[j] is the nominal
                frame index that was actually being displayed during raw
                display-frame slot j. -1 for a leading run of drops with no
                prior valid frame to hold over from.
    """
    dropped = np.zeros(n_raw_frames, dtype=bool)
    idx = np.asarray(dropped_frame_inds, dtype=int)
    idx = idx[(idx >= 0) & (idx < n_raw_frames)]
    dropped[idx] = True

    positions = np.arange(n_raw_frames)
    actual_idx = np.where(dropped, -1, positions)
    np.maximum.accumulate(actual_idx, out=actual_idx)
    return actual_idx
