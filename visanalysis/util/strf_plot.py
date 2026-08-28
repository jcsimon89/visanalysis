"""
Shared plotting helpers for the noizone STRF scripts.

Deliberately dependency-light -- numpy only. analyze_data_strf_noizone.py pulls
in the whole visanalysis/plugin stack (pandas, seaborn, h5py ...), and
average_strf_noizone.py is meant to be a standalone final step that can run
without any of it. Importing these from there rather than from the analysis
script keeps that true.
"""

import numpy as np


def _pole_label(rtz):
    """Axis label for a bar orientation. The position axis is 90 - psi, where
    psi is the angle from that orientation's pole axis a_hat = (cos rtz,
    sin rtz, 0): the anterior-posterior axis at rtz=0, the dorsal-ventral axis
    at rtz=90. Zero is the projector's optical axis (straight lateral);
    positive is anterior at rtz=0 and dorsal at rtz=90.

    Only rtz=90 carries an equivalence in its label, and that asymmetry is
    deliberate. There the pole is the dorsal axis, so constant psi IS constant
    elevation: the axis equals elevation exactly, everywhere on the screen. At
    rtz=0 the pole is the anterior axis and constant psi is a ring about it --
    cos(psi) = cos(el)*cos(az) -- so the axis equals 90 - azimuth only on the
    horizon, departing by ~14 deg at psi=60, el=44. No azimuth equivalence is
    printed because there isn't one to print.
    """
    if abs(rtz % 180) < 1e-6:
        return 'Anterior (+) / posterior (-)'
    if abs(rtz % 180 - 90) < 1e-6:
        return 'Dorsal (+) / ventral (-)  [= elevation]'
    return f'Position from rtz={rtz:.0f}deg pole'

def _panel_size(extent, height_in=1.7):
    """(width, height) inches for one imshow panel, matched to its data aspect.

    extent is [left, right, bottom, top] in degrees. With aspect='equal' a
    panel whose shape does not match the data gets letterboxed, and the leftover
    is whitespace between panels -- which is what made the lag grid look sparse.
    """
    span_x = abs(extent[1] - extent[0])
    span_y = abs(extent[3] - extent[2])
    ratio = (span_x / span_y) if span_y else 1.0
    return height_in * ratio, height_in

def _grid_shape(n, extent):
    """Rows/cols giving a roughly square FIGURE, given the panel's data aspect.

    Prefers a column count that divides n exactly -- a half-empty last row
    wastes as much space as the letterboxing this is meant to fix. Among
    candidates, fewest blanks wins, then closeness to the square-figure ideal.
    """
    span_x = abs(extent[1] - extent[0])
    span_y = abs(extent[3] - extent[2])
    ratio = (span_x / span_y) if span_y else 1.0
    ideal = np.sqrt(max(n, 1) / max(ratio, 1e-6))
    best = min(range(1, int(min(n, 8)) + 1),
               key=lambda c: ((-n) % c, abs(c - ideal)))
    return int(np.ceil(n / best)), int(best)

def _tighten(fig, w_pad=0.01, h_pad=0.01, wspace=0.01, hspace=0.03):
    """Squeeze constrained_layout padding, tolerating matplotlib API drift."""
    try:
        fig.get_layout_engine().set(w_pad=w_pad, h_pad=h_pad,
                                    wspace=wspace, hspace=hspace)
    except Exception:
        pass
