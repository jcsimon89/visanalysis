"""
check_strf_centroids.py

Diagnostic for the registered STRF average: where do the RF centroids sit, and
is the averaged filter's shape trustworthy on both flanks?

    python check_strf_centroids.py --experiment_file_directory "path/to/fly" --tag final

WHY
---
Registration puts every ROI at relative position 0, but an ROI whose RF sits
near the edge of the projected field is a problem twice over:

  1. Its filter is TRUNCATED -- part of the receptive field fell off the
     screen, so the measured flank on that side is cut short.
  2. Its centroid is BIASED inward, because the centroid is computed from the
     visible part only.

Both push the population average toward looking asymmetric in SHAPE, which is
an artefact rather than biology. The tell is a centroid distribution pressed up
against one edge of the field.

Separately, the plotted axis being asymmetric is EXPECTED whenever centroids are
off-centre: an ROI at centroid c covers relative positions [pos0 - c, pos1 - c],
so a population sitting anterior has more screen behind it than in front. That
alone is not a problem -- it just means the posterior flank is characterised
further out than the anterior one.
"""

import argparse
import os

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--experiment_file_directory', required=True)
    ap.add_argument('--tag', default='final', choices=['raw', 'final'])
    ap.add_argument('--edge_margin_frac', type=float, default=1.0,
                    help='flag ROIs whose centroid is within this many RF '
                         'half-widths of the field edge')
    args = ap.parse_args()

    path = os.path.join(args.experiment_file_directory,
                        'fly.hdf5' if args.tag == 'raw' else 'fly_final.hdf5')
    with h5py.File(path, 'r') as f:
        avg = f[f'STRF_{args.tag}/averaged']
        print(f'align_channel = {avg.attrs.get("align_channel")}   '
              f'z_threshold = {avg.attrs.get("z_threshold")}\n')

        for sn in avg['per_roi']:
            for ch in avg['per_roi'][sn]:
                for rtz_key in avg['per_roi'][sn][ch]:
                    for nt_key in avg['per_roi'][sn][ch][rtz_key]:
                        g = avg['per_roi'][sn][ch][rtz_key][nt_key]
                        cen = g['centroid'][()]
                        inc = g['include'][()]
                        c = cen[inc] if cen.ndim == 1 else cen[inc, 0]

                        # the field the bars actually covered, from the 1D group
                        src = f[f'STRF_{args.tag}/{sn}/{ch}/{rtz_key}/{nt_key}']
                        pos = src['bar_pos_deg'][()]
                        lo, hi = float(pos.min()), float(pos.max())

                        print(f'{sn} {ch} {rtz_key} {nt_key}   '
                              f'{int(inc.sum())}/{len(inc)} ROIs included')
                        print(f'   field covered      {lo:+7.1f} .. {hi:+6.1f} deg '
                              f'(centre {0.5*(lo+hi):+.1f})')
                        print(f'   RF centroids       {c.min():+7.1f} .. {c.max():+6.1f}   '
                              f'mean {c.mean():+.1f}  median {np.median(c):+.1f}')

                        # asymmetry of REACH is expected; quantify it
                        reach_neg = lo - c.max()
                        reach_pos = hi - c.min()
                        print(f'   relative reach     {reach_neg:+7.1f} .. {reach_pos:+6.1f} deg'
                              f'   ({"symmetric" if abs(reach_neg+reach_pos) < 10 else "ASYMMETRIC"})')

                        # asymmetry of SHAPE is the thing that would be an artefact
                        half = 0.5 * np.median(np.abs(c - np.median(c))) * 2 or 10.0
                        margin = args.edge_margin_frac * max(half, 10.0)
                        near_lo = int((c - lo < margin).sum())
                        near_hi = int((hi - c < margin).sum())
                        print(f'   within {margin:.0f} deg of an edge: '
                              f'{near_lo} at the {lo:+.0f} end, {near_hi} at the {hi:+.0f} end',
                              end='')
                        if max(near_lo, near_hi) > 0.25 * inc.sum():
                            print('   <- LIKELY TRUNCATED: those RFs run off the')
                            print('      screen, so their measured flank is cut short and their')
                            print('      centroid is pulled inward. Consider excluding them.')
                        else:
                            print('   (few enough to ignore)')
                        print()


if __name__ == '__main__':
    main()
