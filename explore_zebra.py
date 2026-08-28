#!/usr/bin/env python3
"""Play with Zebra noise parameters offline, then take the winners to the rig.

Edit the PARAMS block below, run the file, look at the picture. Nothing here talks to hardware.

    python explore_zebra.py                 one figure at the parameters below
    python explore_zebra.py --fast          the same, but seconds instead of a minute
    python explore_zebra.py --sweep scale 0.8,1.09,1.5
    python explore_zebra.py --defaults      show the protocol's current defaults and exit

The parameters here are the SAME ones JCS_protocol takes, so anything you settle on transfers to
the rig unchanged -- there is no calibration step in between to re-interpret them.

Two rules of thumb for the knobs:

    clusters (deg) ~= 15.9 / scale               larger scale = finer
    cluster time (s) ~= 0.135 / w_per_second     larger rate  = faster
    stripes (deg)  ~= clusters / (n_teeth / 2)

@author: jcsimon
"""
import argparse
import os
import sys

# --------------------------------------------------------------------------
# EDIT ME
# --------------------------------------------------------------------------
PARAMS = dict(
    scale=1.09,          # spatial: clusters ~15 deg, stripes ~16 deg with 4 teeth
    n_teeth=4,           # stripes per cluster; even integer
    w_per_second=0.056,  # temporal: clusters decorrelate in ~2.4 s
    persistence=0.2,     # fine-scale power; higher makes stripe edges more convoluted
    octaves=6,
    seed=0,

    # Screen geometry. Sized to the BrukerJr bowl: it lights -67..+67 azimuth and
    # -49..+38 elevation, and width/height is kept at n_cols/n_rows so texels are
    # square in degrees.
    width=136.0,
    height=102.0,
    n_rows=384,
    n_cols=512,
    gen_n_rows=96,
    gen_n_cols=128,
    fps=120.0,
)

OUT = 'zebra_explore.png'

# Where clandinin_labpack lives, if it is not already importable. Empty means "look next to
# this repo", which is right whenever both are checked out into the same folder -- set it
# explicitly if yours are somewhere else.
LABPACK = ''


def _sibling(name, marker):
    """A repo checked out next to this one, found by a directory it must contain."""
    here = os.path.dirname(os.path.abspath(__file__))
    for candidate in (os.path.join(os.path.dirname(here), name),
                      os.path.expanduser(f'~/Documents/GitHub/{name}')):
        if os.path.isdir(os.path.join(candidate, marker)):
            return candidate
    return ''


def _labpack_dir():
    return LABPACK or _sibling('clandinin_labpack', 'labpack')


def _ensure_imports():
    """Put labpack and visanalysis on the path if the caller has not."""
    here = os.path.dirname(os.path.abspath(__file__))
    # stimpack too: --defaults imports the protocol module, which imports stimpack.
    for path in (here, _labpack_dir(), _sibling('stimpack', 'stimpack')):
        if path and os.path.isdir(path) and path not in sys.path:
            sys.path.insert(0, path)
    try:
        from visanalysis.util import zebra_explorer            # noqa: F401
        from labpack.visual_stim.clandinin import zebra_noise  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            f"could not import the explorer or the generator ({exc}).\n"
            f"Set LABPACK at the top of this file to your clandinin_labpack checkout, "
            f"or put both repos on PYTHONPATH.")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--fast', action='store_true',
                    help='few frames, few seeds, no anisotropy escalation. Seconds rather '
                         'than a minute; the orientation numbers will read UNRESOLVED, which '
                         'is honest rather than a failure')
    ap.add_argument('--sweep', nargs=2, metavar=('PARAM', 'VALUES'),
                    help='compare one parameter across comma-separated values, '
                         'e.g. --sweep n_teeth 2,4,8')
    ap.add_argument('--defaults', action='store_true',
                    help="print JCS_protocol's current defaults and exit, so you can see "
                         "what the rig would do right now")
    ap.add_argument('-o', '--out', default=OUT)
    args = ap.parse_args()

    _ensure_imports()
    from visanalysis.util import zebra_explorer as ze

    if args.defaults:
        import importlib.util
        path = os.path.join(_labpack_dir(), 'labpack', 'protocol', 'JCS_protocol.py')
        spec = importlib.util.spec_from_file_location('jcs_defaults', path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        try:
            spec.loader.exec_module(mod)
        except ImportError as exc:
            raise SystemExit(
                f'reading the protocol defaults needs stimpack importable ({exc}). '
                f'Everything else in this script works without it.')
        proto = mod.ZebraNoiseSphere({'current_rig_name': 'x', 'rig_config': {}})
        print(f'{path}\n')
        for k, v in proto.get_protocol_parameter_defaults().items():
            print(f'  {k:<22} {v!r}')
        return 0

    # Cheap settings while hunting; the full ones once something looks right.
    sampling = (dict(n_frames=256, n_seeds=2, min_snr=0.0, max_frames=256)
                if args.fast else dict(n_frames=1024, n_seeds=8, min_snr=3.0))

    if args.sweep:
        name, raw = args.sweep
        name = name.replace('-', '_')
        if name not in PARAMS:
            raise SystemExit(f'unknown parameter {name!r}; pick from {sorted(PARAMS)}')
        cast = type(PARAMS[name])
        values = [cast(v) for v in raw.split(',')]
        sets = [dict(PARAMS, **{name: v}) for v in values]
        ze.compare(sets, args.out, labels=[f'{name}={v:g}' for v in values],
                   title=f'Zebra noise: {name} sweep',
                   n_frames=sampling['n_frames'], n_seeds=sampling['n_seeds'],
                   min_snr=0.0)
        return 0

    print('rendering...')
    result = ze.analyze(**sampling, **PARAMS)
    ze.plot(result, args.out)

    s = result['summary']
    print()
    print('  what this stimulus is')
    print(f"    stripe width          {s['stripe_deg']:.1f} deg")
    print(f"    spatial correlation   {s['stimulus_corr_deg']:.1f} deg")
    print(f"    temporal correlation  {s['stimulus_corr_s']:.2f} s")
    print(f"    speed                 {s['speed_dps']:.0f} deg/s")
    print(f"    median spatial freq   {s['median_spatial_cpd']:.3f} cycles/deg")
    print(f"    median temporal freq  {s['median_temporal_hz']:.2f} Hz")
    print(f"    white fraction        {s['white_fraction']:.3f}")
    print()
    resolved = s.get('anisotropy_resolved', False)
    print(f"    orientation spread    {s['anisotropy']:.3f} "
          f"(floor {s['anisotropy_floor']:.3f}, "
          f"{'resolved' if resolved else 'UNRESOLVED -- too few samples to judge'})")
    print()
    print('  to run this on the rig, set these in the ZebraNoiseSphere protocol:')
    for k in ('scale', 'n_teeth', 'w_per_second', 'persistence', 'octaves'):
        print(f"    {k:<16} {PARAMS[k]}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
