"""
repack_fly.py

Reclaim the dead space an hdf5 file accumulates when groups are overwritten.

WHY THIS IS NEEDED AT ALL. Deleting an hdf5 object does not return its bytes to
the operating system, and usually not even to the file: `del f[group]` marks the
space free inside the container, the file stays the same size, and the next
write generally appends rather than reusing the hole. The analysis scripts
rewrite their whole results group on every run that is not a reuse -- see the
`if strf_group in f: del f[strf_group]` in analyze_data_strf_one_d.py -- so a
fly that has been reanalysed a few times carries several abandoned copies of a
several-hundred-megabyte group. On fly_008 one such group is 831 MB.

The only way to reclaim it is to write a fresh file, which is what this does.

WHY NOT THE h5repack BINARY. It is the canonical tool, and this would use it if
one were available that was built against the SAME hdf5 the analysis writes
with. On this machine the only h5repack lives in a different conda env and is
hdf5 1.12.1, while h5py here is 1.14.2; running an older tool over a newer file
is a risk worth avoiding on primary data. Copying through h5py is
version-matched by construction. An H5Ocopy preserves compression, chunking and
attributes exactly -- this is a repack, not a re-encode, and no dataset is
decompressed and rewritten.

SAFETY. Rewriting a subject's only data file has to be boring:

* nothing happens unless the file actually has enough dead space to be worth it
  (--min_waste_mb), so this is safe to leave in a shell and run every time
* the free-space check refuses rather than filling the drive, which matters
  here because peak usage is about twice the file size
* the new file is built beside the old one under a temporary name and VERIFIED
  -- every dataset path, shape and dtype, every attribute, and the contents of
  small datasets plus the first and last slice of large ones -- before anything
  replaces anything
* the original is moved aside, not deleted, and is removed only once the
  replacement is in place. --keep_backup leaves it.

Usage:
    python repack_fly.py --experiment_file_directory "path/to/fly_folder"
    python repack_fly.py --experiment_file_directory "..." --dry_run True
"""

import argparse
import io
import os
import shutil
import sys
import time

import h5py
import numpy as np

#: Shared with the analysis scripts on purpose: one file per fly holds the whole
#: history of what was run against it, in order.
ANALYSIS_LOG_NAME = 'analysis_log.txt'

#: Peak disk use is the original plus the copy. The margin covers the copy
#: landing slightly larger than predicted and leaves the drive room to breathe.
SPACE_MARGIN = 2.1

#: Datasets at or below this many elements are compared in full; larger ones on
#: their first and last slice, which catches truncation and corruption without
#: doubling the read of a multi-gigabyte file.
FULL_COMPARE_ELEMS = 1000000


class _Tee:
    """Mirror a stream to a log file, flushed on every write so a run that dies
    partway through still leaves a complete log."""

    def __init__(self, stream, fh):
        self._stream, self._fh = stream, fh

    def write(self, s):
        self._stream.write(s)
        self._fh.write(s)
        self._fh.flush()
        return len(s)

    def flush(self):
        self._stream.flush()
        self._fh.flush()

    def isatty(self):
        return getattr(self._stream, 'isatty', lambda: False)()


def start_log(directory, argv):
    """Tee into the fly's analysis log, or carry on without one.

    The log is a record, not the job. A folder that is read-only to whoever is
    running -- another lab member's data, a share mounted without write access
    -- must not stop a space-reclaiming utility from reporting what it found,
    so a log that cannot be opened is a warning rather than a traceback.
    """
    path = os.path.join(directory, ANALYSIS_LOG_NAME)
    try:
        fh = io.open(path, 'a', encoding='utf-8')
    except (IOError, OSError) as e:
        print('(no analysis log: {}; continuing without one)'.format(e))
        return None
    sys.stdout = _Tee(sys.stdout, fh)
    sys.stderr = _Tee(sys.stderr, fh)
    print('')
    print('=' * 78)
    print('analysis : hdf5 repack (repack_fly.py)')
    print('run at   : ' + time.strftime('%Y-%m-%d %H:%M:%S'))
    print('command  : ' + ' '.join(argv))
    print('=' * 78)
    return path


def storage_bytes(f):
    """Bytes actually occupied by datasets, which is what the file WOULD be if
    it held nothing else. The gap to the file size is the dead space."""
    total = [0]

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            total[0] += obj.id.get_storage_size()

    f.visititems(visit)
    return total[0]


def inventory(f):
    """Every object's identity, for comparing two files structurally."""
    out = {}

    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            out[name] = ('dataset', obj.shape, str(obj.dtype))
        else:
            out[name] = ('group', None, None)

    f.visititems(visit)
    out['/'] = ('root', None, None)
    return out


def _attrs_equal(a, b):
    if sorted(a.keys()) != sorted(b.keys()):
        return False
    for k in a.keys():
        x, y = a[k], b[k]
        if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
            if not np.array_equal(np.asarray(x), np.asarray(y)):
                return False
        elif x != y:
            return False
    return True


def _data_equal(oa, ob):
    try:
        if oa.size <= FULL_COMPARE_ELEMS:
            return bool(np.array_equal(oa[()], ob[()], equal_nan=True))
        return bool(np.array_equal(oa[0], ob[0], equal_nan=True)
                    and np.array_equal(oa[-1], ob[-1], equal_nan=True))
    except TypeError:
        # equal_nan is only defined for floats; strings and opaque types land
        # here, and an exact comparison is the right test for them anyway.
        if oa.size <= FULL_COMPARE_ELEMS:
            return bool(np.array_equal(np.asarray(oa[()]), np.asarray(ob[()])))
        return bool(np.array_equal(np.asarray(oa[0]), np.asarray(ob[0]))
                    and np.array_equal(np.asarray(oa[-1]), np.asarray(ob[-1])))


def verify(src_path, dst_path):
    """Structural and sampled-content comparison. Returns a list of problems."""
    problems = []
    with h5py.File(src_path, 'r') as a, h5py.File(dst_path, 'r') as b:
        ia, ib = inventory(a), inventory(b)
        missing = sorted(set(ia) - set(ib))
        extra = sorted(set(ib) - set(ia))
        if missing:
            problems.append('missing from the copy: {}'.format(missing[:8]))
        if extra:
            problems.append('unexpected in the copy: {}'.format(extra[:8]))
        for name in sorted(set(ia) & set(ib)):
            if ia[name] != ib[name]:
                problems.append('{}: {} -> {}'.format(name, ia[name], ib[name]))
                continue
            oa = a['/'] if name == '/' else a[name]
            ob = b['/'] if name == '/' else b[name]
            if not _attrs_equal(oa.attrs, ob.attrs):
                problems.append('{}: attributes differ'.format(name))
            if isinstance(oa, h5py.Dataset) and oa.size and not _data_equal(oa, ob):
                problems.append('{}: contents differ'.format(name))
    return problems


def repack(path, min_waste_mb=200.0, keep_backup=False, dry_run=False):
    """Rewrite `path` without its dead space. Returns True if it was rewritten."""
    size = os.path.getsize(path)
    with h5py.File(path, 'r') as f:
        used = storage_bytes(f)
        top = sorted(f.keys())
    waste = size - used
    print('{}'.format(path))
    print('  file size      {:9.1f} MB'.format(size / 1e6))
    print('  datasets use   {:9.1f} MB'.format(used / 1e6))
    print('  dead space     {:9.1f} MB  ({:.1f}% of the file)'.format(
        waste / 1e6, 100.0 * waste / max(size, 1)))
    print('  top-level      {}'.format(', '.join(top)))

    if waste / 1e6 < float(min_waste_mb):
        print('  -> below --min_waste_mb ({:g} MB); nothing to do'.format(
            float(min_waste_mb)))
        return False
    if dry_run:
        print('  -> --dry_run: would reclaim about {:.1f} MB'.format(waste / 1e6))
        return False

    free = shutil.disk_usage(os.path.dirname(os.path.abspath(path))).free
    need = int(size * SPACE_MARGIN)
    print('  free on drive  {:9.1f} MB   need about {:.1f} MB'.format(
        free / 1e6, need / 1e6))
    if free < need:
        print('  -> REFUSING: not enough free space to hold the original and '
              'the copy at once. Free {:.1f} MB and retry.'.format(
                  (need - free) / 1e6))
        return False

    tmp = path + '.repack.tmp'
    bak = path + '.bak'
    for stale in (tmp, bak):
        if os.path.exists(stale):
            print('  removing stale {}'.format(os.path.basename(stale)))
            os.remove(stale)

    t0 = time.perf_counter()
    print('  copying...')
    try:
        with h5py.File(path, 'r') as src, h5py.File(tmp, 'w') as dst:
            for k, v in src.attrs.items():
                dst.attrs[k] = v
            for k in src.keys():
                # H5Ocopy: preserves compression, chunking and attributes, and
                # never decompresses. This is a repack, not a re-encode.
                src.copy(k, dst, name=k)
        print('  copied in {:.1f} s; verifying...'.format(
            time.perf_counter() - t0))
        problems = verify(path, tmp)
    except Exception as e:
        print('  -> FAILED during copy/verify: {}: {}'.format(
            type(e).__name__, e))
        if os.path.exists(tmp):
            os.remove(tmp)
        print('  original untouched.')
        raise

    if problems:
        print('  -> VERIFICATION FAILED, original untouched:')
        for p in problems[:12]:
            print('       {}'.format(p))
        os.remove(tmp)
        return False

    new_size = os.path.getsize(tmp)
    # Moved aside rather than deleted, so a failure between the two renames
    # leaves the data on disk under a name that says what it is.
    os.replace(path, bak)
    os.replace(tmp, path)
    print('  verified; replaced.')
    print('  {:.1f} MB -> {:.1f} MB   reclaimed {:.1f} MB ({:.1f}%)'.format(
        size / 1e6, new_size / 1e6, (size - new_size) / 1e6,
        100.0 * (size - new_size) / max(size, 1)))
    if keep_backup:
        print('  backup kept at {}'.format(os.path.basename(bak)))
    else:
        os.remove(bak)
        print('  backup removed.')
    return True


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--experiment_file_directory', nargs='?')
    p.add_argument('--experiment_file_name', nargs='?', default='fly.hdf5')
    p.add_argument('--min_waste_mb', nargs='?', type=float, default=200.0,
                   help='do nothing unless at least this much dead space would '
                        'be reclaimed. Makes the step safe to run every time.')
    p.add_argument('--keep_backup', nargs='?', default='False',
                   help='keep the pre-repack file as fly.hdf5.bak')
    p.add_argument('--dry_run', nargs='?', default='False',
                   help='report what would be reclaimed and stop')
    a = p.parse_args()

    path = os.path.join(a.experiment_file_directory, a.experiment_file_name)
    if not os.path.exists(path):
        print('No such file: ' + path)
        sys.exit(1)

    start_log(a.experiment_file_directory, sys.argv)
    repack(path, min_waste_mb=a.min_waste_mb,
           keep_backup=str(a.keep_backup) == 'True',
           dry_run=str(a.dry_run) == 'True')
