"""The model library: where models live and how they are written and read.

A model is one YAML file plus, when it has arrays, one .npz of the same name:

    <library>/<group>/<name>.yaml      everything scalar and textual
    <library>/<group>/<name>.npz       the arrays the YAML refers to

The YAML is a plain mapping. Any value that is a numpy array in memory is
written to the .npz under a key and replaced in the YAML by the string
"@<key>"; loading reverses it. A model written by hand (a hypothetical one)
usually has no arrays and so no .npz.

What a spec contains (spec_version 1)
-------------------------------------
    spec_version: 1
    name:        the file name, no extension
    kind:        'LN' | 'SwitchedLN'
    output:      {channel, signal: 'imaged' | 'neural', units, indicator}
                 signal says whether the model's output is the measured
                 fluorescence ('imaged', nothing more to apply) or the cell's
                 own signal ('neural', an indicator kernel may be applied).
    provenance:  {source: 'pipeline' | 'group' | 'toy' | 'hypothetical',
                  created, description, ...}  free beyond those three; a
                 pipeline export records fly, series, channel, gate, settings
                 and the stimulus it was fitted to, a group export its members
                 and their hashes, a toy the model it idealised and the fits.
    LN:          filter: {...}, nonlinearity: {...}     (see filters.py,
                 nonlinearities.py for the kinds)
    SwitchedLN:  label: <label name>,
                 states: {<label value as str>: {name, filter, nonlinearity}}

Where the library is
--------------------
resolve_library() looks, in order, at an explicit argument, the environment
variable VISANALYSIS_MODEL_LIBRARY, and `model_library:` in
~/.visanalysis/config.yaml. There is no default: a model written somewhere
nobody chose is a model nobody will find.
"""
import datetime
import hashlib
import io
import os
import pathlib

import numpy as np
import yaml

SPEC_VERSION = 1
ENV_VAR = "VISANALYSIS_MODEL_LIBRARY"
CONFIG_PATH = pathlib.Path("~/.visanalysis/config.yaml")
CONFIG_KEY = "model_library"

__all__ = ["resolve_library", "save_model", "load_model", "list_models", "model_hash",
           "new_spec", "LibraryError"]


class LibraryError(RuntimeError):
    pass


# ── location ────────────────────────────────────────────────────────────────

def resolve_library(path=None, create=False):
    """The library directory as a Path, from the first of: path, $VISANALYSIS_MODEL_LIBRARY,
    ~/.visanalysis/config.yaml. create=True makes the directory."""
    source = None
    if path:
        lib, source = pathlib.Path(path), "argument"
    elif os.environ.get(ENV_VAR):
        lib, source = pathlib.Path(os.environ[ENV_VAR]), "environment variable " + ENV_VAR
    else:
        cfg = CONFIG_PATH.expanduser()
        if cfg.exists():
            with open(cfg, "r", encoding="utf-8") as f:
                data = yaml.safe_load(f) or {}
            if not isinstance(data, dict) or not data.get(CONFIG_KEY):
                raise LibraryError("{} exists but has no '{}:' entry".format(cfg, CONFIG_KEY))
            lib, source = pathlib.Path(str(data[CONFIG_KEY])), str(cfg)
        else:
            raise LibraryError(
                "no model library is set. Give one explicitly, set the environment variable {}, or put\n"
                "    {}: <path>\n"
                "in {}".format(ENV_VAR, CONFIG_KEY, cfg))
    lib = lib.expanduser()
    if create:
        lib.mkdir(parents=True, exist_ok=True)
    elif not lib.is_dir():
        raise LibraryError("model library {} (from {}) is not a directory".format(lib, source))
    return lib


def _target(name, library):
    """(yaml path, npz path) for a model name, which may carry a subdirectory ('contrast_switch/x')."""
    name = str(name).replace("\\", "/")
    if name.endswith((".yaml", ".yml", ".npz")):
        name = name.rsplit(".", 1)[0]
    p = pathlib.Path(library) / name
    return p.with_suffix(".yaml"), p.with_suffix(".npz")


# ── (de)serialisation ───────────────────────────────────────────────────────

def _pack(obj, arrays, path):
    """Replace arrays in a nested structure by '@key' references; collect them."""
    if isinstance(obj, np.ndarray):
        key = ".".join(path) if path else "array"
        base, k = key, 1
        while key in arrays:
            k += 1
            key = "{}_{}".format(base, k)
        arrays[key] = obj
        return "@" + key
    if isinstance(obj, dict):
        return {str(k): _pack(v, arrays, path + [str(k)]) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        if obj and all(isinstance(v, (int, float, np.integer, np.floating)) for v in obj) and len(obj) > 16:
            return _pack(np.asarray(obj), arrays, path)
        return [_pack(v, arrays, path + [str(i)]) for i, v in enumerate(obj)]
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, pathlib.PurePath):
        return str(obj)
    return obj


def _unpack(obj, arrays):
    if isinstance(obj, str) and obj.startswith("@"):
        key = obj[1:]
        if key not in arrays:
            raise LibraryError("spec refers to array '{}' that the .npz does not hold".format(key))
        return np.asarray(arrays[key])
    if isinstance(obj, dict):
        return {k: _unpack(v, arrays) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_unpack(v, arrays) for v in obj]
    return obj


def _yaml_text(spec):
    return yaml.safe_dump(spec, sort_keys=False, default_flow_style=False, allow_unicode=True, width=100)


def new_spec(name, kind, source, description, channel=None, signal="imaged", units="dF/F",
             indicator=None, **provenance):
    """A spec skeleton with the required fields filled in."""
    return dict(spec_version=SPEC_VERSION, name=str(name), kind=str(kind),
                output=dict(channel=channel, signal=signal, units=units, indicator=indicator),
                provenance=dict(source=str(source),
                                created=datetime.datetime.now().isoformat(timespec="seconds"),
                                description=str(description), **provenance))


def save_model(spec, name=None, library=None, overwrite=False):
    """Write a spec to the library. Returns the YAML path.

    spec is a plain dict (see the module docstring); arrays anywhere inside it go to the .npz.
    name defaults to spec['name'] and may include a subdirectory."""
    lib = resolve_library(library, create=True)
    name = spec.get("name") if name is None else name
    if not name:
        raise LibraryError("a model needs a name")
    spec = dict(spec, name=str(name).replace("\\", "/").rsplit("/", 1)[-1])
    spec.setdefault("spec_version", SPEC_VERSION)
    ypath, npath = _target(name, lib)
    if ypath.exists() and not overwrite:
        raise LibraryError("{} exists; pass overwrite=True to replace it".format(ypath))
    arrays = {}
    packed = _pack(spec, arrays, [])
    packed["library_hash"] = _hash_parts(packed, arrays, skip_hash=True)
    ypath.parent.mkdir(parents=True, exist_ok=True)
    with open(ypath, "w", encoding="utf-8", newline="\n") as f:
        f.write(_yaml_text(packed))
    if arrays:
        np.savez(npath, **arrays)
    elif npath.exists():
        npath.unlink()
    return ypath


def load_model(name, library=None):
    """Read a spec back, arrays included. name is a library-relative name or a path to a .yaml."""
    p = pathlib.Path(str(name))
    if p.suffix in (".yaml", ".yml") and p.exists():
        ypath, npath = p, p.with_suffix(".npz")
    else:
        ypath, npath = _target(name, resolve_library(library))
    if not ypath.exists():
        raise LibraryError("no model at {}".format(ypath))
    with open(ypath, "r", encoding="utf-8") as f:
        packed = yaml.safe_load(f)
    if not isinstance(packed, dict):
        raise LibraryError("{} does not hold a mapping".format(ypath))
    arrays = {}
    if npath.exists():
        with np.load(npath, allow_pickle=False) as z:
            arrays = {k: z[k] for k in z.files}
    spec = _unpack(packed, arrays)
    if spec.get("spec_version", 1) > SPEC_VERSION:
        raise LibraryError("{} is spec_version {}, newer than this code ({})".format(
            ypath, spec["spec_version"], SPEC_VERSION))
    spec["_path"] = str(ypath)
    return spec


def list_models(library=None, group=None):
    """Names (library-relative, no extension) of every model, optionally within one subdirectory."""
    lib = resolve_library(library)
    root = lib / group if group else lib
    out = []
    for p in sorted(root.rglob("*.yaml")):
        out.append(p.relative_to(lib).with_suffix("").as_posix())
    return out


# ── identity ────────────────────────────────────────────────────────────────

def _hash_parts(packed, arrays, skip_hash=False):
    h = hashlib.sha256()
    body = {k: v for k, v in packed.items() if not (skip_hash and k == "library_hash")}
    body.pop("_path", None)
    h.update(_yaml_text(body).encode("utf-8"))
    for k in sorted(arrays):
        a = np.ascontiguousarray(arrays[k])
        h.update(k.encode("utf-8"))
        h.update(str(a.dtype).encode("utf-8"))
        h.update(str(a.shape).encode("utf-8"))
        h.update(a.tobytes())
    return h.hexdigest()[:16]


def model_hash(spec):
    """A short content hash of a spec (arrays included), the same whether it has been saved or not.
    A group model records its members' hashes so a re-fitted member is detectable."""
    arrays = {}
    packed = _pack({k: v for k, v in spec.items() if k not in ("_path", "library_hash")}, arrays, [])
    return _hash_parts(packed, arrays)
