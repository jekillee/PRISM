"""
TRANSP data loaders (CDF output + U-File input)

Two parts, kept in one module:

- CDF output   : ``load_transp_cdf`` reads a TRANSP output .CDF (netCDF) file and
                 categorizes variables into profiles (time × radius) and time
                 traces (time). Plus the netCDF backend helpers.
- U-File input : ``scan_run_directory`` / ``get_catalog`` / ``pretty_label`` /
                 ``list_extensions_by_class`` scan a run directory of
                 ``<prefix><shot>.<ext>`` U-Files (via core.ufile_parser) and group
                 them into a run record.
"""

import os
import re
import time

import numpy as np

from core.ufile_parser import parse_ufile, classify_ufile


# =============================================================================
# TRANSP CDF output loader
# =============================================================================

def load_transp_cdf(filepath):
    """Load a TRANSP CDF file and categorize variables.

    Args:
        filepath: Path to .CDF file

    Returns:
        dict with keys:
            filepath, label,
            profiles: {name: {long_name, units, time, x, x_label, data}}
            timetraces: {name: {long_name, units, time, data}}
    """
    nc = _open_netcdf(filepath)
    label = os.path.splitext(os.path.basename(filepath))[0]

    # Identify dimension types
    time_dims = set()
    spatial_dims = set()
    for dim_name in _get_dim_names(nc):
        dn = dim_name.upper()
        if 'TIME' in dn:
            time_dims.add(dim_name)
        if dn in ('X', 'XB', 'X2', 'XB2', 'XRHO', 'XRHOB'):
            spatial_dims.add(dim_name)

    profiles = {}
    timetraces = {}

    for var_name in _get_var_names(nc):
        var = nc.variables[var_name]
        dims = _get_dims(var)

        # Skip coordinate variables and non-numeric
        data = _get_data(var)
        if data.dtype.kind not in ('f', 'i', 'u'):
            continue

        if len(dims) == 1 and dims[0] in time_dims:
            time_arr = _get_coord(nc, dims[0], data.shape[0])
            timetraces[var_name] = {
                'long_name': _get_attr(var, 'long_name', var_name),
                'units': _get_attr(var, 'units', ''),
                'time': time_arr.astype(np.float64),
                'data': data.astype(np.float64),
            }

        elif len(dims) == 2:
            has_time = any(d in time_dims for d in dims)
            has_spatial = any(d in spatial_dims for d in dims)
            if has_time and has_spatial:
                # Orient as (time, x)
                if dims[0] in time_dims:
                    t_dim, x_dim = dims[0], dims[1]
                else:
                    t_dim, x_dim = dims[1], dims[0]
                    data = data.T

                time_arr = _get_coord(nc, t_dim, data.shape[0])
                x = _get_coord(nc, x_dim, data.shape[1])
                if time_arr.ndim == 2:
                    time_arr = time_arr[:, 0]
                if x.ndim == 2:
                    x = x[0, :]

                profiles[var_name] = {
                    'long_name': _get_attr(var, 'long_name', var_name),
                    'units': _get_attr(var, 'units', ''),
                    'time': time_arr.astype(np.float64),
                    'x': x.astype(np.float64),
                    'x_dim': x_dim,
                    'data': data.astype(np.float64),
                }

    # Extract coordinate mapping arrays for x-axis conversion
    coords = {}
    for cname in ('PLFLX', 'RMNMP', 'RMJMP', 'RMAJM'):
        if cname in nc.variables:
            cdata = _get_data(nc.variables[cname]).astype(np.float64)
            coords[cname] = cdata

    _close_netcdf(nc)
    print(f"[TRANSP] Loaded {label}: {len(profiles)} profiles, "
          f"{len(timetraces)} time traces")
    return {
        'filepath': filepath,
        'label': label,
        'profiles': profiles,
        'timetraces': timetraces,
        'coords': coords,
    }


# ---- netCDF backend abstraction ----

def _open_netcdf(filepath):
    try:
        from netCDF4 import Dataset
        return Dataset(filepath, 'r')
    except ImportError:
        from scipy.io import netcdf_file
        return netcdf_file(filepath, 'r', mmap=False)


def _close_netcdf(nc):
    nc.close()


def _get_dim_names(nc):
    return list(nc.dimensions.keys())


def _get_var_names(nc):
    return list(nc.variables.keys())


def _get_dims(var):
    if hasattr(var, 'dimensions'):
        return var.dimensions
    return ()


def _get_data(var):
    try:
        return np.array(var[:])
    except Exception:
        return np.array(var.data)


def _get_coord(nc, dim_name, expected_len):
    if dim_name in nc.variables:
        return _get_data(nc.variables[dim_name]).astype(np.float64)
    return np.arange(expected_len, dtype=np.float64)


def _get_attr(var, attr_name, default=''):
    try:
        val = getattr(var, attr_name, default)
        if isinstance(val, bytes):
            return val.decode('utf-8', errors='replace').strip()
        if isinstance(val, np.ndarray):
            if val.dtype.kind in ('U', 'S', 'O'):
                return str(val.flat[0]).strip()
            return str(val)
        return str(val).strip() if val else default
    except Exception:
        return default


# =============================================================================
# TRANSP U-File directory loader
# =============================================================================

_FNAME_RE = re.compile(r'^([A-Za-z]+)(\d+)\.([A-Za-z0-9]{2,4})$')


def _make_run_key(shot, prefix, label):
    return f"{shot}_{prefix}_{label}"


def _make_label(dir_path, mtime):
    """Derive a human-readable label for a run directory."""
    base = os.path.basename(os.path.normpath(dir_path))
    if base and not base.isdigit():
        return base
    # Fall back to mtime
    return time.strftime("%y%m%d-%H%M", time.localtime(mtime))


def scan_run_directory(dir_path, verbose=True):
    """Scan one directory and return a run record (or None if empty)."""
    if not os.path.isdir(dir_path):
        return None

    candidates = []
    for name in sorted(os.listdir(dir_path)):
        full = os.path.join(dir_path, name)
        if not os.path.isfile(full):
            continue
        m = _FNAME_RE.match(name)
        if not m:
            continue
        candidates.append((m.group(1), int(m.group(2)), m.group(3).upper(), full))

    if not candidates:
        return None

    # Pick dominant (prefix, shot) pair (most files share it)
    pairs = {}
    for prefix, shot, ext, _ in candidates:
        pairs[(prefix, shot)] = pairs.get((prefix, shot), 0) + 1
    (prefix, shot) = max(pairs, key=pairs.get)

    ufiles = {}
    classes = {}
    skipped = []
    for fprefix, fshot, ext, full in candidates:
        if fprefix != prefix or fshot != shot:
            continue
        if ext in ufiles:
            continue  # already loaded
        try:
            uf = parse_ufile(full)
            if uf is None:
                continue
            ufiles[ext] = uf
            classes[ext] = classify_ufile(uf)
        except Exception as e:
            skipped.append((ext, str(e)))
            if verbose:
                print(f"[UFILE] Failed to parse {os.path.basename(full)}: {e}")

    if not ufiles:
        return None

    mtime = max(os.path.getmtime(p) for p in
                (os.path.join(dir_path, f"{prefix}{shot}.{e.upper()}")
                 for e in ufiles)
                if os.path.isfile(p))

    label = _make_label(dir_path, mtime)
    run_key = _make_run_key(shot, prefix, label)

    if verbose:
        n_classes = {}
        for cls in classes.values():
            n_classes[cls] = n_classes.get(cls, 0) + 1
        cls_summary = ", ".join(f"{k}={v}" for k, v in n_classes.items())
        print(f"[UFILE] {prefix}{shot} · {label}: {len(ufiles)} files "
              f"({cls_summary})")
        for ext in sorted(ufiles):
            uf = ufiles[ext]
            x0 = uf.get('X0', {})
            x1 = uf.get('X1', {})
            shape = uf.get('F', {}).get('data')
            shape = shape.shape if hasattr(shape, 'shape') else ()
            print(f"  {ext:>4}  NDIM={uf.get('NDIM')}  shape={shape}  "
                  f"X0={x0.get('name','')!r}/{x0.get('units','')!r}  "
                  f"X1={x1.get('name','')!r}/{x1.get('units','')!r}  "
                  f"-> {classes[ext]}")

    return {
        'run_key': run_key,
        'shot': shot,
        'prefix': prefix,
        'label': label,
        'dir': dir_path,
        'mtime': mtime,
        'ufiles': ufiles,
        'classes': classes,
        'skipped': skipped,
    }


def list_extensions_by_class(run, target_classes):
    """Return sorted list of extensions in `run` matching given classes."""
    if isinstance(target_classes, str):
        target_classes = (target_classes,)
    target = set(target_classes)
    return sorted([ext for ext, cls in run['classes'].items() if cls in target])


# ----- Catalog (optional pretty labels) ---------------------------------------

_CATALOG_CACHE = None


def get_catalog():
    """Load (and cache) the optional `config/transp_ufiles.json` catalog."""
    global _CATALOG_CACHE
    if _CATALOG_CACHE is not None:
        return _CATALOG_CACHE
    import json
    cfg_path = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'config', 'transp_ufiles.json'
    )
    try:
        with open(cfg_path, 'r') as f:
            _CATALOG_CACHE = json.load(f)
    except Exception:
        _CATALOG_CACHE = {}
    return _CATALOG_CACHE


def pretty_label(ext, uf):
    """Build a display label for a UFILE entry: '{ext} {name} [{units}]'."""
    cat = get_catalog().get(ext.upper(), {})
    name = cat.get('label') or (uf.get('F', {}).get('name') or '').strip()
    units = (uf.get('F', {}).get('units') or '').strip()
    parts = [ext]
    if name and name != ext:
        parts.append(name)
    if units:
        parts.append(f"[{units}]")
    return '  '.join(parts)
