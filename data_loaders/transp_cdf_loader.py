"""
TRANSP CDF file loader

Reads TRANSP output CDF (netCDF) files and categorizes variables into
profile (2D: time × radius) and time trace (1D: time) data.
"""

import os
import numpy as np


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
        if dn in ('X', 'XB', 'X2', 'XB2'):
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
            time = _get_coord(nc, dims[0], data.shape[0])
            timetraces[var_name] = {
                'long_name': _get_attr(var, 'long_name', var_name),
                'units': _get_attr(var, 'units', ''),
                'time': time.astype(np.float64),
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

                time = _get_coord(nc, t_dim, data.shape[0])
                x = _get_coord(nc, x_dim, data.shape[1])
                if time.ndim == 2:
                    time = time[:, 0]
                if x.ndim == 2:
                    x = x[0, :]

                profiles[var_name] = {
                    'long_name': _get_attr(var, 'long_name', var_name),
                    'units': _get_attr(var, 'units', ''),
                    'time': time.astype(np.float64),
                    'x': x.astype(np.float64),
                    'x_label': x_dim,
                    'data': data.astype(np.float64),
                }

    _close_netcdf(nc)
    print(f"[TRANSP] Loaded {label}: {len(profiles)} profiles, "
          f"{len(timetraces)} time traces")
    return {
        'filepath': filepath,
        'label': label,
        'profiles': profiles,
        'timetraces': timetraces,
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
