"""
TRANSP U-File parser (text-based 1D / 2D / 3D format).

Ported from OMFIT's OMFITuFile.load() — minimal, dependency-free.

Parsed object is a plain dict with keys:
    SHOT, DEVICE, NDIM, COMPRESSION, LSHOT, DATE,
    SCALARS  : list of {'name', 'desc', 'data'}
    X0, X1, X2 : {'name', 'units', 'data': ndarray}  (only those that exist)
    F        : {'name', 'units', 'data': ndarray}    (shape == tuple of len(Xi))
    PROCESSING_CODE, COMMENTS, BEAMS (optional)
    SOURCE   : {'path', 'prefix', 'shot', 'ext'}
"""

import os
import re
import numpy as np

_HEADER0_RE = re.compile(
    r'\s*(?P<shot>\d+)\s*(?P<device>\w+)\s+(?P<ndim>\d)\s+(?P<comp>\d)\s+(?P<lshot>\d*)'
)


def parse_ufile(path):
    """Parse a TRANSP U-File and return a dict (see module docstring)."""
    with open(path, 'r', errors='replace') as f:
        lines = f.readlines()

    if not lines:
        return None

    out = {}
    m = _HEADER0_RE.match(lines[0])
    if not m:
        raise ValueError(f"Not a U-File (bad header): {os.path.basename(path)}")
    out['SHOT'] = int(m.group('shot'))
    out['DEVICE'] = m.group('device').strip()
    out['NDIM'] = int(m.group('ndim'))
    out['COMPRESSION'] = int(m.group('comp'))
    try:
        out['LSHOT'] = max(6, int(m.group('lshot') or 0))
    except Exception:
        out['LSHOT'] = 0

    out['DATE'] = lines[1].split(';')[0].strip()
    nscal = int(lines[2].split(';')[0].strip())

    out['SCALARS'] = []
    n = 2
    end_scalars = 3
    for k in range(nscal):
        try:
            value = float(lines[3 + 2 * k].split(';')[0].strip())
            name_desc = lines[3 + 2 * k + 1].split(';')[0]
            if ':' in name_desc:
                name, desc = name_desc.split(':', 1)
            else:
                name, desc = name_desc, ''
            out['SCALARS'].append({
                'name': name.strip(),
                'desc': desc.strip(),
                'data': value,
            })
            end_scalars = 3 + 2 * k + 2
        except Exception:
            break

    # Determine number of independent dimensions by looking ahead until we
    # find a line whose first ';' field is a plain integer (the first dim
    # length).
    ndim = -1
    for line in lines[end_scalars:]:
        try:
            int(line.split(';')[0].strip())
            break
        except Exception:
            ndim += 1

    if ndim < 0:
        raise ValueError(f"Could not determine NDIM in {os.path.basename(path)}")

    def _setvar(line):
        """Parse a variable label line.

        Format (OMFIT writer): '{:<20}{:<10};-LABEL TEXT-'.
        We split on the trailing ';' first, then on multi-space gaps so
        small column-padding mismatches don't merge name and units.
        """
        head = line.split(';', 1)[0].rstrip()
        # split on runs of 2+ spaces — name and units are typically
        # separated by enough padding even when name is short
        parts = re.split(r'\s{2,}', head.strip())
        if len(parts) >= 2:
            return {'name': parts[0].strip(), 'units': parts[1].strip()}
        # Fallback: OMFIT writes ' ' + name(20 chars) + units(10 chars) + ';-..'
        # so the data starts at column 1 (one leading space). Column-split on
        # the lstripped data part to avoid eating units' first char.
        body = line.split(';', 1)[0].lstrip()
        return {
            'name': body[:20].strip(),
            'units': body[20:30].strip(),
        }

    dims = []
    for i in range(ndim):
        out['X{}'.format(i)] = _setvar(lines[end_scalars + i])
        try:
            dims.append(int(lines[end_scalars + ndim + 2 + i].split(';')[0]))
        except Exception:
            dims.append(0)
    out['F'] = _setvar(lines[end_scalars + ndim])

    try:
        out['PROCESSING_CODE'] = int(lines[end_scalars + ndim + 1].split(';')[0])
    except Exception:
        out['PROCESSING_CODE'] = 0

    end_header = end_scalars + ndim * 2 + 2

    # Find trailer
    start_trailer = len(lines)
    for i, line in enumerate(lines[end_header:]):
        if 'END-OF-DATA' in line:
            start_trailer = end_header + i
            break

    # Read numeric data. TRANSP scruncher sometimes leaves no space between
    # negatives in fixed-format output, so insert spaces before stray '-'.
    flat = ' '.join(lines[end_header:start_trailer])
    flat = flat.replace('E-', 'E~').replace('e-', 'e~')
    flat = flat.replace('-', ' -')
    flat = flat.replace('E~', 'E-').replace('e~', 'e-')
    arr = np.fromstring(flat, sep=' ', dtype=float)

    cursor = 0
    for i in range(ndim):
        n = dims[i] if dims[i] > 0 else 0
        out['X{}'.format(i)]['data'] = arr[cursor:cursor + n]
        cursor += n

    nF = int(np.prod(dims)) if dims else 0
    if nF > 0 and len(arr) >= cursor + nF:
        F = arr[cursor:cursor + nF]
        # Fortran ordering: reshape with reversed dims, then transpose.
        try:
            F = F.reshape(dims[::-1]).T
        except Exception:
            F = F.reshape(dims)
        out['F']['data'] = F
    else:
        out['F']['data'] = np.array([])

    # Comments and beams trailer
    out['COMMENTS'] = []
    if start_trailer + 1 < len(lines):
        for i, line in enumerate(lines[start_trailer + 1:]):
            if 'Beams used:' in line:
                segs = (line + (lines[start_trailer + 2 + i]
                                if start_trailer + 2 + i < len(lines) else ''))
                tokens = [t.strip(', \n') for t in segs.split(' ')[3:]]
                tokens = [t for t in tokens if t]
                out['BEAMS'] = tokens
            stripped = line.rstrip('\n')
            if stripped.strip():
                out['COMMENTS'].append(stripped)

    base = os.path.basename(path)
    pm = re.match(r'^([A-Za-z]+)(\d+)\.([A-Za-z0-9]+)$', base)
    out['SOURCE'] = {
        'path': path,
        'prefix': pm.group(1) if pm else '',
        'shot': pm.group(2) if pm else str(out.get('SHOT', '')),
        'ext': pm.group(3) if pm else '',
    }
    return out


# ----- Classification ---------------------------------------------------------

# Names that indicate the axis is a radial/spatial profile coordinate.
# Compared after upper() / split-on-spaces; first token is checked too.
_SPATIAL_NAMES = {
    'R', 'RMAJ', 'RMAJOR', 'RMAJM', 'RMIN', 'RMINOR',
    'RHO', 'RHOPOL', 'RHOTOR', 'RHO_POL', 'RHO_TOR',
    'PSI', 'PSIN', 'PSI_N', 'PSI_NORM',
    'X', 'XB', 'XI', 'A', 'R/A', 'RA',
    'ZONE', 'ZONEBND', 'ZONEBNDY', 'ZONE_BOUNDARY',
    'CHORD', 'CHORDS', 'LOS',
}

# Names that indicate the axis is a discrete channel/source/beam index.
_CHANNEL_NAMES = {
    'CHANNEL', 'CHANNELS', 'CH', 'INDEX', 'IDX',
    'BEAM', 'BEAMS', 'BEAM_INDEX', 'BEAM_NUMBER', 'BEAM_NO', 'BEAMNO',
    'SOURCE', 'SOURCES', 'SOURCE_NO', 'SOURCENO',
    'SYSTEM', 'SYSTEMS', 'SYSTEM_NO', 'SYSTEMNO',
    'NUMBER', 'NO', 'ID', 'GROUP',
}

_LENGTH_UNITS = {'M', 'CM', 'MM', 'METER', 'METERS',
                 'CENTIMETER', 'CENTIMETERS'}
_TIME_UNITS = {'S', 'SEC', 'SECS', 'SECOND', 'SECONDS', 'MS', 'MSEC'}


def _norm_name(axis_dict):
    return (axis_dict.get('name') or '').strip().upper() if axis_dict else ''


def _norm_units(axis_dict):
    return (axis_dict.get('units') or '').strip().upper() if axis_dict else ''


def _data_looks_like_time(data):
    """Heuristic: monotonic ascending values consistent with seconds."""
    if data is None or not hasattr(data, '__len__') or len(data) < 3:
        return False
    try:
        d = np.asarray(data, dtype=float).ravel()
    except Exception:
        return False
    if d.size < 3 or not np.all(np.isfinite(d)):
        return False
    diffs = np.diff(d)
    if not np.all(diffs > 0):
        return False
    # Plasma-shot timeline — wide range to also accept ms-scale axes
    return 0.0 <= float(d[0]) and float(d[-1]) < 1e6


def _is_time_axis_by_name(axis_dict):
    """Strong time signal: name contains TIME or units are time units."""
    if not axis_dict:
        return False
    name = _norm_name(axis_dict)
    units = _norm_units(axis_dict)
    if 'TIME' in name or name == 'T':
        return True
    if units in _TIME_UNITS:
        return True
    return False


def _is_time_axis(axis_dict):
    """Full time-axis test (name/units OR data-pattern fallback)."""
    if _is_time_axis_by_name(axis_dict):
        return True
    # Data fallback: only if name doesn't look spatial.
    if not _is_spatial_axis(axis_dict) and _data_looks_like_time(
            axis_dict.get('data')):
        return True
    return False


def _is_channel_axis(axis_dict):
    """Heuristic: discrete-channel axis (NBI beam, ECH source, etc.)."""
    if not axis_dict:
        return False
    name = _norm_name(axis_dict)
    if not name:
        return False
    if name in _CHANNEL_NAMES:
        return True
    # Match first token: "BEAM NUMBER", "SYSTEM NO", "SOURCE INDEX", ...
    head = name.split()[0] if name.split() else ''
    if head in _CHANNEL_NAMES:
        return True
    return False


def _is_spatial_axis(axis_dict):
    """Heuristic: spatial/radial-profile axis."""
    if not axis_dict:
        return False
    name = _norm_name(axis_dict)
    units = _norm_units(axis_dict)
    if name in _SPATIAL_NAMES:
        return True
    head = name.split()[0] if name.split() else ''
    if head in _SPATIAL_NAMES:
        return True
    if units in _LENGTH_UNITS:
        return True
    return False


def classify_ufile(uf):
    """Classify a U-File. Returns one of:
        'time_trace', 'time_trace_multi', 'profile_time', 'geometry'.

    Side effect: sets `uf['_time_axis']` to the X-index (0/1) of the time
    axis when present (None for geometry).

    Rules (2D):
        TIME on X0 or X1 + other axis spatial  -> profile_time
        TIME on X0 or X1 + other axis channel  -> time_trace_multi
        TIME present + ambiguous other axis:
            other length <= 16 -> time_trace_multi
            other length >  16 -> profile_time
        no TIME axis -> geometry
    """
    ndim = uf.get('NDIM', 0)
    F = uf.get('F', {}).get('data')
    if F is None or F.size == 0:
        uf['_time_axis'] = None
        return 'geometry'
    sq = np.squeeze(F).shape if F.ndim > 0 else ()
    eff_ndim = len(sq) if sq else ndim

    x0 = uf.get('X0')
    x1 = uf.get('X1')

    if eff_ndim <= 1:
        if _is_time_axis(x0):
            uf['_time_axis'] = 0
            return 'time_trace'
        uf['_time_axis'] = None
        return 'geometry'

    if eff_ndim == 2:
        # Strong (name/units-based) signal wins first; only fall back to
        # the data-pattern heuristic if neither axis has a time-like name.
        n0 = _is_time_axis_by_name(x0)
        n1 = _is_time_axis_by_name(x1)
        if n0 or n1:
            if n1 and not n0:
                time_idx, other_idx, other_ax = 1, 0, x0
            else:
                time_idx, other_idx, other_ax = 0, 1, x1
        else:
            t0 = _is_time_axis(x0)
            t1 = _is_time_axis(x1)
            if t1 and not t0:
                time_idx, other_idx, other_ax = 1, 0, x0
            elif t0 and not t1:
                time_idx, other_idx, other_ax = 0, 1, x1
            elif t0 and t1:
                # Both look time-like via data — prefer X0
                time_idx, other_idx, other_ax = 0, 1, x1
            else:
                uf['_time_axis'] = None
                return 'geometry'

        uf['_time_axis'] = time_idx

        if 'BEAMS' in uf or _is_channel_axis(other_ax):
            return 'time_trace_multi'
        if _is_spatial_axis(other_ax):
            return 'profile_time'

        n_other = sq[other_idx] if len(sq) > other_idx else 0
        if 0 < n_other <= 16:
            return 'time_trace_multi'
        return 'profile_time'

    uf['_time_axis'] = None
    return 'geometry'  # 3D and beyond → bucket as geometry
