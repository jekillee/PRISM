"""
EFIT Viewer data loader

Loads EFIT equilibrium data from three sources:
  1. MDS+ trees (efit01, efitrt1, efit02, efitrt2, efit04)
  2. g-files (GEQDSK format - profiles + 2D psi(R,Z))
  3. a-files (AEQDSK format - scalar time traces)

MDS+ node structure (KSTAR):
  AEQDSK scalars : \\{TREE}::TOP.RESULTS.AEQDSK:{NODE}  shape=(N_time,)
  GEQDSK profiles: \\{TREE}::TOP.RESULTS.GEQDSK:{NODE}  shape=(N_time, 129)
  GEQDSK 2D      : PSIRZ shape=(N_time, 129, 129)
  Time           : ATIME for scalars, GTIME for profiles/2D
"""

import os
from typing import Dict, Any, Optional, List

import numpy as np


# ---------------------------------------------------------------------------
# MDS+ node definitions
# ---------------------------------------------------------------------------

AEQDSK_SCALAR_NODES = [
    'BETAN', 'BETAP', 'BETAT', 'LI', 'Q0', 'Q95', 'QMIN',
    'AMINOR', 'KAPPA', 'TRITOP', 'TRIBOT', 'DRSEP',
    'IPMHD', 'WMHD', 'VOLUME', 'AREA', 'AREAO',
    'BPOLAV', 'BT0', 'ROUT', 'AOUT', 'EOUT',
    'WBDOT', 'QMERCI', 'CHISQ',
]

GEQDSK_PROFILE_NODES = [
    'PRES', 'PPRIME', 'FPOL', 'FFPRIM', 'QPSI', 'RHOVN',
]

GEQDSK_2D_NODES = [
    'PSIRZ', 'RGRID', 'Z',
    'RBBBS', 'ZBBBS', 'RMAXIS', 'ZMAXIS',
    'SSIMAG', 'SSIBRY', 'PSIN', 'BCENTR', 'RZERO',
]

# Human-readable labels for AEQDSK scalars
AEQDSK_LABELS = {
    'BETAN':  r'$\beta_N$',
    'BETAP':  r'$\beta_p$',
    'BETAT':  r'$\beta_T$',
    'LI':     r'$l_i$',
    'Q0':     r'$q_0$',
    'Q95':    r'$q_{95}$',
    'QMIN':   r'$q_{min}$',
    'QSTAR':  r'$q^*$ (diamag)',
    'QOUT':   r'$q_{edge}$',
    'IPMEAS': r'$I_p$ (meas) [A]',
    'RCEN':   r'$R_{cen}$ [cm]',
    'RMAXIS': r'$R_{axis}$ [cm]',
    'ZMAXIS': r'$Z_{axis}$ [cm]',
    'RCUR':   r'$R_{cur}$ [cm]',
    'ZCUR':   r'$Z_{cur}$ [cm]',
    'ROUT':   r'$R_{out}$ [cm]',
    'ZOUT':   r'$Z_{out}$ [cm]',
    'ELONGM': r'$\kappa_{mid}$',
    'BETAPD': r'$\beta_p$ (diamag)',
    'BETATD': r'$\beta_T$ (diamag)',
    'WDIA':   r'$W_{dia}$ [J]',
    'SHEARB': r'Shear$_b$',
    'GAPIN':  r'Gap$_{in}$ [cm]',
    'GAPOUT': r'Gap$_{out}$ [cm]',
    'GAPTOP': r'Gap$_{top}$ [cm]',
    'GAPBOT': r'Gap$_{bot}$ [cm]',
    'VERTN':  r'Vert. stab. index',
    'BTM':    r'$B_{T,mid}$ [T]',
    'BTVAC':  r'$B_{T,vac}$ [T]',
    'SIMAGX': r'$\psi_{axis}$',
    'SIBDRY': r'$\psi_{bdy}$',
    'AMINOR': r'$a$ [m]',
    'KAPPA':  r'$\kappa$',
    'TRITOP': r'$\delta_{top}$',
    'TRIBOT': r'$\delta_{bot}$',
    'DRSEP':  r'$\Delta R_{sep}$ [m]',
    'IPMHD':  r'$I_p$ (MHD) [A]',
    'WMHD':   r'$W_{MHD}$ [J]',
    'VOLUME': r'Volume [$m^3$]',
    'AREA':   r'Area [$m^2$]',
    'AREAO':  r'Area$_O$ [$m^2$]',
    'BPOLAV': r'$\langle B_p \rangle$ [T]',
    'BT0':    r'$B_{T0}$ [T]',
    'ROUT':   r'$R_{out}$ [m]',
    'AOUT':   r'$a_{out}$ [m]',
    'EOUT':   r'$\epsilon_{out}$',
    'WBDOT':  r'$dW/dt$ [W]',
    'QMERCI': r'$q_{Mercier}$',
    'CHISQ':  r'$\chi^2$',
}

# Human-readable labels for GEQDSK profiles
GEQDSK_PROFILE_LABELS = {
    'PRES':   r'Pressure [Pa]',
    'PPRIME': r"$P'$ [Pa/Wb]",
    'FPOL':   r'$F_{pol}$ [T$\cdot$m]',
    'FFPRIM': r"$FF'$ [T$^2$$\cdot$m$^2$/Wb]",
    'QPSI':   r'$q$',
    'RHOVN':  r'$\rho_N$ (vol)',
}

# Available EFIT trees on KSTAR
EFIT_TREES = ['efit01', 'efitrt1', 'efit02', 'efitrt2', 'efit04']


# ---------------------------------------------------------------------------
# Helper: auto-detect time units and convert to seconds
# ---------------------------------------------------------------------------

def _convert_time_to_seconds(time_arr: np.ndarray) -> np.ndarray:
    """Auto-detect whether time is in ms or s, return array in seconds.

    Heuristic: if max(|time|) > 100 the data is assumed to be in ms.
    """
    if time_arr is None or len(time_arr) == 0:
        return time_arr
    if np.max(np.abs(time_arr)) > 100.0:
        return time_arr / 1e3
    return time_arr.copy()


# ===================================================================
# 1. MDS+ loader
# ===================================================================

def load_efit_mds(shot: int, tree: str = 'efit01',
                  mds_ip: str = 'mdsr.kstar.kfe.re.kr:8005',
                  load_scalars: bool = True,
                  load_profiles: bool = True,
                  load_2d: bool = True) -> Dict[str, Any]:
    """Load EFIT data from an MDS+ tree.

    Args:
        shot: KSTAR shot number.
        tree: EFIT tree name (e.g. 'efit01', 'efitrt1').
        mds_ip: MDS+ server address.
        load_scalars: Load AEQDSK scalar time traces.
        load_profiles: Load GEQDSK 1D profiles.
        load_2d: Load GEQDSK 2D (PSIRZ, boundary).

    Returns:
        Dictionary with keys:
            'source'  : str  (e.g. '#40848 efit01')
            'shot'    : int
            'tree'    : str
            'scalars' : {node: {'time': ndarray[s], 'data': ndarray, 'label': str}}
            'profiles': {node: {'time': ndarray[s], 'psin': ndarray, 'data': ndarray(time,psin), 'label': str}}
            'eq2d'    : {'time': ndarray[s], 'rgrid': ndarray, 'zgrid': ndarray, ...}
    """
    from MDSplus import Connection

    print(f"[EFIT] Loading #{shot} from MDS+ tree '{tree}'...")

    mds = Connection(mds_ip)
    mds.openTree(tree, shot)

    result = {
        'source': f'#{shot} {tree}',
        'shot': shot,
        'tree': tree,
        'scalars': {},
        'profiles': {},
        'eq2d': {},
    }

    # ---- AEQDSK scalar time base ----
    atime = None
    if load_scalars:
        atime = _mds_get_safe(mds, tree, 'ATIME', section='AEQDSK')
        if atime is not None:
            atime = _convert_time_to_seconds(atime)

    # ---- GEQDSK time base ----
    gtime = None
    if load_profiles or load_2d:
        gtime = _mds_get_safe(mds, tree, 'GTIME', section='GEQDSK')
        if gtime is not None:
            gtime = _convert_time_to_seconds(gtime)

    # ---- Load AEQDSK scalars ----
    if not load_scalars:
        print("[EFIT]   Skipping AEQDSK scalars")
    for node in (AEQDSK_SCALAR_NODES if load_scalars else []):
        data = _mds_get_safe(mds, tree, node, section='AEQDSK')
        if data is not None and atime is not None:
            result['scalars'][node] = {
                'time': atime,
                'data': np.asarray(data, dtype=np.float64),
                'label': AEQDSK_LABELS.get(node, node),
            }

    if result['scalars']:
        print(f"[EFIT]   Loaded {len(result['scalars'])} AEQDSK scalars "
              f"({len(atime)} time points)")
    else:
        print("[EFIT]   Warning: No AEQDSK scalars loaded")

    # ---- GEQDSK psi_n grid ----
    psin = None
    if not load_profiles:
        print("[EFIT]   Skipping GEQDSK profiles")
    if load_profiles:
        psin = _mds_get_safe(mds, tree, 'PSIN', section='GEQDSK')
    if psin is not None:
        # PSIN may be 1D (129,) or 2D (N_time, 129). Use first slice if 2D.
        if psin.ndim == 2:
            psin = psin[0]
        psin = np.asarray(psin, dtype=np.float64)

    # ---- Load GEQDSK profiles ----
    for node in (GEQDSK_PROFILE_NODES if load_profiles else []):
        data = _mds_get_safe(mds, tree, node, section='GEQDSK')
        if data is not None and gtime is not None:
            data = np.asarray(data, dtype=np.float64)
            # Ensure 2D: (N_time, N_psin)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            psin_axis = psin if psin is not None else np.linspace(0, 1, data.shape[1])
            result['profiles'][node] = {
                'time': gtime,
                'psin': psin_axis,
                'data': data,
                'label': GEQDSK_PROFILE_LABELS.get(node, node),
            }

    if result['profiles']:
        nprof = list(result['profiles'].values())[0]['data'].shape[1]
        print(f"[EFIT]   Loaded {len(result['profiles'])} GEQDSK profiles "
              f"({len(gtime)} time points, {nprof} psin grid)")
    else:
        print("[EFIT]   Warning: No GEQDSK profiles loaded")

    # ---- Load GEQDSK 2D / boundary data ----
    eq2d = {}
    if not load_2d:
        print("[EFIT]   Skipping GEQDSK 2D/boundary")
        result['eq2d'] = eq2d
        mds.closeTree(tree, shot)
        print(f"[EFIT] Done loading #{shot} {tree}")
        return result

    eq2d['time'] = gtime

    rgrid = _mds_get_safe(mds, tree, 'RGRID', section='GEQDSK')
    zgrid = _mds_get_safe(mds, tree, 'Z', section='GEQDSK')
    psirz = _mds_get_safe(mds, tree, 'PSIRZ', section='GEQDSK')

    if rgrid is not None:
        eq2d['rgrid'] = np.asarray(rgrid, dtype=np.float64)
    if zgrid is not None:
        eq2d['zgrid'] = np.asarray(zgrid, dtype=np.float64)
    if psirz is not None:
        eq2d['psirz'] = np.asarray(psirz, dtype=np.float64)

    # Boundary
    rbbbs = _mds_get_safe(mds, tree, 'RBBBS', section='GEQDSK')
    zbbbs = _mds_get_safe(mds, tree, 'ZBBBS', section='GEQDSK')
    if rbbbs is not None:
        eq2d['rbbbs'] = np.asarray(rbbbs, dtype=np.float64)
    if zbbbs is not None:
        eq2d['zbbbs'] = np.asarray(zbbbs, dtype=np.float64)

    # Limiter
    xlim = _mds_get_safe(mds, tree, 'XLIM', section='GEQDSK')
    ylim = _mds_get_safe(mds, tree, 'YLIM', section='GEQDSK')
    if xlim is not None and ylim is not None:
        eq2d['rlim'] = np.asarray(xlim, dtype=np.float64)
        eq2d['zlim'] = np.asarray(ylim, dtype=np.float64)

    # Magnetic axis
    for node in ['RMAXIS', 'ZMAXIS', 'SSIMAG', 'SSIBRY', 'BCENTR', 'RZERO']:
        val = _mds_get_safe(mds, tree, node, section='GEQDSK')
        if val is not None:
            eq2d[node.lower()] = np.asarray(val, dtype=np.float64)

    # X-points (from AEQDSK, time-dependent)
    for node in ['RXPT1', 'ZXPT1', 'RXPT2', 'ZXPT2']:
        val = _mds_get_safe(mds, tree, node, section='AEQDSK')
        if val is not None:
            eq2d[node.lower()] = np.asarray(val, dtype=np.float64)

    # RHOVN for rho_tor contour option (from profiles, but store in eq2d too)
    rhovn = _mds_get_safe(mds, tree, 'RHOVN', section='GEQDSK')
    if rhovn is not None:
        eq2d['rhovn'] = np.asarray(rhovn, dtype=np.float64)

    if psin is not None:
        eq2d['psin_grid'] = psin

    if 'psirz' in eq2d:
        print(f"[EFIT]   Loaded 2D equilibrium: PSIRZ shape={eq2d['psirz'].shape}")
    else:
        print("[EFIT]   Warning: PSIRZ not loaded")

    result['eq2d'] = eq2d

    mds.closeTree(tree, shot)
    print(f"[EFIT] Done loading #{shot} {tree}")

    return result


def _mds_get_safe(mds, tree: str, node: str,
                  section: str = 'GEQDSK') -> Optional[np.ndarray]:
    """Safely fetch a single MDS+ node, returning None on failure.

    Args:
        mds: MDSplus Connection object (tree already open).
        tree: EFIT tree name (for TDI path construction).
        node: Node name (e.g. 'PSIRZ', 'BETAN').
        section: 'AEQDSK' or 'GEQDSK'.

    Returns:
        numpy array or None if the node does not exist / fails.
    """
    tdi = f'\\{tree}::TOP.RESULTS.{section}:{node}'
    try:
        data = mds.get(tdi).data()
        return np.asarray(data)
    except Exception as e:
        print(f"[EFIT]   Warning: Could not read {tdi} - {e}")
        return None


# ===================================================================
# 2. g-file (GEQDSK) parser
# ===================================================================

def load_gfile(filepath: str) -> Dict[str, Any]:
    """Parse a standard GEQDSK (g-file) and return an EFIT data dict.

    The returned dictionary has the same structure as load_efit_mds() but
    contains a single time slice.  The 'scalars' dict is empty for g-files;
    scalar-like quantities (rmaxis, bcentr, etc.) are stored in 'eq2d'.

    Standard GEQDSK format:
        Line 1: case string(s), idum, nw, nh
        Then 5-per-line floats: header scalars + 1D/2D arrays
        Then nbbbs, limitr  followed by boundary/limiter points

    Args:
        filepath: Path to the GEQDSK g-file.

    Returns:
        EFIT data dictionary (single time slice).

    Raises:
        RuntimeError: If parsing fails critically.
    """
    print(f"[EFIT] Loading g-file: {os.path.basename(filepath)}")

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        raise RuntimeError(f"Failed to open g-file: {e}")

    # ---- Parse header line ----
    header = lines[0]
    # Last 3 integers on header line: idum, nw, nh
    tokens = header.split()
    try:
        nh = int(tokens[-1])
        nw = int(tokens[-2])
    except (ValueError, IndexError):
        raise RuntimeError(f"Cannot parse nw, nh from g-file header: {header.strip()}")

    case_string = header.strip()
    print(f"[EFIT]   Header: {case_string}")
    print(f"[EFIT]   Grid: nw={nw}, nh={nh}")

    # ---- Read all remaining numbers in 5-per-line format ----
    all_values = []  # type: List[float]
    for line in lines[1:]:
        # Try to parse floats; skip lines that fail entirely
        parts = _split_fortran_floats(line)
        all_values.extend(parts)

    idx = 0  # reading cursor into all_values

    def _read(n: int = 1) -> np.ndarray:
        nonlocal idx
        arr = np.array(all_values[idx:idx + n], dtype=np.float64)
        idx += n
        return arr

    def _read1() -> float:
        nonlocal idx
        val = all_values[idx]
        idx += 1
        return val

    # ---- Header scalars (20 values) ----
    rdim   = _read1()
    zdim   = _read1()
    rcentr = _read1()
    rleft  = _read1()
    zmid   = _read1()
    rmaxis = _read1()
    zmaxis = _read1()
    simag  = _read1()
    sibry  = _read1()
    bcentr = _read1()
    cpasma = _read1()
    _read1()  # simag (duplicate)
    _read1()  # xdum
    _read1()  # rmaxis (duplicate)
    _read1()  # xdum
    _read1()  # zmaxis (duplicate)
    _read1()  # xdum
    _read1()  # sibry (duplicate)
    _read1()  # xdum
    _read1()  # xdum

    # ---- 1D arrays ----
    fpol   = _read(nw)
    pres   = _read(nw)
    ffprim = _read(nw)
    pprime = _read(nw)

    # ---- 2D psi(R,Z) stored row-major: psirz[nh][nw] ----
    psirz_flat = _read(nw * nh)
    psirz = psirz_flat.reshape((nh, nw))

    qpsi = _read(nw)

    # ---- Boundary and limiter counts ----
    nbbbs  = int(_read1())
    limitr = int(_read1())

    # ---- Boundary points ----
    rbbbs = np.zeros(nbbbs)
    zbbbs = np.zeros(nbbbs)
    if nbbbs > 0:
        bdry = _read(2 * nbbbs)
        rbbbs = bdry[0::2]
        zbbbs = bdry[1::2]

    # ---- Limiter points ----
    rlim = np.zeros(limitr)
    zlim = np.zeros(limitr)
    if limitr > 0:
        lim = _read(2 * limitr)
        rlim = lim[0::2]
        zlim = lim[1::2]

    # ---- Build grids ----
    rgrid = np.linspace(rleft, rleft + rdim, nw)
    zgrid = np.linspace(zmid - zdim / 2.0, zmid + zdim / 2.0, nh)
    psin  = np.linspace(0.0, 1.0, nw)

    # ---- Construct result ----
    # Try to extract shot number, time, suffix from filename
    source = os.path.basename(filepath)
    shot_from_file, time_from_file, suffix = _parse_gfile_name(filepath)

    time_arr = np.array([time_from_file]) if time_from_file is not None else np.array([0.0])

    result = {
        'source': source,
        'shot': shot_from_file,
        'tree': suffix,  # e.g. 'kin_0' from g040848.02000_kin_0
        'scalars': {},
        'profiles': {},
        'eq2d': {},
    }

    # Profiles (single time slice, add leading time dimension)
    for name, data, label in [
        ('PRES',   pres,   GEQDSK_PROFILE_LABELS.get('PRES', 'PRES')),
        ('PPRIME', pprime, GEQDSK_PROFILE_LABELS.get('PPRIME', 'PPRIME')),
        ('FPOL',   fpol,   GEQDSK_PROFILE_LABELS.get('FPOL', 'FPOL')),
        ('FFPRIM', ffprim, GEQDSK_PROFILE_LABELS.get('FFPRIM', 'FFPRIM')),
        ('QPSI',   qpsi,   GEQDSK_PROFILE_LABELS.get('QPSI', 'QPSI')),
    ]:
        result['profiles'][name] = {
            'time': time_arr,
            'psin': psin,
            'data': data.reshape(1, -1),
            'label': label,
        }

    # 2D equilibrium
    result['eq2d'] = {
        'time': time_arr,
        'rgrid': rgrid,
        'zgrid': zgrid,
        'psirz': psirz.reshape(1, nh, nw),
        'rbbbs': rbbbs.reshape(1, -1),
        'zbbbs': zbbbs.reshape(1, -1),
        'rmaxis': np.array([rmaxis]),
        'zmaxis': np.array([zmaxis]),
        'ssimag': np.array([simag]),
        'ssibry': np.array([sibry]),
        'psin_grid': psin,
        'bcentr': np.array([bcentr]),
        'rzero': np.array([rcentr]),
        'rlim': rlim,
        'zlim': zlim,
        'cpasma': cpasma,
    }

    print(f"[EFIT]   Parsed g-file: nw={nw}, nh={nh}, nbbbs={nbbbs}, "
          f"Rmag={rmaxis:.3f}, Zmag={zmaxis:.3f}")
    print(f"[EFIT] Done loading g-file")

    return result


def _split_fortran_floats(line: str) -> List[float]:
    """Parse a line of Fortran-formatted floats.

    Handles cases where negative signs act as delimiters (no space between
    numbers), e.g. " 1.234567E+00-2.345678E-01".
    """
    values = []
    # First try simple split
    tokens = line.split()
    for tok in tokens:
        try:
            values.append(float(tok))
        except ValueError:
            # Handle concatenated Fortran floats like "1.23E+00-4.56E-01"
            sub_vals = _split_concatenated_fortran(tok)
            values.extend(sub_vals)
    return values


def _split_concatenated_fortran(tok: str) -> List[float]:
    """Split a string that may contain concatenated Fortran floats.

    Example: "1.234E+00-2.345E-01" -> [1.234e+00, -2.345e-01]
    """
    import re
    # Pattern: optional sign, digits, dot, digits, optional exponent
    pattern = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[eEdD][+-]?\d+)?'
    matches = re.findall(pattern, tok)
    result = []
    for m in matches:
        try:
            result.append(float(m.replace('D', 'E').replace('d', 'e')))
        except ValueError:
            pass
    return result


def _parse_gfile_name(filepath: str):
    """Try to extract shot number, time, and suffix from g-file name.

    Common naming: g{SHOT}.{TIME_ms}_{SUFFIX}
    e.g. g040848.02000_kin_0 -> shot=40848, t=2.0s, suffix='kin_0'

    Returns:
        (shot: int or None, time_s: float or None, suffix: str or None)
    """
    import re
    basename = os.path.basename(filepath)
    match = re.match(r'[ga](\d+)\.(\d+)(?:_(.+))?', basename)
    if match:
        shot = int(match.group(1))
        time_ms = int(match.group(2))
        time_s = time_ms / 1000.0
        suffix = match.group(3)  # None if no underscore suffix
        return shot, time_s, suffix
    return None, None, None


# ===================================================================
# 3. a-file (AEQDSK) parser
# ===================================================================

def load_afile(filepath: str) -> Dict[str, Any]:
    """Parse a standard AEQDSK (a-file) and return an EFIT data dict.

    The AEQDSK file contains scalar time-trace quantities.  The format is
    simpler than GEQDSK: a header followed by blocks of scalar values at
    each time step.

    Standard AEQDSK format:
        Line 1: header (shot, time info)
        Then blocks of scalar values.  The exact layout varies by EFIT version,
        but the common KSTAR layout has a fixed number of values per time step.

    This parser reads the file and attempts to extract the most commonly
    used scalar quantities.

    Args:
        filepath: Path to the AEQDSK a-file.

    Returns:
        EFIT data dictionary with 'scalars' populated.
    """
    print(f"[EFIT] Loading a-file: {os.path.basename(filepath)}")

    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        raise RuntimeError(f"Failed to open a-file: {e}")

    source = os.path.basename(filepath)
    shot_from_file, _, suffix = _parse_afile_name(filepath)

    # ---- Parse header ----
    # First line typically contains: shot, time, etc.
    header = lines[0].strip()
    print(f"[EFIT]   Header: {header}")

    # ---- Read all numeric values (skip first 4 header lines) ----
    all_values = []  # type: List[float]
    for line in lines[4:]:
        parts = _split_fortran_floats(line)
        all_values.extend(parts)

    # ---- FreeQDSK-standard AEQDSK layout ----
    # Reference: https://github.com/freegs-plasma/FreeQDSK
    #
    # Block 1: 24 general scalars
    # Laser block: rco2v(mco2v) + dco2v(mco2v) + rco2r(mco2r) + dco2r(mco2r)
    # Block 2: 41 more scalars (shearb, bpolav, qout, simagx, sibdry, ...)
    # Extended section: optional arrays + scalars

    # Parse mco2v, mco2r from metadata line (line 3)
    import re
    mco2v, mco2r = 0, 0
    meta_line = lines[3] if len(lines) > 3 else ''
    meta_ints = re.findall(r'\d+', meta_line)
    # Format: * time jflag lflag limloc mco2v mco2r ...
    # The numbers after "SNB" are mco2v, mco2r
    if 'SNB' in meta_line:
        snb_idx = meta_line.index('SNB')
        after_snb = re.findall(r'\d+', meta_line[snb_idx:])
        if len(after_snb) >= 2:
            mco2v = int(after_snb[0])
            mco2r = int(after_snb[1])

    # EFITVIEWER_KFE standard AEQDSK layout (KSTAR)
    # Reference: /vmefittools/source/efit/reada.pro
    #
    # Block 1: 24 general scalars (0-23)
    BLOCK1_MAP = {
        0:  'CHISQ',       # chi-squared
        1:  'RCEN',        # R geometric center [cm]
        2:  'BT0',         # B toroidal at center [T]
        3:  'IPMEAS',      # measured Ip [A]
        4:  'IPMHD',       # MHD Ip [A]
        5:  'ROUT',        # R outer surface [cm]
        6:  'ZOUT',        # Z outer surface [cm]
        7:  'AMINOR',      # minor radius [cm]
        8:  'KAPPA',       # elongation
        9:  'TRITOP',      # upper triangularity
        10: 'TRIBOT',      # lower triangularity
        11: 'VOLUME',      # plasma volume [cm³]
        12: 'RCUR',        # R current centroid [cm]
        13: 'ZCUR',        # Z current centroid [cm]
        14: 'QSTAR',       # q-star (diamagnetic q, not same as Q0)
        15: 'BETAT',       # beta toroidal
        16: 'BETAP',       # beta poloidal
        17: 'LI',          # internal inductance
        18: 'GAPIN',       # gap inboard [cm]
        19: 'GAPOUT',      # gap outboard [cm]
        20: 'GAPTOP',      # gap top [cm]
        21: 'GAPBOT',      # gap bottom [cm]
        22: 'Q95',         # q at 95% flux
        23: 'VERTN',       # vertical stability index
    }

    # Laser block (variable length)
    laser_len = 2 * mco2v + 2 * mco2r

    # Block 2: 41 more scalars (EFITVIEWER positions 33+)
    b2_off = 24 + laser_len
    BLOCK2_MAP = {
        b2_off + 0:  'SHEARB',    # shear at boundary
        b2_off + 1:  'BPOLAV',    # avg poloidal field [T]
        b2_off + 2:  'QOUT',      # q at edge
        b2_off + 3:  'SIMAGX',    # psi at magnetic axis
        b2_off + 4:  'SIBDRY',    # psi at boundary
        b2_off + 10: 'AREA',      # cross-section area [cm²]
        b2_off + 11: 'WMHD',      # MHD stored energy [J]
        b2_off + 13: 'ELONGM',    # elongation at midplane
        b2_off + 14: 'QM',        # q at midplane
        b2_off + 20: 'RXPT1',     # R lower X-point [cm]
        b2_off + 21: 'ZXPT1',     # Z lower X-point [cm]
        b2_off + 22: 'RXPT2',     # R upper X-point [cm]
        b2_off + 23: 'ZXPT2',     # Z upper X-point [cm]
        b2_off + 26: 'BTM',       # B at midplane [T]
        b2_off + 27: 'BTVAC',     # B vacuum [T]
        b2_off + 32: 'RMAXIS',    # R magnetic axis [cm]
        b2_off + 33: 'ZMAXIS',    # Z magnetic axis [cm]
        b2_off + 36: 'BETAPD',    # betap (diamagnetic)
        b2_off + 37: 'BETATD',    # betat (diamagnetic)
        b2_off + 38: 'WDIA',      # diamagnetic energy [J]
    }

    # BETAN is NOT directly stored — calculated per EFITVIEWER:
    # BETAN = BETAT / IN, where IN = Ip[MA] / (a[m] * BT[T])

    # Extended section: positions after block 2
    ext_off = b2_off + 41
    # 4 size integers, then arrays, then more scalars
    # Positions 98+ in EFITVIEWER:
    EXT_MAP = {}
    # QMIN at EFITVIEWER pos 98 → ext_off + 4 + nsilop + magpri + nfcoil + nesum + ...
    # Too complex to map without knowing exact sizes. Skip for now.

    AFILE_POS_MAP = {**BLOCK1_MAP, **BLOCK2_MAP}

    n_total = len(all_values)

    result = {
        'source': source,
        'shot': shot_from_file,
        'tree': suffix,
        'scalars': {},
        'profiles': {},
        'eq2d': {},
    }

    # Use time from filename
    _, time_from_name, _ = _parse_afile_name(filepath)
    time_arr = np.array([time_from_name if time_from_name is not None else 0.0])

    for pos, node in AFILE_POS_MAP.items():
        if pos < n_total:
            result['scalars'][node] = {
                'time': time_arr,
                'data': np.array([all_values[pos]]),
                'label': AEQDSK_LABELS.get(node, node),
            }

    # Compute BETAN = BETAT / (Ip[MA] / (a[m] * BT[T]))
    if all(p < n_total for p in [2, 4, 7, 15]):
        betat = all_values[15]
        ip_a = all_values[4]
        aminor_cm = all_values[7]
        bt_t = all_values[2]
        if abs(ip_a) > 0 and abs(aminor_cm) > 0 and abs(bt_t) > 0:
            i_n = abs(ip_a) * 1e-6 / (aminor_cm * 0.01 * abs(bt_t))
            betan = betat / i_n if i_n > 0 else 0
            result['scalars']['BETAN'] = {
                'time': time_arr,
                'data': np.array([betan]),
                'label': AEQDSK_LABELS.get('BETAN', r'$\beta_N$'),
            }

    # Also add AOUT = AMINOR (same value per KSTAR convention)
    if 'AMINOR' in result['scalars']:
        result['scalars']['AOUT'] = {
            'time': time_arr,
            'data': result['scalars']['AMINOR']['data'].copy(),
            'label': AEQDSK_LABELS.get('AOUT', r'$a_{out}$'),
        }

    n_loaded = len(result['scalars'])
    print(f"[EFIT]   Loaded {n_loaded} scalars from a-file")
    print(f"[EFIT] Done loading a-file")

    return result


def _parse_afile_name(filepath: str):
    """Try to extract shot number, time, suffix from a-file name.

    Common naming: a{SHOT}.{TIME_ms}_{SUFFIX} e.g. a040848.02000_kin_0

    Returns:
        (shot: int or None, time_s: float or None, suffix: str or None)
    """
    return _parse_gfile_name(filepath)  # same format as g-file
