"""
BIPROFILE MDS+ data loader

MDS+ Tree: 'biprofile'
Prefixes: BITI(Ti), BIVT(vT), BITE(Te), BINE(ne)

Structure:
    \\PREFIX:PSIN, :TIME, :EFIT_USED, :FIT_FUNC
    \\PREFIX.PREFIX###:MEAN, :UNC, :LFSR
    \\PREFIX.PED:WIDTH, :HEIGHT

DIAG_PARAMS (in 'biprofile' tree):
    \\BITS.BITS_CORE/EDGE{XX}: SCALE, NE_USE, TE_USE, NE_OUTLIER_P, TE_OUTLIER_P
    \\BICES.BICES_TI/VT{XX}: TI_USE, TI_OUTLIER_P, VT_USE, VT_OUTLIER_P
    \\BINE_TCI.BINE_TCI{XX}: NE_USE, NE_OUTLIER_P
"""

import json
from pathlib import Path

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader

PARAM_TREES = {'Ti': 'BITI', 'vT': 'BIVT', 'Te': 'BITE', 'ne': 'BINE'}

# Thomson channel R positions: read from the shared config so biprofile_loader
# and thomson_loader stay in lockstep (single source of truth).
_TS_CONFIG_PATH = Path(__file__).parent.parent / 'config' / 'thomson_positions.json'
with open(_TS_CONFIG_PATH, 'r') as _f:
    _THOMSON_CONFIG = json.load(_f)


def _try_get(mds, path, default=None):
    try:
        return mds.get(path).data()
    except Exception:
        return default


def _get_thomson_positions(shot):
    """Return (core_m, edge_m) channel R positions for the given shot.
    Reads from config/thomson_positions.json — same data thomson_loader uses."""
    for entry in _THOMSON_CONFIG['positions']:
        lo, hi = entry['shot_range']
        if hi is None:
            hi = float('inf')
        if lo <= shot <= hi:
            return ([r / 1000.0 for r in entry['core']],
                    [r / 1000.0 for r in entry['edge'] if r is not None])
    last = _THOMSON_CONFIG['positions'][-1]
    return ([r / 1000.0 for r in last['core']],
            [r / 1000.0 for r in last['edge'] if r is not None])


class BiProfileLoader(BaseDiagnosticLoader):
    """Loader for BIPROFILE Bayesian inference fitting results"""

    def load_data(self, shot_number, analysis_type=None):
        pass  # Not used; use the module-level functions instead


def load_biprofile(shot_number, param='Ti'):
    """Load BIPROFILE fitting results for one parameter."""
    from MDSplus import Connection
    from config.app_config import AppConfig

    prefix = PARAM_TREES[param]
    mds = Connection(AppConfig().MDS_IP)
    mds.openTree('biprofile', shot_number)

    try:
        psin = mds.get(f'\\{prefix}:PSIN').data()
        time = mds.get(f'\\{prefix}:TIME').data()
        n_psi, n_time = len(psin), len(time)

        efit_used = str(_try_get(mds, f'\\{prefix}:EFIT_USED', 'unknown'))
        fit_func = str(_try_get(mds, f'\\{prefix}:FIT_FUNC', 'unknown'))

        mean = np.full((n_psi, n_time), np.nan)
        unc = np.full((n_psi, n_time), np.nan)

        for i in range(n_psi):
            node = f'\\{prefix}.{prefix}{i+1:03d}'
            d = _try_get(mds, f'{node}:MEAN')
            if d is not None:
                mean[i] = d
            d = _try_get(mds, f'{node}:UNC')
            if d is not None:
                unc[i] = d

        ped_width = _try_get(mds, f'\\{prefix}.PED:WIDTH')
        ped_height = _try_get(mds, f'\\{prefix}.PED:HEIGHT')

        return {
            'psin': psin, 'time': time,
            'mean': mean, 'unc': unc,
            'efit_used': efit_used, 'fit_func': fit_func,
            'ped_width': ped_width, 'ped_height': ped_height,
        }
    finally:
        mds.closeTree('biprofile', shot_number)


def load_diag_params(shot_number):
    """Load DIAG_PARAMS: channel usage and outlier probabilities."""
    from MDSplus import Connection
    from config.app_config import AppConfig

    mds = Connection(AppConfig().MDS_IP)
    mds.openTree('biprofile', shot_number)

    try:
        result = {}

        thomson = {'core': [], 'edge': []}
        for i in range(1, 15):
            base = f'\\BITS.BITS_CORE{i:02d}'
            thomson['core'].append({
                'scale': _try_get(mds, f'{base}:SCALE'),
                'ne_use': _try_get(mds, f'{base}:NE_USE'),
                'te_use': _try_get(mds, f'{base}:TE_USE'),
                'ne_outlier_p': _try_get(mds, f'{base}:NE_OUTLIER_P'),
                'te_outlier_p': _try_get(mds, f'{base}:TE_OUTLIER_P'),
            })
        for i in range(1, 18):
            base = f'\\BITS.BITS_EDGE{i:02d}'
            thomson['edge'].append({
                'scale': _try_get(mds, f'{base}:SCALE'),
                'ne_use': _try_get(mds, f'{base}:NE_USE'),
                'te_use': _try_get(mds, f'{base}:TE_USE'),
                'ne_outlier_p': _try_get(mds, f'{base}:NE_OUTLIER_P'),
                'te_outlier_p': _try_get(mds, f'{base}:TE_OUTLIER_P'),
            })
        result['thomson'] = thomson

        ces = {'ti': [], 'vt': []}
        for i in range(1, 33):
            xx = f'{i:02d}'
            ces['ti'].append({
                'ti_use': _try_get(mds, f'\\BICES.BICES_TI{xx}:TI_USE'),
                'ti_outlier_p': _try_get(mds, f'\\BICES.BICES_TI{xx}:TI_OUTLIER_P'),
            })
            ces['vt'].append({
                'vt_use': _try_get(mds, f'\\BICES.BICES_VT{xx}:VT_USE'),
                'vt_outlier_p': _try_get(mds, f'\\BICES.BICES_VT{xx}:VT_OUTLIER_P'),
            })
        result['ces'] = ces

        tci = {'channels': []}
        for i in range(1, 6):
            xx = f'{i:02d}'
            tci['channels'].append({
                'ne_use': _try_get(mds, f'\\BINE_TCI.BINE_TCI{xx}:NE_USE'),
                'ne_outlier_p': _try_get(mds, f'\\BINE_TCI.BINE_TCI{xx}:NE_OUTLIER_P'),
            })
        result['tci'] = tci

        return result
    finally:
        mds.closeTree('biprofile', shot_number)


def load_efit_psin(shot_number, efit_tree='EFIT01'):
    """Load EFIT psi_N(R) at midplane (z=0) plus 1D equilibrium profiles
    needed by derived quantities (q, fpol, B0, R_axis).

    Backwards compatible: returns the original {'time', 'r_grid', 'psin'} keys
    and adds 'q_psin_grid', 'q_t_psi', 'fpol_t_psi', 'bcentr', 'rmaxis'.

    1D EFIT profiles (qpsi, fpol) are stored vs an internal ψ_N grid built from
    \\mw (linspace 0→1). Scalars are interpolated to EFIT time.
    """
    from MDSplus import Connection
    from scipy.interpolate import interp2d
    from config.app_config import AppConfig

    mds = Connection(AppConfig().MDS_IP)
    mds.openTree(efit_tree, shot_number)

    try:
        time_arr = mds.get('\\gtime').data()
        if efit_tree not in ['efitrt1', 'efitrt2']:
            time_arr = time_arr / 1e3

        r_grid = mds.get('\\r').data()
        z_grid = mds.get('\\z').data()
        simag = mds.get('\\ssimag').data()
        sibry = mds.get('\\ssibry').data()
        psirz = mds.get('\\psirz').data()

        n_time = len(time_arr)
        n_r = len(r_grid)
        psin_all = np.zeros((n_time, n_r))

        for t in range(n_time):
            f = interp2d(r_grid, z_grid, psirz[t])
            psi_z0 = f(r_grid, 0).flatten()
            psin_all[t] = (psi_z0 - simag[t]) / (sibry[t] - simag[t])

        # Equilibrium 1D profiles + scalars for derived quantities
        try:
            nw = int(np.asarray(mds.get('\\mw').data()).flatten()[0])
        except Exception:
            nw = 65
        q_psin_grid = np.linspace(0.0, 1.0, nw)

        def _safe_get_2d(node, fallback_shape):
            try:
                return np.asarray(mds.get(node).data())
            except Exception:
                return np.full(fallback_shape, np.nan)

        def _safe_get_1d(node, fallback_len):
            try:
                return np.asarray(mds.get(node).data()).flatten()
            except Exception:
                return np.full(fallback_len, np.nan)

        def _check_and_get_2d(node, fallback_shape):
            """Coerce to (n_time, nw). 1D-shape (nw,) is broadcast across time."""
            try:
                arr = np.asarray(mds.get(node).data())
                if arr.ndim == 2:
                    return arr, True
                if arr.ndim == 1 and arr.size == fallback_shape[1]:
                    # profile-only (no time dim) — broadcast
                    return np.tile(arr, (fallback_shape[0], 1)), True
                print(f"[EFIT {efit_tree}] {node} shape {arr.shape} unexpected "
                      f"(expected ({fallback_shape[0]},{fallback_shape[1]})) — using NaN")
                return np.full(fallback_shape, np.nan), False
            except Exception as e:
                print(f"[EFIT {efit_tree}] {node} not available: {e}")
                return np.full(fallback_shape, np.nan), False

        def _check_and_get_1d(node, fallback_len):
            """Coerce to (n_time,). Scalars are broadcast; small length mismatch
            (off-by-one, common in KSTAR realtime EFIT) is resized."""
            try:
                arr = np.asarray(mds.get(node).data()).flatten()
                if arr.size == fallback_len:
                    return arr, True
                if arr.size == 1:
                    # scalar (e.g. constant vacuum BT) — broadcast
                    return np.full(fallback_len, float(arr[0])), True
                # Small off-by-N mismatch (e.g. 1524 vs 1525): linearly resample
                # along an implicit normalized index so the value at each target
                # time is the nearest-corresponding source sample.
                if 0 < arr.size < fallback_len * 10:   # accept any reasonable size
                    src_x = np.linspace(0.0, 1.0, arr.size)
                    tgt_x = np.linspace(0.0, 1.0, fallback_len)
                    resampled = np.interp(tgt_x, src_x, arr)
                    print(f"[EFIT {efit_tree}] {node} size {arr.size} ≠ {fallback_len} "
                          f"(time frames) — linearly resampled to {fallback_len}")
                    return resampled, True
                print(f"[EFIT {efit_tree}] {node} size {arr.size} unexpected "
                      f"(expected {fallback_len}) — using NaN")
                return np.full(fallback_len, np.nan), False
            except Exception as e:
                print(f"[EFIT {efit_tree}] {node} not available: {e}")
                return np.full(fallback_len, np.nan), False

        q_t_psi, q_ok = _check_and_get_2d('\\qpsi', (n_time, nw))
        fpol_t_psi, fpol_ok = _check_and_get_2d('\\fpol', (n_time, nw))

        # bcentr: normally one value per time frame. If the size happens to
        # disagree with n_time (e.g. a missing afile frame at shot start —
        # 1524 vs 1525 etc.), fall back to broadcasting the median across the
        # available valid samples. This is justified for KSTAR because the
        # vacuum BT is set by the toroidal coil current and is constant
        # during a pulse.
        try:
            _bc_raw = np.asarray(mds.get('\\bcentr').data()).flatten()
            if _bc_raw.size == n_time:
                bcentr = _bc_raw
                bcentr_ok = True
            else:
                _bc_valid = np.isfinite(_bc_raw) & (_bc_raw != 0)
                if _bc_valid.any():
                    _bc_val = float(np.nanmedian(_bc_raw[_bc_valid]))
                    bcentr = np.full(n_time, _bc_val)
                    bcentr_ok = True
                    print(f"[EFIT {efit_tree}] bcentr size {_bc_raw.size} ≠ {n_time}: "
                          f"using shot median {_bc_val:.3f} T "
                          f"({_bc_valid.sum()} valid samples; KSTAR BT is constant)")
                else:
                    bcentr = np.full(n_time, np.nan)
                    bcentr_ok = False
                    print(f"[EFIT {efit_tree}] bcentr: no valid samples — using NaN")
        except Exception as e:
            bcentr = np.full(n_time, np.nan)
            bcentr_ok = False
            print(f"[EFIT {efit_tree}] bcentr not available: {e}")

        rmaxis, rmaxis_ok = _check_and_get_1d('\\rmaxis', n_time)
        print(f"[EFIT {efit_tree}] extended geometry: "
              f"qpsi={'OK' if q_ok else 'NaN'}, "
              f"fpol={'OK' if fpol_ok else 'NaN'}, "
              f"bcentr={'OK' if bcentr_ok else 'NaN'}, "
              f"rmaxis={'OK' if rmaxis_ok else 'NaN'}")
        # bcentr and rmaxis may be stored as constants; broadcast
        if bcentr.size == 1:
            bcentr = np.full(n_time, float(bcentr[0]))
        if rmaxis.size == 1:
            rmaxis = np.full(n_time, float(rmaxis[0]))
        if bcentr.size != n_time:
            bcentr = np.full(n_time, np.nan)
        if rmaxis.size != n_time:
            rmaxis = np.full(n_time, np.nan)

        # Midplane B_p(R) at z=0 from numerical derivative of ψ:
        #   B_p ≈ B_z(z=0) = (1/R) · dψ/dR
        #   dψ/dR = (sibry - simag) · d(ψ_N)/dR
        Bp_z0 = np.zeros((n_time, n_r))
        for t in range(n_time):
            dpsin_dR = np.gradient(psin_all[t], r_grid)
            dpsi_dR = (sibry[t] - simag[t]) * dpsin_dR
            with np.errstate(divide='ignore', invalid='ignore'):
                Bp_z0[t] = dpsi_dR / r_grid

        return {
            'time': time_arr, 'r_grid': r_grid, 'psin': psin_all,
            'q_psin_grid': q_psin_grid,
            'q_t_psi': q_t_psi,
            'fpol_t_psi': fpol_t_psi,
            'bcentr': bcentr,
            'rmaxis': rmaxis,
            # Δψ = sibry - simag is needed for derivatives w.r.t. ψ (poloidal flux)
            # required by Sauter bootstrap / neoclassical conductivity.
            'simag': simag,
            'sibry': sibry,
            # B_p(R) at z=0 midplane (n_efit_time, n_r) — used by full neoclassical poloidal velocity
            'Bp_z0': Bp_z0,
        }
    finally:
        mds.closeTree(efit_tree, shot_number)


def load_ces_raw(shot_number, analysis_type='nn'):
    """Load CES raw data (Ti, vT) from kstar tree."""
    from MDSplus import Connection
    from config.app_config import AppConfig

    prefix = '\\CESNN_' if analysis_type == 'nn' else '\\CES_'
    mds = Connection(AppConfig().MDS_IP)
    mds.openTree('kstar', shot_number)

    try:
        time_arr = mds.get(f'dim_of({prefix}TI01)').data()
        R, Ti, Ti_err, vT, vT_err = [], [], [], [], []
        for ch in range(1, 33):
            try:
                R.append(mds.get(f'{prefix}RT{ch:02d}').data()[0] * 1e-3)
                Ti.append(mds.get(f'{prefix}TI{ch:02d}').data() * 1e-3)
                Ti_err.append(mds.get(f'{prefix}TI{ch:02d}:err_bar').data() * 1e-3)
                vT.append(mds.get(f'{prefix}VT{ch:02d}').data())
                vT_err.append(mds.get(f'{prefix}VT{ch:02d}:err_bar').data())
            except Exception:
                pass
        mds.closeTree('kstar', shot_number)
        if not R:
            raise RuntimeError("No CES channels loaded")
        return {
            'time': time_arr, 'R': np.array(R),
            'Ti': np.array(Ti), 'Ti_err': np.array(Ti_err),
            'vT': np.array(vT), 'vT_err': np.array(vT_err),
        }
    except Exception as e:
        try:
            mds.closeTree('kstar', shot_number)
        except Exception:
            pass
        raise RuntimeError(f"CES load failed: {e}")


def load_thomson_raw(shot_number):
    """Load Thomson raw data (Te, ne) from kstar tree."""
    from MDSplus import Connection
    from config.app_config import AppConfig

    mds = Connection(AppConfig().MDS_IP)
    mds.openTree('kstar', shot_number)
    core_pos, edge_pos = _get_thomson_positions(shot_number)

    try:
        R, Te, Te_err, ne, ne_err = [], [], [], [], []
        time_arr = None

        for i in range(14):
            try:
                te = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_TE').data().flatten()
                ne_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_NE').data().flatten()
                te_errh = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_TERRH').data().flatten()
                ne_errh = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_NERRH').data().flatten()
                if time_arr is None:
                    time_arr = mds.get(f'dim_of(\\TS_CORE{i+1}:CORE{i+1}_TE)').data().flatten()
                if np.nanmax(te) > 0 and np.nanmax(ne_raw) > 0 and i < len(core_pos):
                    R.append(core_pos[i])
                    Te.append(te * 1e-3)
                    Te_err.append(te_errh * 1e-3)
                    ne.append(ne_raw * 1e-19)
                    ne_err.append(ne_errh * 1e-19)
            except Exception:
                pass

        for i in range(17):
            try:
                te = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_TE').data().flatten()
                ne_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_NE').data().flatten()
                te_errh = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_TERRH').data().flatten()
                ne_errh = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_NERRH').data().flatten()
                if time_arr is None:
                    time_arr = mds.get(f'dim_of(\\TS_EDGE{i+1}:EDGE{i+1}_TE)').data().flatten()
                if np.nanmax(te) > 0 and np.nanmax(ne_raw) > 0 and i < len(edge_pos):
                    R.append(edge_pos[i])
                    Te.append(te * 1e-3)
                    Te_err.append(te_errh * 1e-3)
                    ne.append(ne_raw * 1e-19)
                    ne_err.append(ne_errh * 1e-19)
            except Exception:
                pass

        mds.closeTree('kstar', shot_number)
        if not R:
            raise RuntimeError("No Thomson channels loaded")
        idx = np.argsort(R)
        return {
            'time': time_arr, 'R': np.array(R)[idx],
            'Te': np.array(Te)[idx], 'Te_err': np.array(Te_err)[idx],
            'ne': np.array(ne)[idx], 'ne_err': np.array(ne_err)[idx],
        }
    except Exception as e:
        try:
            mds.closeTree('kstar', shot_number)
        except Exception:
            pass
        raise RuntimeError(f"Thomson load failed: {e}")


def map_R_to_psin(R_positions, efit_data, t_target):
    """Map R positions to psi_N using EFIT data at given time."""
    from scipy.interpolate import interp1d
    t_idx = np.argmin(np.abs(efit_data['time'] - t_target))
    f = interp1d(efit_data['r_grid'], efit_data['psin'][t_idx],
                 fill_value='extrapolate')
    return np.abs(f(R_positions))
