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

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader

PARAM_TREES = {'Ti': 'BITI', 'vT': 'BIVT', 'Te': 'BITE', 'ne': 'BINE'}

# Thomson channel R positions [mm] by shot range
_THOMSON_POSITIONS = [
    {'range': (0, 21778),
     'core': [1806,1826,1848,1871,1894,1917,1942,1966,1991,2016,2041,2068,2093,2119],
     'edge': [2124,2137,2143,2149,2156,2162,2177,2191,2202,2216,2229,2242,2257,2271,2285,2297,2311]},
    {'range': (24082, 27400),
     'core': [1802,1823,1845,1866,1889,1912,1937,1961,1986,2011,2036,2061,2086,2114],
     'edge': [2114,2126,2138,2150,2164,2177,2189,2196,2202,2208,2214,2221,2235,2248,2261,2275,2289]},
    {'range': (27401, 30445),
     'core': [1797,1818,1841,1862,1884,1908,1931,1954,1979,2004,2030,2056,2082,2108],
     'edge': [2108,2120,2133,2146,2153,2171,2183,2190,2197,2203,2209,2216,2229,2243,2254,2269,2282]},
    {'range': (30446, 32768),
     'core': [1800,1820,1842,1864,1887,1909,1933,1956,1983,2008,2033,2058,2083,2109],
     'edge': [2111,2124,2136,2149,2161,2174,2187,2194,2200,2207,2212,2219,2232,2245,2259,2272,2285]},
    {'range': (32769, 34836),
     'core': [1800,1820,1842,1864,1887,1911,1933,1958,1981,2006,2030,2058,2083,2110],
     'edge': [2111,2124,2137,2147,2162,2175,2187,2195,2200,2205,2210,2217,2233,2245,2259,2271,2290]},
    {'range': (34837, 37741),
     'core': [1795,1817,1839,1861,1884,1908,1930,1954,1978,2003,2029,2054,2081,2108],
     'edge': [2108,2116,2125,2136,2145,2155,2163,2171,2180,2191,2200,2213,2222,2232,2253,2274,2296]},
    {'range': (37742, 999999),
     'core': [1785,1807,1828,1851,1873,1896,1920,1944,1968,1993,2018,2043,2070,2095],
     'edge': [2096,2105,2113,2123,2131,2141,2149,2158,2168,2178,2197,2217,2237,2258,2280]},
]


def _try_get(mds, path, default=None):
    try:
        return mds.get(path).data()
    except Exception:
        return default


def _get_thomson_positions(shot):
    for entry in _THOMSON_POSITIONS:
        lo, hi = entry['range']
        if lo <= shot <= hi:
            return ([r / 1000.0 for r in entry['core']],
                    [r / 1000.0 for r in entry['edge'] if r is not None])
    e = _THOMSON_POSITIONS[-1]
    return ([r/1000 for r in e['core']], [r/1000 for r in e['edge'] if r])


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
    """Load EFIT psi_N(R) at midplane (z=0)."""
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

        return {'time': time_arr, 'r_grid': r_grid, 'psin': psin_all}
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
