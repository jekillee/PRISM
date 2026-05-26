"""
Derived plasma quantities from BIPROFILE Bayesian fits + EFIT geometry.

Inputs (per shot, loaded by data_loaders/biprofile_loader.py):
    bi:   dict {param: {'psin', 'time', 'mean', 'unc', ...}} for param in {Ti, vT, Te, ne}
          - psin: 1D array (n_psi,), normalized poloidal flux grid
          - time: 1D array (n_time,), in seconds
          - mean, unc: 2D (n_psi, n_time)
          - Ti, Te in [keV]; ne in [1e19 m^-3]; vT in [km/s]
    efit: dict with geometry — see biprofile_loader.load_efit_psin
          - r_grid (1D), psin (n_efit_time, n_r) on z=0 midplane (LFS+HFS combined)
          - q_psin_grid (1D, n_psi_eq), q (n_efit_time, n_psi_eq)
          - fpol (n_efit_time, n_psi_eq), bcentr (n_efit_time,), rmaxis (n_efit_time,)

Quantities are returned as a flat dict keyed by short name:
    {name: {'psin': psin_1d, 'time': time_1d, 'value': 2D[psin, time],
            'label': latex_str, 'unit': str}}

All physics in SI internally; conversions on output handled per quantity.

Constants follow scipy.constants.
"""

import numpy as np
from scipy.interpolate import interp1d


# Physical constants (SI)
E_CHARGE = 1.602176634e-19      # C
M_E = 9.1093837015e-31          # kg
M_P = 1.67262192369e-27         # kg
EPS0 = 8.8541878128e-12         # F/m
MU0 = 1.25663706212e-6          # T·m/A
K_BOLTZ = 1.380649e-23          # J/K  (not used; energies in eV/J directly)

# Ion species: (mass_amu, charge)
ION_SPECIES = {
    'H': (1.00784, 1),
    'D': (2.01410, 1),
    'T': (3.01605, 1),
}
# Impurity (CES is typically C6+ on KSTAR)
IMP_SPECIES = {
    'C6': (12.011, 6),
    'B5': (10.81, 5),
}


# ------------------------------------------------------------------
# Catalog of available derived quantities (for UI population)
# ------------------------------------------------------------------

PROFILE_QUANTITIES = [
    # (key, name, label, unit, requires)
    ('p_e',           'Electron pressure',                 r'$p_e$',                'kPa',                ('Te', 'ne')),
    ('p_i',           'Ion pressure',                      r'$p_i$',                'kPa',                ('Ti', 'ne')),
    ('p_tot',         'Total pressure',                    r'$p_\mathrm{tot}$',     'kPa',                ('Te', 'Ti', 'ne')),
    ('dTe_dR',        'Te gradient (LFS)',                 r'$dT_e/dR$',            'keV/m',              ('Te',)),
    ('dTi_dR',        'Ti gradient (LFS)',                 r'$dT_i/dR$',            'keV/m',              ('Ti',)),
    ('dne_dR',        'ne gradient (LFS)',                 r'$dn_e/dR$',            r'$10^{19}$/m$^4$',   ('ne',)),
    ('R_over_LTe',    'Inverse Te scale length',           r'$R/L_{T_e}$',          '-',                  ('Te',)),
    ('R_over_LTi',    'Inverse Ti scale length',           r'$R/L_{T_i}$',          '-',                  ('Ti',)),
    ('R_over_Lne',    'Inverse ne scale length',           r'$R/L_{n_e}$',          '-',                  ('ne',)),
    ('beta_e',        'Electron toroidal beta',            r'$\beta_e$',            '%',                  ('Te', 'ne')),
    ('beta_i',        'Ion toroidal beta',                 r'$\beta_i$',            '%',                  ('Ti', 'ne')),
    ('beta_tot',      'Total toroidal beta',               r'$\beta_\mathrm{tot}$', '%',                  ('Te', 'Ti', 'ne')),
    ('nu_star_e',     'Electron collisionality',           r'$\nu^\ast_e$',         '-',                  ('Te', 'ne')),
    ('nu_star_i',     'Ion collisionality',                r'$\nu^\ast_i$',         '-',                  ('Ti', 'ne')),
    ('omega_pe',      'Electron plasma frequency',         r'$\omega_{pe}$',        'rad/s',              ('ne',)),
    ('omega_ci',      'Ion cyclotron frequency',           r'$\omega_{ci}$',        'rad/s',              ()),
    ('c_s',           'Ion sound speed',                   r'$c_s$',                'm/s',                ('Te',)),
    ('rho_i',         'Ion gyroradius (Larmor radius)',    r'$\rho_i$',             'm',                  ('Ti',)),
    ('E_r_vp0',       'Radial electric field (vθ=0)',      r'$E_r$ (vθ=0)',         'kV/m',               ('Ti', 'ne', 'vT')),
    ('E_r_neo',       'Radial electric field (neoclassical vθ)', r'$E_r$ (neo)',    'kV/m',               ('Ti', 'ne', 'vT')),
    ('vpol_neo_i',    'Neoclassical poloidal velocity (main ion D)', r'$v_{\theta,i}^\mathrm{neo}$',  'km/s',     ('Ti', 'ne')),
    ('vpol_neo_imp',  'Neoclassical poloidal velocity (impurity C⁶⁺)', r'$v_{\theta,\mathrm{imp}}^\mathrm{neo}$', 'km/s', ('Ti', 'ne')),
    ('omega_neo_diff','Neoclassical ω_imp − ω_main',       r'$\omega_I - \omega_i$', 'krad/s',          ('Ti', 'ne')),
    ('sigma_neo',     'Neoclassical conductivity',         r'$\sigma_\mathrm{neo}$', '1/(Ω·m)',          ('Te', 'ne')),
    ('eta_neo',       'Neoclassical resistivity',          r'$\eta_\mathrm{neo}$',  'Ω·m',                ('Te', 'ne')),
    ('j_BS',          'Bootstrap current density',         r'$\langle j_\mathrm{BS}\cdot B\rangle/B_0$', 'A/m²', ('Te', 'Ti', 'ne')),
]

PROFILE_QUANTITY_KEYS = [q[0] for q in PROFILE_QUANTITIES]


def get_quantity_meta(key):
    """Return (name, label, unit, requires) for a quantity key, or None."""
    for k, name, lab, unit, req in PROFILE_QUANTITIES:
        if k == key:
            return name, lab, unit, req
    return None


# ------------------------------------------------------------------
# Geometry helpers
# ------------------------------------------------------------------

def _interp_efit_to_bi_time(efit_arr, efit_time, bi_time):
    """Interpolate an EFIT-time-dependent array to biprofile time axis.

    efit_arr: shape (n_efit_time, ...) or (n_efit_time,)
    Returns: shape (n_bi_time, ...) or (n_bi_time,)
    """
    efit_arr = np.asarray(efit_arr)
    if efit_arr.ndim == 1:
        f = interp1d(efit_time, efit_arr, bounds_error=False,
                     fill_value=(efit_arr[0], efit_arr[-1]))
        return f(bi_time)
    # >=2D: interp along axis 0
    out = np.empty((len(bi_time),) + efit_arr.shape[1:], dtype=float)
    for idx in np.ndindex(efit_arr.shape[1:]):
        col = efit_arr[(slice(None),) + idx]
        f = interp1d(efit_time, col, bounds_error=False,
                     fill_value=(col[0], col[-1]))
        out[(slice(None),) + idx] = f(bi_time)
    return out


def _build_R_of_psin(efit, bi_time, psin_grid):
    """Build R_LFS(psin, time) on biprofile (time, psin) grid.

    Returns R [m] of shape (n_bi_time, n_psi).
    """
    e_time = efit['time']
    r_grid = efit['r_grid']
    psin_2d = efit['psin']   # (n_efit_time, n_r) — psi_N along z=0 horizontal line

    n_bi = len(bi_time)
    n_psi = len(psin_grid)
    R_out = np.full((n_bi, n_psi), np.nan)

    for it, t in enumerate(bi_time):
        et_idx = int(np.argmin(np.abs(e_time - t)))
        psin_line = psin_2d[et_idx]
        # Find magnetic-axis index = minimum of psi_N (works for both single-null/limited)
        i_axis = int(np.nanargmin(psin_line))
        # LFS branch only
        R_lfs = r_grid[i_axis:]
        psin_lfs = psin_line[i_axis:]
        # Ensure monotonic
        order = np.argsort(psin_lfs)
        try:
            f = interp1d(psin_lfs[order], R_lfs[order],
                         bounds_error=False, fill_value=np.nan)
            R_out[it] = f(psin_grid)
        except Exception:
            pass
    return R_out


def _interp_midplane_to_psin(efit, bi_time, psin_grid, key):
    """Interpolate a midplane R-quantity (e.g. Bp_z0(R)) onto (bi_time, psin)
    by inverting psi_N(R) on the LFS branch. Shape: (n_bi_time, n_psi)."""
    e_time = efit.get('time')
    r_grid = efit.get('r_grid')
    psin_2d = efit.get('psin')
    arr = efit.get(key)
    if arr is None or e_time is None or r_grid is None or psin_2d is None:
        return np.full((len(bi_time), len(psin_grid)), np.nan)
    n_bi = len(bi_time); n_psi = len(psin_grid)
    out = np.full((n_bi, n_psi), np.nan)
    for it, t in enumerate(bi_time):
        et_idx = int(np.argmin(np.abs(e_time - t)))
        psin_line = psin_2d[et_idx]
        Y_line = arr[et_idx]
        i_axis = int(np.nanargmin(psin_line))
        psin_lfs = psin_line[i_axis:]
        Y_lfs = Y_line[i_axis:]
        order = np.argsort(psin_lfs)
        try:
            f = interp1d(psin_lfs[order], Y_lfs[order],
                         bounds_error=False, fill_value=np.nan)
            out[it] = f(psin_grid)
        except Exception:
            pass
    return out


def _interp_eq_profile_to_bi(efit, bi_time, psin_grid, name):
    """Interp EFIT-1D-profile (e.g. qpsi, fpol) onto (bi_time, psin_grid).

    Returns shape (n_bi_time, n_psi); NaN where ψ_N out of [0,1] OR when the
    requested key (or its grid 'q_psin_grid') is missing from `efit` (graceful
    fallback so quantities not depending on `name` still compute).
    """
    e_time = efit.get('time')
    psin_eq = efit.get('q_psin_grid')
    prof = efit.get(name)
    n_bi = len(bi_time)
    n_psi = len(psin_grid)
    if psin_eq is None or prof is None or e_time is None:
        return np.full((n_bi, n_psi), np.nan)

    out = np.full((n_bi, n_psi), np.nan)
    for it, t in enumerate(bi_time):
        et_idx = int(np.argmin(np.abs(e_time - t)))
        try:
            f = interp1d(psin_eq, prof[et_idx], bounds_error=False, fill_value=np.nan)
            valid_mask = (psin_grid >= 0) & (psin_grid <= 1.0)
            out[it, valid_mask] = f(psin_grid[valid_mask])
        except Exception:
            pass
    return out


# ------------------------------------------------------------------
# Core computation
# ------------------------------------------------------------------

class DerivedQuantities:
    """Computes derived plasma quantities from biprofile + EFIT.

    Convention: outputs stored as `value[psi, time]` to match biprofile layout.
    """

    def __init__(self, bi, efit, zeff=2.0, ion_species='D',
                 imp_species='C6', er_mode='neo'):
        """
        bi:   dict from biprofile_loader.load_biprofile output (keyed by param)
        efit: dict from biprofile_loader.load_efit_psin output (extended with q, fpol, bcentr)
        zeff: effective charge (default 2)
        ion_species: 'H' | 'D' | 'T' (main ion for collisionality, gyrofrequency, c_s)
        imp_species: 'C6' | 'B5' (impurity for Er rotation-balance calculation)
        er_mode: 'vp_zero' | 'neo' (poloidal velocity assumption for E_r);
                 'neo' uses the full neoclassical (Kim-Diamond-Groebner) v_θ_imp.
        """
        self.bi = bi
        self.efit = efit
        self.zeff = float(zeff)
        self.ion = ion_species
        self.imp = imp_species
        self.er_mode = er_mode

        # Reference grid: pick the first available param's (psin, time).
        # All other params will be 2D-interpolated onto this grid so derived
        # ops like ne * Ti, vT / R, etc. always broadcast cleanly.
        ref_key = next((k for k in ('Te', 'ne', 'Ti', 'vT') if k in bi), None)
        if ref_key is None:
            raise ValueError("No biprofile parameters available")
        self.psin = np.asarray(bi[ref_key]['psin'])
        self.time = np.asarray(bi[ref_key]['time'])

        # Pre-interpolate each available biprofile param onto the reference
        # (psin, time) grid. Output shape is (n_time_ref, n_psi_ref) so it
        # matches all EFIT-derived arrays (self.R, self.q, self.Bt, ...).
        self._fields = {}
        for k in ('Ti', 'vT', 'Te', 'ne'):
            if k not in bi:
                continue
            self._fields[k] = self._interp_param_to_ref(bi[k])

        # Pre-interpolate EFIT geom onto bi (time, psin) grid
        # R_LFS(psin, t) for R/L_x and gradient conversion
        self.R = _build_R_of_psin(efit, self.time, self.psin)  # (n_t, n_psi)
        # q(psin, t) — collisionality
        self.q = _interp_eq_profile_to_bi(efit, self.time, self.psin, 'q_t_psi')
        # F = R*Bt profile -- gives Bt(R) = F/R
        self.fpol = _interp_eq_profile_to_bi(efit, self.time, self.psin, 'fpol_t_psi')
        # Scalars vs time — fall back to NaN if not provided by EFIT loader
        n_t = len(self.time)
        if 'bcentr' in efit and 'time' in efit:
            self.bcentr = _interp_efit_to_bi_time(efit['bcentr'], efit['time'], self.time)
        else:
            self.bcentr = np.full(n_t, np.nan)
        if 'rmaxis' in efit and 'time' in efit:
            self.rmaxis = _interp_efit_to_bi_time(efit['rmaxis'], efit['time'], self.time)
        else:
            # crude fallback: use min of midplane ψ_N → R(ψ_min) as axis R
            r_grid = efit.get('r_grid')
            psin_2d = efit.get('psin')
            if r_grid is not None and psin_2d is not None:
                self.rmaxis = np.array([r_grid[int(np.nanargmin(psin_2d[
                    int(np.argmin(np.abs(efit['time'] - t)))]))]
                    for t in self.time])
            else:
                self.rmaxis = np.full(n_t, np.nan)

        # B_t(psin, t) at LFS: B_t = F/R  (where R is LFS major radius at that ψ_N)
        with np.errstate(invalid='ignore', divide='ignore'):
            self.Bt = self.fpol / self.R                          # (n_t, n_psi), T
        # B_p(ψ_N, t) at LFS midplane — numerical derivative of ψ_N at z=0
        # (loaded by load_efit_psin as 'Bp_z0'). Fallback to q-based rough
        # estimate only if missing.
        if 'Bp_z0' in efit:
            self.Bp = _interp_midplane_to_psin(efit, self.time, self.psin, 'Bp_z0')
        else:
            with np.errstate(invalid='ignore', divide='ignore'):
                self.Bp = self.Bt / np.where(np.abs(self.q) > 1e-6, self.q, np.nan)

        # ε = a/R, take a ≈ R_LCFS - R_axis from R(psin=1) - R_axis
        # rmaxis is per-time; take a(t) from R at psin closest to 1.
        a_idx = int(np.argmin(np.abs(self.psin - 1.0)))
        self.a_minor = self.R[:, a_idx] - self.rmaxis           # (n_t,) m
        # eps profile (only meaningful for ψ_N≤1)
        with np.errstate(invalid='ignore', divide='ignore'):
            self.eps = (self.R - self.rmaxis[:, None]) / self.rmaxis[:, None]
            self.eps = np.clip(self.eps, 1e-6, None)

        # Δψ = ψ_LCFS - ψ_axis (per time), Wb. Used for d/dψ conversions in Sauter.
        if 'simag' in efit and 'sibry' in efit and 'time' in efit:
            simag = _interp_efit_to_bi_time(efit['simag'], efit['time'], self.time)
            sibry = _interp_efit_to_bi_time(efit['sibry'], efit['time'], self.time)
            self.delta_psi = sibry - simag
        else:
            self.delta_psi = np.full(n_t, np.nan)

        # Trapped particle fraction f_T(ψ_N, t), Sauter Eq. 14
        eps_c = np.clip(self.eps, 1e-6, 0.999)
        self.fT = np.clip(
            1.0 - (1.0 - eps_c)**2 /
            (np.sqrt(1.0 - eps_c**2) * (1.0 + 1.46 * np.sqrt(eps_c))),
            0.0, 1.0)

        # Total field magnitude at midplane LFS
        self.B_tot = np.sqrt(self.Bt**2 + self.Bp**2)

        # Flux-surface-averaged ⟨B²⟩ — analytic approximation for circular
        # surfaces with finite ε (good to a few % for KSTAR ε ~ 0.2-0.3):
        #   ⟨B²⟩ ≈ B_t,axis² · (1 + 1.5 ε²) + B_p,midplane²
        # The first term comes from ⟨R₀²/R²⟩_FSA = 1 + 1.5ε² for circular ψ.
        # True flux-surface contour integration is left as a future upgrade.
        B0 = self.bcentr[:, None]
        self.B2_FSA = B0**2 * (1.0 + 1.5 * eps_c**2) + self.Bp**2

        # Cache
        self._cache = {}
        self._sauter_cache = None
        self._neo_cache = None

    # -- helpers ----------------------------------------------------

    def _interp_param_to_ref(self, param_data):
        """Interpolate a biprofile parameter (psin × time) onto the reference
        (self.psin, self.time) grid, returning a (n_time_ref, n_psi_ref) array.

        Skips interpolation if grids are already identical (common case).
        """
        from scipy.interpolate import RegularGridInterpolator
        src_psin = np.asarray(param_data['psin'])
        src_time = np.asarray(param_data['time'])
        # biprofile 'mean' is stored as (n_psi, n_time)
        src = np.asarray(param_data['mean'])
        # Fast path: identical grids
        if (src_psin.shape == self.psin.shape
                and src_time.shape == self.time.shape
                and np.allclose(src_psin, self.psin)
                and np.allclose(src_time, self.time)):
            return src.T   # → (n_time, n_psi)
        # General path: 2D interpolation (psin, time)
        try:
            f = RegularGridInterpolator(
                (src_psin, src_time), src,
                bounds_error=False, fill_value=np.nan)
            P, T = np.meshgrid(self.psin, self.time, indexing='ij')
            pts = np.stack([P.ravel(), T.ravel()], axis=-1)
            out = f(pts).reshape(P.shape)   # (n_psi, n_time)
            return out.T                    # (n_time, n_psi)
        except Exception as e:
            print(f"[Derived] Param interp failed: {e}")
            return np.full((len(self.time), len(self.psin)), np.nan)

    def _T(self, key):
        """Legacy: return field as (n_psi, n_time)."""
        f = self._fields.get(key) if hasattr(self, '_fields') else None
        return None if f is None else f.T

    def _T_T(self, key):
        """Return field as (n_time, n_psi) on the reference grid — directly
        compatible with EFIT-derived arrays (self.R, self.q, self.Bt, ...)."""
        return self._fields.get(key) if hasattr(self, '_fields') else None

    def _output(self, value_TT, label, unit):
        """Wrap a (n_time, n_psi) array as standard output (psin, time, value[psin,time])."""
        return {
            'psin': self.psin, 'time': self.time,
            'value': value_TT.T,    # back to (n_psi, n_time)
            'label': label, 'unit': unit,
        }

    def _grad_R(self, y_TT):
        """∂y/∂R at LFS midplane.  y_TT shape (n_t, n_psi). Returns same shape.

        Uses non-uniform finite differences along psin axis, divided by dR/dpsin.
        """
        n_t, n_psi = y_TT.shape
        out = np.full_like(y_TT, np.nan, dtype=float)
        for it in range(n_t):
            R_row = self.R[it]
            y_row = y_TT[it]
            valid = np.isfinite(R_row) & np.isfinite(y_row)
            if valid.sum() < 3:
                continue
            R_v = R_row[valid]
            y_v = y_row[valid]
            order = np.argsort(R_v)
            R_s, y_s = R_v[order], y_v[order]
            dy_dR = np.gradient(y_s, R_s)
            # write back in original ordering
            out_v = np.full(R_v.shape, np.nan)
            out_v[order] = dy_dR
            out_row = np.full(n_psi, np.nan)
            out_row[valid] = out_v
            out[it] = out_row
        return out

    # -- quantities -------------------------------------------------

    def _pressure_e(self):
        """p_e = n_e · T_e in kPa.  ne[1e19], Te[keV] → p[kPa] = ne_19·Te_keV·1.602."""
        ne, Te = self._T_T('ne'), self._T_T('Te')
        if ne is None or Te is None:
            return None
        p = ne * Te * 1.602          # kPa
        return self._output(p, r'$p_e$', 'kPa')

    def _pressure_i(self):
        """p_i ≈ (n_e/Z_eff) · T_i  (single-ion equivalent, scaled by 1/Zeff).

        Quasineutrality with one impurity species: n_i = n_e·(Z_imp - Zeff)/(Z_imp·(Z_imp - 1)·...) — complex.
        For first cut: p_i ≈ n_e · T_i / Z_eff is the conventional approximation when Z_main=1.
        """
        ne, Ti = self._T_T('ne'), self._T_T('Ti')
        if ne is None or Ti is None:
            return None
        # main ion density n_i ≈ n_e·(Z_imp - Zeff)/(Z_imp - Z_main); Z_main=1, Z_imp from imp_species
        Z_imp = IMP_SPECIES[self.imp][1]
        ni_over_ne = (Z_imp - self.zeff) / (Z_imp - 1.0)
        ni_over_ne = max(ni_over_ne, 0.0)
        p = ne * ni_over_ne * Ti * 1.602
        return self._output(p, r'$p_i$', 'kPa')

    def _pressure_tot(self):
        pe = self._pressure_e()
        pi = self._pressure_i()
        if pe is None or pi is None:
            return None
        v = pe['value'] + pi['value']
        return {'psin': self.psin, 'time': self.time, 'value': v,
                'label': r'$p_\mathrm{tot}$', 'unit': 'kPa'}

    def _dT_dR(self, key, label):
        T = self._T_T(key)
        if T is None:
            return None
        dT_dR = self._grad_R(T)  # keV/m
        return self._output(dT_dR, label, 'keV/m')

    def _dne_dR(self):
        ne = self._T_T('ne')
        if ne is None:
            return None
        return self._output(self._grad_R(ne), r'$dn_e/dR$', r'$10^{19}$/m$^4$')

    def _R_over_LT(self, key, label):
        """R/L_T = -R · (1/T) · dT/dR  evaluated at LFS midplane."""
        T = self._T_T(key)
        if T is None:
            return None
        dT_dR = self._grad_R(T)
        with np.errstate(invalid='ignore', divide='ignore'):
            R_over_LT = -self.R * dT_dR / T
        return self._output(R_over_LT, label, '-')

    def _R_over_Ln(self):
        ne = self._T_T('ne')
        if ne is None:
            return None
        dn_dR = self._grad_R(ne)
        with np.errstate(invalid='ignore', divide='ignore'):
            R_over_Ln = -self.R * dn_dR / ne
        return self._output(R_over_Ln, r'$R/L_{n_e}$', '-')

    def _beta(self, kind):
        """beta = 2·μ0·p / B² (toroidal beta, B = B_t at axis).  Output [%]."""
        # pressure in Pa
        if kind == 'e':
            ne, Te = self._T_T('ne'), self._T_T('Te')
            if ne is None or Te is None:
                return None
            p_Pa = ne * Te * 1.602e3
            label = r'$\beta_e$'
        elif kind == 'i':
            ne, Ti = self._T_T('ne'), self._T_T('Ti')
            if ne is None or Ti is None:
                return None
            Z_imp = IMP_SPECIES[self.imp][1]
            ni_over_ne = max((Z_imp - self.zeff) / (Z_imp - 1.0), 0.0)
            p_Pa = ne * ni_over_ne * Ti * 1.602e3
            label = r'$\beta_i$'
        else:  # tot
            pe = self._beta('e'); pi = self._beta('i')
            if pe is None or pi is None:
                return None
            return {'psin': self.psin, 'time': self.time,
                    'value': pe['value'] + pi['value'],
                    'label': r'$\beta_\mathrm{tot}$', 'unit': '%'}
        with np.errstate(invalid='ignore', divide='ignore'):
            B0 = self.bcentr[:, None]
            beta = 2 * MU0 * p_Pa / B0**2 * 100.0   # percent
        return self._output(beta, label, '%')

    def _log_lambda(self, Te_TT, ne_TT):
        """Coulomb logarithm (Wesson, electron-ion): ln Λ = 31.3 - ln(√n_e/T_e) for T_e in eV, n_e in m^-3.

        Te in keV, ne in 1e19 m^-3 here.
        """
        with np.errstate(invalid='ignore', divide='ignore'):
            Te_eV = Te_TT * 1e3
            ne_m3 = ne_TT * 1e19
            lnL = 31.3 - np.log(np.sqrt(ne_m3) / np.maximum(Te_eV, 1.0))
        return np.clip(lnL, 10.0, 25.0)

    def _nu_star(self, kind):
        """Collisionality ν* (Hinton-Hazeltine).

        ν*_e = ν_ei · q · R / (v_th_e · ε^{3/2})
        ν_ei = (4/3)√(2π) · n_e · e⁴ · Zeff · lnΛ / [(4πε0)² · √m_e · T_e^{3/2}]
        """
        ne_TT = self._T_T('ne')
        if ne_TT is None:
            return None
        if kind == 'e':
            T_TT = self._T_T('Te')
            if T_TT is None:
                return None
            m_s = M_E
            Z_factor = self.zeff      # ν_ee + ν_ei effective
            label = r'$\nu^\ast_e$'
        else:  # 'i'
            T_TT = self._T_T('Ti')
            if T_TT is None:
                return None
            m_s = ION_SPECIES[self.ion][0] * M_P
            Z_factor = self.zeff
            label = r'$\nu^\ast_i$'

        lnL = self._log_lambda(self._T_T('Te'), ne_TT)
        ne_m3 = ne_TT * 1e19
        T_J = T_TT * 1e3 * E_CHARGE

        with np.errstate(invalid='ignore', divide='ignore'):
            # collision frequency for species s with background ions of charge Zeff
            nu_s = ((4.0/3.0) * np.sqrt(2.0 * np.pi)
                    * ne_m3 * (E_CHARGE**4) * Z_factor * lnL
                    / ((4.0 * np.pi * EPS0)**2 * np.sqrt(m_s) * T_J**1.5))
            v_th = np.sqrt(T_J / m_s)
            nu_star = nu_s * self.q * self.rmaxis[:, None] / (v_th * self.eps**1.5)
        return self._output(nu_star, label, '-')

    def _omega_pe(self):
        ne = self._T_T('ne')
        if ne is None:
            return None
        ne_m3 = ne * 1e19
        with np.errstate(invalid='ignore'):
            wpe = np.sqrt(ne_m3 * E_CHARGE**2 / (EPS0 * M_E))
        return self._output(wpe, r'$\omega_{pe}$', 'rad/s')

    def _omega_ci(self):
        m_i = ION_SPECIES[self.ion][0] * M_P
        wci_2d = E_CHARGE * np.abs(self.Bt) / m_i   # (n_t, n_psi)
        return self._output(wci_2d, r'$\omega_{ci}$', 'rad/s')

    def _c_s(self):
        Te = self._T_T('Te')
        if Te is None:
            return None
        m_i = ION_SPECIES[self.ion][0] * M_P
        cs = np.sqrt(Te * 1e3 * E_CHARGE / m_i)
        return self._output(cs, r'$c_s$', 'm/s')

    def _rho_i(self):
        Ti = self._T_T('Ti')
        if Ti is None:
            return None
        m_i = ION_SPECIES[self.ion][0] * M_P
        v_th = np.sqrt(2.0 * Ti * 1e3 * E_CHARGE / m_i)
        with np.errstate(invalid='ignore', divide='ignore'):
            wci = E_CHARGE * np.abs(self.Bt) / m_i
            rho = v_th / wci
        return self._output(rho, r'$\rho_i$', 'm')

    def _Er_kvm(self, er_mode):
        """Compute radial electric field [kV/m] from impurity force balance.

        E_r = (1/(Z_imp·e·n_imp))·∂p_imp/∂R - v_θ·B_φ + v_φ·B_θ

        - n_imp/n_e ≈ const, so ∂p_imp/∂R/n_imp = T_i·∂lnn_e/∂R + ∂T_i/∂R
        - v_φ from CES toroidal rotation (vT, in km/s)
        - v_θ: 0 ('vp_zero') OR full neoclassical analytic ('neo')
              v_θ ≈ k1 · (T_i_J / (Z_imp·e·B_t)) · ∂lnT_i/∂R

        Returns (n_time, n_psi) array in kV/m, or None if required data missing.
        """
        Ti = self._T_T('Ti'); ne = self._T_T('ne'); vT = self._T_T('vT')
        if Ti is None or ne is None or vT is None:
            return None
        Z_imp = IMP_SPECIES[self.imp][1]

        with np.errstate(invalid='ignore', divide='ignore'):
            dTi_dR = self._grad_R(Ti)            # keV/m
            dlnne_dR = self._grad_R(ne) / ne     # 1/m
            term_p = (Ti * 1e3 * dlnne_dR + dTi_dR * 1e3) / Z_imp     # V/m
            term_t = (vT * 1e3) * self.Bp                              # V/m
            if er_mode == 'vp_zero':
                v_pol = np.zeros_like(Ti)
            else:  # 'neo' — full neoclassical (Kim-Diamond-Groebner) regime-correct v_θ
                v_pol_dict = self._vpol_neo_imp()
                if v_pol_dict is None:
                    v_pol = np.zeros_like(Ti)
                else:
                    # vpol returned as (n_psi, n_time) km/s; convert back to
                    # (n_t, n_psi) m/s for the local computation
                    v_pol = v_pol_dict['value'].T * 1e3
            term_pol = -v_pol * self.Bt
            Er_kVm = (term_p + term_t + term_pol) / 1e3
        return Er_kVm

    def _E_r_vp0(self):
        v = self._Er_kvm('vp_zero')
        return None if v is None else self._output(v, r'$E_r$ (vθ=0)', 'kV/m')

    def _E_r_neo(self):
        v = self._Er_kvm('neo')
        return None if v is None else self._output(v, r'$E_r$ (neo)', 'kV/m')

    # ------------------------------------------------------------------
    # Full neoclassical poloidal velocity (Kim-Diamond-Groebner,
    # Phys. Fluids B 3 (1991) 2050). Regime coverage: banana / plateau /
    # Pfirsch-Schlüter — via viscosity coefficients μ_xx that absorb the
    # collisionality dependence through g = f_T/f_C and the impurity
    # strength α. The simple "k₁ = −1.17" is the asymptotic banana limit
    # of K₁ in this framework.
    # ------------------------------------------------------------------

    def _compute_neo(self):
        """Compute and cache the full neoclassical (Kim-Diamond-Groebner)
        framework.

        Returns dict with K1, K2, vti, rhoi, n_i, n_imp, fT/fc, plus all
        intermediate quantities (mu_00, mu_01, mu_11, β, α, g, D).
        Returns None if required Ti / ne / B data is missing.
        """
        if self._neo_cache is not None:
            return self._neo_cache
        Ti = self._T_T('Ti'); ne = self._T_T('ne')
        if Ti is None or ne is None:
            return None
        Bt = self.Bt; Bp = self.Bp
        if Bt is None or Bp is None:
            return None

        # Species
        Z_imp = IMP_SPECIES[self.imp][1]            # 6 for C
        m_imp = IMP_SPECIES[self.imp][0] * M_P
        Z_i = 1                                     # main ion D
        m_i = ION_SPECIES[self.ion][0] * M_P

        # Quasineutrality split of n_e into n_i and n_imp (constant Zeff)
        ni_over_ne = max((Z_imp - self.zeff) / (Z_imp - 1.0), 0.0)
        nimp_over_ne = max((self.zeff - 1.0) / (Z_imp * (Z_imp - 1.0)), 0.0)

        ne_m3 = ne * 1e19
        n_i = ne_m3 * ni_over_ne
        n_imp = ne_m3 * nimp_over_ne

        # Impurity strength parameter α (after Eq. 9)
        with np.errstate(invalid='ignore', divide='ignore'):
            alpha = (n_imp * Z_imp**2) / np.maximum(n_i * Z_i**2, 1e-30)

        # Trapped / circulating fractions per Kim-Diamond-Groebner App. C
        eps_c = np.clip(self.eps, 1e-6, 0.999)
        fc = 1.0 - 1.46 * np.sqrt(eps_c) + 0.46 * eps_c * np.sqrt(eps_c)
        ft = 1.0 - fc
        with np.errstate(invalid='ignore', divide='ignore'):
            g = ft / fc

        # Viscosity coefficients (Eq. C15) for the main ion
        sqrt2 = np.sqrt(2.0)
        ln1psqrt2 = np.log(1.0 + sqrt2)
        with np.errstate(invalid='ignore'):
            mu_00 = g * (alpha + sqrt2 - ln1psqrt2)
            mu_01 = g * (1.5 * alpha + 4.0 / sqrt2 - 2.5 * ln1psqrt2)
            mu_11 = g * (3.25 * alpha + (9.75) / sqrt2 - 6.25 * ln1psqrt2)

        # Thermal speeds (Eq. 34)
        Ti_J = Ti * 1e3 * E_CHARGE
        vti = np.sqrt(2.0 * Ti_J / m_i)             # m/s
        vtI = np.sqrt(2.0 * Ti_J / m_imp)           # m/s (impurity same T as ion)

        # Larmor radii (m)
        with np.errstate(invalid='ignore', divide='ignore'):
            rho_i = (m_i * vti) / (Z_i * E_CHARGE * self.B_tot)
            rho_ip = (m_i * vti) / (Z_i * E_CHARGE * np.abs(self.Bp))

        # β coupling (Eq. 32)
        with np.errstate(invalid='ignore', divide='ignore'):
            beta = ((27.0 / 4.0) ** 2 * (m_i / m_imp) ** 2
                    / (7.5 + np.sqrt(2.0 * alpha) * (vti / vtI)))

        # D, K1, K2 (Eq. 29-31)
        with np.errstate(invalid='ignore', divide='ignore'):
            D = mu_00 * (mu_11 + sqrt2 + alpha - alpha * beta) - mu_01**2
            D = np.where(np.abs(D) > 1e-30, D, np.nan)
            K1 = mu_01 * (sqrt2 + alpha - alpha * beta) / D
            K2 = (mu_00 * mu_11 - mu_01**2) / D

        # Pressure gradients needed for V_pol_imp (Eq. 34)
        # p_i = n_i T_i, p_imp = n_imp T_imp  (T_imp ≈ T_i — CES measures C6+)
        p_i = n_i * Ti_J
        p_imp = n_imp * Ti_J
        T_imp = Ti                                  # same temperature

        with np.errstate(invalid='ignore', divide='ignore'):
            dTi_dR = self._grad_R(Ti)               # keV/m
            dp_i_dR = self._grad_R(p_i)             # Pa/m
            dpimp_dR = self._grad_R(p_imp)
            LTi_inv = dTi_dR / Ti                   # 1/m
            Lpi_inv = dp_i_dR / np.maximum(p_i, 1e-30)
            Lpimp_inv = dpimp_dR / np.maximum(p_imp, 1e-30)

        self._neo_cache = {
            'alpha': alpha, 'fc': fc, 'ft': ft, 'g': g,
            'mu_00': mu_00, 'mu_01': mu_01, 'mu_11': mu_11,
            'beta': beta, 'D': D, 'K1': K1, 'K2': K2,
            'vti': vti, 'vtI': vtI, 'rho_i': rho_i, 'rho_ip': rho_ip,
            'LTi_inv': LTi_inv, 'Lpi_inv': Lpi_inv, 'Lpimp_inv': Lpimp_inv,
            'Z_i': Z_i, 'Z_imp': Z_imp, 'T_imp': T_imp, 'Ti': Ti,
        }
        return self._neo_cache

    def _vpol_neo_i(self):
        """Main-ion (D) poloidal velocity per Eq. 33:

            v_θ,i = (1/2)·v_thi·ρ_i · K₁ · (1/L_Ti) · B_tot·B_t / ⟨B²⟩
        """
        neo = self._compute_neo()
        if neo is None:
            return None
        with np.errstate(invalid='ignore', divide='ignore'):
            geom = self.B_tot * self.Bt / self.B2_FSA
            v_pol_ms = (0.5 * neo['vti'] * neo['rho_i'] * neo['K1']
                        * neo['LTi_inv'] * geom)
        return self._output(v_pol_ms / 1e3, r'$v_{\theta,i}^\mathrm{neo}$', 'km/s')

    def _vpol_neo_imp(self):
        """Impurity (C⁶⁺) poloidal velocity per Eq. 34:

            v_θ,I = (1/2)·v_thi·ρ_i · [(K₁ + (3/2)K₂)/L_Ti − 1/L_pi
                                       + (Z_i/Z_I)·(T_I/T_i)·1/L_pI] · B_tot·B_t/⟨B²⟩
        """
        neo = self._compute_neo()
        if neo is None:
            return None
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = (neo['Z_i'] / neo['Z_imp']) * (neo['T_imp'] / neo['Ti'])
            bracket = ((neo['K1'] + 1.5 * neo['K2']) * neo['LTi_inv']
                       - neo['Lpi_inv']
                       + ratio * neo['Lpimp_inv'])
            geom = self.B_tot * self.Bt / self.B2_FSA
            v_pol_ms = 0.5 * neo['vti'] * neo['rho_i'] * bracket * geom
        return self._output(v_pol_ms / 1e3, r'$v_{\theta,\mathrm{imp}}^\mathrm{neo}$', 'km/s')

    def _omega_neo_diff(self):
        """Neoclassical impurity-main-ion toroidal rotation difference (Eq. 40):

            ω_I − ω_i = (3/4)·K₂·(v_thi/R)·ρ_ip · (1/L_Ti)

        Useful for converting CES-measured ω_C6 → main-ion ω_D in TRANSP input.
        Output in krad/s.
        """
        neo = self._compute_neo()
        if neo is None:
            return None
        with np.errstate(invalid='ignore', divide='ignore'):
            R_safe = np.where(self.R > 0, self.R, np.nan)
            domega = (0.75 * neo['K2'] * (neo['vti'] / R_safe)
                      * neo['rho_ip'] * neo['LTi_inv'])   # rad/s
            domega_krads = domega / 1e3
        return self._output(domega_krads, r'$\omega_I - \omega_i$', 'krad/s')

    # ------------------------------------------------------------------
    # Sauter neoclassical conductivity / bootstrap current
    # References:
    #   O. Sauter et al., Phys. Plasmas 6 (1999) 2834
    #   O. Sauter et al., Phys. Plasmas 9 (2002) 5140 [erratum]
    # Implementation follows the original 1999/2002 formulation (not Redl 2021),
    # ------------------------------------------------------------------

    def _compute_sauter(self):
        """Compute and cache the full Sauter coefficient set (L31, L32, L34, α,
        σ_neo, F33, etc.). Called once per (zeff, bi, efit) configuration.
        Returns None if required inputs are missing."""
        if self._sauter_cache is not None:
            return self._sauter_cache

        Te_TT = self._T_T('Te')
        ne_TT = self._T_T('ne')
        Ti_TT = self._T_T('Ti')
        if Te_TT is None or ne_TT is None or Ti_TT is None:
            return None

        Te = Te_TT * 1e3          # eV
        ne = ne_TT * 1e19         # m^-3
        Ti = Ti_TT * 1e3          # eV
        Z = self.zeff             # scalar
        fT = self.fT              # (n_t, n_psi)
        q = np.abs(self.q)        # |q|
        # Use R_axis for the Sauter collisionality (eqn 18b/18c convention)
        R0 = self.rmaxis[:, None]
        eps = np.maximum(self.eps, 1e-6)

        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            # Coulomb log (NRL formula, more accurate at low Te)
            lnL_e = (23.5 - np.log(np.sqrt(ne / 1e6) * Te**(-5.0 / 4.0))
                     - np.sqrt(1e-5 + (np.log(Te) - 2.0)**2 / 16.0))
            lnL_e = np.clip(lnL_e, 10.0, 25.0)
            # For ions (Sauter 18e): Z=Zeff approximation
            lnL_i = 30.0 - np.log(Z**3 * np.sqrt(ne / Z) / (Ti**1.5))
            lnL_i = np.clip(lnL_i, 10.0, 25.0)

            # Collisionality (Sauter 18b, 18c)
            nuestar = (6.921e-18 * q * R0 * ne * Z * lnL_e
                       / (Te**2 * eps**1.5))
            nuistar = (4.90e-18 * q * R0 * (ne / Z) * Z**4 * lnL_i
                       / (Ti**2 * eps**1.5))

            # Spitzer conductivity (Sauter 18a)
            NZ = 0.58 + 0.74 / (0.76 + Z)
            sigma_sp = 1.9012e4 * Te**1.5 / (Z * NZ * lnL_e)   # 1/(Ω·m)

            # Neoclassical conductivity (Sauter 13a, 13b)
            f33teff = fT / (1.0 + (0.55 - 0.1 * fT) * np.sqrt(nuestar)
                            + 0.45 * (1.0 - fT) * nuestar * Z**(-1.5))
            F33 = (1.0 - (1.0 + 0.36 / Z) * f33teff
                   + 0.59 / Z * f33teff**2
                   - 0.23 / Z * f33teff**3)
            sigma_neo = sigma_sp * F33

            # L31 (Sauter 14)
            f31teff = fT / (1.0 + (1.0 - 0.1 * fT) * np.sqrt(nuestar)
                            + 0.5 * (1.0 - fT) * nuestar / Z)

            def _F31(X):
                return ((1.0 + 1.4 / (Z + 1.0)) * X
                        - 1.9 / (Z + 1.0) * X**2
                        + 0.3 / (Z + 1.0) * X**3
                        + 0.2 / (Z + 1.0) * X**4)

            L31 = _F31(f31teff)

            # L32 (Sauter 15)
            f32eeteff = fT / (1.0 + 0.26 * (1.0 - fT) * np.sqrt(nuestar)
                              + 0.18 * (1.0 - 0.37 * fT) * nuestar / np.sqrt(Z))
            f32eiteff = fT / (1.0 + (1.0 + 0.6 * fT) * np.sqrt(nuestar)
                              + 0.85 * (1.0 - 0.37 * fT) * nuestar * (1.0 + Z))
            X = f32eeteff
            F32ee = ((0.05 + 0.62 * Z) / (Z * (1.0 + 0.44 * Z)) * (X - X**4)
                     + 1.0 / (1.0 + 0.22 * Z) * (X**2 - X**4 - 1.2 * (X**3 - X**4))
                     + 1.2 / (1.0 + 0.5 * Z) * X**4)
            Y = f32eiteff
            F32ei = (-(0.56 + 1.93 * Z) / (Z * (1.0 + 0.44 * Z)) * (Y - Y**4)
                     + 4.95 / (1.0 + 2.48 * Z) * (Y**2 - Y**4 - 0.55 * (Y**3 - Y**4))
                     - 1.2 / (1.0 + 0.5 * Z) * Y**4)
            L32 = F32ee + F32ei

            # L34 (Sauter 16)
            f34teff = fT / (1.0 + (1.0 - 0.1 * fT) * np.sqrt(nuestar)
                            + 0.5 * (1.0 - 0.5 * fT) * nuestar / Z)
            L34 = _F31(f34teff)

            # α (Sauter 17a + 2002 erratum)
            alpha0 = -1.17 * (1.0 - fT) / (1.0 - 0.22 * fT - 0.19 * fT**2)
            alpha = ((alpha0 + 0.25 * (1.0 - fT**2) * np.sqrt(nuistar))
                     / (1.0 + 0.5 * np.sqrt(nuistar))
                     + 0.315 * nuistar**2 * fT**6) / (1.0 + 0.15 * nuistar**2 * fT**6)

        self._sauter_cache = {
            'sigma_sp': sigma_sp, 'sigma_neo': sigma_neo, 'F33': F33,
            'L31': L31, 'L32': L32, 'L34': L34, 'alpha': alpha,
            'nuestar': nuestar, 'nuistar': nuistar,
            'lnL_e': lnL_e, 'lnL_i': lnL_i,
        }
        return self._sauter_cache

    def _sigma_neo(self):
        s = self._compute_sauter()
        if s is None:
            return None
        return self._output(s['sigma_neo'], r'$\sigma_\mathrm{neo}$', '1/(Ω·m)')

    def _eta_neo(self):
        s = self._compute_sauter()
        if s is None:
            return None
        with np.errstate(invalid='ignore', divide='ignore'):
            eta = 1.0 / s['sigma_neo']
        return self._output(eta, r'$\eta_\mathrm{neo}$', 'Ω·m')

    def _grad_psin(self, y_TT):
        """∂y/∂ψ_N at each (time, ψ_N). y_TT shape (n_t, n_psi)."""
        n_t, n_psi = y_TT.shape
        out = np.full_like(y_TT, np.nan, dtype=float)
        for it in range(n_t):
            y_row = y_TT[it]
            valid = np.isfinite(y_row)
            if valid.sum() < 3:
                continue
            psin_v = self.psin[valid]
            y_v = y_row[valid]
            dy_dpsi = np.gradient(y_v, psin_v)
            out_row = np.full(n_psi, np.nan)
            out_row[valid] = dy_dpsi
            out[it] = out_row
        return out

    def _j_BS(self):
        """Bootstrap current density ⟨j_BS·B⟩/B₀ from Sauter 1999 jB form
        (the smoother of the two expressions in the paper's conclusion).

            ⟨j_BS·B⟩ = −I_ψ · p · [L31·∂lnn_e/∂ψ + R_pe·(L31+L32)·∂lnT_e/∂ψ
                                    + (1−R_pe)·(L31+α·L34)·∂lnT_i/∂ψ] · sgn(q)

        We normalize by the vacuum BT₀ (self.bcentr) so the output is in A/m²
        and directly comparable with j_ohmic = σ_neo · E_∥.
        """
        s = self._compute_sauter()
        if s is None:
            return None
        Te_TT = self._T_T('Te'); ne_TT = self._T_T('ne'); Ti_TT = self._T_T('Ti')
        if Te_TT is None or ne_TT is None or Ti_TT is None:
            return None
        # SI conversions
        Te_eV = Te_TT * 1e3
        Ti_eV = Ti_TT * 1e3
        ne_m3 = ne_TT * 1e19
        # Pressure (Pa)
        p_e = ne_m3 * Te_eV * E_CHARGE
        # ion density: quasineutrality split (n_i / n_e = (Z_imp - Zeff)/(Z_imp - 1))
        Z_imp = IMP_SPECIES[self.imp][1]
        ni_over_ne = max((Z_imp - self.zeff) / (Z_imp - 1.0), 0.0)
        p_i = ne_m3 * ni_over_ne * Ti_eV * E_CHARGE
        p_tot = p_e + p_i
        with np.errstate(invalid='ignore', divide='ignore'):
            R_pe = p_e / p_tot
            # Gradients in ψ_N space, converted to ψ-space via 1/Δψ
            # (we use d/dψ = (1/Δψ) · d/dψ_N, so all terms get the same factor)
            d_psi = self.delta_psi[:, None]
            dlnTe_dpsi = self._grad_psin(Te_eV) / Te_eV / d_psi
            dlnTi_dpsi = self._grad_psin(Ti_eV) / Ti_eV / d_psi
            dlnne_dpsi = self._grad_psin(ne_m3) / ne_m3 / d_psi
            # I_ψ ≡ F = fpol (T·m)
            I_psi = self.fpol
            sign_q = np.sign(self.q)
            bracket = (s['L31'] * dlnne_dpsi
                       + R_pe * (s['L31'] + s['L32']) * dlnTe_dpsi
                       + (1.0 - R_pe) * (s['L31'] + s['alpha'] * s['L34']) * dlnTi_dpsi)
            # ⟨j_BS·B⟩ has units N/m³ = T·A/m². Divide by BT₀ to get A/m².
            jB = -I_psi * p_tot * bracket * sign_q
            B0 = self.bcentr[:, None]
            j_BS = jB / B0
        return self._output(j_BS, r'$\langle j_\mathrm{BS}\cdot B\rangle/B_0$', 'A/m²')

    # ------------------------------------------------------------------
    # Public dispatcher
    # ------------------------------------------------------------------

    _DISPATCH = {
        'p_e':        '_pressure_e',
        'p_i':        '_pressure_i',
        'p_tot':      '_pressure_tot',
        'dTe_dR':     None,  # special: use lambda
        'dTi_dR':     None,
        'dne_dR':     '_dne_dR',
        'R_over_LTe': None,
        'R_over_LTi': None,
        'R_over_Lne': '_R_over_Ln',
        'beta_e':     None,
        'beta_i':     None,
        'beta_tot':   None,
        'nu_star_e':  None,
        'nu_star_i':  None,
        'omega_pe':   '_omega_pe',
        'omega_ci':   '_omega_ci',
        'c_s':        '_c_s',
        'rho_i':      '_rho_i',
        'E_r_vp0':       '_E_r_vp0',
        'E_r_neo':       '_E_r_neo',
        'vpol_neo_i':    '_vpol_neo_i',
        'vpol_neo_imp':  '_vpol_neo_imp',
        'omega_neo_diff':'_omega_neo_diff',
        'sigma_neo':     '_sigma_neo',
        'eta_neo':       '_eta_neo',
        'j_BS':          '_j_BS',
    }

    def compute(self, key):
        """Compute one derived quantity by key, with caching."""
        if key in self._cache:
            return self._cache[key]
        try:
            if key == 'dTe_dR':
                out = self._dT_dR('Te', r'$dT_e/dR$')
            elif key == 'dTi_dR':
                out = self._dT_dR('Ti', r'$dT_i/dR$')
            elif key == 'R_over_LTe':
                out = self._R_over_LT('Te', r'$R/L_{T_e}$')
            elif key == 'R_over_LTi':
                out = self._R_over_LT('Ti', r'$R/L_{T_i}$')
            elif key in ('beta_e', 'beta_i', 'beta_tot'):
                out = self._beta(key.split('_')[1])
            elif key == 'nu_star_e':
                out = self._nu_star('e')
            elif key == 'nu_star_i':
                out = self._nu_star('i')
            else:
                m = self._DISPATCH.get(key)
                if m is None:
                    raise ValueError(f"Unknown derived quantity: {key}")
                out = getattr(self, m)()
        except Exception as e:
            print(f"[Derived] Failed {key}: {e}")
            out = None
        self._cache[key] = out
        return out
