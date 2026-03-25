"""
Profile fitting functions for PRISM

Provides:
- FitResult dataclass for storing fit results
- Fitting functions: mtanh_prof, tanh_prof (ptanh), eped_prof, spline
- fit_profile driver function
- Default parameter sets for Ti, Te, ne, vT
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any
import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit


@dataclass
class FitResult:
    """Container for profile fitting results

    Attributes:
        x_fit: Fine x-grid for fitted curve (typically 200 points)
        y_fit: Fitted values on fine grid
        y_fit_std: Standard deviation of fitted values (GPR only, None otherwise)
        params: Fitted parameter values dict {name: value}
        param_errors: Parameter uncertainties dict {name: error}
        func_name: Name of fitting function used
        chi_squared: Reduced chi-squared of fit
        success: Whether fit converged
        message: Status message
    """
    x_fit: np.ndarray
    y_fit: np.ndarray
    params: Dict[str, float]
    param_errors: Dict[str, float]
    func_name: str
    chi_squared: float = 0.0
    success: bool = True
    message: str = ""
    y_fit_std: Optional[np.ndarray] = None


# ===== Fitting functions =====

def _mtanh(x, b):
    """Modified tanh basis function (numerically stable)

    Algebraically equivalent to:
        ((1 + b*x) * exp(x) - exp(-x)) / (exp(x) + exp(-x))
    but rewritten using np.tanh to avoid exp overflow for large |x|.
    """
    t = np.tanh(x)
    return t + b * x * (1.0 + t) / 2.0


def mtanh_prof(x, a1, a2, a3, a4, a5, a6, a7, a8):
    """Modified tanh profile

    8-parameter modified tanh with exponential core modification.

    Parameters:
        a1: edge value (pedestal foot)
        a2: pedestal height
        a3: pedestal width
        a4: pedestal position (symmetry point)
        a5: slope parameter for mtanh
        a6: core peak value
        a7: core width
        a8: core exponent
    """
    y = (a2 - a1) / 2.0 * (_mtanh((a4 - x) / 2 / (a3 / 4 + 1.e-7), a5) + 1.0) + a1
    # Use abs(x) to avoid NaN from fractional power of negative numbers
    x_abs = np.abs(x)
    y = y + (a6 - y) * np.exp(-1.0 * (x_abs / (np.abs(a7) + 0.2)) ** np.abs(a8))
    return y


def tanh_prof(x, a1, a2, a3, a4, a5, a6, a7):
    """Standard tanh profile (ptanh)

    7-parameter tanh with polynomial core.

    Parameters:
        a1: edge/SOL value
        a2: pedestal height
        a3: linear core coefficient
        a4: quadratic core coefficient
        a5: cubic core coefficient
        a6: pedestal width
        a7: pedestal position
    """
    y = (a2 - a1) * (1.0 + a3 * x + a4 * x**2 + a5 * x**3) * \
        0.5 * (1.0 - np.tanh((x - a7) / a6 * 2.0)) + a1
    return y


def eped_prof(x, a1, a2, a3, a4, a5, a6):
    """EPED-style profile

    6-parameter profile with EPED pedestal shape.

    Parameters:
        a1: edge/SOL value
        a2: pedestal height
        a3: pedestal width
        a4: core height
        a5: core exponent 1
        a6: core exponent 2
    """
    y = a1
    y = y + a2 * (np.tanh(1) - np.tanh((x - 1.0 + 0.5 * a3) / a3 * 2.0))
    cf = 0.5 * (1.0 + np.tanh((a5 - 0.2) / 1.e-2)) * a5 + 0.1
    # Use abs to avoid fractional power of negative numbers
    yt = (np.abs(x) / (1.0 - a3 + 1.e-10)) ** cf
    cf2 = 0.5 * (1.0 + np.tanh((a6 - 0.2) / 1.e-2)) * a6 + 0.1
    npstep = 1.0 + np.tanh((1.0 - a3 - x) / 1.e-6)
    y = y + a4 * ((np.abs(1 - yt)) ** cf2) * 0.5 * npstep
    return y


# ===== Function registry =====

FIT_FUNCTIONS = {
    'mtanh': {
        'func': mtanh_prof, 'n_params': 8,
        'param_names': ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'a8'],
        'description': 'Modified tanh — 8 parameters. Pedestal + core shape with asymmetric tanh. Best for profiles with distinct pedestal and core peaking.',
        'param_descriptions': {
            'a1': 'Edge value (pedestal foot)',
            'a2': 'Pedestal height (top)',
            'a3': 'Pedestal width',
            'a4': 'Pedestal position (symmetry point)',
            'a5': 'Slope parameter for mtanh',
            'a6': 'Core peak value',
            'a7': 'Core width',
            'a8': 'Core exponent',
        },
        'formula_latex': (
            r'$y = a_1 + \frac{a_2 - a_1}{2}\left[\mathrm{mtanh}\!\left(\frac{a_4 - x}{a_3/2},\, a_5\right)'
            r' + 1\right] + (a_6 - y)\,\exp\!\left[-\left(\frac{|x|}{a_7}\right)^{a_8}\right]$'
            '\n'
            r'$\mathrm{where,}\ \ \mathrm{mtanh}(z,b) = \tanh(z) + \frac{bz(1+\tanh z)}{2}$'
        ),
    },
    'ptanh': {
        'func': tanh_prof, 'n_params': 7,
        'param_names': ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7'],
        'description': 'Polynomial tanh — 7 parameters. Pedestal tanh × polynomial core shape. Flexible core profile with polynomial expansion.',
        'param_descriptions': {
            'a1': 'Edge/SOL value',
            'a2': 'Pedestal height',
            'a3': 'Linear core coefficient',
            'a4': 'Quadratic core coefficient',
            'a5': 'Cubic core coefficient',
            'a6': 'Pedestal width',
            'a7': 'Pedestal position',
        },
        'formula_latex': (
            r'$y = a_1 + (a_2 - a_1)(1 + a_3 x + a_4 x^2 + a_5 x^3)'
            r'\times\frac{1}{2}\left[1 - \tanh\!\left(\frac{x - a_7}{a_6}\cdot 2\right)\right]$'
        ),
    },
    'EPED': {
        'func': eped_prof, 'n_params': 6,
        'param_names': ['a1', 'a2', 'a3', 'a4', 'a5', 'a6'],
        'description': 'EPED profile — 6 parameters. Pedestal tanh + power-law core. Used in EPED-type pedestal stability analysis.',
        'param_descriptions': {
            'a1': 'Edge/SOL value',
            'a2': 'Pedestal height',
            'a3': 'Pedestal width',
            'a4': 'Core height',
            'a5': 'Core exponent 1',
            'a6': 'Core exponent 2',
        },
        'formula_latex': (
            r'$y = a_1 + a_2\left[\tanh(1) - \tanh\!\left(\frac{x - 1 + a_3/2}{a_3}\cdot 2\right)\right]$'
            '\n'
            r'$\quad + a_4\left(1 - \left(\frac{|x|}{1 - a_3}\right)^{a_5}\right)^{a_6}\ \ (x < 1 - a_3)$'
        ),
    },
    'spline': {
        'func': None, 'n_params': 0, 'param_names': [],
        'description': 'Smoothing spline — automatic knot selection. Non-parametric. Smoothness controlled by data noise. No pedestal parameters extracted.',
        'param_descriptions': {},
        'formula_latex': '',
    },
    'RBF': {
        'func': None, 'n_params': 0, 'param_names': [],
        'description': 'Radial Basis Function (RBF) interpolation. Non-parametric. Places N Gaussian basis functions evenly spaced in the fitting range, then solves for weights via least-squares against all data points. Fewer bases = smoother, more bases = more detailed. No pedestal parameters.',
        'param_descriptions': {},
        'formula_latex': (
            r'$f(x) = \sum_i w_i\,\exp\!\left(-\left(\frac{\|x - c_i\|}{\epsilon}\right)^2\right) + p(x)$'
        ),
    },
    'GPR': {
        'func': None, 'n_params': 0, 'param_names': [],
        'description': 'Gaussian Process Regression (GPR) with uncertainty. Non-parametric Bayesian method. Provides mean prediction with confidence bands (±2σ). Uses squared-exponential kernel with automatic hyperparameter optimization (marginal likelihood maximization). No pedestal parameters.',
        'param_descriptions': {},
        'formula_latex': (
            r'$f(x) \sim \mathcal{GP}\!\left(0,\, k(x, x^\prime)\right)$'
            '\n'
            r'$k(x, x^\prime) = \sigma_f^2 \exp\!\left(-\frac{(x-x^\prime)^2}{2\ell^2}\right) + \sigma_n^2\delta_{xx^\prime}$'
        ),
    },
}


# ===== Default parameters =====
# Format: {param_name: (default_value, min_bound, max_bound)}

_DEFAULT_PARAMS = {
    'mtanh': {
        'Ti': {
            'a1': (0.15, 0.05, 0.3),
            'a2': (1.0, 0.0, np.inf),
            'a3': (0.05, 0.04, 0.15),
            'a4': (0.96, 0.85, 1.02),
            'a5': (1.0, 0.0, 2.0),
            'a6': (3.0, 0.0, 10.0),
            'a7': (0.1, 0.0, 0.5),
            'a8': (1.5, 1.01, 10.0),
        },
        'Te': {
            'a1': (0.15, 0.05, 0.3),
            'a2': (1.0, 0.0, np.inf),
            'a3': (0.05, 0.04, 0.15),
            'a4': (0.96, 0.85, 1.02),
            'a5': (1.0, 0.0, 2.0),
            'a6': (3.0, 0.0, 10.0),
            'a7': (0.1, 0.0, 0.5),
            'a8': (1.5, 1.01, 10.0),
        },
        'ne': {
            'a1': (0.5, 0.3, 1.0),
            'a2': (1.0, 0.0, np.inf),
            'a3': (0.05, 0.035, 0.1),
            'a4': (0.96, 0.85, 1.02),
            'a5': (1.0, 0.0, 2.0),
            'a6': (3.0, 0.0, 10.0),
            'a7': (0.1, 0.0, 0.5),
            'a8': (1.5, 1.01, 10.0),
        },
        'vT': {
            'a1': (50.0, 0.0, np.inf),
            'a2': (1.0, 0.0, np.inf),
            'a3': (0.05, 0.04, 0.15),
            'a4': (0.96, 0.85, 1.02),
            'a5': (1.0, 0.0, 2.0),
            'a6': (100.0, 0.0, np.inf),
            'a7': (0.1, 0.0, 0.5),
            'a8': (1.5, 1.01, 10.0),
        },
    },
    'ptanh': {
        'Ti': {
            'a1': (0.15, 0.05, 0.2),
            'a2': (1.0, 0.3, 20.0),
            'a3': (-0.1, -np.inf, 0.0),
            'a4': (1.0, -np.inf, np.inf),
            'a5': (1.0, -np.inf, np.inf),
            'a6': (0.05, 0.04, 0.15),
            'a7': (0.98, 0.88, 1.04),
        },
        'Te': {
            'a1': (0.15, 0.05, 0.2),
            'a2': (1.0, 0.3, 20.0),
            'a3': (-0.1, -np.inf, 0.0),
            'a4': (1.0, -np.inf, np.inf),
            'a5': (1.0, -np.inf, np.inf),
            'a6': (0.05, 0.04, 0.15),
            'a7': (0.98, 0.88, 1.04),
        },
        'ne': {
            'a1': (0.5, 0.2, 1.0),
            'a2': (1.0, 0.5, 10.0),
            'a3': (-0.1, -np.inf, 0.0),
            'a4': (1.0, -np.inf, np.inf),
            'a5': (1.0, -np.inf, np.inf),
            'a6': (0.05, 0.03, 0.1),
            'a7': (0.98, 0.88, 1.04),
        },
        'vT': {
            'a1': (50.0, 0.0, 100.0),
            'a2': (100.0, 0.0, 500.0),
            'a3': (-0.1, -np.inf, 0.0),
            'a4': (1.0, -np.inf, np.inf),
            'a5': (1.0, -np.inf, np.inf),
            'a6': (0.05, 0.04, 0.15),
            'a7': (0.98, 0.88, 1.04),
        },
    },
    'EPED': {
        'Ti': {
            'a1': (0.15, 0.05, 0.3),
            'a2': (0.5, 0.15, 5.0),
            'a3': (0.05, 0.03, 0.1),
            'a4': (3.0, 0.2, 10.0),
            'a5': (1.3, 1.01, 4.0),
            'a6': (2.0, 1.01, 4.0),
        },
        'Te': {
            'a1': (0.15, 0.05, 0.3),
            'a2': (0.5, 0.15, 5.0),
            'a3': (0.05, 0.03, 0.1),
            'a4': (3.0, 0.2, 10.0),
            'a5': (1.3, 1.01, 4.0),
            'a6': (2.0, 1.01, 4.0),
        },
        'ne': {
            'a1': (0.5, 0.3, 1.0),
            'a2': (1.0, 0.15, 5.0),
            'a3': (0.05, 0.03, 0.1),
            'a4': (3.0, 0.2, 10.0),
            'a5': (1.3, 1.01, 2.5),
            'a6': (2.0, 1.01, 2.5),
        },
        'vT': {
            'a1': (50.0, 0.0, 100.0),
            'a2': (100.0, 0.0, np.inf),
            'a3': (0.05, 0.04, 0.15),
            'a4': (200.0, 0.0, np.inf),
            'a5': (1.3, 1.01, 4.0),
            'a6': (2.0, 1.01, 4.0),
        },
    },
}


def get_default_params(func_name: str, param_type: str) -> Dict[str, Tuple[float, float, float]]:
    """Get default parameters for a fitting function and parameter type.

    Args:
        func_name: 'mtanh', 'ptanh', or 'EPED'
        param_type: 'Ti', 'Te', 'ne', or 'vT'

    Returns:
        Dict of {param_name: (default, min, max)}
    """
    if func_name not in _DEFAULT_PARAMS:
        return {}

    func_defaults = _DEFAULT_PARAMS[func_name]

    # Fall back to Ti defaults if specific type not found
    if param_type not in func_defaults:
        param_type = 'Ti'

    return func_defaults.get(param_type, {})


def fit_profile(x_data: np.ndarray, y_data: np.ndarray,
                func_name: str = 'mtanh',
                param_type: str = 'Ti',
                user_params: Optional[Dict[str, Tuple[float, float, float, bool]]] = None,
                sigma: Optional[np.ndarray] = None,
                x_fit_range: Tuple[float, float] = (0.0, 1.2),
                n_fit_points: int = 200,
                n_bases: Optional[int] = None) -> FitResult:
    """Fit a profile using the specified function.

    Args:
        x_data: x coordinates (e.g., psi_N, rho_pol)
        y_data: y values to fit
        func_name: Fitting function name ('mtanh', 'ptanh', 'EPED', 'spline')
        param_type: Parameter type for defaults ('Ti', 'Te', 'ne', 'vT')
        user_params: User-specified parameters {name: (value, min, max, fixed)}
                     If None, uses defaults.
        sigma: Uncertainties for weighted fitting (same length as y_data)
        x_fit_range: Range for fitted curve output
        n_fit_points: Number of points in fitted curve

    Returns:
        FitResult object
    """
    # Remove NaN/Inf values and limit to fitting range
    valid = np.isfinite(x_data) & np.isfinite(y_data)
    valid &= (x_data >= x_fit_range[0]) & (x_data <= x_fit_range[1])

    # Handle sigma: replace NaN/invalid with median of valid sigmas
    # so that data without error bars (e.g. ECE) still participates in fitting
    if sigma is not None:
        sig_valid = np.isfinite(sigma) & (sigma > 0)
        if sig_valid.any() and not sig_valid[valid].all():
            # Some valid data points lack sigma — fill with median
            median_sig = np.nanmedian(sigma[sig_valid & valid])
            sigma = sigma.copy()
            sigma[~sig_valid] = median_sig
        elif not sig_valid.any():
            sigma = None  # no valid sigma at all — use unweighted

    x = x_data[valid]
    y = y_data[valid]
    sig = sigma[valid] if sigma is not None else None

    if len(x) < 3:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name=func_name, success=False,
            message="Not enough valid data points for fitting"
        )

    # Sort by x
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]
    if sig is not None:
        sig = sig[sort_idx]

    x_fine = np.linspace(x_fit_range[0], x_fit_range[1], n_fit_points)

    # Non-parametric fitting
    if func_name == 'spline':
        return _fit_spline(x, y, sig, x_fine)
    if func_name == 'RBF':
        return _fit_rbf(x, y, sig, x_fine, n_bases=n_bases)
    if func_name == 'GPR':
        return _fit_gpr(x, y, sig, x_fine)

    # Parametric fitting
    func_info = FIT_FUNCTIONS.get(func_name)
    if func_info is None or func_info['func'] is None:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name=func_name, success=False,
            message=f"Unknown fitting function: {func_name}"
        )

    func = func_info['func']
    param_names = func_info['param_names']

    # Get initial values and bounds
    defaults = get_default_params(func_name, param_type)

    p0 = []
    lower_bounds = []
    upper_bounds = []
    fixed_params = {}

    for pname in param_names:
        if user_params and pname in user_params:
            val, lb, ub, fixed = user_params[pname]
            if fixed:
                fixed_params[pname] = val
            else:
                p0.append(val)
                lower_bounds.append(lb)
                upper_bounds.append(ub)
        elif pname in defaults:
            val, lb, ub = defaults[pname]
            p0.append(val)
            lower_bounds.append(lb)
            upper_bounds.append(ub)
        else:
            p0.append(1.0)
            lower_bounds.append(-np.inf)
            upper_bounds.append(np.inf)

    # Build wrapper function that handles fixed params
    free_names = [n for n in param_names if n not in fixed_params]

    def fit_func(x_val, *free_vals):
        all_vals = []
        free_idx = 0
        for pname in param_names:
            if pname in fixed_params:
                all_vals.append(fixed_params[pname])
            else:
                all_vals.append(free_vals[free_idx])
                free_idx += 1
        return func(x_val, *all_vals)

    try:
        popt, pcov = curve_fit(
            fit_func, x, y,
            p0=p0,
            bounds=(lower_bounds, upper_bounds),
            sigma=sig,
            absolute_sigma=True if sig is not None else False,
            maxfev=10000
        )

        perr = np.sqrt(np.diag(pcov))

        # Build full param dicts
        params = {}
        param_errors = {}
        free_idx = 0
        for pname in param_names:
            if pname in fixed_params:
                params[pname] = fixed_params[pname]
                param_errors[pname] = 0.0
            else:
                params[pname] = popt[free_idx]
                param_errors[pname] = perr[free_idx]
                free_idx += 1

        # Evaluate on fine grid
        y_fine = func(x_fine, *[params[n] for n in param_names])

        # Compute reduced chi-squared
        y_pred = func(x, *[params[n] for n in param_names])
        residuals = y - y_pred
        if sig is not None:
            chi2 = np.sum((residuals / sig) ** 2)
        else:
            chi2 = np.sum(residuals ** 2)
        n_free = len(free_names)
        dof = len(x) - n_free
        reduced_chi2 = chi2 / dof if dof > 0 else chi2

        return FitResult(
            x_fit=x_fine, y_fit=y_fine,
            params=params, param_errors=param_errors,
            func_name=func_name,
            chi_squared=reduced_chi2,
            success=True,
            message="Fit converged"
        )

    except Exception as e:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name=func_name, success=False,
            message=f"Fit failed: {str(e)}"
        )


def _fit_spline(x, y, sigma, x_fine):
    """Fit using smoothing spline"""
    try:
        # Use error as weight (1/sigma)
        weights = 1.0 / sigma if sigma is not None else None

        # Determine smoothing factor based on data
        s = len(x) if sigma is not None else None

        spline = UnivariateSpline(x, y, w=weights, s=s, k=3)
        y_fine = spline(x_fine)

        # Chi-squared estimate
        y_pred = spline(x)
        residuals = y - y_pred
        if sigma is not None:
            chi2 = np.sum((residuals / sigma) ** 2)
        else:
            chi2 = np.sum(residuals ** 2)
        dof = max(1, len(x) - spline.get_knots().size)
        reduced_chi2 = chi2 / dof

        return FitResult(
            x_fit=x_fine, y_fit=y_fine,
            params={'knots': spline.get_knots().size},
            param_errors={},
            func_name='spline',
            chi_squared=reduced_chi2,
            success=True,
            message=f"Spline fit with {spline.get_knots().size} knots"
        )

    except Exception as e:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name='spline', success=False,
            message=f"Spline fit failed: {str(e)}"
        )


def _fit_rbf(x, y, sigma, x_fine, n_bases=None):
    """Fit using Radial Basis Function interpolation.

    Places n_bases Gaussian basis functions evenly spaced in the data range,
    then solves for weights via least-squares against ALL data points.

    Args:
        n_bases: Number of basis function centers (default 5).
                 More bases → more detailed, fewer → smoother.
    """
    try:
        if n_bases is None or n_bases <= 0:
            n_bases = 5
        n_bases = min(n_bases, len(x))

        # Place centers evenly in data range
        centers = np.linspace(x.min(), x.max(), n_bases)

        # Epsilon (width): scale with center spacing so bases overlap
        spacing = (x.max() - x.min()) / max(n_bases - 1, 1)
        epsilon = spacing  # each basis covers ~1 spacing

        # Build design matrix: Phi[i,j] = exp(-(x_i - c_j)^2 / epsilon^2)
        Phi = np.exp(-((x[:, None] - centers[None, :]) / epsilon) ** 2)

        # Weighted least squares: minimize sum((y - Phi @ w)^2 / sigma^2)
        if sigma is not None:
            W = 1.0 / sigma
            Phi_w = Phi * W[:, None]
            y_w = y * W
        else:
            Phi_w = Phi
            y_w = y

        # Solve with Tikhonov regularization (small lambda for stability)
        lam = 1e-6
        ATA = Phi_w.T @ Phi_w + lam * np.eye(n_bases)
        ATy = Phi_w.T @ y_w
        weights = np.linalg.solve(ATA, ATy)

        # Evaluate on fine grid
        Phi_fine = np.exp(-((x_fine[:, None] - centers[None, :]) / epsilon) ** 2)
        y_fine = Phi_fine @ weights

        # Chi-squared
        y_pred = Phi @ weights
        residuals = y - y_pred
        if sigma is not None:
            chi2 = np.sum((residuals / sigma) ** 2)
        else:
            chi2 = np.sum(residuals ** 2)
        dof = max(1, len(x) - n_bases)
        reduced_chi2 = chi2 / dof

        return FitResult(
            x_fit=x_fine, y_fit=y_fine,
            params={'epsilon': float(epsilon), 'n_bases': n_bases},
            param_errors={},
            func_name='RBF',
            chi_squared=reduced_chi2,
            success=True,
            message=f"RBF fit with {n_bases} bases, epsilon={epsilon:.4f}"
        )

    except Exception as e:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name='RBF', success=False,
            message=f"RBF fit failed: {str(e)}"
        )


def _fit_gpr(x, y, sigma, x_fine):
    """Fit using Gaussian Process Regression with uncertainty bands"""
    try:
        from sklearn.gaussian_process import GaussianProcessRegressor
        from sklearn.gaussian_process.kernels import (
            ConstantKernel, RBF as GPR_RBF, WhiteKernel
        )

        x_2d = x.reshape(-1, 1)
        x_fine_2d = x_fine.reshape(-1, 1)

        # Normalize y for numerical stability
        y_mean = np.mean(y)
        y_std = np.std(y)
        if y_std < 1e-10:
            y_std = 1.0
        y_norm = (y - y_mean) / y_std

        # Estimate noise level from sigma
        if sigma is not None:
            alpha = (sigma / y_std) ** 2
        else:
            alpha = 1e-4

        # Kernel: constant * RBF + white noise
        length_scale = (x.max() - x.min()) / 5.0
        kernel = (ConstantKernel(1.0, (0.01, 100.0))
                  * GPR_RBF(length_scale, (length_scale * 0.1, length_scale * 10.0))
                  + WhiteKernel(0.01, (1e-6, 1.0)))

        gpr = GaussianProcessRegressor(
            kernel=kernel, alpha=alpha,
            n_restarts_optimizer=5, normalize_y=False
        )
        gpr.fit(x_2d, y_norm)

        y_pred_norm, y_std_norm = gpr.predict(x_fine_2d, return_std=True)

        # De-normalize
        y_fine = y_pred_norm * y_std + y_mean
        y_fine_std = y_std_norm * y_std

        # Chi-squared estimate
        y_pred_data = gpr.predict(x_2d) * y_std + y_mean
        residuals = y - y_pred_data
        if sigma is not None:
            chi2 = np.sum((residuals / sigma) ** 2)
        else:
            chi2 = np.sum(residuals ** 2)
        dof = max(1, len(x) - 3)
        reduced_chi2 = chi2 / dof

        # Extract optimized kernel parameters
        opt_kernel = gpr.kernel_
        params_dict = {}
        for pname, pval in opt_kernel.get_params().items():
            if isinstance(pval, (int, float)):
                params_dict[pname] = float(pval)

        return FitResult(
            x_fit=x_fine, y_fit=y_fine,
            params=params_dict,
            param_errors={},
            func_name='GPR',
            chi_squared=reduced_chi2,
            success=True,
            message=f"GPR fit with optimized kernel",
            y_fit_std=y_fine_std
        )

    except ImportError:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name='GPR', success=False,
            message="GPR requires scikit-learn: pip install scikit-learn"
        )
    except Exception as e:
        return FitResult(
            x_fit=np.array([]), y_fit=np.array([]),
            params={}, param_errors={},
            func_name='GPR', success=False,
            message=f"GPR fit failed: {str(e)}"
        )
