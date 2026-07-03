"""
Compute orchestrators for the PRISM batch SDK.

These wrap the Qt-free compute cores (core/nmode.py) with the parameter handling
that used to live in the GUI tab's _run_calculation / _save_data, returning a
typed result with full provenance metadata. No PySide6 / matplotlib here.
"""

from datetime import datetime, timezone

import numpy as np

from core.nmode import (
    NModeConfig,
    load_mirnov_data,
    calculate_fft,
    calculate_mode_numbers,
    compute_amp_evolution,
)
from batch.results import NModeResult, IRVBResult

try:
    from config.app_config import VERSION
except Exception:  # pragma: no cover - app_config should always import (Qt-free)
    VERSION = "unknown"


_SIGN_STR = {0: "abs", 1: "pos", -1: "neg", 2: "all"}


def _filter_mode(mode, msign):
    """Sign-filter the mode-number map exactly as NModeSpectrumTab._save_data."""
    if msign == 1:        # positive modes only
        return np.where(mode > 0, mode, 0)
    if msign == -1:       # negative modes only
        return np.where(mode < 0, mode, 0)
    if msign == 0:        # absolute value
        return np.abs(mode)
    return mode.copy()    # msign == 2 (all): untouched


def _slice_amp(amp_evolution, msign, nmodes):
    """Slice an amplitude-evolution array to the sign-selected mode rows.

    Rows are [n=+1..+nmodes, n=-1..-nmodes]; pos keeps the first half, neg the
    second, abs/all keep both. Applied identically to the Max and Sum traces so
    amplitude_max[i] / amplitude_sum[i] both map to n_modes_list[i].
    """
    amp_save = amp_evolution.copy()
    if msign == 1:        # positive modes only
        return amp_save[:nmodes, :]
    if msign == -1:       # negative modes only
        return amp_save[nmodes:, :]
    return amp_save       # abs / all keep both halves


def _n_modes_list(msign, nmodes):
    """Mode-number list stored in metadata, ORDERED to match the amplitude rows.

    amp_evolution rows are physically [n=+1..+nmodes, n=-1..-nmodes] (ascending |n|
    within each sign half), so the list must follow that same order for the saved-NPZ
    contract amplitude[i] <-> n_modes_list[i] to hold. Mirrors the GUI _save_data fix.
    """
    if msign == -1:
        return list(range(-1, -nmodes - 1, -1))                       # [-1, -2, ..., -nmodes]
    if msign == 2:
        return list(range(1, nmodes + 1)) + list(range(-1, -nmodes - 1, -1))
    return list(range(1, nmodes + 1))   # abs / pos


def compute_nmode(spec, *, mds_server=None):
    """Run the full n-mode pipeline for one NModeJobSpec and return an NModeResult.

    load_mirnov_data -> calculate_fft -> calculate_mode_numbers, then sign-filter
    to match the GUI's saved arrays and attach provenance metadata.
    """
    if mds_server is None:
        mds_server = NModeConfig.MDS_SERVER

    # tmin/tmax = None means "full shot": load wide (like the GUI's -1..100) and
    # let the actual data range decide. Otherwise load the requested window.
    load_tmin = -1.0 if spec.tmin is None else spec.tmin
    load_tmax = 100.0 if spec.tmax is None else spec.tmax
    mirnov = load_mirnov_data(spec.shot, load_tmin, load_tmax, mds_server=mds_server)

    # Effective window, clamped to the loaded data range exactly as the GUI does
    # (NModeSpectrumTab._run_calculation) so time_range metadata is accurate. For a
    # None bound: full-shot tmin starts at max(0, data_start) and tmax = data_end.
    # Arrays are unaffected either way (calculate_fft clips via searchsorted).
    data_tmin = float(mirnov.time[0])
    data_tmax = float(mirnov.time[-1])
    eff_tmin = max(0.0, data_tmin) if spec.tmin is None else max(float(spec.tmin), data_tmin)
    eff_tmax = data_tmax if spec.tmax is None else min(float(spec.tmax), data_tmax)
    fft = calculate_fft(
        mirnov, spec.t_interval, spec.integrate, spec.detrend, eff_tmin, eff_tmax
    )
    # The mode map is independent of the amplitude reduction, so compute it once
    # ('max' trace) and derive the 'sum' trace from the same mode map. Both are
    # saved (matching the GUI's NPZ).
    mode, freq, amp_evolution_max = calculate_mode_numbers(
        fft, spec.fmin, spec.fmax, spec.tol, spec.nmodes, spec.frac,
        spec.msign, spec.integrate, 'max',
    )
    amp_evolution_sum = compute_amp_evolution(
        fft.amp, freq, mode, spec.fmin, spec.fmax, spec.nmodes,
        spec.integrate, 'sum',
    )

    mode_save = _filter_mode(mode, spec.msign)
    amp_max_save = _slice_amp(amp_evolution_max, spec.msign, spec.nmodes)
    amp_sum_save = _slice_amp(amp_evolution_sum, spec.msign, spec.nmodes)

    meta = {
        "shot": spec.shot,
        "subsystem": spec.SUBSYSTEM,
        "time_range": [eff_tmin, eff_tmax],
        "freq_range": [spec.fmin, spec.fmax],
        "time_interval": spec.t_interval,
        "n_modes": spec.nmodes,
        "n_modes_list": _n_modes_list(spec.msign, spec.nmodes),
        "tolerance": spec.tol,
        "fraction": spec.frac,
        "sign": _SIGN_STR.get(spec.msign, "pos"),
        "integrate": bool(spec.integrate),
        "detrend": bool(spec.detrend),
        "prism_version": VERSION,
        "created_iso": datetime.now(timezone.utc).isoformat(),
        "mds_server": mds_server,
        "params_json": spec.canonical(),
    }

    # Keep the GUI's native dtypes so the .npz is byte-for-byte the same shape as
    # what NModeSpectrumTab._save_data writes (time/freq float32, mode int, amp float64).
    return NModeResult(
        shot=spec.shot,
        time=fft.time,
        frequency=freq,
        mode_spectrum=mode_save,
        amplitude=amp_max_save,
        amplitude_max=amp_max_save,
        amplitude_sum=amp_sum_save,
        meta=meta,
    )


def compute_irvb(spec, *, mds_server=None):
    """Run the full IRVB pipeline for one IRVBJobSpec and return an IRVBResult.

    Mirrors NModeSpectrumTab / IRVBTab exactly (Qt-free): load the IRVB
    reconstruction, load 2D EFIT, slice IRVB frames to the EFIT time range, then
    for every frame compute the regional Prad and resample psi_N onto the IRVB
    grid. `mds_server` (if given) overrides the EFIT MDS+ server.
    """
    from scipy.interpolate import RectBivariateSpline
    from config.app_config import AppConfig
    from data_loaders.irvb_loader import IRVBLoader
    from data_loaders.efit_loader import EFITLoader

    cfg = AppConfig()
    if mds_server:
        cfg.MDS_IP = mds_server

    irvb = IRVBLoader(cfg).load_data(spec.shot)
    efit2d = EFITLoader(cfg).load_efit_2d(spec.shot, spec.efit_tree)

    # Slice IRVB to the valid EFIT time range (0 .. EFIT end), like the GUI's
    # _slice_by_efit_time.
    m = (irvb.time >= 0) & (irvb.time <= efit2d.time[-1])
    if not np.any(m):
        raise RuntimeError(f"shot {spec.shot}: no IRVB frames within the EFIT time range")
    time = irvb.time[m]
    recon = irvb.recon[m]
    ptot = np.asarray(irvb.ptot)[:len(irvb.time)][m]
    n_frames = len(time)

    x_grid, y_grid = irvb.x_grid, irvb.y_grid
    X, _ = np.meshgrid(x_grid, y_grid)
    vol_factor = 2 * np.pi * X * (x_grid[1] - x_grid[0]) * (y_grid[1] - y_grid[0])

    psi_boundaries = list(spec.psi_boundaries)
    edges = [0] + psi_boundaries + [np.inf]
    n_regions = len(edges) - 1
    region_prad = np.zeros((n_regions, n_frames))

    max_bdry = efit2d.bdry_r.shape[1]
    efit_bdry_r = np.zeros((n_frames, max_bdry))
    efit_bdry_z = np.zeros((n_frames, max_bdry))
    efit_nbdry = np.zeros(n_frames, dtype=int)
    efit_rmaxis = np.zeros(n_frames)
    efit_zmaxis = np.zeros(n_frames)
    efit_psi_n = np.zeros((n_frames, len(efit2d.z_grid), len(efit2d.r_grid)))
    psi_n = np.zeros((n_frames, len(y_grid), len(x_grid)))
    psi_irvb_cache = {}

    for i, t in enumerate(time):
        efit_idx = efit2d.find_time_index(t)
        bdry_r, bdry_z = efit2d.get_boundary(efit_idx)
        nbdry = len(bdry_r)
        efit_bdry_r[i, :nbdry] = bdry_r
        efit_bdry_z[i, :nbdry] = bdry_z
        efit_nbdry[i] = nbdry
        efit_rmaxis[i], efit_zmaxis[i] = efit2d.get_magnetic_axis(efit_idx)
        psi_n_2d = efit2d.get_psi_normalized(efit_idx)
        efit_psi_n[i] = psi_n_2d

        if efit_idx not in psi_irvb_cache:
            spline = RectBivariateSpline(efit2d.z_grid, efit2d.r_grid, psi_n_2d)
            psi_irvb_cache[efit_idx] = spline(y_grid, x_grid)
        psi_n[i] = psi_irvb_cache[efit_idx]

        for r in range(n_regions):
            mask = (psi_n[i] >= edges[r]) & (psi_n[i] < edges[r + 1])
            region_prad[r, i] = np.sum(recon[i] * mask * vol_factor)

    region_labels = [
        f"psi_N > {edges[r]}" if edges[r + 1] == np.inf
        else f"{edges[r]} < psi_N < {edges[r + 1]}"
        for r in range(n_regions)
    ]

    # Limiter is optional in EFIT; store empty arrays (not None) so the saved NPZ
    # and the example reader stay well-typed.
    lim_r = efit2d.limiter_r if efit2d.limiter_r is not None else np.array([])
    lim_z = efit2d.limiter_z if efit2d.limiter_z is not None else np.array([])

    meta = {
        "shot": spec.shot,
        "subsystem": spec.SUBSYSTEM,
        "psi_boundaries": psi_boundaries,
        "efit_tree": spec.efit_tree,
        "time_range": [float(time[0]), float(time[-1])],
        "region_labels": region_labels,
        "prism_version": VERSION,
        "created_iso": datetime.now(timezone.utc).isoformat(),
        "mds_server": cfg.MDS_IP,
        "params_json": spec.canonical(),
    }

    return IRVBResult(
        shot=spec.shot,
        time=time, R=x_grid, Z=y_grid,
        prad_2d=recon, psi_n=psi_n, region_prad=region_prad, ptot=ptot,
        efit_r_grid=efit2d.r_grid, efit_z_grid=efit2d.z_grid, efit_psi_n=efit_psi_n,
        efit_bdry_r=efit_bdry_r, efit_bdry_z=efit_bdry_z, efit_nbdry=efit_nbdry,
        efit_rmaxis=efit_rmaxis, efit_zmaxis=efit_zmaxis,
        efit_limiter_r=lim_r, efit_limiter_z=lim_z,
        meta=meta,
    )
