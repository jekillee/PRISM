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
)
from batch.results import NModeResult

try:
    from config.app_config import VERSION
except Exception:  # pragma: no cover - app_config should always import (Qt-free)
    VERSION = "unknown"


_SIGN_STR = {0: "abs", 1: "pos", -1: "neg", 2: "all"}


def _apply_msign(mode, amp_evolution, msign, nmodes):
    """Sign-filter the mode/amplitude arrays exactly as NModeSpectrumTab._save_data.

    pos/neg zero out the other sign and slice the matching half of amp_evolution;
    abs takes |mode|; all keeps the full signed array and both amplitude halves.
    """
    mode_save = mode.copy()
    amp_save = amp_evolution.copy()
    if msign == 1:        # positive modes only
        mode_save = np.where(mode_save > 0, mode_save, 0)
        amp_save = amp_save[:nmodes, :]
    elif msign == -1:     # negative modes only
        mode_save = np.where(mode_save < 0, mode_save, 0)
        amp_save = amp_save[nmodes:, :]
    elif msign == 0:      # absolute value
        mode_save = np.abs(mode_save)
    # msign == 2 (all): leave both arrays untouched
    return mode_save, amp_save


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
    mode, freq, amp_evolution = calculate_mode_numbers(
        fft, spec.fmin, spec.fmax, spec.tol, spec.nmodes, spec.frac,
        spec.msign, spec.integrate,
    )

    mode_save, amp_save = _apply_msign(mode, amp_evolution, spec.msign, spec.nmodes)

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
        amplitude=amp_save,
        meta=meta,
    )
