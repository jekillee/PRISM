"""
N-mode (toroidal mode number) compute core for PRISM.

Qt-free numerical pipeline for Mirnov-coil toroidal mode-number analysis,
extracted verbatim from ui/tabs/diagnostics/spectral/nmode_spectrum_tab.py so
it can be imported and run headlessly (batch / CLI) without pulling in PySide6.

Public pipeline:
    mirnov = load_mirnov_data(shot, tmin, tmax)          # MDS+ fetch + time align
    fft    = calculate_fft(mirnov, t_interval, ...)      # windowed FFT
    mode, freq, amp_evolution = calculate_mode_numbers(fft, fmin, fmax, ...)

The only behavioural change vs. the GUI version is that the MDS+ server is now a
parameter (`mds_server`) instead of the hardcoded NModeConfig.MDS_SERVER, with the
same value as the default so existing callers are unaffected.
"""

import json
from pathlib import Path
import time as tclock
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
from scipy.signal import detrend
from scipy.fft import fft
# NOTE: MDSplus is imported lazily inside _load_single_channel so this module can be
# imported (and its FFT / mode-number math used) on a host without MDSplus.

# Load Mirnov coil configuration (core/ is one level below the version root)
_MIRNOV_CONFIG_PATH = Path(__file__).parents[1] / 'config' / 'mirnov_config.json'
with open(_MIRNOV_CONFIG_PATH, 'r') as f:
    _MIRNOV_CONFIG = json.load(f)


# =============================================================================
# Configuration
# =============================================================================

class NModeConfig:
    """N-Mode spectrum configuration"""
    MDS_SERVER = 'mdsr.kstar.kfe.re.kr:8005'

    DEFAULT_SHOT = 0
    DEFAULT_TMIN = 0.0
    DEFAULT_TMAX = 10.0
    DEFAULT_FMIN = 0
    DEFAULT_FMAX = 100
    DEFAULT_TINTERVAL = 0.01
    DEFAULT_NMODES = 5
    DEFAULT_TOL = 0.8
    DEFAULT_FRAC = 1e-2
    DEFAULT_INTEGRATE = False
    DEFAULT_DETREND = True
    DEFAULT_MSIGN = 1
    DEFAULT_PLOT_TYPE = 'contour'
    DEFAULT_NUMC = 50

    FIGURE_SIZE = (10, 6)


# =============================================================================
# Data Structures
# =============================================================================

class MirnovData:
    """Container for Mirnov coil data"""
    def __init__(self):
        self.time = None
        self.signals = None
        self.positions = None
        self.channel_names = None
        self.shot = None


class FFTResult:
    """Container for FFT results"""
    def __init__(self):
        self.time = None
        self.freq = None
        self.amp = None
        self.phase = None
        self.positions = None
        self.shot = None
        self.array = 'TOROIDAL'
        self.p = None
        self.dt = None


# =============================================================================
# Mirnov Coil Configuration
# =============================================================================

def get_shot_year(shot):
    """Get year from shot number using config file"""
    for year_str, shot_range in _MIRNOV_CONFIG['shot_year_ranges'].items():
        if year_str.startswith('_'):
            continue
        if not isinstance(shot_range, list) or len(shot_range) != 2:
            continue
        shot_min, shot_max = shot_range
        if shot_max is None:
            shot_max = float('inf')
        if shot_min <= shot <= shot_max:
            return int(year_str)
    return 2025  # Default to latest year


def get_mirnov_config(year):
    """Get Mirnov coil configuration for given year from config file"""
    year_str = str(year)
    yearly_signs = _MIRNOV_CONFIG['yearly_signs']

    # Get signs for the year, fallback to default
    signs = yearly_signs.get(year_str, yearly_signs['default'])

    # Get base channel names and angles
    names = _MIRNOV_CONFIG['channel_names'].copy()
    angles = _MIRNOV_CONFIG['toroidal_angles']

    # Apply channel name overrides for 2017+
    if year >= 2017:
        names[9] = 'MC1P03'

    # Filter to valid channels only (signs != 0)
    valid_idx = [i for i, s in enumerate(signs) if s != 0]

    return {
        'names': [names[i] for i in valid_idx],
        'angles': [angles[i] for i in valid_idx],
        'signs': [signs[i] for i in valid_idx]
    }


# =============================================================================
# Data Loading
# =============================================================================

def _load_single_channel(args):
    """Load single Mirnov channel"""
    from MDSplus import Connection   # lazy: keep core.nmode importable without MDSplus
    shot, name, angle, sign, tmin, tmax, mds_server = args

    try:
        mds = Connection(mds_server)
        mds.openTree('kstar', shot)
        mds.get(f'SetTimeContext({tmin},{tmax},)').data()
        data = mds.get(f'\\{name}').data()
        time_arr = mds.get(f'dim_of(\\{name})').data()
        mds.get('SetTimeContext(,,)').data()
        mds.closeTree('kstar', shot)

        return {
            'name': name, 'angle': angle, 'sign': sign,
            'data': np.array(data, dtype=np.float32) * sign,
            'time': np.array(time_arr, dtype=np.float32),
            'n_points': len(data), 'error': None
        }
    except Exception as e:
        return {
            'name': name, 'angle': angle, 'sign': sign,
            'data': None, 'time': None, 'n_points': 0, 'error': str(e)
        }


def load_mirnov_data(shot, tmin, tmax, *, mds_server=None):
    """Load Mirnov coil data from MDS+ with parallel processing"""
    if mds_server is None:
        mds_server = NModeConfig.MDS_SERVER

    year = get_shot_year(shot)
    config = get_mirnov_config(year)
    n_channels = len(config['names'])

    print(f"[n-Mode] Loading Mirnov data for shot #{shot} (year {year})...")
    print(f"[n-Mode]   Time range: {tmin:.2f} - {tmax:.2f} s")
    print(f"[n-Mode]   Channels ({n_channels}):")
    for i, (name, angle) in enumerate(zip(config['names'], config['angles']), 1):
        print(f"[n-Mode]     {i:>2d}. {name}  {angle:>6.2f}°")
    print(f"[n-Mode]   Threads: {n_channels}")

    args_list = [
        (shot, name, angle, sign, tmin, tmax, mds_server)
        for name, angle, sign in zip(config['names'], config['angles'], config['signs'])
    ]

    results = []
    completed = 0
    t0 = tclock.time()

    with ThreadPoolExecutor(max_workers=n_channels) as executor:
        futures = {executor.submit(_load_single_channel, args): args[1] for args in args_list}
        for future in as_completed(futures):
            results.append(future.result())
            completed += 1
            progress = completed / n_channels
            bar = '█' * int(40 * progress) + '░' * (40 - int(40 * progress))
            print(f'\r[n-Mode]   [{bar}] {completed}/{n_channels}', end='', flush=True)

    print(f'\n[n-Mode]   Completed in {tclock.time() - t0:.2f} sec')

    valid_results = []
    results.sort(key=lambda x: x['angle'])

    # Collect valid results
    for result in results:
        if result['error'] is None and result['n_points'] > 0:
            valid_results.append(result)

    if len(valid_results) == 0:
        raise RuntimeError("No valid Mirnov channels found")

    # Find common time range across all channels
    common_tmin = max(r['time'][0] for r in valid_results)
    common_tmax = min(r['time'][-1] for r in valid_results)

    if common_tmin >= common_tmax:
        raise RuntimeError(f"No overlapping time range found across channels")

    print(f"[n-Mode]   Common time range: {common_tmin:.4f} - {common_tmax:.4f} s")

    # Use the channel with finest sampling as reference
    ref_result = min(valid_results, key=lambda r: (r['time'][-1] - r['time'][0]) / r['n_points'])
    ref_time = ref_result['time']

    # Trim reference time to common range
    ref_mask = (ref_time >= common_tmin) & (ref_time <= common_tmax)
    time_arr = ref_time[ref_mask]
    n_points = len(time_arr)

    print(f"[n-Mode]   Reference channel: {ref_result['name']} ({n_points} points)")

    # Interpolate all channels to common time grid
    signals, valid_names, valid_angles = [], [], []

    for result in valid_results:
        ch_time = result['time']
        ch_data = result['data']

        # Check if interpolation is needed
        if len(ch_time) == n_points and np.allclose(ch_time[ref_mask] if len(ch_time) == len(ref_time) else ch_time, time_arr, rtol=1e-6):
            # Same time grid, just trim
            ch_mask = (ch_time >= common_tmin) & (ch_time <= common_tmax)
            interp_data = ch_data[ch_mask]
        else:
            # Interpolate to common time grid
            interp_data = np.interp(time_arr, ch_time, ch_data).astype(np.float32)

        signals.append(interp_data)
        valid_names.append(result['name'])
        valid_angles.append(result['angle'])

    print(f"[n-Mode]   Aligned {len(signals)} channels to common time grid")

    mirnov = MirnovData()
    mirnov.time = np.array(time_arr, dtype=np.float32)
    mirnov.signals = np.array(signals, dtype=np.float32)
    mirnov.positions = np.array(valid_angles)
    mirnov.channel_names = valid_names
    mirnov.shot = shot

    print(f"[n-Mode]   Loaded {len(valid_names)} channels, {len(time_arr)} points")
    return mirnov


# =============================================================================
# FFT Calculation
# =============================================================================

def calculate_fft(mirnov, t_interval, integrate=False, run_detrend=True, tmin=None, tmax=None):
    """Calculate FFT for all channels"""
    t0 = tclock.time()

    time = mirnov.time
    signals = mirnov.signals.astype(np.float32)
    n_channels = signals.shape[0]

    if tmin is not None and tmax is not None:
        idx_start = np.searchsorted(time, tmin)
        idx_end = np.searchsorted(time, tmax)
        time = time[idx_start:idx_end]
        signals = signals[:, idx_start:idx_end]

    n_points = signals.shape[1]
    dt = (time[100] - time[0]) / 100.0
    fs = 1.0 / dt

    q = int(n_points * dt / t_interval)
    p = int(n_points / q)
    n_windows, points_per_window = q, p

    print(f"[n-Mode] FFT: fs={fs/1e6:.2f}MHz, windows={n_windows}, pts/win={points_per_window}")

    # Step 1: Preprocessing
    print('[n-Mode]   1/4 Preprocessing...')
    if integrate:
        for i in range(n_channels):
            signals[i] = np.cumsum(signals[i]) * dt * 1e4

    n_use = n_windows * points_per_window
    signals = signals[:, :n_use]
    time_trimmed = time[:n_use]
    signals_reshaped = signals.reshape(n_channels, n_windows, points_per_window)

    # Step 2: Detrend
    print('[n-Mode]   2/4 Detrending...')
    if run_detrend:
        n_segments = 8
        denom = points_per_window % n_segments
        main_size = points_per_window - denom
        seg_size = main_size // n_segments

        main_part = signals_reshaped[:, :, :main_size].reshape(
            n_channels, n_windows, n_segments, seg_size)
        main_part = detrend(main_part, axis=3, type='linear')
        signals_reshaped[:, :, :main_size] = main_part.reshape(
            n_channels, n_windows, main_size)

        # Handle remaining points (denom)
        if denom > 0:
            remainder = signals_reshaped[:, :, -denom:].reshape(
                n_channels, n_windows, denom, 1)
            remainder = detrend(remainder, axis=2, type='linear')
            signals_reshaped[:, :, -denom:] = remainder.reshape(
                n_channels, n_windows, denom)

    # Step 3: FFT calculation
    print('[n-Mode]   3/4 FFT...')
    n_freq_half = points_per_window // 2 + 1
    overp = 1.0 / float(points_per_window)

    amp = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)
    phase = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)

    fft_result_raw = fft(signals_reshaped, axis=2, workers=-1)
    fft_half = fft_result_raw[:, :, :n_freq_half]
    amp = np.transpose(np.abs(fft_half) * overp, (1, 2, 0)).astype(np.float32)
    phase = np.transpose(np.angle(fft_half), (1, 2, 0)).astype(np.float32)

    # Step 4: Post-processing
    print('[n-Mode]   4/4 Finalizing...')
    window_start_indices = np.arange(n_windows) * points_per_window
    time_centers = (time_trimmed[window_start_indices] +
                   time_trimmed[window_start_indices + points_per_window - 1]) / 2.0

    idt = 1.0 / (points_per_window * dt)
    freq = np.arange(n_freq_half, dtype='f') * idt * 1e-3

    result = FFTResult()
    result.time = time_centers.astype(np.float32)
    result.freq = freq
    result.amp = amp
    result.phase = phase
    result.positions = mirnov.positions
    result.shot = mirnov.shot
    result.p = points_per_window
    result.dt = dt

    print(f'\n[n-Mode]   Completed in {tclock.time() - t0:.2f} sec')
    return result


# =============================================================================
# Mode Number Calculation
# =============================================================================

def calculate_mode_numbers(fft_result, fmin, fmax, tol, nmodes, frac, msign, integrate=False):
    """Calculate toroidal mode numbers"""
    t0 = tclock.time()

    amp = fft_result.amp
    phase = fft_result.phase
    freq = fft_result.freq
    positions = fft_result.positions

    n_time, freq_range, n_channels = amp.shape
    nch = 1

    fmax_threshold = frac * np.max(amp[:, :, nch])
    mode = np.zeros((n_time, freq_range), dtype='int')

    wh = np.where(amp[:, :, nch] >= fmax_threshold)
    wh_i, wh_j = wh[0], wh[1]
    valid = wh_j != 0
    wh_i, wh_j = wh_i[valid], wh_j[valid]

    if len(wh_i) == 0:
        return mode, freq, np.zeros((2 * nmodes, n_time))

    pos_rad = np.array(positions) * np.pi / 180.0
    nmeas = n_channels
    nvec = np.arange(nmodes) + 1.0

    pharr = nvec * pos_rad[:, np.newaxis] * 1j
    exparr = np.exp(pharr) / float(nmeas)

    bamp = amp[wh_i, wh_j, :]
    barg = phase[wh_i, wh_j, :]
    bvec = bamp * np.exp(1j * barg)
    b2avg = np.sum(bamp ** 2, axis=1) / nmeas

    bmodep = np.dot(bvec, exparr)
    bmoden = np.dot(bvec, np.conj(exparr))

    if bmodep.ndim == 1:
        bmodep = bmodep.reshape(1, -1)
        bmoden = bmoden.reshape(1, -1)
        b2avg = np.array([b2avg]) if np.isscalar(b2avg) else b2avg

    b2avg_vec = np.tile(b2avg[:, np.newaxis], (1, nmodes))

    std_pos = np.sqrt(np.abs(1.0 - np.abs(bmodep) ** 2 / b2avg_vec))
    std_neg = np.sqrt(np.abs(1.0 - np.abs(bmoden) ** 2 / b2avg_vec))

    kmn_pos = std_pos.argmin(axis=1)
    mn_pos = std_pos[np.arange(len(kmn_pos)), kmn_pos]
    kmn_neg = std_neg.argmin(axis=1)
    mn_neg = std_neg[np.arange(len(kmn_neg)), kmn_neg]

    pos_assign = (mn_pos <= mn_neg) & (mn_pos <= tol)
    neg_assign = (mn_pos > mn_neg) & (mn_neg <= tol)

    mode[wh_i[pos_assign], wh_j[pos_assign]] = kmn_pos[pos_assign] + 1
    mode[wh_i[neg_assign], wh_j[neg_assign]] = -kmn_neg[neg_assign] - 1

    amp_for_evolution = np.abs(amp[:, :, 0]).copy()
    if not integrate:
        amp_for_evolution = 2.0 * amp_for_evolution
        freq_hz = freq * 1e3
        freq_hz_safe = np.where(freq_hz > 0, freq_hz, 1.0)
        amp_for_evolution = amp_for_evolution / (2.0 * np.pi * freq_hz_safe)
        amp_for_evolution[:, freq_hz == 0] = 0.0
        amp_for_evolution = amp_for_evolution * 1e4

    freq_idx = np.where((freq >= fmin) & (freq <= fmax))[0]
    amp_evolution = np.zeros((2 * nmodes, n_time))
    temp = np.zeros_like(amp_for_evolution)

    for i in range(nmodes):
        for sign_mult, offset in [(1, 0), (-1, nmodes)]:
            ind1 = np.where(mode == sign_mult * (i + 1))
            if len(ind1[0]) > 0:
                temp[ind1] = amp_for_evolution[ind1]
                amp_evolution[i + offset, :] = temp[:, freq_idx].max(axis=1)
                temp[ind1] = 0.0

    print(f"[n-Mode] Mode calculation completed in {tclock.time() - t0:.2f} sec")
    return mode, freq, amp_evolution
