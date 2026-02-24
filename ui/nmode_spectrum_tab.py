"""
N-Mode Spectrum tab for toroidal mode number analysis
Based on KSTARtMirnovAnalysis_module.py, integrated into PRISM
"""

import json
from pathlib import Path
import numpy as np
import gc
import os
from scipy.signal import detrend
from scipy.fft import fft
import scipy.ndimage as ndimage
import time as tclock
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from MDSplus import Connection

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QComboBox,
    QCheckBox, QRadioButton, QButtonGroup, QSlider, QFrame,
    QMessageBox, QFileDialog, QApplication, QSplitter,
    QDialog, QTextEdit, QStyle, QScrollArea,
    QSpinBox, QDialogButtonBox,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QColor, QTextCharFormat, QSyntaxHighlighter, QGuiApplication

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from config.user_settings import get_tab_settings, set_tab_settings

# Load Mirnov coil configuration
_MIRNOV_CONFIG_PATH = Path(__file__).parent.parent / 'config' / 'mirnov_config.json'
with open(_MIRNOV_CONFIG_PATH, 'r') as f:
    _MIRNOV_CONFIG = json.load(f)


# Example script for loading n-mode spectrum NPZ file
NMODE_EXAMPLE_SCRIPT = '''"""
Example script for loading and plotting PRISM n-mode spectrum NPZ file
"""

import numpy as np
import matplotlib.pyplot as plt


def get_mode_color(n):
    """Get color for mode number"""
    colors = ['#000000', '#FF0000', '#00FF00', '#0000FF', '#FF9100',
              '#1f78b4', '#b2df8a', '#33a02c', '#ffffb3', '#b3bada',
              '#1b9e77', '#b15928', '#CCEBC5', '#ff3d6f', '#e6ab02']
    return colors[abs(n) % len(colors)]


# Load NPZ file
filepath = 'nmode_XXXXXX_3.0-5.0s.npz'
data = np.load(filepath, allow_pickle=True)

# Extract data
time = data['time']                   # Time array [s]
frequency = data['frequency']         # Frequency array [kHz]
mode_spectrum = data['mode_spectrum'] # Mode number array (time, freq)
amplitude = data['amplitude']         # Amplitude evolution (2*nmodes, time)

# Get metadata
metadata = data['metadata'].item()
print(f"Shot: {metadata['shot']}")
print(f"Time range: {metadata['time_range']} s")
print(f"Freq range: {metadata['freq_range']} kHz")
print(f"n-modes: {metadata['n_modes']}")
print(f"n_modes_list: {metadata['n_modes_list']}")
print(f"Sign: {metadata['sign']}")

# Plot settings
n_modes_list = metadata['n_modes_list']
fmin, fmax = metadata['freq_range']
tmin, tmax = metadata['time_range']

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

# Plot n-mode spectrum with colors for each mode
for n in n_modes_list:
    if n == 0:
        continue
    color = get_mode_color(abs(n))

    # Create masked array for this mode
    mode_mask = (mode_spectrum == n)
    if not np.any(mode_mask):
        continue

    # Get positions where this mode exists
    t_idx, f_idx = np.where(mode_mask)
    if len(t_idx) == 0:
        continue

    ax1.scatter(time[t_idx], frequency[f_idx], c=color, s=1, label=f'n={n}')

ax1.set_ylabel('Frequency [kHz]')
ax1.set_xlim(tmin, tmax)
ax1.set_ylim(fmin, fmax)
ax1.set_title(f"Shot #{metadata['shot']} n-mode spectrum ({metadata['sign']})")
ax1.legend(loc='upper right', fontsize=8, markerscale=5)
ax1.grid(True, alpha=0.3)

# Plot amplitude evolution with matching colors
for i, n in enumerate(n_modes_list):
    color = get_mode_color(abs(n))
    ax2.plot(time, amplitude[i], color=color, label=f'n={n}')

ax2.set_xlabel('Time [s]')
ax2.set_ylabel('Amplitude [a.u.]')
ax2.set_xlim(tmin, tmax)
ax2.legend(loc='upper right', fontsize=8)
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
'''


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
    shot, name, angle, sign, tmin, tmax = args

    try:
        mds = Connection(NModeConfig.MDS_SERVER)
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


def load_mirnov_data(shot, tmin, tmax):
    """Load Mirnov coil data from MDS+ with parallel processing"""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    year = get_shot_year(shot)
    config = get_mirnov_config(year)
    n_channels = len(config['names'])

    print(f"[n-Mode] Loading Mirnov data for shot #{shot} (year {year})...")
    print(f"[n-Mode]   Time range: {tmin:.2f} - {tmax:.2f} s")
    print(f"[n-Mode]   Channels ({n_channels}):")
    for i, (name, angle) in enumerate(zip(config['names'], config['angles']), 1):
        print(f"[n-Mode]     {i:>2d}. {name}  {angle:>6.2f}\u00b0")
    print(f"[n-Mode]   Threads: {n_channels}")

    args_list = [
        (shot, name, angle, sign, tmin, tmax)
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
            bar = '\u2588' * int(40 * progress) + '\u2591' * (40 - int(40 * progress))
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


# =============================================================================
# Plotting
# =============================================================================

def get_mode_color(n):
    """Get color for mode number"""
    colors = ['#000000', '#FF0000', '#00FF00', '#0000FF', '#FF9100',
              '#1f78b4', '#b2df8a', '#33a02c', '#ffffb3', '#b3bada',
              '#1b9e77', '#b15928', '#CCEBC5', '#ff3d6f', '#e6ab02']
    return colors[abs(n) % len(colors)]


def plot_mode_spectrum(ax, fft_result, mode, freq_use, amp_evolution,
                       tmin, tmax, fmin, fmax, nmodes, msign, integrate,
                       plot_type='contour', numc=50, colors=None,
                       legend_fontsize=8, label_fontsize=None, tick_fontsize=None):
    """Plot n-mode spectrum"""
    ax.clear()

    time = fft_result.time
    amp = fft_result.amp

    sz = amp.shape
    n_time = sz[0]

    j1 = np.argmin(np.abs(time - tmin))
    j2 = np.argmin(np.abs(time - tmax))
    k1 = np.argmin(np.abs(freq_use - fmin))
    k2 = np.argmin(np.abs(freq_use - fmax))

    mode_plot = mode[j1:j2, k1:k2]

    if msign == 1:
        mode_filtered = (mode_plot > 0) * mode_plot
    elif msign == -1:
        mode_filtered = np.abs((mode_plot < 0) * mode_plot)
    elif msign == 0:
        mode_filtered = np.abs(mode_plot)
    else:
        mode_filtered = mode_plot

    amp_plot = np.abs(amp[j1:j2, k1:k2, 0])

    if not integrate:
        amp_plot = 2.0 * amp_plot
        freq_hz = freq_use[k1:k2] * 1e3
        freq_hz_safe = np.where(freq_hz > 0, freq_hz, 1.0)
        amp_plot = amp_plot / (2.0 * np.pi * freq_hz_safe)
        amp_plot[:, freq_hz == 0] = 0
        amp_plot = amp_plot * 1e4

    x = np.arange(mode_plot.shape[0]) * (tmax - tmin) / (mode_plot.shape[0] - 1) + tmin
    y = np.arange(k2 - k1) * (fmax - fmin) / (k2 - k1) + fmin

    if msign == 2:
        mxmode = np.max(np.abs((mode_filtered > 0) * mode_filtered))
        mnmode = -np.max(np.abs((mode_filtered < 0) * mode_filtered))
    else:
        mxmode = np.max(mode_filtered) if np.max(mode_filtered) > 0 else nmodes
        mnmode = 1

    lv = np.arange(int(mnmode), int(mxmode) + 1)
    lv = lv[lv != 0]

    from matplotlib.colors import LogNorm

    for i, n_val in enumerate(lv):
        if n_val == 0:
            continue

        k = np.where(mode_filtered == n_val)
        if len(k[0]) == 0:
            continue

        thismode = np.zeros_like(amp_plot)
        thismode[k] = amp_plot[k]

        color = colors[abs(n_val) - 1] if colors and abs(n_val) - 1 < len(colors) else get_mode_color(abs(n_val))

        if plot_type == 'contour':
            valid_data = thismode[thismode > 0]
            if len(valid_data) == 0:
                continue

            try:
                levels = np.logspace(
                    np.log10(valid_data.min()),
                    np.log10(valid_data.max()),
                    numc
                )
                ax.contour(x, y, thismode.T, levels, colors=[color], linewidths=0.8)
            except Exception:
                ax.contour(x, y, thismode.T, numc, colors=[color], linewidths=0.8)
        else:
            thismode_masked = np.ma.masked_where(thismode == 0, thismode)
            ax.pcolormesh(x, y, thismode_masked.T, shading='auto', alpha=0.7)

    if label_fontsize:
        ax.set_ylabel('Frequency [kHz]', fontsize=label_fontsize)
    else:
        ax.set_ylabel('Frequency [kHz]')
    if tick_fontsize:
        ax.tick_params(labelsize=tick_fontsize)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(fmin, fmax)
    ax.set_axisbelow(False)
    import matplotlib as mpl
    ax.grid(True, linestyle='--', linewidth=0.3, color=mpl.rcParams.get('grid.color', '#444444'))

    ax.set_title(f'#{fft_result.shot}')

    # Always show all modes in legend based on nmodes setting
    legend_elements = []
    from matplotlib.lines import Line2D
    for n in range(1, nmodes + 1):
        color = colors[n - 1] if colors and n - 1 < len(colors) else get_mode_color(n)
        if msign == 2:
            legend_elements.append(Line2D([0], [0], color=color, label=f'n=+{n}'))
            legend_elements.append(Line2D([0], [0], color=color, linestyle='--', label=f'n=-{n}'))
        elif msign == 1:
            legend_elements.append(Line2D([0], [0], color=color, label=f'n=+{n}'))
        elif msign == -1:
            legend_elements.append(Line2D([0], [0], color=color, label=f'n=-{n}'))
        else:  # msign == 0 (abs)
            legend_elements.append(Line2D([0], [0], color=color, label=f'n={n}'))


    if legend_elements:
        ax.legend(handles=legend_elements, loc='upper right', fontsize=legend_fontsize,
                  framealpha=0.8, ncol=min(len(legend_elements), 5))

    return lv


def plot_amplitude_evolution(ax, fft_result, amp_evolution, tmin, tmax, nmodes, msign,
                              colors=None, legend_fontsize=8, label_fontsize=None, tick_fontsize=None):
    """Plot amplitude evolution for each mode"""
    ax.clear()

    time = fft_result.time
    j1 = np.argmin(np.abs(time - tmin))
    j2 = np.argmin(np.abs(time - tmax))
    time_plot = time[j1:j2]

    for i in range(nmodes):
        n = i + 1
        if msign in [0, 1, 2]:
            amp_pos = amp_evolution[i, j1:j2]
            c = colors[n - 1] if colors and n - 1 < len(colors) else get_mode_color(n)
            ax.plot(time_plot, amp_pos, color=c,
                   linewidth=1.5, label=f'n=+{n}')
        if msign in [0, -1, 2]:
            amp_neg = amp_evolution[i + nmodes, j1:j2]
            ls = '--' if msign == 2 else '-'
            c = colors[n - 1] if colors and n - 1 < len(colors) else get_mode_color(n)
            ax.plot(time_plot, amp_neg, color=c,
                   linewidth=1.5, linestyle=ls, label=f'n=-{n}')

    if label_fontsize:
        ax.set_xlabel('Time [s]', fontsize=label_fontsize)
        ax.set_ylabel('Amplitude [Gauss]', fontsize=label_fontsize)
    else:
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Amplitude [Gauss]')
    if tick_fontsize:
        ax.tick_params(labelsize=tick_fontsize)
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=legend_fontsize, framealpha=0.8, ncol=5)
    import matplotlib as mpl
    ax.grid(True, linestyle='--', linewidth=0.3, color=mpl.rcParams.get('grid.color', '#444444'))


# =============================================================================
# N-Mode Spectrum Tab
# =============================================================================

class NModeSpectrumTab:
    """N-Mode Spectrum analysis tab for PRISM"""

    # Shared column 0 width for consistent alignment across panels
    LABEL_COLUMN_WIDTH = 105

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = QWidget()
        self.toolbar = None

        self.mirnov_data = None
        self.fft_result = None
        self.mode_result = None
        self.freq_use = None
        self.amp_evolution = None
        self.current_shot = None

        # Plot options
        self.color_mode = "Default"
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

    def create_widgets(self):
        """Create n-mode spectrum tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)
        self.figure.subplots_adjust(left=0.10, right=0.97, top=0.95, bottom=0.10, hspace=0.15)
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212, sharex=self.ax1)
        self.ax1.set_ylabel('Frequency [kHz]')
        self.ax2.set_xlabel('Time [s]')
        self.ax2.set_ylabel('Amplitude')
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        # Left side: canvas
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.addWidget(self.canvas)

        # Right side: scrollable control panel with fixed width
        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QVBoxLayout.SetNoConstraint)

        scroll_area.setWidget(control_frame)
        scroll_area.viewport().setAutoFillBackground(False)
        control_frame.setAutoFillBackground(False)

        # Splitter for left/right layout
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(canvas_widget)
        splitter.addWidget(scroll_area)

        # Set the main layout of self.frame
        main_layout = QVBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(splitter)

        self._create_parameters_panel(control_frame)
        self._create_run_plot_panel(control_frame)
        self._create_save_controls(control_frame)

        # Add stretch at the bottom of control panel
        control_layout.addStretch()

        # Load saved settings
        self.load_settings()

    def _create_parameters_panel(self, parent):
        """Create parameters panel"""
        frame = QGroupBox("1. Parameters", parent)
        grid = QGridLayout(frame)

        grid.setColumnMinimumWidth(0, self.LABEL_COLUMN_WIDTH)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)

        row = 0

        # Shot with up/down buttons in same row
        grid.addWidget(QLabel('Shot'), row, 0)

        shot_widget = QWidget()
        shot_layout = QHBoxLayout(shot_widget)
        shot_layout.setContentsMargins(0, 0, 0, 0)

        self.shot_entry = QLineEdit()
        self.shot_entry.setText(str(NModeConfig.DEFAULT_SHOT))
        self.shot_entry.setFixedWidth(80)
        self.shot_entry.returnPressed.connect(self._run_calculation)
        shot_layout.addWidget(self.shot_entry)

        btn_updown = QWidget()
        btn_updown_layout = QVBoxLayout(btn_updown)
        btn_updown_layout.setContentsMargins(0, 0, 0, 0)
        btn_updown_layout.setSpacing(0)
        mini_btn_style = "padding: 0px; border-radius: 2px;"
        up_btn = QPushButton()
        up_btn.setIcon(get_icon(QStyle.SP_ArrowUp))
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_updown_layout.addWidget(up_btn)
        down_btn = QPushButton()
        down_btn.setIcon(get_icon(QStyle.SP_ArrowDown))
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_updown_layout.addWidget(down_btn)
        shot_layout.addWidget(btn_updown)

        shot_layout.addStretch()

        grid.addWidget(shot_widget, row, 1)

        row += 1

        # Time range
        grid.addWidget(QLabel('Time [s]'), row, 0)

        time_widget = QWidget()
        time_grid = QGridLayout(time_widget)
        time_grid.setContentsMargins(0, 0, 0, 0)

        self.tmin_entry = QLineEdit()
        self.tmin_entry.setText(str(NModeConfig.DEFAULT_TMIN))
        self.tmin_entry.setFixedWidth(80)
        time_grid.addWidget(self.tmin_entry, 0, 0)

        time_grid.addWidget(QLabel('-'), 0, 1)

        self.tmax_entry = QLineEdit()
        self.tmax_entry.setText(str(NModeConfig.DEFAULT_TMAX))
        self.tmax_entry.setFixedWidth(80)
        time_grid.addWidget(self.tmax_entry, 0, 2)

        # Use full shot checkbox (below tmin entry)
        self.use_full_shot_checkbox = QCheckBox('Use full shot length')
        self.use_full_shot_checkbox.setChecked(False)
        self.use_full_shot_checkbox.stateChanged.connect(self._toggle_time_entries)
        time_grid.addWidget(self.use_full_shot_checkbox, 1, 0, 1, 3)

        grid.addWidget(time_widget, row, 1, 1, 3)

        row += 1

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        grid.addWidget(separator, row, 0, 1, 4)
        row += 1

        # Time interval
        grid.addWidget(QLabel('dt [s]'), row, 0)
        self.tinterval_entry = QLineEdit()
        self.tinterval_entry.setText('0.01')
        self.tinterval_entry.setFixedWidth(80)
        grid.addWidget(self.tinterval_entry, row, 1)
        row += 1

        # Frequency range
        grid.addWidget(QLabel('Freq [kHz]'), row, 0)

        freq_widget = QWidget()
        freq_layout = QHBoxLayout(freq_widget)
        freq_layout.setContentsMargins(0, 0, 0, 0)

        self.fmin_entry = QLineEdit()
        self.fmin_entry.setText(str(NModeConfig.DEFAULT_FMIN))
        self.fmin_entry.setFixedWidth(80)
        freq_layout.addWidget(self.fmin_entry)

        freq_layout.addWidget(QLabel('-'))

        self.fmax_entry = QLineEdit()
        self.fmax_entry.setText(str(NModeConfig.DEFAULT_FMAX))
        self.fmax_entry.setFixedWidth(80)
        freq_layout.addWidget(self.fmax_entry)
        freq_layout.addStretch()

        grid.addWidget(freq_widget, row, 1, 1, 3)
        row += 1

        # n-modes slider
        grid.addWidget(QLabel('n-modes'), row, 0)

        nmodes_widget = QWidget()
        nmodes_layout = QHBoxLayout(nmodes_widget)
        nmodes_layout.setContentsMargins(0, 0, 0, 0)

        self._nmodes_value = NModeConfig.DEFAULT_NMODES
        self.nmodes_slider = QSlider(Qt.Horizontal)
        self.nmodes_slider.setMinimum(1)
        self.nmodes_slider.setMaximum(8)
        self.nmodes_slider.setValue(NModeConfig.DEFAULT_NMODES)
        self.nmodes_slider.valueChanged.connect(self._on_nmodes_changed)
        nmodes_layout.addWidget(self.nmodes_slider, 1)

        self.nmodes_label = QLabel(str(NModeConfig.DEFAULT_NMODES))
        self.nmodes_label.setFixedWidth(25)
        nmodes_layout.addWidget(self.nmodes_label)

        grid.addWidget(nmodes_widget, row, 1, 1, 3)
        row += 1

        # Tolerance
        grid.addWidget(QLabel('Tolerance'), row, 0)
        self.tol_entry = QLineEdit()
        self.tol_entry.setText(str(NModeConfig.DEFAULT_TOL))
        self.tol_entry.setFixedWidth(80)
        grid.addWidget(self.tol_entry, row, 1)
        row += 1

        # Fraction
        grid.addWidget(QLabel('Fraction'), row, 0)
        self.frac_entry = QLineEdit()
        self.frac_entry.setText(str(NModeConfig.DEFAULT_FRAC))
        self.frac_entry.setFixedWidth(80)
        grid.addWidget(self.frac_entry, row, 1)
        row += 1

        # Sign (radio buttons)
        grid.addWidget(QLabel('Sign'), row, 0)

        sign_widget = QWidget()
        sign_layout = QHBoxLayout(sign_widget)
        sign_layout.setContentsMargins(0, 0, 0, 0)

        self.msign_button_group = QButtonGroup()
        self._msign_value = NModeConfig.DEFAULT_MSIGN

        radio_abs = QRadioButton('abs')
        radio_abs.setProperty('sign_value', 0)
        sign_layout.addWidget(radio_abs)
        self.msign_button_group.addButton(radio_abs)

        radio_pos = QRadioButton('pos')
        radio_pos.setProperty('sign_value', 1)
        sign_layout.addWidget(radio_pos)
        self.msign_button_group.addButton(radio_pos)

        radio_neg = QRadioButton('neg')
        radio_neg.setProperty('sign_value', -1)
        sign_layout.addWidget(radio_neg)
        self.msign_button_group.addButton(radio_neg)

        radio_all = QRadioButton('all')
        radio_all.setProperty('sign_value', 2)
        sign_layout.addWidget(radio_all)
        self.msign_button_group.addButton(radio_all)

        # Set default selection
        radio_pos.setChecked(True)

        self.msign_button_group.buttonClicked.connect(self._on_msign_changed)

        grid.addWidget(sign_widget, row, 1, 1, 3)
        row += 1

        # Options (integrate + detrend checkboxes)
        option_widget = QWidget()
        option_layout = QHBoxLayout(option_widget)
        option_layout.setContentsMargins(0, 0, 0, 0)

        self.integrate_checkbox = QCheckBox('Integrate (dB/dt -> B)')
        self.integrate_checkbox.setChecked(NModeConfig.DEFAULT_INTEGRATE)
        option_layout.addWidget(self.integrate_checkbox)

        self.detrend_checkbox = QCheckBox('Detrend')
        self.detrend_checkbox.setChecked(NModeConfig.DEFAULT_DETREND)
        option_layout.addWidget(self.detrend_checkbox)

        grid.addWidget(option_widget, row, 0, 1, 4)

        parent.layout().addWidget(frame)

    def _create_run_plot_panel(self, parent):
        """Create run and plot options panel"""
        frame = QGroupBox("2. Plot", parent)
        grid = QGridLayout(frame)

        # Match Parameters panel column structure
        grid.setColumnMinimumWidth(0, self.LABEL_COLUMN_WIDTH)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)

        row = 0

        # Calculate and Plot button + Option button
        self.run_button = QPushButton('Calculate and Plot')
        self.run_button.clicked.connect(self._run_calculation)
        grid.addWidget(self.run_button, row, 0, 1, 3)

        plot_options_btn = QPushButton('Option')
        plot_options_btn.clicked.connect(self._show_plot_options_dialog)
        grid.addWidget(plot_options_btn, row, 3)
        row += 1

        # Status
        self.status_label = QLabel('Ready')
        self.status_label.setStyleSheet('color: gray; font-weight: bold;')
        grid.addWidget(self.status_label, row, 0, 1, 4)
        row += 1

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        grid.addWidget(separator, row, 0, 1, 4)
        row += 1

        # Plot type (radio buttons)
        grid.addWidget(QLabel('Plot type'), row, 0)

        plot_type_widget = QWidget()
        plot_type_layout = QHBoxLayout(plot_type_widget)
        plot_type_layout.setContentsMargins(0, 0, 0, 0)

        self.plot_type_button_group = QButtonGroup()
        self._plot_type_value = NModeConfig.DEFAULT_PLOT_TYPE

        radio_contour = QRadioButton('contour')
        radio_contour.setProperty('plot_type_value', 'contour')
        radio_contour.setChecked(True)
        plot_type_layout.addWidget(radio_contour)
        self.plot_type_button_group.addButton(radio_contour)

        radio_imshow = QRadioButton('imshow')
        radio_imshow.setProperty('plot_type_value', 'imshow')
        plot_type_layout.addWidget(radio_imshow)
        self.plot_type_button_group.addButton(radio_imshow)

        self.plot_type_button_group.buttonClicked.connect(self._on_plot_type_changed)

        grid.addWidget(plot_type_widget, row, 1, 1, 3)
        row += 1

        # Contour levels
        self.contour_levels_label = QLabel('Contour levels')
        grid.addWidget(self.contour_levels_label, row, 0)

        self.contour_levels_entry = QLineEdit()
        self.contour_levels_entry.setText(str(NModeConfig.DEFAULT_NUMC))
        self.contour_levels_entry.setFixedWidth(80)
        grid.addWidget(self.contour_levels_entry, row, 1)
        row += 1

        # Update Plot button
        update_btn = QPushButton('Update Plot')
        update_btn.clicked.connect(self._update_plot)
        grid.addWidget(update_btn, row, 0, 1, 4)

        parent.layout().addWidget(frame)

    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = QGroupBox("3. Save Data", parent)
        layout = QVBoxLayout(frame)

        btn_widget = QWidget()
        btn_layout = QHBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        self.save_button = QPushButton('Save as .npz')
        self.save_button.clicked.connect(self._save_data)
        self.save_button.setEnabled(False)
        btn_layout.addWidget(self.save_button, 1)

        example_btn = QPushButton('Example Script')
        example_btn.clicked.connect(self._show_example_script)
        btn_layout.addWidget(example_btn, 1)

        layout.addWidget(btn_widget)

        parent.layout().addWidget(frame)

    def _show_example_script(self):
        """Show example script for loading NPZ file with syntax highlighting"""
        script = NMODE_EXAMPLE_SCRIPT

        popup = QDialog(self.frame)
        popup.setWindowTitle("Example Script - N-Mode Spectrum NPZ")
        popup.resize(700, 550)

        dialog_layout = QVBoxLayout(popup)

        text_widget = QTextEdit()
        text_widget.setReadOnly(True)
        text_widget.setFont(QFont('Courier', 10))
        text_widget.setStyleSheet(
            'background-color: #19232d; color: #ffffff;'
        )
        text_widget.setPlainText(script)

        # Apply syntax highlighting
        self._highlighter = NModePythonHighlighter(text_widget.document())

        dialog_layout.addWidget(text_widget)

        btn_widget = QWidget()
        btn_layout = QHBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        def copy_to_clipboard():
            clipboard = QGuiApplication.clipboard()
            clipboard.setText(script)
            QMessageBox.information(popup, "Copied", "Script copied to clipboard")

        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(copy_to_clipboard)
        btn_layout.addWidget(copy_btn)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(popup.accept)
        btn_layout.addWidget(close_btn)

        dialog_layout.addWidget(btn_widget)
        popup.exec()

    def _save_data(self):
        """Save n-mode spectrum data to NPZ file"""
        if self.fft_result is None or self.mode_result is None:
            QMessageBox.warning(self.frame, "Warning", "No data to save. Run calculation first.")
            return

        # Get parameters
        shot = self.current_shot
        tmin = float(self.tmin_entry.text())
        tmax = float(self.tmax_entry.text())
        fmin = float(self.fmin_entry.text())
        fmax = float(self.fmax_entry.text())
        nmodes = self._nmodes_value

        # Default filename
        default_name = f"nmode_{shot}_{tmin:.1f}-{tmax:.1f}s.npz"

        # File dialog
        filepath, _ = QFileDialog.getSaveFileName(
            self.frame,
            'Save N-Mode Spectrum Data',
            os.path.join(os.path.expanduser("~"), default_name),
            'NumPy NPZ (*.npz)'
        )

        if not filepath:
            return

        try:
            # Build n_modes list based on sign setting
            msign = self._msign_value
            sign_str = {0: 'abs', 1: 'pos', -1: 'neg', 2: 'all'}

            if msign == 0:  # abs
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == 1:  # pos
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == -1:  # neg
                n_modes_list = list(range(-nmodes, 0))
            else:  # all
                n_modes_list = list(range(-nmodes, 0)) + list(range(1, nmodes + 1))

            # Filter mode_spectrum and amplitude by sign setting
            mode_save = self.mode_result.copy()
            amp_save = self.amp_evolution.copy()

            if msign == 1:  # pos only
                mode_save = np.where(mode_save > 0, mode_save, 0)
                amp_save = amp_save[:nmodes, :]
            elif msign == -1:  # neg only
                mode_save = np.where(mode_save < 0, mode_save, 0)
                amp_save = amp_save[nmodes:, :]
            elif msign == 0:  # abs
                mode_save = np.abs(mode_save)

            # Prepare metadata
            metadata = {
                'shot': shot,
                'time_range': [tmin, tmax],
                'freq_range': [fmin, fmax],
                'time_interval': float(self.tinterval_entry.text()),
                'n_modes': nmodes,
                'n_modes_list': n_modes_list,
                'tolerance': float(self.tol_entry.text()),
                'fraction': float(self.frac_entry.text()),
                'sign': sign_str.get(msign, 'pos'),
                'integrate': self.integrate_checkbox.isChecked(),
                'detrend': self.detrend_checkbox.isChecked(),
            }

            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=self.fft_result.time,
                     frequency=self.freq_use,
                     mode_spectrum=mode_save,
                     amplitude=amp_save
            )

            print(f"[n-Mode] Data saved to: {filepath}")
            QMessageBox.information(self.frame, "Saved", f"Data saved to:\n{filepath}")

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to save data: {str(e)}")

    def _show_plot_options_dialog(self):
        """Show Plot Options dialog for N-Mode spectrum"""
        WIDGET_WIDTH = 150

        dialog = QDialog(self.frame)
        dialog.setWindowTitle("Plot Options")
        dialog.setMinimumWidth(280)
        dlg_layout = QVBoxLayout(dialog)

        # Color mode
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        color_combo = QComboBox()
        color_combo.setFixedWidth(WIDGET_WIDTH)
        color_combo.addItems([
            "Default",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
        ])
        color_combo.setCurrentText(self.color_mode)
        color_row.addWidget(color_combo)
        dlg_layout.addLayout(color_row)

        # Label font size
        label_row = QHBoxLayout()
        label_row.addWidget(QLabel("Label font size"))
        label_spin = QSpinBox()
        label_spin.setFixedWidth(WIDGET_WIDTH)
        label_spin.setRange(6, 24)
        label_spin.setValue(self.label_fontsize)
        label_row.addWidget(label_spin)
        dlg_layout.addLayout(label_row)

        # Legend font size
        legend_row = QHBoxLayout()
        legend_row.addWidget(QLabel("Legend font size"))
        legend_spin = QSpinBox()
        legend_spin.setFixedWidth(WIDGET_WIDTH)
        legend_spin.setRange(4, 20)
        legend_spin.setValue(self.legend_fontsize)
        legend_row.addWidget(legend_spin)
        dlg_layout.addLayout(legend_row)

        # Tick font size
        tick_row = QHBoxLayout()
        tick_row.addWidget(QLabel("Tick font size"))
        tick_spin = QSpinBox()
        tick_spin.setFixedWidth(WIDGET_WIDTH)
        tick_spin.setRange(6, 20)
        tick_spin.setValue(self.tick_fontsize)
        tick_row.addWidget(tick_spin)
        dlg_layout.addLayout(tick_row)

        # Default / OK / Cancel
        btn_box = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)

        def reset_defaults():
            color_combo.setCurrentText("Default")
            label_spin.setValue(12)
            legend_spin.setValue(8)
            tick_spin.setValue(10)
        btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(reset_defaults)
        dlg_layout.addWidget(btn_box)

        if dialog.exec() == QDialog.Accepted:
            self.color_mode = color_combo.currentText()
            self.label_fontsize = label_spin.value()
            self.legend_fontsize = legend_spin.value()
            self.tick_fontsize = tick_spin.value()
            # Auto-replot if data exists
            if self.fft_result is not None and self.mode_result is not None:
                self._update_plot()

    def _get_mode_colors(self, nmodes):
        """Get colors for mode numbers using selected color mode"""
        if self.color_mode == "Default":
            # Use built-in n-mode color palette
            return [get_mode_color(i + 1) for i in range(nmodes)]

        from plotting.plot_manager import ColorManager
        start = self.color_mode.find('(')
        end = self.color_mode.find(')')
        if start != -1 and end != -1:
            cmap_name = self.color_mode[start + 1:end]
        else:
            cmap_name = 'tab10'
        cm = ColorManager()
        entries = list(range(nmodes))
        return cm.get_colors_for_entries(entries, colormap=cmap_name)

    def _toggle_time_entries(self):
        """Enable/disable time entry fields based on checkbox state"""
        if self.use_full_shot_checkbox.isChecked():
            self.tmin_entry.setEnabled(False)
            self.tmax_entry.setEnabled(False)
        else:
            self.tmin_entry.setEnabled(True)
            self.tmax_entry.setEnabled(True)

    def _on_plot_type_changed(self, button):
        """Handle plot type radio button change and toggle contour levels"""
        self._plot_type_value = button.property('plot_type_value')
        self._toggle_contour_levels()

    def _toggle_contour_levels(self):
        """Enable/disable contour levels entry based on plot type"""
        if self._plot_type_value == 'contour':
            self.contour_levels_label.setEnabled(True)
            self.contour_levels_entry.setEnabled(True)
        else:
            self.contour_levels_label.setEnabled(False)
            self.contour_levels_entry.setEnabled(False)

    def _on_nmodes_changed(self, value):
        """Update nmodes label"""
        self._nmodes_value = value
        self.nmodes_label.setText(str(value))

    def _on_msign_changed(self, button):
        """Handle sign radio button change"""
        self._msign_value = button.property('sign_value')

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _run_calculation(self):
        """Run complete workflow"""
        t0 = tclock.time()

        try:
            shot = int(self.shot_entry.text())
            use_full_shot = self.use_full_shot_checkbox.isChecked()

            # When using full shot, use wide time range for initial load
            if use_full_shot:
                load_tmin = -1.0
                load_tmax = 100.0
            else:
                load_tmin = float(self.tmin_entry.text())
                load_tmax = float(self.tmax_entry.text())

            t_interval = float(self.tinterval_entry.text())
            fmin = float(self.fmin_entry.text())
            fmax = float(self.fmax_entry.text())
            tol = float(self.tol_entry.text())
            nmodes = self._nmodes_value
            frac = float(self.frac_entry.text())
            msign = self._msign_value
            integrate = self.integrate_checkbox.isChecked()
            run_detrend = self.detrend_checkbox.isChecked()

            self.run_button.setEnabled(False)
            self.status_label.setText('Running...')
            self.status_label.setStyleSheet('color: blue; font-weight: bold;')
            QApplication.processEvents()

            need_reload = (
                self.mirnov_data is None or
                self.current_shot != shot or
                getattr(self, '_last_load_tmin', None) != load_tmin or
                getattr(self, '_last_load_tmax', None) != load_tmax
            )

            if need_reload:
                self.mirnov_data = None
                self.fft_result = None
                gc.collect()

                self.status_label.setText('Loading data...')
                self.status_label.setStyleSheet('color: blue; font-weight: bold;')
                QApplication.processEvents()

                self.mirnov_data = load_mirnov_data(shot, load_tmin, load_tmax)
                self.current_shot = shot
                self._last_load_tmin = load_tmin
                self._last_load_tmax = load_tmax

            # Get actual data time range
            actual_tmin = self.mirnov_data.time[0]
            actual_tmax = self.mirnov_data.time[-1]

            # Determine tmin/tmax to use for calculation
            if use_full_shot:
                # Start from 0 if actual_tmin is negative
                tmin = max(0.0, actual_tmin)
                tmax = actual_tmax
                # Update display (entry fields are disabled but show the values)
                self.tmin_entry.setText(f"{tmin:.3f}")
                self.tmax_entry.setText(f"{tmax:.3f}")
                print(f"[n-Mode] Using full shot: {tmin:.3f} - {tmax:.3f} s")
            else:
                tmin = float(self.tmin_entry.text())
                tmax = float(self.tmax_entry.text())
                # Adjust if exceeds actual range
                if tmin < actual_tmin:
                    print(f"[n-Mode] tmin adjusted: {tmin:.3f} -> {actual_tmin:.3f} s")
                    tmin = actual_tmin
                    self.tmin_entry.setText(f"{tmin:.3f}")

                if tmax > actual_tmax:
                    print(f"[n-Mode] tmax adjusted: {tmax:.3f} -> {actual_tmax:.3f} s")
                    tmax = actual_tmax
                    self.tmax_entry.setText(f"{tmax:.3f}")

            self.status_label.setText('Calculating FFT...')
            self.status_label.setStyleSheet('color: blue; font-weight: bold;')
            QApplication.processEvents()

            self.fft_result = None
            gc.collect()

            self.fft_result = calculate_fft(self.mirnov_data, t_interval, integrate, run_detrend, tmin, tmax)
            gc.collect()

            self.status_label.setText('Calculating modes...')
            self.status_label.setStyleSheet('color: blue; font-weight: bold;')
            QApplication.processEvents()

            self.mode_result, self.freq_use, self.amp_evolution = calculate_mode_numbers(
                self.fft_result, fmin, fmax, tol, nmodes, frac, msign, integrate
            )
            gc.collect()

            self._update_plot()

            elapsed = tclock.time() - t0
            self.status_label.setText(f'Done ({elapsed:.1f}s)')
            self.status_label.setStyleSheet('color: green; font-weight: bold;')
            self.run_button.setEnabled(True)
            print(f"[n-Mode] Total calculation completed in {elapsed:.2f}s")

            # Enable save button
            self.save_button.setEnabled(True)

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed: {str(e)}")
            self.status_label.setText('Error')
            self.status_label.setStyleSheet('color: red; font-weight: bold;')
            self.run_button.setEnabled(True)

    def _update_plot(self):
        """Update plot with current settings"""
        if self.fft_result is None or self.mode_result is None:
            return

        try:
            tmin = float(self.tmin_entry.text())
            tmax = float(self.tmax_entry.text())
            fmin = float(self.fmin_entry.text())
            fmax = float(self.fmax_entry.text())
            nmodes = self._nmodes_value
            msign = self._msign_value
            plot_type = self._plot_type_value
            numc = int(self.contour_levels_entry.text())
            integrate = self.integrate_checkbox.isChecked()

            # Adjust tmin/tmax to FFT result time range
            actual_tmin = self.fft_result.time[0]
            actual_tmax = self.fft_result.time[-1]

            if tmin < actual_tmin:
                tmin = actual_tmin
            if tmax > actual_tmax:
                tmax = actual_tmax

            mode_colors = self._get_mode_colors(nmodes)

            plot_mode_spectrum(self.ax1, self.fft_result, self.mode_result, self.freq_use,
                              self.amp_evolution, tmin, tmax, fmin, fmax, nmodes, msign,
                              integrate, plot_type, numc, colors=mode_colors,
                              legend_fontsize=self.legend_fontsize,
                              label_fontsize=self.label_fontsize,
                              tick_fontsize=self.tick_fontsize)
            plot_amplitude_evolution(self.ax2, self.fft_result, self.amp_evolution,
                                    tmin, tmax, nmodes, msign, colors=mode_colors,
                                    legend_fontsize=self.legend_fontsize,
                                    label_fontsize=self.label_fontsize,
                                    tick_fontsize=self.tick_fontsize)

            self.canvas.draw()
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Plot failed: {str(e)}")

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.text(),
            "tmin": self.tmin_entry.text(),
            "tmax": self.tmax_entry.text(),
            "use_full_shot": self.use_full_shot_checkbox.isChecked(),
            "tinterval": self.tinterval_entry.text(),
            "fmin": self.fmin_entry.text(),
            "fmax": self.fmax_entry.text(),
            "nmodes": self._nmodes_value,
            "tolerance": self.tol_entry.text(),
            "fraction": self.frac_entry.text(),
            "sign": self._msign_value,
            "integrate": self.integrate_checkbox.isChecked(),
            "detrend": self.detrend_checkbox.isChecked(),
            "plot_type": self._plot_type_value,
            "contour_levels": self.contour_levels_entry.text(),
            "color_mode": self.color_mode,
            "label_fontsize": self.label_fontsize,
            "legend_fontsize": self.legend_fontsize,
            "tick_fontsize": self.tick_fontsize,
        }
        set_tab_settings("nmode", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("nmode")

        if settings.get("shot"):
            self.shot_entry.setText(str(settings["shot"]))

        if settings.get("tmin"):
            self.tmin_entry.setText(str(settings["tmin"]))

        if settings.get("tmax"):
            self.tmax_entry.setText(str(settings["tmax"]))

        if settings.get("use_full_shot") is not None:
            self.use_full_shot_checkbox.setChecked(bool(settings["use_full_shot"]))
            self._toggle_time_entries()

        if settings.get("tinterval"):
            self.tinterval_entry.setText(str(settings["tinterval"]))

        if settings.get("fmin"):
            self.fmin_entry.setText(str(settings["fmin"]))

        if settings.get("fmax"):
            self.fmax_entry.setText(str(settings["fmax"]))

        if settings.get("nmodes") is not None:
            val = int(settings["nmodes"])
            self._nmodes_value = val
            self.nmodes_slider.setValue(val)
            self.nmodes_label.setText(str(val))

        if settings.get("tolerance"):
            self.tol_entry.setText(str(settings["tolerance"]))

        if settings.get("fraction"):
            self.frac_entry.setText(str(settings["fraction"]))

        if settings.get("sign") is not None:
            self._msign_value = int(settings["sign"])
            for button in self.msign_button_group.buttons():
                if button.property('sign_value') == self._msign_value:
                    button.setChecked(True)
                    break

        if settings.get("integrate") is not None:
            self.integrate_checkbox.setChecked(bool(settings["integrate"]))

        if settings.get("detrend") is not None:
            self.detrend_checkbox.setChecked(bool(settings["detrend"]))

        if settings.get("plot_type"):
            self._plot_type_value = settings["plot_type"]
            for button in self.plot_type_button_group.buttons():
                if button.property('plot_type_value') == self._plot_type_value:
                    button.setChecked(True)
                    break
            self._toggle_contour_levels()

        if settings.get("contour_levels"):
            self.contour_levels_entry.setText(str(settings["contour_levels"]))

        if settings.get("color_mode"):
            self.color_mode = settings["color_mode"]
        if settings.get("label_fontsize"):
            self.label_fontsize = int(settings["label_fontsize"])
        if settings.get("legend_fontsize"):
            self.legend_fontsize = int(settings["legend_fontsize"])
        if settings.get("tick_fontsize"):
            self.tick_fontsize = int(settings["tick_fontsize"])


class NModePythonHighlighter(QSyntaxHighlighter):
    """Simple Python syntax highlighter for the n-mode example script dialog"""

    def __init__(self, document):
        super().__init__(document)
        import re
        self._rules = []

        # Keywords (purple)
        kw_format = QTextCharFormat()
        kw_format.setForeground(QColor('#c670e0'))
        keywords = r'\b(import|from|as|def|class|if|elif|else|for|while|try|except|with|return|yield|lambda|and|or|not|in|is|None|True|False|print|continue)\b'
        self._rules.append((re.compile(keywords), kw_format))

        # Builtins (orange)
        builtin_format = QTextCharFormat()
        builtin_format.setForeground(QColor('#fab16c'))
        builtins = r'\b(np|plt|data|metadata|time|frequency|mode_spectrum|amplitude|fig|ax1|ax2|im)\b'
        self._rules.append((re.compile(builtins), builtin_format))

        # Numbers (yellow)
        num_format = QTextCharFormat()
        num_format.setForeground(QColor('#faed5c'))
        numbers = r'\b(\d+\.?\d*)\b'
        self._rules.append((re.compile(numbers), num_format))

        # Strings (green)
        str_format = QTextCharFormat()
        str_format.setForeground(QColor('#b0e686'))
        strings = r'(\"[^\"]*\"|\'[^\']*\'|f\"[^\"]*\"|f\'[^\']*\')'
        self._rules.append((re.compile(strings), str_format))

        # Comments (gray)
        comment_format = QTextCharFormat()
        comment_format.setForeground(QColor('#999999'))
        comments = r'(#[^\n]*)'
        self._rules.append((re.compile(comments), comment_format))

    def highlightBlock(self, text):
        for pattern, fmt in self._rules:
            for match in pattern.finditer(text):
                start = match.start()
                length = match.end() - match.start()
                self.setFormat(start, length, fmt)
