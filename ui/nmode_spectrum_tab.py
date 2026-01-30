#!/usr/bin/python3.8

"""
N-Mode Spectrum tab for toroidal mode number analysis
Based on KSTARtMirnovAnalysis_module.py, integrated into PRISM
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import gc
import os
from scipy.signal import detrend
from scipy.fft import fft
import scipy.ndimage as ndimage
import time as tclock
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from MDSplus import Connection

from ui.ui_constants import CONTROL_PANEL_WIDTH, PAD_X, PAD_Y
from config.user_settings import get_tab_settings, set_tab_settings


# Example script for loading n-mode spectrum NPZ file
NMODE_EXAMPLE_SCRIPT = '''"""
Example script for loading and plotting PRISM n-mode spectrum NPZ file
"""

import numpy as np
import matplotlib.pyplot as plt

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

# Define mode colors (matching PRISM)
def get_mode_color(n):
    colors = ['#000000', '#FF0000', '#00FF00', '#0000FF', '#FF9100',
              '#1f78b4', '#b2df8a', '#33a02c', '#ffffb3', '#b3bada']
    return colors[abs(n) % len(colors)]

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
    """Get year from shot number (KSTAR shot ranges)"""
    if shot > 37741:
        year = 2025
    elif shot > 34836:
        year = 2024
    elif shot > 32768:
        year = 2023
    elif shot > 30445:
        year = 2022
    elif shot > 27400:
        year = 2021
    elif shot > 24081:
        year = 2020
    elif shot > 21758:
        year = 2019
    elif shot > 19396:
        year = 2018
    elif shot > 17376:
        year = 2017
    elif shot > 14407:
        year = 2016
    elif shot > 11724:
        year = 2015
    elif shot > 9427:
        year = 2014
    elif shot > 8354:
        year = 2013
    elif shot > 6470:
        year = 2012
    elif shot > 4468:
        year = 2011
    elif shot > 2342:
        year = 2010
    elif shot > 1283:
        year = 2009
    else:
        year = 2008
    return year


def get_mirnov_config(year):
    """Get Mirnov coil configuration for given year"""
    if year == 2010:
        signs = [ 1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1,  0,
                  1,  1,  1,  1,  1,  1,
                  1,  0]
    elif year == 2011:
        signs = [ 1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1, -1,
                 -1, -1, -1, -1, -1, -1,
                 -1,  0]
    elif year == 2012:
        signs = [ 1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1, -1,
                 -1, -1, -1, -1, -1, -1,
                 -1,  1]
    elif year == 2013:
        signs = [ 1,  1,  1,  1,  1,  1,
                  1,  1,  1,  1,  1, -1,
                 -1, -1, -1, -1, -1, -1,
                 -1,  0]
    elif year == 2014:
        signs = [ 0,  1,  1,  1,  1,  1,
                  0,  0,  0,  1,  0, -1,
                 -1, -1, -1, -1,  0, -1,
                 -1,  0]
    elif year == 2015:
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  1, -1,  1, -1,
                 -1, -1, -1, -1,  0,  0,
                 -1, -1]
    elif year == 2016:
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  1,  0,  1, -1,
                 -1, -1, -1, -1,  0,  1,
                 -1, -1]
    elif year == 2017:
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  0, -1,  1, -1,
                 -1, -1, -1, -1,  0,  0,
                 -1,  0]
    elif year == 2023 or year == 2024:
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  0, -1,  0, -1,
                 -1, -1, -1, -1,  0,  0,
                 -1,  0]
    elif year == 2025:
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  0, -1,  0, -1,
                 -1, -1, -1,  0,  0,  0,
                 -1,  0]
    else:
        # Default for 2018-2022 and other years
        signs = [ 0,  1,  1,  1,  1,  1,
                  1,  1,  0, -1,  1, -1,
                 -1, -1, -1, -1,  0,  0,
                 -1,  0]
    
    names = ['MC1T01', 'MC1T02', 'MC1T03', 'MC1T04', 'MC1T05', 'MC1T06',
             'MC1T07', 'MC1T08', 'MC1T09', 'MC1T10', 'MC1T11', 'MC1T12',
             'MC1T13', 'MC1T14', 'MC1T15', 'MC1T16', 'MC1T17', 'MC1T18',
             'MC1T19', 'MC1T20']
    
    angles = [  1.60,  20.35,  35.35,  49.30,  70.50,  91.60,
              110.35, 132.50, 142.70, 160.50, 181.60, 200.35,
              215.35, 229.30, 257.30, 271.60, 290.35, 312.50,
              319.30, 343.90]
    
    if year >= 2017:
        names[9] = 'MC1P03'
    
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
    
    print(f"Loading Mirnov data for shot #{shot} (year {year})...")
    print(f"  Time range: {tmin:.2f} - {tmax:.2f} s")
    print(f"  Channels ({n_channels}): {config['names']}")
    print(f"  Threads: {n_channels}")
    
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
            bar = '█' * int(40 * progress) + '░' * (40 - int(40 * progress))
            print(f'\r  [{bar}] {completed}/{n_channels}', end='', flush=True)
    
    print(f'\n  Completed in {tclock.time() - t0:.2f} sec')
    
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

    print(f"  Common time range: {common_tmin:.4f} - {common_tmax:.4f} s")

    # Use the channel with finest sampling as reference
    ref_result = min(valid_results, key=lambda r: (r['time'][-1] - r['time'][0]) / r['n_points'])
    ref_time = ref_result['time']

    # Trim reference time to common range
    ref_mask = (ref_time >= common_tmin) & (ref_time <= common_tmax)
    time_arr = ref_time[ref_mask]
    n_points = len(time_arr)

    print(f"  Reference channel: {ref_result['name']} ({n_points} points in common range)")

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

    print(f"  Aligned {len(signals)} channels to common time grid")
    
    mirnov = MirnovData()
    mirnov.time = np.array(time_arr, dtype=np.float32)
    mirnov.signals = np.array(signals, dtype=np.float32)
    mirnov.positions = np.array(valid_angles)
    mirnov.channel_names = valid_names
    mirnov.shot = shot
    
    print(f"  Loaded {len(valid_names)} channels, {len(time_arr)} points")
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
    
    print(f"FFT: fs={fs/1e6:.2f}MHz, windows={n_windows}, pts/win={points_per_window}")
    
    # Step 1: Preprocessing
    print('  1/4 Preprocessing...')
    if integrate:
        for i in range(n_channels):
            signals[i] = np.cumsum(signals[i]) * dt * 1e4
    
    n_use = n_windows * points_per_window
    signals = signals[:, :n_use]
    time_trimmed = time[:n_use]
    signals_reshaped = signals.reshape(n_channels, n_windows, points_per_window)
    
    # Step 2: Detrend
    print('  2/4 Detrending...')
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
    print('  3/4 FFT...')
    n_freq_half = points_per_window // 2 + 1
    overp = 1.0 / float(points_per_window)
    
    amp = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)
    phase = np.zeros((n_windows, n_freq_half, n_channels), dtype=np.float32)
    
    fft_result_raw = fft(signals_reshaped, axis=2, workers=-1)
    fft_half = fft_result_raw[:, :, :n_freq_half]
    amp = np.transpose(np.abs(fft_half) * overp, (1, 2, 0)).astype(np.float32)
    phase = np.transpose(np.angle(fft_half), (1, 2, 0)).astype(np.float32)
    
    # Step 4: Post-processing
    print('  4/4 Finalizing...')
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
    
    print(f'\n  Completed in {tclock.time() - t0:.2f} sec')
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
    
    print(f"Mode calculation completed in {tclock.time() - t0:.2f} sec")
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
                       plot_type='contour', numc=50):
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
        
        color = get_mode_color(abs(n_val))
        
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
                ax.contour(x, y, thismode.T, levels, colors=color, linewidths=0.8)
            except Exception:
                ax.contour(x, y, thismode.T, numc, colors=color, linewidths=0.8)
        else:
            thismode_masked = np.ma.masked_where(thismode == 0, thismode)
            ax.pcolormesh(x, y, thismode_masked.T, shading='auto', alpha=0.7)
    
    ax.set_ylabel('Frequency [kHz]')
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(fmin, fmax)
    ax.set_axisbelow(False)
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)
    
    ax.set_title(f'#{fft_result.shot}')
    
    # Always show all modes in legend based on nmodes setting
    legend_elements = []
    from matplotlib.lines import Line2D
    for n in range(1, nmodes + 1):
        color = get_mode_color(n)
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
        ax.legend(handles=legend_elements, loc='upper right', fontsize=8, 
                  framealpha=0.8, ncol=min(len(legend_elements), 5))
    
    return lv


def plot_amplitude_evolution(ax, fft_result, amp_evolution, tmin, tmax, nmodes, msign):
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
            ax.plot(time_plot, amp_pos, color=get_mode_color(n), 
                   linewidth=1.5, label=f'n=+{n}')
        if msign in [0, -1, 2]:
            amp_neg = amp_evolution[i + nmodes, j1:j2]
            ls = '--' if msign == 2 else '-'
            ax.plot(time_plot, amp_neg, color=get_mode_color(n),
                   linewidth=1.5, linestyle=ls, label=f'n=-{n}')
    
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Amplitude [Gauss]')
    ax.set_xlim(tmin, tmax)
    ax.set_ylim(bottom=0)
    ax.legend(loc='upper right', fontsize=8, framealpha=0.8, ncol=5)
    ax.grid(True, linestyle='--', linewidth=0.5, color='gray', alpha=0.5)


# =============================================================================
# N-Mode Spectrum Tab
# =============================================================================

class NModeSpectrumTab:
    """N-Mode Spectrum analysis tab for PRISM"""
    
    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        
        self.frame = ttk.Frame(parent)
        self.toolbar = None
        
        self.mirnov_data = None
        self.fft_result = None
        self.mode_result = None
        self.freq_use = None
        self.amp_evolution = None
        self.current_shot = None
    
    def create_widgets(self):
        """Create n-mode spectrum tab widgets"""
        self.figure = Figure(NModeConfig.FIGURE_SIZE, tight_layout=True)
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212, sharex=self.ax1)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True)
        
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)
        
        self._create_parameters_panel(control_frame)
        self._create_run_plot_panel(control_frame)
        self._create_save_controls(control_frame)
        
        # Load saved settings
        self.load_settings()
    
    # Shared column 0 width for consistent alignment across panels
    LABEL_COLUMN_WIDTH = 105

    def _create_parameters_panel(self, parent):
        """Create parameters panel"""
        frame = ttk.LabelFrame(parent, text="1. Parameters", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(3, weight=1)
        
        row = 0
        
        # Shot with up/down buttons in same row
        tk.Label(frame, text='Shot').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        
        shot_frame = ttk.Frame(frame)
        shot_frame.grid(row=row, column=1, padx=5, pady=5, sticky='w')
        
        self.shot_var = tk.StringVar(value=str(NModeConfig.DEFAULT_SHOT))
        self.shot_entry = ttk.Entry(shot_frame, textvariable=self.shot_var, width=10)
        self.shot_entry.pack(side=tk.LEFT)
        self.shot_entry.bind('<Return>', lambda e: self._run_calculation())
        
        ttk.Button(shot_frame, text='\u25B2', width=2, 
                   command=lambda: self._adjust_shot(1)).pack(side=tk.LEFT, padx=(2, 0))
        ttk.Button(shot_frame, text='\u25BC', width=2, 
                   command=lambda: self._adjust_shot(-1)).pack(side=tk.LEFT)
        
        row += 1
        
        # Time range
        tk.Label(frame, text='Time [s]').grid(row=row, column=0, padx=5, pady=5, sticky='nw')
        time_frame = tk.Frame(frame)
        time_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='w')

        # Time inputs with grid layout
        self.tmin_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TMIN))
        self.tmin_entry = tk.Entry(time_frame, textvariable=self.tmin_var, width=10)
        self.tmin_entry.grid(row=0, column=0, padx=(0, 2))
        tk.Label(time_frame, text='-').grid(row=0, column=1)
        self.tmax_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TMAX))
        self.tmax_entry = tk.Entry(time_frame, textvariable=self.tmax_var, width=10)
        self.tmax_entry.grid(row=0, column=2, padx=(2, 0))

        # Use full shot checkbox (below tmin entry)
        self.use_full_shot_var = tk.BooleanVar(value=False)
        tk.Checkbutton(time_frame, text='Use full shot length',
                       variable=self.use_full_shot_var,
                       command=self._toggle_time_entries).grid(
            row=1, column=0, columnspan=3, sticky='w')
        row += 1
        
        # Separator
        ttk.Separator(frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=5)
        row += 1
        
        # Time interval
        tk.Label(frame, text='Time interval [s]').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.tinterval_var = tk.StringVar(value='0.01')
        tk.Entry(frame, textvariable=self.tinterval_var, width=10).grid(
            row=row, column=1, padx=5, pady=5, sticky='w')
        row += 1
        
        # Frequency range
        tk.Label(frame, text='Freq [kHz]').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        freq_frame = tk.Frame(frame)
        freq_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='w')
        self.fmin_var = tk.StringVar(value=str(NModeConfig.DEFAULT_FMIN))
        tk.Entry(freq_frame, textvariable=self.fmin_var, width=10).pack(side=tk.LEFT, padx=(0, 2))
        tk.Label(freq_frame, text='-').pack(side=tk.LEFT)
        self.fmax_var = tk.StringVar(value=str(NModeConfig.DEFAULT_FMAX))
        tk.Entry(freq_frame, textvariable=self.fmax_var, width=10).pack(side=tk.LEFT, padx=(2, 0))
        row += 1
        
        # n-modes
        tk.Label(frame, text='n-modes').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        nmodes_frame = tk.Frame(frame)
        nmodes_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='ew')
        self.nmodes_var = tk.IntVar(value=NModeConfig.DEFAULT_NMODES)
        ttk.Scale(nmodes_frame, from_=1, to=8, variable=self.nmodes_var,
                  orient='horizontal', command=self._on_nmodes_changed).pack(
                      side=tk.LEFT, fill='x', expand=True)
        self.nmodes_label = tk.Label(nmodes_frame, text=str(NModeConfig.DEFAULT_NMODES), width=3)
        self.nmodes_label.pack(side=tk.RIGHT, padx=(5, 0))
        row += 1
        
        # Tolerance
        tk.Label(frame, text='Tolerance').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.tol_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TOL))
        tk.Entry(frame, textvariable=self.tol_var, width=10).grid(row=row, column=1, padx=5, pady=5, sticky='w')
        row += 1
        
        # Fraction
        tk.Label(frame, text='Fraction').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.frac_var = tk.StringVar(value=str(NModeConfig.DEFAULT_FRAC))
        tk.Entry(frame, textvariable=self.frac_var, width=10).grid(row=row, column=1, padx=5, pady=5, sticky='w')
        row += 1
        
        # Sign
        tk.Label(frame, text='Sign').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        sign_frame = tk.Frame(frame)
        sign_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='w')
        self.msign_var = tk.IntVar(value=NModeConfig.DEFAULT_MSIGN)
        tk.Radiobutton(sign_frame, text='abs', variable=self.msign_var, value=0).pack(side=tk.LEFT, padx=3)
        tk.Radiobutton(sign_frame, text='pos', variable=self.msign_var, value=1).pack(side=tk.LEFT, padx=3)
        tk.Radiobutton(sign_frame, text='neg', variable=self.msign_var, value=-1).pack(side=tk.LEFT, padx=3)
        tk.Radiobutton(sign_frame, text='all', variable=self.msign_var, value=2).pack(side=tk.LEFT, padx=3)
        row += 1
        
        # Options
        option_frame = tk.Frame(frame)
        option_frame.grid(row=row, column=0, columnspan=4, padx=5, pady=5, sticky='w')
        self.integrate_var = tk.BooleanVar(value=NModeConfig.DEFAULT_INTEGRATE)
        tk.Checkbutton(option_frame, text='Integrate (dB/dt -> B)',
                       variable=self.integrate_var).pack(side=tk.LEFT, padx=5)
        self.detrend_var = tk.BooleanVar(value=NModeConfig.DEFAULT_DETREND)
        tk.Checkbutton(option_frame, text='Detrend',
                       variable=self.detrend_var).pack(side=tk.LEFT, padx=5)
    
    def _create_run_plot_panel(self, parent):
        """Create run and plot options panel"""
        frame = ttk.LabelFrame(parent, text="2. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        # Match Parameters panel column structure
        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(3, weight=1)

        row = 0

        # Calculate and Plot button
        self.run_button = ttk.Button(frame, text='Calculate and Plot', command=self._run_calculation)
        self.run_button.grid(row=row, column=0, columnspan=4, padx=5, pady=5, sticky='ew')
        row += 1

        # Status
        self.status_label = tk.Label(frame, text='Ready', fg='gray', font=('TkDefaultFont', 9, 'bold'))
        self.status_label.grid(row=row, column=0, columnspan=4, padx=5, pady=2, sticky='w')
        row += 1

        # Separator
        ttk.Separator(frame, orient='horizontal').grid(row=row, column=0, columnspan=4, sticky='ew', pady=5)
        row += 1

        # Plot type (radio buttons)
        tk.Label(frame, text='Plot type').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        plot_type_frame = tk.Frame(frame)
        plot_type_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='w')
        self.plot_type_var = tk.StringVar(value=NModeConfig.DEFAULT_PLOT_TYPE)
        tk.Radiobutton(plot_type_frame, text='contour', variable=self.plot_type_var,
                       value='contour', command=self._toggle_contour_levels).pack(side=tk.LEFT)
        tk.Radiobutton(plot_type_frame, text='imshow', variable=self.plot_type_var,
                       value='imshow', command=self._toggle_contour_levels).pack(side=tk.LEFT, padx=(5, 0))
        row += 1

        # Contour levels
        self.contour_levels_label = tk.Label(frame, text='Contour levels')
        self.contour_levels_label.grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.numc_var = tk.StringVar(value=str(NModeConfig.DEFAULT_NUMC))
        self.contour_levels_entry = tk.Entry(frame, textvariable=self.numc_var, width=10)
        self.contour_levels_entry.grid(row=row, column=1, padx=5, pady=5, sticky='w')
        row += 1

        # Update Plot button
        ttk.Button(frame, text='Update Plot', command=self._update_plot).grid(
            row=row, column=0, columnspan=4, padx=5, pady=5, sticky='ew')
    
    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = ttk.LabelFrame(parent, text="3. Save Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        self.save_button = ttk.Button(btn_frame, text='Save as NPZ',
                                       command=self._save_data, state='disabled')
        self.save_button.pack(side=tk.LEFT, expand=True, fill='x', padx=(0, 2))

        ttk.Button(btn_frame, text='Example Script',
                   command=self._show_example_script).pack(side=tk.LEFT, expand=True, fill='x', padx=(2, 0))
    
    def _show_example_script(self):
        """Show example script for loading NPZ file with syntax highlighting"""
        script = NMODE_EXAMPLE_SCRIPT
        
        popup = tk.Toplevel(self.frame)
        popup.title("Example Script - N-Mode Spectrum NPZ")
        popup.geometry("700x550")
        
        text_frame = tk.Frame(popup)
        text_frame.pack(fill='both', expand=True, padx=10, pady=10)
        
        scrollbar_y = tk.Scrollbar(text_frame, orient='vertical')
        scrollbar_y.pack(side=tk.RIGHT, fill='y')
        
        scrollbar_x = tk.Scrollbar(text_frame, orient='horizontal')
        scrollbar_x.pack(side=tk.BOTTOM, fill='x')
        
        text_widget = tk.Text(text_frame, wrap='none', font=('Courier', 10),
                              yscrollcommand=scrollbar_y.set,
                              xscrollcommand=scrollbar_x.set,
                              bg='#19232d', fg='#ffffff',
                              insertbackground='white')
        text_widget.pack(side=tk.LEFT, fill='both', expand=True)
        
        scrollbar_y.config(command=text_widget.yview)
        scrollbar_x.config(command=text_widget.xview)
        
        # Define syntax highlighting tags (Spyder dark theme)
        text_widget.tag_configure('keyword', foreground='#c670e0')
        text_widget.tag_configure('builtin', foreground='#fab16c')
        text_widget.tag_configure('string', foreground='#b0e686')
        text_widget.tag_configure('comment', foreground='#999999')
        text_widget.tag_configure('number', foreground='#faed5c')
        
        text_widget.insert('1.0', script)
        
        # Apply syntax highlighting
        self._apply_syntax_highlighting(text_widget)
        
        text_widget.config(state='disabled')
        
        btn_frame = tk.Frame(popup)
        btn_frame.pack(fill='x', padx=10, pady=(0, 10))
        
        def copy_to_clipboard():
            popup.clipboard_clear()
            popup.clipboard_append(script)
            messagebox.showinfo("Copied", "Script copied to clipboard")
        
        tk.Button(btn_frame, text="Copy to Clipboard", command=copy_to_clipboard).pack(side=tk.LEFT)
        tk.Button(btn_frame, text="Close", command=popup.destroy).pack(side=tk.RIGHT)
    
    def _apply_syntax_highlighting(self, text_widget):
        """Apply Python syntax highlighting to text widget"""
        import re
        
        content = text_widget.get('1.0', 'end')
        
        # Keywords
        keywords = r'\b(import|from|as|def|class|if|elif|else|for|while|try|except|with|return|yield|lambda|and|or|not|in|is|None|True|False|print|continue)\b'
        # Builtins
        builtins = r'\b(np|plt|data|metadata|time|frequency|mode_spectrum|amplitude|fig|ax1|ax2|im)\b'
        # Strings
        strings = r'(\"\"\"[\s\S]*?\"\"\"|\'\'\'[\s\S]*?\'\'\'|\"[^\"]*\"|\'[^\']*\'|f\"[^\"]*\"|f\'[^\']*\')'
        # Comments
        comments = r'(#[^\n]*)'
        # Numbers
        numbers = r'\b(\d+\.?\d*)\b'
        
        patterns = [
            (comments, 'comment'),
            (strings, 'string'),
            (keywords, 'keyword'),
            (builtins, 'builtin'),
            (numbers, 'number'),
        ]
        
        for pattern, tag in patterns:
            for match in re.finditer(pattern, content):
                start_idx = f"1.0+{match.start()}c"
                end_idx = f"1.0+{match.end()}c"
                text_widget.tag_add(tag, start_idx, end_idx)
    
    def _save_data(self):
        """Save n-mode spectrum data to NPZ file"""
        if self.fft_result is None or self.mode_result is None:
            messagebox.showwarning("Warning", "No data to save. Run calculation first.")
            return
        
        # Get parameters
        shot = self.current_shot
        tmin = float(self.tmin_var.get())
        tmax = float(self.tmax_var.get())
        fmin = float(self.fmin_var.get())
        fmax = float(self.fmax_var.get())
        nmodes = int(self.nmodes_var.get())
        
        # Default filename
        default_name = f"nmode_{shot}_{tmin:.1f}-{tmax:.1f}s.npz"
        
        # Hide hidden files in file dialog (Linux Tk)
        try:
            self.frame.tk.call('tk_getOpenFile', '-foption')
        except:
            pass
        self.frame.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')

        # File dialog
        filepath = filedialog.asksaveasfilename(
            initialdir=os.path.expanduser("~"),
            defaultextension='.npz',
            filetypes=[('NumPy NPZ', '*.npz')],
            initialfile=default_name,
            title='Save N-Mode Spectrum Data'
        )
        
        if not filepath:
            return
        
        try:
            # Build n_modes list based on sign setting
            msign = self.msign_var.get()
            sign_str = {0: 'abs', 1: 'pos', -1: 'neg', 2: 'all'}
            
            if msign == 0:  # abs
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == 1:  # pos
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == -1:  # neg
                n_modes_list = list(range(-nmodes, 0))
            else:  # all
                n_modes_list = list(range(-nmodes, 0)) + list(range(1, nmodes + 1))
            
            # Prepare metadata
            metadata = {
                'shot': shot,
                'time_range': [tmin, tmax],
                'freq_range': [fmin, fmax],
                'time_interval': float(self.tinterval_var.get()),
                'n_modes': nmodes,
                'n_modes_list': n_modes_list,
                'tolerance': float(self.tol_var.get()),
                'fraction': float(self.frac_var.get()),
                'sign': sign_str.get(msign, 'pos'),
                'integrate': self.integrate_var.get(),
                'detrend': self.detrend_var.get(),
            }
            
            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=self.fft_result.time,
                     frequency=self.freq_use,
                     mode_spectrum=self.mode_result,
                     amplitude=self.amp_evolution
            )
            
            print(f"N-mode spectrum data saved to: {filepath}")
            messagebox.showinfo("Saved", f"Data saved to:\n{filepath}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {str(e)}")
    
    def _toggle_time_entries(self):
        """Enable/disable time entry fields based on checkbox state"""
        if self.use_full_shot_var.get():
            self.tmin_entry.config(state='disabled')
            self.tmax_entry.config(state='disabled')
        else:
            self.tmin_entry.config(state='normal')
            self.tmax_entry.config(state='normal')

    def _toggle_contour_levels(self):
        """Enable/disable contour levels entry based on plot type"""
        if self.plot_type_var.get() == 'contour':
            self.contour_levels_label.config(state='normal')
            self.contour_levels_entry.config(state='normal')
        else:
            self.contour_levels_label.config(state='disabled')
            self.contour_levels_entry.config(state='disabled')

    def _on_nmodes_changed(self, value):
        """Update nmodes label"""
        self.nmodes_label.config(text=str(int(float(value))))
    
    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_var.get())
            new_shot = max(1, current + delta)
            self.shot_var.set(str(new_shot))
        except ValueError:
            pass
    
    def _run_calculation(self):
        """Run complete workflow"""
        try:
            shot = int(self.shot_var.get())
            use_full_shot = self.use_full_shot_var.get()

            # When using full shot, use wide time range for initial load
            if use_full_shot:
                load_tmin = -1.0
                load_tmax = 100.0
            else:
                load_tmin = float(self.tmin_var.get())
                load_tmax = float(self.tmax_var.get())

            t_interval = float(self.tinterval_var.get())
            fmin = float(self.fmin_var.get())
            fmax = float(self.fmax_var.get())
            tol = float(self.tol_var.get())
            nmodes = int(self.nmodes_var.get())
            frac = float(self.frac_var.get())
            msign = self.msign_var.get()
            integrate = self.integrate_var.get()
            run_detrend = self.detrend_var.get()

            self.run_button.config(state='disabled')
            self.status_label.config(text='Running...', fg='blue')
            self.frame.update()

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

                self.status_label.config(text='Loading data...', fg='blue')
                self.frame.update()

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
                self.tmin_var.set(f"{tmin:.3f}")
                self.tmax_var.set(f"{tmax:.3f}")
                print(f"  Using full shot: {tmin:.3f} - {tmax:.3f} s")
            else:
                tmin = float(self.tmin_var.get())
                tmax = float(self.tmax_var.get())
                # Adjust if exceeds actual range
                if tmin < actual_tmin:
                    print(f"  tmin adjusted: {tmin:.3f} -> {actual_tmin:.3f} s")
                    tmin = actual_tmin
                    self.tmin_var.set(f"{tmin:.3f}")

                if tmax > actual_tmax:
                    print(f"  tmax adjusted: {tmax:.3f} -> {actual_tmax:.3f} s")
                    tmax = actual_tmax
                    self.tmax_var.set(f"{tmax:.3f}")
            
            self.status_label.config(text='Calculating FFT...', fg='blue')
            self.frame.update()
            
            self.fft_result = None
            gc.collect()
            
            self.fft_result = calculate_fft(self.mirnov_data, t_interval, integrate, run_detrend, tmin, tmax)
            gc.collect()
            
            self.status_label.config(text='Calculating modes...', fg='blue')
            self.frame.update()
            
            self.mode_result, self.freq_use, self.amp_evolution = calculate_mode_numbers(
                self.fft_result, fmin, fmax, tol, nmodes, frac, msign, integrate
            )
            gc.collect()
            
            self._update_plot()
            
            self.status_label.config(text='Done', fg='green')
            self.run_button.config(state='normal')
            
            # Enable save button
            self.save_button.config(state='normal')
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed: {str(e)}")
            self.status_label.config(text='Error', fg='red')
            self.run_button.config(state='normal')
    
    def _update_plot(self):
        """Update plot with current settings"""
        if self.fft_result is None or self.mode_result is None:
            return
        
        try:
            tmin = float(self.tmin_var.get())
            tmax = float(self.tmax_var.get())
            fmin = float(self.fmin_var.get())
            fmax = float(self.fmax_var.get())
            nmodes = int(self.nmodes_var.get())
            msign = self.msign_var.get()
            plot_type = self.plot_type_var.get()
            numc = int(self.numc_var.get())
            integrate = self.integrate_var.get()
            
            # Adjust tmin/tmax to FFT result time range
            actual_tmin = self.fft_result.time[0]
            actual_tmax = self.fft_result.time[-1]
            
            if tmin < actual_tmin:
                tmin = actual_tmin
            if tmax > actual_tmax:
                tmax = actual_tmax
            
            plot_mode_spectrum(self.ax1, self.fft_result, self.mode_result, self.freq_use,
                              self.amp_evolution, tmin, tmax, fmin, fmax, nmodes, msign,
                              integrate, plot_type, numc)
            plot_amplitude_evolution(self.ax2, self.fft_result, self.amp_evolution,
                                    tmin, tmax, nmodes, msign)
            
            self.figure.tight_layout()
            self.canvas.draw()
        except Exception as e:
            messagebox.showerror("Error", f"Plot failed: {str(e)}")
    
    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_var.get(),
            "tmin": self.tmin_var.get(),
            "tmax": self.tmax_var.get(),
            "use_full_shot": self.use_full_shot_var.get(),
            "tinterval": self.tinterval_var.get(),
            "fmin": self.fmin_var.get(),
            "fmax": self.fmax_var.get(),
            "nmodes": self.nmodes_var.get(),
            "tolerance": self.tol_var.get(),
            "fraction": self.frac_var.get(),
            "sign": self.msign_var.get(),
            "integrate": self.integrate_var.get(),
            "detrend": self.detrend_var.get(),
            "plot_type": self.plot_type_var.get(),
            "contour_levels": self.numc_var.get()
        }
        set_tab_settings("nmode", settings)
    
    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("nmode")
        
        if settings.get("shot"):
            self.shot_var.set(settings["shot"])
        
        if settings.get("tmin"):
            self.tmin_var.set(settings["tmin"])
        
        if settings.get("tmax"):
            self.tmax_var.set(settings["tmax"])

        if settings.get("use_full_shot") is not None:
            self.use_full_shot_var.set(settings["use_full_shot"])
            self._toggle_time_entries()

        if settings.get("tinterval"):
            self.tinterval_var.set(settings["tinterval"])
        
        if settings.get("fmin"):
            self.fmin_var.set(settings["fmin"])
        
        if settings.get("fmax"):
            self.fmax_var.set(settings["fmax"])
        
        if settings.get("nmodes") is not None:
            self.nmodes_var.set(settings["nmodes"])
            self.nmodes_label.config(text=str(settings["nmodes"]))
        
        if settings.get("tolerance"):
            self.tol_var.set(settings["tolerance"])
        
        if settings.get("fraction"):
            self.frac_var.set(settings["fraction"])
        
        if settings.get("sign") is not None:
            self.msign_var.set(settings["sign"])
        
        if settings.get("integrate") is not None:
            self.integrate_var.set(settings["integrate"])
        
        if settings.get("detrend") is not None:
            self.detrend_var.set(settings["detrend"])
        
        if settings.get("plot_type"):
            self.plot_type_var.set(settings["plot_type"])
            self._toggle_contour_levels()

        if settings.get("contour_levels"):
            self.numc_var.set(settings["contour_levels"])