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
    
    signals, valid_names, valid_angles = [], [], []
    time_arr = None
    results.sort(key=lambda x: x['angle'])
    
    for result in results:
        if result['error'] is None:
            signals.append(result['data'])
            valid_names.append(result['name'])
            valid_angles.append(result['angle'])
            if time_arr is None:
                time_arr = result['time']
    
    if len(signals) == 0:
        raise RuntimeError("No valid Mirnov channels found")
    
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
    
    sign_str = {0: 'abs', 1: 'pos', -1: 'neg', 2: 'all'}
    ax.set_title(f'KSTAR #{fft_result.shot} n-mode spectrum ({sign_str.get(msign, "pos")})')
    
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
    ax.set_title('Mode Amplitude Evolution (max)')
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
        self._create_plot_options_panel(control_frame)
        self._create_save_controls(control_frame)
    
    def _create_parameters_panel(self, parent):
        """Create parameters panel"""
        frame = tk.LabelFrame(parent, text="Parameters", font=('TkDefaultFont', 9, 'bold'))
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(3, weight=1)
        
        row = 0
        
        # Shot
        tk.Label(frame, text='Shot').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.shot_var = tk.StringVar(value=str(NModeConfig.DEFAULT_SHOT))
        self.shot_entry = tk.Entry(frame, textvariable=self.shot_var, width=12)
        self.shot_entry.grid(row=row, column=1, padx=5, pady=5, sticky='w')
        self.shot_entry.bind('<Return>', lambda e: self._run_calculation())
        row += 1
        
        # Time range
        tk.Label(frame, text='Time [s]').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        time_frame = tk.Frame(frame)
        time_frame.grid(row=row, column=1, columnspan=3, padx=5, pady=5, sticky='w')
        self.tmin_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TMIN))
        tk.Entry(time_frame, textvariable=self.tmin_var, width=8).pack(side=tk.LEFT, padx=(0, 2))
        tk.Label(time_frame, text='-').pack(side=tk.LEFT)
        self.tmax_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TMAX))
        tk.Entry(time_frame, textvariable=self.tmax_var, width=8).pack(side=tk.LEFT, padx=(2, 0))
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
        tk.Entry(freq_frame, textvariable=self.fmin_var, width=8).pack(side=tk.LEFT, padx=(0, 2))
        tk.Label(freq_frame, text='-').pack(side=tk.LEFT)
        self.fmax_var = tk.StringVar(value=str(NModeConfig.DEFAULT_FMAX))
        tk.Entry(freq_frame, textvariable=self.fmax_var, width=8).pack(side=tk.LEFT, padx=(2, 0))
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
        
        # Tolerance and Fraction
        tk.Label(frame, text='Tolerance').grid(row=row, column=0, padx=5, pady=5, sticky='w')
        self.tol_var = tk.StringVar(value=str(NModeConfig.DEFAULT_TOL))
        tk.Entry(frame, textvariable=self.tol_var, width=8).grid(row=row, column=1, padx=5, pady=5, sticky='w')
        tk.Label(frame, text='Fraction').grid(row=row, column=2, padx=5, pady=5, sticky='w')
        self.frac_var = tk.StringVar(value=str(NModeConfig.DEFAULT_FRAC))
        tk.Entry(frame, textvariable=self.frac_var, width=8).grid(row=row, column=3, padx=5, pady=5, sticky='w')
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
        row += 1
        
        # Run button
        self.run_button = ttk.Button(frame, text='Run', command=self._run_calculation)
        self.run_button.grid(row=row, column=0, columnspan=4, padx=5, pady=10, sticky='ew')
        row += 1
        
        # Status
        self.status_label = tk.Label(frame, text='Ready', fg='gray', font=('TkDefaultFont', 9, 'bold'))
        self.status_label.grid(row=row, column=0, columnspan=4, padx=5, pady=2, sticky='w')
    
    def _create_plot_options_panel(self, parent):
        """Create plot options panel"""
        frame = tk.LabelFrame(parent, text="Plot Options", font=('TkDefaultFont', 9, 'bold'))
        frame.pack(fill='x', padx=5, pady=5)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Plot type
        tk.Label(frame, text='Plot type:').grid(row=0, column=0, padx=5, pady=5, sticky='w')
        type_frame = tk.Frame(frame)
        type_frame.grid(row=0, column=1, padx=5, pady=5, sticky='w')
        self.plot_type_var = tk.StringVar(value=NModeConfig.DEFAULT_PLOT_TYPE)
        tk.Radiobutton(type_frame, text='contour', variable=self.plot_type_var, 
                       value='contour').pack(side=tk.LEFT, padx=5)
        tk.Radiobutton(type_frame, text='imshow', variable=self.plot_type_var, 
                       value='imshow').pack(side=tk.LEFT, padx=5)
        
        # Contour levels
        tk.Label(frame, text='Contour levels:').grid(row=1, column=0, padx=5, pady=5, sticky='w')
        self.numc_var = tk.StringVar(value=str(NModeConfig.DEFAULT_NUMC))
        tk.Entry(frame, textvariable=self.numc_var, width=10).grid(
            row=1, column=1, padx=5, pady=5, sticky='w')
        
        # Update button
        ttk.Button(frame, text='Update Plot', command=self._update_plot).grid(
            row=2, column=0, columnspan=2, padx=5, pady=10, sticky='ew')
    
    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = tk.LabelFrame(parent, text="Save Data", font=('TkDefaultFont', 9, 'bold'))
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        self.save_button = ttk.Button(frame, text='Save as NPZ', 
                                       command=self._save_data, state='disabled')
        self.save_button.pack(fill='x', padx=PAD_X, pady=PAD_Y)
    
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
    
    def _on_nmodes_changed(self, value):
        """Update nmodes label"""
        self.nmodes_label.config(text=str(int(float(value))))
    
    def _run_calculation(self):
        """Run complete workflow"""
        try:
            shot = int(self.shot_var.get())
            tmin = float(self.tmin_var.get())
            tmax = float(self.tmax_var.get())
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
                getattr(self, '_last_tmin', None) != tmin or
                getattr(self, '_last_tmax', None) != tmax
            )
            
            if need_reload:
                self.mirnov_data = None
                self.fft_result = None
                gc.collect()
                
                self.status_label.config(text='Loading data...', fg='blue')
                self.frame.update()
                
                self.mirnov_data = load_mirnov_data(shot, tmin, tmax)
                self.current_shot = shot
                self._last_tmin = tmin
                self._last_tmax = tmax
            
            # Adjust tmin/tmax if they exceed actual data range
            actual_tmin = self.mirnov_data.time[0]
            actual_tmax = self.mirnov_data.time[-1]
            
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