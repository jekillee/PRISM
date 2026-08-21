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

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, apply_dark_figure_style, apply_shot_arrow_icons,
    save_file_async,
)
from config.user_settings import get_tab_settings, set_tab_settings
from core.nmode import (
    NModeConfig, MirnovData, FFTResult,
    load_mirnov_data, calculate_fft, calculate_mode_numbers,
    compute_amp_evolution,
)



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

# Amplitude evolution (2*nmodes, time) [Gauss], both reductions:
#   amplitude_max = peak amplitude bin in the [fmin, fmax] band (raw)
#   amplitude_sum = band sum with a size-3 moving average over time
# Older files (before both were saved) only have 'amplitude' (= max).
amplitude_max = data['amplitude_max'] if 'amplitude_max' in data.files else data['amplitude']
amplitude_sum = data['amplitude_sum'] if 'amplitude_sum' in data.files else None

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

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

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

# ax2 = Max reduction, ax3 = Sum reduction (matching colors)
for i, n in enumerate(n_modes_list):
    color = get_mode_color(abs(n))
    ax2.plot(time, amplitude_max[i], color=color, label=f'n={n}')
    if amplitude_sum is not None:
        ax3.plot(time, amplitude_sum[i], color=color, label=f'n={n}')

ax2.set_ylabel('Amplitude (Max) [Gauss]')
ax2.set_xlim(tmin, tmax)
ax2.legend(loc='upper right', fontsize=8)
ax2.grid(True, alpha=0.3)

ax3.set_xlabel('Time [s]')
ax3.set_ylabel('Amplitude (Sum) [Gauss]')
ax3.set_xlim(tmin, tmax)
ax3.legend(loc='upper right', fontsize=8)
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
'''



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
                       legend_fontsize=8, label_fontsize=None, tick_fontsize=None,
                       title_fontsize=None, contour_linewidth=0.8, selected_modes=None):
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
        if selected_modes and abs(n_val) not in selected_modes:
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
                ax.contour(x, y, thismode.T, levels, colors=[color], linewidths=contour_linewidth)
            except Exception:
                ax.contour(x, y, thismode.T, numc, colors=[color], linewidths=contour_linewidth)
        else:
            from matplotlib.colors import ListedColormap
            thismode_masked = np.ma.masked_where(thismode == 0, thismode)
            # Single flat color for this mode (masked zeros are transparent)
            mono_cmap = ListedColormap([color])
            ax.pcolormesh(x, y, thismode_masked.T, shading='auto',
                         cmap=mono_cmap, alpha=1.0)

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

    ax.set_title(f'#{fft_result.shot}', fontsize=title_fontsize)

    # Show legend for selected modes
    legend_elements = []
    from matplotlib.lines import Line2D
    show_modes = selected_modes if selected_modes else list(range(1, nmodes + 1))
    for n in show_modes:
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
                              colors=None, legend_fontsize=8, label_fontsize=None, tick_fontsize=None,
                              amp_linewidth=1.5, selected_modes=None, amp_mode='max'):
    """Plot amplitude evolution for each mode"""
    ax.clear()

    # y-label reflects the reduction: peak bin (Max) vs band-sum (Sum)
    ylabel = f"Amplitude ({'Sum' if amp_mode == 'sum' else 'Max'}) [Gauss]"

    time = fft_result.time
    j1 = np.argmin(np.abs(time - tmin))
    j2 = np.argmin(np.abs(time - tmax))
    time_plot = time[j1:j2]

    show_modes = selected_modes if selected_modes else list(range(1, nmodes + 1))

    for i in range(nmodes):
        n = i + 1
        if n not in show_modes:
            continue
        if msign in [0, 1, 2]:
            amp_pos = amp_evolution[i, j1:j2]
            c = colors[n - 1] if colors and n - 1 < len(colors) else get_mode_color(n)
            ax.plot(time_plot, amp_pos, color=c,
                   linewidth=amp_linewidth, label=f'n=+{n}')
        if msign in [0, -1, 2]:
            amp_neg = amp_evolution[i + nmodes, j1:j2]
            ls = '--' if msign == 2 else '-'
            c = colors[n - 1] if colors and n - 1 < len(colors) else get_mode_color(n)
            ax.plot(time_plot, amp_neg, color=c,
                   linewidth=amp_linewidth, linestyle=ls, label=f'n=-{n}')

    if label_fontsize:
        ax.set_xlabel('Time [s]', fontsize=label_fontsize)
        ax.set_ylabel(ylabel, fontsize=label_fontsize)
    else:
        ax.set_xlabel('Time [s]')
        ax.set_ylabel(ylabel)
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
        self.title_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10
        self.contour_linewidth = 0.8
        self.amp_linewidth = 1.5

    def create_widgets(self):
        """Create n-mode spectrum tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)
        self.figure.subplots_adjust(left=0.10, right=0.97, top=0.92, bottom=0.10, hspace=0.15)
        self.ax1 = self.figure.add_subplot(211)
        self.ax2 = self.figure.add_subplot(212, sharex=self.ax1)
        self.ax1.set_ylabel('Frequency [kHz]', fontsize=self.label_fontsize)
        self.ax1.set_label('Frequency [kHz]')
        self.ax2.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        self.ax2.set_ylabel('Amplitude (Max) [Gauss]', fontsize=self.label_fontsize)
        self.ax2.set_label('Amplitude (Max) [Gauss]')
        for ax in [self.ax1, self.ax2]:
            ax.tick_params(labelsize=self.tick_fontsize)
            import matplotlib as mpl
            ax.grid(ls='--', lw=0.3, c=mpl.rcParams.get('grid.color', '#444444'))
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
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_updown_layout.addWidget(up_btn)
        down_btn = QPushButton()
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_updown_layout.addWidget(down_btn)
        apply_shot_arrow_icons(up_btn, down_btn)
        shot_layout.addWidget(btn_updown)

        shot_layout.addStretch()

        grid.addWidget(shot_widget, row, 1)

        row += 1

        # Time range (same layout as Freq)
        grid.addWidget(QLabel('Time [s]'), row, 0)

        time_widget = QWidget()
        time_layout = QHBoxLayout(time_widget)
        time_layout.setContentsMargins(0, 0, 0, 0)

        self.tmin_entry = QLineEdit()
        self.tmin_entry.setText(str(NModeConfig.DEFAULT_TMIN))
        self.tmin_entry.setFixedWidth(80)
        time_layout.addWidget(self.tmin_entry)

        time_layout.addWidget(QLabel('-'))

        self.tmax_entry = QLineEdit()
        self.tmax_entry.setText(str(NModeConfig.DEFAULT_TMAX))
        self.tmax_entry.setFixedWidth(80)
        time_layout.addWidget(self.tmax_entry)
        time_layout.addStretch()

        grid.addWidget(time_widget, row, 1, 1, 3)

        row += 1

        # Use full shot checkbox
        self.use_full_shot_checkbox = QCheckBox('Use full shot length')
        self.use_full_shot_checkbox.setChecked(False)
        self.use_full_shot_checkbox.stateChanged.connect(self._toggle_time_entries)
        grid.addWidget(self.use_full_shot_checkbox, row, 1, 1, 3)

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

        plot_options_btn = QPushButton('Style')
        plot_options_btn.clicked.connect(self._show_plot_options_dialog)
        grid.addWidget(plot_options_btn, row, 3)
        row += 1

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        grid.addWidget(separator, row, 0, 1, 4)
        row += 1

        # n-modes checkboxes
        grid.addWidget(QLabel('n-modes'), row, 0)

        nmodes_widget = QWidget()
        nmodes_layout = QHBoxLayout(nmodes_widget)
        nmodes_layout.setContentsMargins(0, 0, 0, 0)

        self._selected_modes = [1, 2, 3, 4, 5]
        self.mode_checkboxes = {}
        for n in range(1, 6):
            cb = QCheckBox(str(n))
            cb.setChecked(True)
            cb.stateChanged.connect(self._on_modes_changed)
            nmodes_layout.addWidget(cb)
            self.mode_checkboxes[n] = cb
        nmodes_layout.addStretch()

        grid.addWidget(nmodes_widget, row, 1, 1, 3)
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

        # Amplitude reduction (Max / Sum across the frequency band). This toggle
        # controls only what the plot shows; the saved NPZ always stores BOTH
        # reductions (amplitude_max and amplitude_sum).
        grid.addWidget(QLabel('Amplitude'), row, 0)

        amp_mode_widget = QWidget()
        amp_mode_layout = QHBoxLayout(amp_mode_widget)
        amp_mode_layout.setContentsMargins(0, 0, 0, 0)

        self.amp_mode_button_group = QButtonGroup()
        self._amp_mode_value = 'max'

        radio_amp_max = QRadioButton('Max')
        radio_amp_max.setProperty('amp_mode_value', 'max')
        radio_amp_max.setChecked(True)
        amp_mode_layout.addWidget(radio_amp_max)
        self.amp_mode_button_group.addButton(radio_amp_max)

        radio_amp_sum = QRadioButton('Sum')
        radio_amp_sum.setProperty('amp_mode_value', 'sum')
        amp_mode_layout.addWidget(radio_amp_sum)
        self.amp_mode_button_group.addButton(radio_amp_sum)
        amp_mode_layout.addStretch()

        self.amp_mode_button_group.buttonClicked.connect(self._on_amp_mode_changed)

        grid.addWidget(amp_mode_widget, row, 1, 1, 3)
        row += 1

        # Equalize control-row heights so vertical spacing is uniform: line edits
        # are naturally taller than radio/checkbox rows, which otherwise makes the
        # gaps between n-modes / Plot type / Contour levels / Amplitude uneven.
        row_h = self.contour_levels_entry.sizeHint().height()
        for w in (nmodes_widget, plot_type_widget, self.contour_levels_entry, amp_mode_widget):
            w.setFixedHeight(row_h)

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
        self._example_script_popup = popup
        popup.show()

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
        nmodes = 5  # Always export all 5 modes

        # Default filename
        default_name = f"nmode_{shot}_{tmin:.1f}-{tmax:.1f}s.npz"

        # Non-modal save dialog; the write runs in the callback so the main
        # window stays movable while the dialog is open.
        save_file_async(
            self.frame, 'Save N-Mode Spectrum Data',
            os.path.join(os.path.expanduser("~"), default_name),
            'NumPy NPZ (*.npz)',
            lambda filepath: self._write_npz(filepath, shot, tmin, tmax, fmin, fmax, nmodes),
        )

    def _write_npz(self, filepath, shot, tmin, tmax, fmin, fmax, nmodes):
        """Write the n-mode NPZ to `filepath` (non-modal save-dialog callback)."""
        if not filepath:
            return

        try:
            # Build n_modes list based on sign setting
            msign = self._msign_value
            sign_str = {0: 'abs', 1: 'pos', -1: 'neg', 2: 'all'}

            # n_modes_list must follow the physical amp_evolution row order
            # ([+1..+nmodes, -1..-nmodes]) so amplitude[i] <-> n_modes_list[i] holds
            # for the saved NPZ (and the bundled example reader). Mirror in batch._n_modes_list.
            if msign == 0:  # abs
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == 1:  # pos
                n_modes_list = list(range(1, nmodes + 1))
            elif msign == -1:  # neg
                n_modes_list = list(range(-1, -nmodes - 1, -1))
            else:  # all
                n_modes_list = list(range(1, nmodes + 1)) + list(range(-1, -nmodes - 1, -1))

            # Filter mode_spectrum and amplitude by sign setting. Both the Max and
            # Sum amplitude traces are sliced identically (same row layout
            # [+1..+n, -1..-n]) so amplitude_max[i] / amplitude_sum[i] both map to
            # n_modes_list[i].
            def _slice_amp(amp):
                if msign == 1:      # pos only
                    return amp[:nmodes, :]
                if msign == -1:     # neg only
                    return amp[nmodes:, :]
                return amp          # abs / all keep both halves

            # The plot shows only the currently-selected reduction, but the NPZ
            # stores BOTH: compute Max and Sum from the current FFT/mode result at
            # save time (the mode map is independent of the reduction).
            integrate = self.integrate_checkbox.isChecked()
            amp_max_full = compute_amp_evolution(
                self.fft_result.amp, self.freq_use, self.mode_result,
                fmin, fmax, nmodes, integrate, 'max')
            amp_sum_full = compute_amp_evolution(
                self.fft_result.amp, self.freq_use, self.mode_result,
                fmin, fmax, nmodes, integrate, 'sum')

            mode_save = self.mode_result.copy()
            amp_max_save = _slice_amp(amp_max_full)
            amp_sum_save = _slice_amp(amp_sum_full)

            if msign == 1:  # pos only
                mode_save = np.where(mode_save > 0, mode_save, 0)
            elif msign == -1:  # neg only
                mode_save = np.where(mode_save < 0, mode_save, 0)
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

            # Save to NPZ. Both amplitude reductions are stored (amplitude_max =
            # peak bin, amplitude_sum = band sum). `amplitude` is kept as an alias
            # of amplitude_max for backward compatibility with older readers.
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=self.fft_result.time,
                     frequency=self.freq_use,
                     mode_spectrum=mode_save,
                     amplitude=amp_max_save,
                     amplitude_max=amp_max_save,
                     amplitude_sum=amp_sum_save
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

        # Title font size
        title_row = QHBoxLayout()
        title_row.addWidget(QLabel("Title font size"))
        title_spin = QSpinBox()
        title_spin.setFixedWidth(WIDGET_WIDTH)
        title_spin.setRange(6, 24)
        title_spin.setValue(self.title_fontsize)
        title_row.addWidget(title_spin)
        dlg_layout.addLayout(title_row)

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

        # Separator
        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        dlg_layout.addWidget(sep)

        # Contour line width
        lw_row = QHBoxLayout()
        lw_row.addWidget(QLabel("Contour line width"))
        from PySide6.QtWidgets import QDoubleSpinBox
        lw_spin = QDoubleSpinBox()
        lw_spin.setFixedWidth(WIDGET_WIDTH)
        lw_spin.setRange(0.1, 5.0)
        lw_spin.setSingleStep(0.1)
        lw_spin.setDecimals(1)
        lw_spin.setValue(self.contour_linewidth)
        lw_row.addWidget(lw_spin)
        dlg_layout.addLayout(lw_row)

        # Amplitude line width
        alw_row = QHBoxLayout()
        alw_row.addWidget(QLabel("Amplitude line width"))
        alw_spin = QDoubleSpinBox()
        alw_spin.setFixedWidth(WIDGET_WIDTH)
        alw_spin.setRange(0.1, 5.0)
        alw_spin.setSingleStep(0.1)
        alw_spin.setDecimals(1)
        alw_spin.setValue(self.amp_linewidth)
        alw_row.addWidget(alw_spin)
        dlg_layout.addLayout(alw_row)

        # Default / OK / Cancel
        btn_box = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)

        def reset_defaults():
            color_combo.setCurrentText("Default")
            label_spin.setValue(12)
            title_spin.setValue(12)
            legend_spin.setValue(8)
            tick_spin.setValue(10)
            lw_spin.setValue(0.8)
            alw_spin.setValue(1.5)
        btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(reset_defaults)
        dlg_layout.addWidget(btn_box)

        def _apply():
            self.color_mode = color_combo.currentText()
            self.label_fontsize = label_spin.value()
            self.title_fontsize = title_spin.value()
            self.legend_fontsize = legend_spin.value()
            self.tick_fontsize = tick_spin.value()
            self.contour_linewidth = lw_spin.value()
            self.amp_linewidth = alw_spin.value()
            # Auto-replot if data exists
            if self.fft_result is not None and self.mode_result is not None:
                self._update_plot()
        dialog.accepted.connect(_apply)
        self._style_dialog = dialog
        dialog.show()

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
        is_contour = self._plot_type_value == 'contour'
        self.contour_levels_label.setEnabled(is_contour)
        self.contour_levels_entry.setEnabled(is_contour)
        self.contour_levels_entry.setStyleSheet(
            "" if is_contour else "background-color: rgba(128,128,128,0.3);"
        )

    def _on_modes_changed(self):
        """Update selected modes list from checkboxes"""
        self._selected_modes = [n for n in range(1, 6) if self.mode_checkboxes[n].isChecked()]

    def _on_msign_changed(self, button):
        """Handle sign radio button change"""
        self._msign_value = button.property('sign_value')

    def _on_amp_mode_changed(self, button):
        """Handle Max/Sum amplitude reduction change: recompute and replot.

        Only affects the displayed panel; the saved NPZ always contains both
        reductions (see _save_data).
        """
        self._amp_mode_value = button.property('amp_mode_value')

        # Recompute the amplitude evolution with the selected reduction without
        # re-running the FFT / mode fit, then replot.
        if self.fft_result is None or self.mode_result is None:
            return
        try:
            fmin = float(self.fmin_entry.text())
            fmax = float(self.fmax_entry.text())
            integrate = self.integrate_checkbox.isChecked()
            self.amp_evolution = compute_amp_evolution(
                self.fft_result.amp, self.freq_use, self.mode_result,
                fmin, fmax, 5, integrate, self._amp_mode_value
            )
            self._update_plot()
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed: {str(e)}")

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
            nmodes = 5  # Always calculate all 5 modes
            frac = float(self.frac_entry.text())
            msign = self._msign_value
            integrate = self.integrate_checkbox.isChecked()
            run_detrend = self.detrend_checkbox.isChecked()

            self.run_button.setEnabled(False)
            self._set_status('Running...', 'blue')

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

                self._set_status('Loading data...', 'blue')

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

            self._set_status('Calculating FFT...', 'blue')

            self.fft_result = None
            gc.collect()

            self.fft_result = calculate_fft(self.mirnov_data, t_interval, integrate, run_detrend, tmin, tmax)
            gc.collect()

            self._set_status('Calculating modes...', 'blue')

            self.mode_result, self.freq_use, self.amp_evolution = calculate_mode_numbers(
                self.fft_result, fmin, fmax, tol, nmodes, frac, msign, integrate,
                self._amp_mode_value
            )
            gc.collect()

            self._update_plot()

            elapsed = tclock.time() - t0
            self._set_status(f'Done ({elapsed:.1f}s)', 'green')
            self.run_button.setEnabled(True)
            print(f"[n-Mode] Total calculation completed in {elapsed:.2f}s")
            print("[n-Mode] " + "=" * 60 + "\n")

            # Enable save button
            self.save_button.setEnabled(True)

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed: {str(e)}")
            self._set_status('Error', 'red')
            self.run_button.setEnabled(True)

    def _set_status(self, text, color='blue', timeout=0):
        from ui.ui_constants import show_status
        show_status(self.frame, "n-Mode", text, color, timeout)

    def _update_plot(self):
        """Update plot with current settings"""
        if self.fft_result is None or self.mode_result is None:
            return

        try:
            tmin = float(self.tmin_entry.text())
            tmax = float(self.tmax_entry.text())
            fmin = float(self.fmin_entry.text())
            fmax = float(self.fmax_entry.text())
            nmodes = 5  # Always use all 5 modes for data access
            selected_modes = self._selected_modes
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
                              tick_fontsize=self.tick_fontsize,
                              title_fontsize=self.title_fontsize,
                              contour_linewidth=self.contour_linewidth,
                              selected_modes=selected_modes)
            plot_amplitude_evolution(self.ax2, self.fft_result, self.amp_evolution,
                                    tmin, tmax, nmodes, msign, colors=mode_colors,
                                    legend_fontsize=self.legend_fontsize,
                                    label_fontsize=self.label_fontsize,
                                    tick_fontsize=self.tick_fontsize,
                                    amp_linewidth=self.amp_linewidth,
                                    selected_modes=selected_modes,
                                    amp_mode=self._amp_mode_value)

            self.canvas.draw()

            # Reset the navigation history so the toolbar's Home button returns
            # to the view just drawn. Without this, after replotting (e.g.
            # switching Max <-> Sum, which rescales the amplitude axis) Home
            # would jump back to the stale view captured for the previous data.
            if self.toolbar is not None:
                try:
                    self.toolbar.update()
                    self.toolbar.push_current()
                except Exception:
                    pass
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
            "selected_modes": self._selected_modes,
            "tolerance": self.tol_entry.text(),
            "fraction": self.frac_entry.text(),
            "sign": self._msign_value,
            "integrate": self.integrate_checkbox.isChecked(),
            "detrend": self.detrend_checkbox.isChecked(),
            "plot_type": self._plot_type_value,
            "contour_levels": self.contour_levels_entry.text(),
            "amp_mode": self._amp_mode_value,
            "color_mode": self.color_mode,
            "label_fontsize": self.label_fontsize,
            "title_fontsize": self.title_fontsize,
            "legend_fontsize": self.legend_fontsize,
            "tick_fontsize": self.tick_fontsize,
            "contour_linewidth": self.contour_linewidth,
            "amp_linewidth": self.amp_linewidth,
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

        if settings.get("selected_modes") is not None:
            self._selected_modes = [int(m) for m in settings["selected_modes"]]
            for n in range(1, 6):
                self.mode_checkboxes[n].setChecked(n in self._selected_modes)

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

        if settings.get("amp_mode"):
            self._amp_mode_value = settings["amp_mode"]
            for button in self.amp_mode_button_group.buttons():
                if button.property('amp_mode_value') == self._amp_mode_value:
                    button.setChecked(True)
                    break

        if settings.get("color_mode"):
            self.color_mode = settings["color_mode"]
        if settings.get("label_fontsize"):
            self.label_fontsize = int(settings["label_fontsize"])
        if settings.get("title_fontsize"):
            self.title_fontsize = int(settings["title_fontsize"])
        if settings.get("legend_fontsize"):
            self.legend_fontsize = int(settings["legend_fontsize"])
        if settings.get("tick_fontsize"):
            self.tick_fontsize = int(settings["tick_fontsize"])
        if settings.get("contour_linewidth"):
            self.contour_linewidth = float(settings["contour_linewidth"])
        if settings.get("amp_linewidth"):
            self.amp_linewidth = float(settings["amp_linewidth"])


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
