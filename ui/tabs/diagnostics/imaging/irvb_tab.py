"""
IRVB (Infra-Red Video Bolometer) tab
2D radiation profile with EFIT overlay and regional Prad time traces
"""

import os
import re
import time as tclock
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.gridspec import GridSpec
from scipy.interpolate import RectBivariateSpline

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QGroupBox,
    QCheckBox, QSlider, QSplitter, QMessageBox, QFileDialog,
    QApplication, QDialog, QTextEdit, QScrollBar, QStyle,
    QScrollArea, QSpinBox, QDialogButtonBox, QFrame,
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QFont, QColor, QTextCharFormat, QSyntaxHighlighter

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon,
    apply_shot_arrow_icons,
)
from config.user_settings import get_tab_settings, set_tab_settings


# Example script for loading IRVB NPZ file
IRVB_EXAMPLE_SCRIPT = '''"""
Example script for loading and plotting PRISM IRVB NPZ file
"""

import numpy as np
import matplotlib.pyplot as plt

# Load NPZ file
filepath = 'irvb_XXXXXX_0.0-10.0s.npz'
data = np.load(filepath, allow_pickle=True)

# Extract data
time = data['time']               # Time array [s]
R = data['R']                     # R coordinates [m]
Z = data['Z']                     # Z coordinates [m]
prad_2d = data['prad_2d']         # 2D Prad array (time, Z, R) [MW/m^3]
region_prad = data['region_prad'] # Regional Prad (n_regions, time) [MW]
ptot = data['ptot']               # Total radiated power [MW]

# EFIT data
efit_r_grid = data['efit_r_grid']       # EFIT R grid [m]
efit_z_grid = data['efit_z_grid']       # EFIT Z grid [m]
efit_psi_n = data['efit_psi_n']         # Normalized psi (time, Z, R)
efit_bdry_r = data['efit_bdry_r']       # LCFS R coordinates (time, max_points) [m]
efit_bdry_z = data['efit_bdry_z']       # LCFS Z coordinates (time, max_points) [m]
efit_nbdry = data['efit_nbdry']         # Number of boundary points per frame
efit_rmaxis = data['efit_rmaxis']       # Magnetic axis R (time) [m]
efit_zmaxis = data['efit_zmaxis']       # Magnetic axis Z (time) [m]
efit_limiter_r = data['efit_limiter_r'] # Limiter R [m]
efit_limiter_z = data['efit_limiter_z'] # Limiter Z [m]

# Get metadata
metadata = data['metadata'].item()
print(f"Shot: {metadata['shot']}")
print(f"EFIT tree: {metadata['efit_tree']}")
print(f"Psi boundaries: {metadata['psi_boundaries']}")
print(f"Region labels: {metadata['region_labels']}")
print(f"Total frames: {len(time)}")

# Plot 2D Prad at middle frame
frame_idx = len(time) // 2  # Middle frame
fig, ax = plt.subplots(figsize=(6, 8))

# Use same colormap settings as PRISM
vmin, vmax = 0.0, 1.5
levels = np.linspace(vmin, vmax, 101)

im = ax.contourf(R, Z, prad_2d[frame_idx], levels=levels)
plt.colorbar(im, ax=ax, label='Prad [MW/m$^3$]')

# Plot LCFS
nbdry = efit_nbdry[frame_idx]
ax.plot(efit_bdry_r[frame_idx, :nbdry], efit_bdry_z[frame_idx, :nbdry], 'k-', lw=2)

# Plot limiter
ax.plot(efit_limiter_r, efit_limiter_z, 'k--', lw=1)

# Plot magnetic axis
ax.plot(efit_rmaxis[frame_idx], efit_zmaxis[frame_idx], 'k+', markersize=10, mew=2)

ax.set_xlabel('R [m]')
ax.set_ylabel('Z [m]')
ax.set_title(f"Shot #{metadata['shot']} t={time[frame_idx]:.3f}s")
ax.set_aspect('equal')
plt.tight_layout()
plt.show()

# Plot regional Prad time traces with Ptot
# Note: region_prad shape is (n_regions, time)
fig, ax = plt.subplots(figsize=(10, 5))
region_labels = metadata['region_labels']
for i, label in enumerate(region_labels):
    ax.plot(time, region_prad[i, :], label=label)
ax.plot(time, ptot, 'k--', lw=1.5, label='Ptot')
ax.set_xlabel('Time [s]')
ax.set_ylabel('Prad [MW]')
ax.set_title(f"Shot #{metadata['shot']} Regional Prad")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
'''


class IRVBTab:
    """IRVB visualization tab"""

    # Default psi boundaries for region separation
    DEFAULT_PSI_BOUNDARIES = "0.7, 1.0"
    MAX_BOUNDARIES = 5

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90

    # 2D plot color settings
    PRAD_VMIN = 0.0
    PRAD_VMAX = 1.5
    PRAD_LEVELS = 101

    # Background psi contour levels (9 levels from 0.1 to 0.9)
    PSI_BG_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    # Region colors (tab10 colormap)
    REGION_COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']

    def __init__(self, parent, app_config, diagnostic_config, efit_loader):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        self.efit_loader = efit_loader

        self.frame = QWidget()
        self.toolbar = None

        # Data storage
        self.irvb_data = None
        self.efit_2d = None
        self.region_prad = None
        self.psi_boundaries = []
        self.shot_number = None
        self.ip_fault_time = None

        # Frame control
        self.current_frame = 0
        self.total_frames = 0

        # Playback state
        self.is_playing = False
        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_next_frame)
        self._last_frame_time = 0
        self._fps_history = []

        # Plot references for fast update
        self.im_2d = None
        self.contour_psi_bg = None
        self.contour_psi_bounds = []
        self.contour_bdry = None
        self.limiter_line = None
        self.maxis_marker = None
        self.time_vlines = []
        self.colorbar = None

        # Slider value tracking
        self._slider_value = 0

        # Plot Options settings - Time Trace
        self.trace_color_mode = "Fixed(tab10)"
        self.trace_label_fontsize = 12
        self.trace_legend_fontsize = 8
        self.trace_tick_fontsize = 10

        # Plot Options settings - 2D Plot
        self.plot2d_colormap = "viridis"
        self.plot2d_lcfs_color = "black"
        self.plot2d_maxis_color = "black"
        self.plot2d_limiter_color = "white"
        self.plot2d_flux_color = "gray"
        self.plot2d_label_fontsize = 12
        self.plot2d_tick_fontsize = 10

    def create_widgets(self):
        """Create IRVB tab widgets"""
        main_layout = QHBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)

        # Left: canvas area
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)

        # Main figure
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=False)
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()
        canvas_layout.addWidget(self.canvas)

        # Bind mouse wheel for frame navigation
        self.canvas.mpl_connect('scroll_event', self._on_mouse_wheel)

        splitter.addWidget(canvas_widget)

        # Right: Scrollable control panel
        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QVBoxLayout.SetNoConstraint)

        self._create_load_controls(control_layout)
        self._create_efit_controls(control_layout)
        self._create_plot_controls(control_layout)
        self._create_frame_controls(control_layout)
        self._create_playback_controls(control_layout)
        self._create_save_controls(control_layout)
        control_layout.addStretch()

        scroll_area.setWidget(control_widget)
        scroll_area.viewport().setAutoFillBackground(False)
        control_widget.setAutoFillBackground(False)
        splitter.addWidget(scroll_area)

        # Load saved settings
        self.load_settings()

    def _create_load_controls(self, parent_layout):
        """Create data loading section"""
        group = QGroupBox("1. Load IRVB Data")
        grid = QGridLayout(group)

        # Shot number with up/down buttons in same row
        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_layout = QHBoxLayout()

        self.shot_entry = QLineEdit()
        self.shot_entry.setMinimumWidth(80)
        self.shot_entry.returnPressed.connect(self._load_data)
        shot_layout.addWidget(self.shot_entry, 1)

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

        grid.addLayout(shot_layout, 0, 1)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self._load_data)
        grid.addWidget(self.fetch_button, 0, 2)

        # Status label (left-aligned)

        parent_layout.addWidget(group)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _create_efit_controls(self, parent_layout):
        """Create EFIT settings section with psi boundaries"""
        group = QGroupBox("2. EFIT Settings")
        grid = QGridLayout(group)

        # EFIT Tree selection
        grid.addWidget(QLabel("EFIT Tree"), 0, 0)

        efit_options = list(self.app_config.EFIT_TREES.keys())
        self.efit_tree_combo = QComboBox()
        self.efit_tree_combo.addItems(efit_options)
        self.efit_tree_combo.setCurrentIndex(0)
        grid.addWidget(self.efit_tree_combo, 0, 1)

        # Psi boundaries
        grid.addWidget(QLabel('psi bounds'), 1, 0)
        self.psi_entry = QLineEdit(self.DEFAULT_PSI_BOUNDARIES)
        grid.addWidget(self.psi_entry, 1, 1)

        # Hint label
        hint = QLabel("(comma-separated, max 5 values)")
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        grid.addWidget(hint, 2, 0, 1, 2)

        parent_layout.addWidget(group)

    def _create_plot_controls(self, parent_layout):
        """Create plot button section"""
        group = QGroupBox("3. Plot")
        layout = QVBoxLayout(group)

        # Plot + Option buttons in same row
        plot_row = QHBoxLayout()

        plot_btn = QPushButton('Plot')
        plot_btn.clicked.connect(self._plot_data)
        plot_row.addWidget(plot_btn, 3)

        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_plot_options_dialog)
        plot_row.addWidget(style_btn, 1)

        layout.addLayout(plot_row)

        parent_layout.addWidget(group)

    def _create_frame_controls(self, parent_layout):
        """Create frame navigation controls"""
        group = QGroupBox("4. Frame Control")
        layout = QVBoxLayout(group)

        # Frame slider with < > buttons
        slider_layout = QHBoxLayout()

        prev_btn = QPushButton()
        prev_btn.setIcon(get_icon(QStyle.SP_ArrowBack))
        prev_btn.setFixedWidth(30)
        prev_btn.clicked.connect(lambda: self._step_frame(-1))
        slider_layout.addWidget(prev_btn)

        self.frame_slider = QSlider(Qt.Horizontal)
        self.frame_slider.setMinimum(0)
        self.frame_slider.setMaximum(0)
        self.frame_slider.valueChanged.connect(self._on_slider_change)
        slider_layout.addWidget(self.frame_slider)

        next_btn = QPushButton()
        next_btn.setIcon(get_icon(QStyle.SP_ArrowForward))
        next_btn.setFixedWidth(30)
        next_btn.clicked.connect(lambda: self._step_frame(1))
        slider_layout.addWidget(next_btn)

        layout.addLayout(slider_layout)

        # Frame number input
        frame_input_layout = QHBoxLayout()
        frame_input_layout.addWidget(QLabel('Frame'))

        self.frame_entry = QLineEdit('1')
        self.frame_entry.setFixedWidth(60)
        self.frame_entry.returnPressed.connect(self._goto_frame)
        frame_input_layout.addWidget(self.frame_entry)

        frame_input_layout.addWidget(QLabel('/'))

        self.frame_total_entry = QLineEdit()
        self.frame_total_entry.setFixedWidth(60)
        self.frame_total_entry.setReadOnly(True)
        frame_input_layout.addWidget(self.frame_total_entry)

        go_btn = QPushButton()
        go_btn.setIcon(self.frame.style().standardIcon(QStyle.SP_DialogOkButton))
        go_btn.setFixedSize(24, 24)
        go_btn.setToolTip("Go to frame")
        go_btn.clicked.connect(self._goto_frame)
        frame_input_layout.addWidget(go_btn)
        frame_input_layout.addStretch()

        layout.addLayout(frame_input_layout)

        # Time input
        time_input_layout = QHBoxLayout()
        time_input_layout.addWidget(QLabel('Time [s]'))

        self.time_entry = QLineEdit('0.0')
        self.time_entry.setFixedWidth(75)
        self.time_entry.returnPressed.connect(self._goto_time)
        time_input_layout.addWidget(self.time_entry)

        go_time_btn = QPushButton()
        go_time_btn.setIcon(self.frame.style().standardIcon(QStyle.SP_DialogOkButton))
        go_time_btn.setFixedSize(24, 24)
        go_time_btn.setToolTip("Go to time")
        go_time_btn.clicked.connect(self._goto_time)
        time_input_layout.addStretch()

        layout.addLayout(time_input_layout)

        # Mouse wheel hint
        hint = QLabel("(Mouse wheel: navigate frames)")
        hint.setStyleSheet("color: gray; font-size: 8pt;")
        hint.setAlignment(Qt.AlignCenter)
        layout.addWidget(hint)

        parent_layout.addWidget(group)

    def _create_playback_controls(self, parent_layout):
        """Create playback control section"""
        group = QGroupBox("5. Playback Control")
        layout = QVBoxLayout(group)

        # Row 1: Speed dropdown | Play | Stop
        row1 = QHBoxLayout()

        self.speed_combo = QComboBox()
        self.speed_combo.addItems(['0.5x', '1x', 'Max'])
        self.speed_combo.setCurrentText('1x')
        row1.addWidget(self.speed_combo)

        self.play_btn = QPushButton()
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self.play_btn.clicked.connect(self._toggle_play)
        row1.addWidget(self.play_btn)

        stop_btn = QPushButton()
        stop_btn.setIcon(get_icon(QStyle.SP_MediaStop))
        stop_btn.clicked.connect(self._stop_play)
        row1.addWidget(stop_btn)

        layout.addLayout(row1)

        # Row 2: FPS + Loop
        row2 = QHBoxLayout()
        self.actual_fps_label = QLabel('Actual: -- FPS')
        self.actual_fps_label.setStyleSheet("color: gray;")
        row2.addWidget(self.actual_fps_label)
        self.loop_checkbox = QCheckBox('Loop')
        self.loop_checkbox.setChecked(True)
        row2.addWidget(self.loop_checkbox)
        row2.addStretch()
        layout.addLayout(row2)

        parent_layout.addWidget(group)

    def _create_save_controls(self, parent_layout):
        """Create save data section"""
        group = QGroupBox("6. Save Data")
        layout = QVBoxLayout(group)

        btn_layout = QHBoxLayout()

        self.save_button = QPushButton('Save as .npz')
        self.save_button.setEnabled(False)
        self.save_button.clicked.connect(self._save_data)
        btn_layout.addWidget(self.save_button)

        example_btn = QPushButton('Example Script')
        example_btn.clicked.connect(self._show_example_script)
        btn_layout.addWidget(example_btn)

        layout.addLayout(btn_layout)

        parent_layout.addWidget(group)

    def _show_plot_options_dialog(self):
        """Show Plot Options dialog with Time Trace and 2D Plot sections"""
        WIDGET_WIDTH = 150

        dialog = QDialog(self.frame)
        dialog.setWindowTitle("Plot Options")
        dialog.setMinimumWidth(320)
        dlg_layout = QVBoxLayout(dialog)

        # === Time Trace Section ===
        trace_group = QGroupBox("Time Trace")
        trace_layout = QVBoxLayout(trace_group)

        # Color mode
        color_row = QHBoxLayout()
        color_row.addWidget(QLabel("Color"))
        trace_color_combo = QComboBox()
        trace_color_combo.setFixedWidth(WIDGET_WIDTH)
        trace_color_combo.addItems([
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
        ])
        trace_color_combo.setCurrentText(self.trace_color_mode)
        color_row.addWidget(trace_color_combo)
        trace_layout.addLayout(color_row)

        # Label font size
        label_row = QHBoxLayout()
        label_row.addWidget(QLabel("Label font size"))
        trace_label_spin = QSpinBox()
        trace_label_spin.setFixedWidth(WIDGET_WIDTH)
        trace_label_spin.setRange(6, 24)
        trace_label_spin.setValue(self.trace_label_fontsize)
        label_row.addWidget(trace_label_spin)
        trace_layout.addLayout(label_row)

        # Legend font size
        legend_row = QHBoxLayout()
        legend_row.addWidget(QLabel("Legend font size"))
        trace_legend_spin = QSpinBox()
        trace_legend_spin.setFixedWidth(WIDGET_WIDTH)
        trace_legend_spin.setRange(4, 20)
        trace_legend_spin.setValue(self.trace_legend_fontsize)
        legend_row.addWidget(trace_legend_spin)
        trace_layout.addLayout(legend_row)

        # Tick font size
        tick_row = QHBoxLayout()
        tick_row.addWidget(QLabel("Tick font size"))
        trace_tick_spin = QSpinBox()
        trace_tick_spin.setFixedWidth(WIDGET_WIDTH)
        trace_tick_spin.setRange(6, 20)
        trace_tick_spin.setValue(self.trace_tick_fontsize)
        tick_row.addWidget(trace_tick_spin)
        trace_layout.addLayout(tick_row)

        dlg_layout.addWidget(trace_group)

        # === 2D Plot Section ===
        plot2d_group = QGroupBox("2D Plot")
        plot2d_layout = QVBoxLayout(plot2d_group)

        # Colorbar
        cbar_row = QHBoxLayout()
        cbar_row.addWidget(QLabel("Colorbar"))
        plot2d_cmap_combo = QComboBox()
        plot2d_cmap_combo.setFixedWidth(WIDGET_WIDTH)
        plot2d_cmap_combo.addItems([
            "viridis", "hot", "jet", "coolwarm", "inferno", "plasma", "magma", "cividis"
        ])
        plot2d_cmap_combo.setCurrentText(self.plot2d_colormap)
        cbar_row.addWidget(plot2d_cmap_combo)
        plot2d_layout.addLayout(cbar_row)

        # LCFS color
        lcfs_row = QHBoxLayout()
        lcfs_row.addWidget(QLabel("LCFS color"))
        lcfs_combo = QComboBox()
        lcfs_combo.setFixedWidth(WIDGET_WIDTH)
        lcfs_combo.addItems(["black", "white", "red", "blue", "green", "yellow", "cyan", "magenta"])
        lcfs_combo.setCurrentText(self.plot2d_lcfs_color)
        lcfs_row.addWidget(lcfs_combo)
        plot2d_layout.addLayout(lcfs_row)

        # Magnetic axis color
        maxis_row = QHBoxLayout()
        maxis_row.addWidget(QLabel("Mag. axis color"))
        maxis_combo = QComboBox()
        maxis_combo.setFixedWidth(WIDGET_WIDTH)
        maxis_combo.addItems(["black", "white", "red", "blue", "green", "yellow", "cyan", "magenta"])
        maxis_combo.setCurrentText(self.plot2d_maxis_color)
        maxis_row.addWidget(maxis_combo)
        plot2d_layout.addLayout(maxis_row)

        # Limiter color
        limiter_row = QHBoxLayout()
        limiter_row.addWidget(QLabel("Limiter color"))
        limiter_combo = QComboBox()
        limiter_combo.setFixedWidth(WIDGET_WIDTH)
        limiter_combo.addItems(["white", "black", "red", "blue", "green", "yellow", "cyan", "magenta", "gray"])
        limiter_combo.setCurrentText(self.plot2d_limiter_color)
        limiter_row.addWidget(limiter_combo)
        plot2d_layout.addLayout(limiter_row)

        # Flux contour color
        flux_row = QHBoxLayout()
        flux_row.addWidget(QLabel("Flux contour color"))
        flux_combo = QComboBox()
        flux_combo.setFixedWidth(WIDGET_WIDTH)
        flux_combo.addItems(["gray", "white", "black", "silver", "lightblue", "lightgreen"])
        flux_combo.setCurrentText(self.plot2d_flux_color)
        flux_row.addWidget(flux_combo)
        plot2d_layout.addLayout(flux_row)

        # Label font size
        p2d_label_row = QHBoxLayout()
        p2d_label_row.addWidget(QLabel("Label font size"))
        plot2d_label_spin = QSpinBox()
        plot2d_label_spin.setFixedWidth(WIDGET_WIDTH)
        plot2d_label_spin.setRange(6, 24)
        plot2d_label_spin.setValue(self.plot2d_label_fontsize)
        p2d_label_row.addWidget(plot2d_label_spin)
        plot2d_layout.addLayout(p2d_label_row)

        # Tick font size
        p2d_tick_row = QHBoxLayout()
        p2d_tick_row.addWidget(QLabel("Tick font size"))
        plot2d_tick_spin = QSpinBox()
        plot2d_tick_spin.setFixedWidth(WIDGET_WIDTH)
        plot2d_tick_spin.setRange(6, 20)
        plot2d_tick_spin.setValue(self.plot2d_tick_fontsize)
        p2d_tick_row.addWidget(plot2d_tick_spin)
        plot2d_layout.addLayout(p2d_tick_row)

        dlg_layout.addWidget(plot2d_group)

        # Default / OK / Cancel
        btn_box = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)

        def reset_defaults():
            trace_color_combo.setCurrentText("Fixed(tab10)")
            trace_label_spin.setValue(12)
            trace_legend_spin.setValue(8)
            trace_tick_spin.setValue(10)
            plot2d_cmap_combo.setCurrentText("viridis")
            lcfs_combo.setCurrentText("black")
            maxis_combo.setCurrentText("black")
            limiter_combo.setCurrentText("white")
            flux_combo.setCurrentText("gray")
            plot2d_label_spin.setValue(12)
            plot2d_tick_spin.setValue(10)
        btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(reset_defaults)
        dlg_layout.addWidget(btn_box)

        def _apply():
            self.trace_color_mode = trace_color_combo.currentText()
            self.trace_label_fontsize = trace_label_spin.value()
            self.trace_legend_fontsize = trace_legend_spin.value()
            self.trace_tick_fontsize = trace_tick_spin.value()
            self.plot2d_colormap = plot2d_cmap_combo.currentText()
            self.plot2d_lcfs_color = lcfs_combo.currentText()
            self.plot2d_maxis_color = maxis_combo.currentText()
            self.plot2d_limiter_color = limiter_combo.currentText()
            self.plot2d_flux_color = flux_combo.currentText()
            self.plot2d_label_fontsize = plot2d_label_spin.value()
            self.plot2d_tick_fontsize = plot2d_tick_spin.value()
            # Auto-replot if data exists
            if self.region_prad is not None:
                self._setup_figure()
                self._update_plot(self.current_frame)
        dialog.accepted.connect(_apply)
        self._style_dialog = dialog
        dialog.show()

    @staticmethod
    def _parse_color_mode(text):
        """Extract colormap name from dropdown text like 'Fixed(tab10)' → 'tab10'"""
        start = text.find('(')
        end = text.find(')')
        if start != -1 and end != -1:
            return text[start + 1:end]
        return 'tab10'

    def _get_trace_colors(self, n_regions):
        """Get colors for time trace regions based on color mode setting"""
        cmap_name = self._parse_color_mode(self.trace_color_mode)
        is_gradient = self.trace_color_mode.startswith("Gradient")

        try:
            cmap = plt.cm.get_cmap(cmap_name)
        except ValueError:
            cmap = plt.cm.get_cmap('tab10')
            is_gradient = False

        colors = []
        for i in range(n_regions):
            if is_gradient:
                colors.append(cmap(i / max(1, n_regions - 1)))
            else:
                if hasattr(cmap, 'N'):
                    colors.append(cmap(i % cmap.N))
                else:
                    colors.append(cmap(i / 10))

        return colors

    def _show_example_script(self):
        """Show example script for loading NPZ file with syntax highlighting"""
        script = IRVB_EXAMPLE_SCRIPT

        popup = QDialog(self.frame)
        popup.setWindowTitle("Example Script - IRVB NPZ")
        popup.resize(700, 600)

        main_layout = QVBoxLayout(popup)

        text_widget = QTextEdit()
        text_widget.setFont(QFont('Courier', 10))
        text_widget.setStyleSheet(
            "background-color: #19232d; color: #ffffff;"
        )
        text_widget.setPlainText(script)

        # Apply syntax highlighting
        self._apply_syntax_highlighting(text_widget)

        text_widget.setReadOnly(True)
        main_layout.addWidget(text_widget)

        btn_layout = QHBoxLayout()

        def copy_to_clipboard():
            clipboard = QApplication.clipboard()
            clipboard.setText(script)
            QMessageBox.information(popup, "Copied", "Script copied to clipboard")

        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(copy_to_clipboard)
        btn_layout.addWidget(copy_btn)

        btn_layout.addStretch()

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(popup.close)
        btn_layout.addWidget(close_btn)

        main_layout.addLayout(btn_layout)

        self._example_script_popup = popup
        popup.show()

    def _apply_syntax_highlighting(self, text_widget):
        """Apply Python syntax highlighting to QTextEdit"""
        content = text_widget.toPlainText()

        # Keywords
        keywords = r'\b(import|from|as|def|class|if|elif|else|for|while|try|except|with|return|yield|lambda|and|or|not|in|is|None|True|False|print)\b'
        # Builtins
        builtins = r'\b(np|plt|data|metadata|time|R|Z|prad_2d|region_prad|ptot|fig|ax|im)\b'
        # Strings
        strings = r'(\"\"\"[\s\S]*?\"\"\"|\'\'\'[\s\S]*?\'\'\'|\"[^\"]*\"|\'[^\']*\'|f\"[^\"]*\"|f\'[^\']*\')'
        # Comments
        comments = r'(#[^\n]*)'
        # Numbers
        numbers = r'\b(\d+\.?\d*)\b'

        patterns = [
            (comments, QColor('#999999')),
            (strings, QColor('#b0e686')),
            (keywords, QColor('#c670e0')),
            (builtins, QColor('#fab16c')),
            (numbers, QColor('#faed5c')),
        ]

        cursor = text_widget.textCursor()
        for pattern, color in patterns:
            fmt = QTextCharFormat()
            fmt.setForeground(color)
            for match in re.finditer(pattern, content):
                cursor.setPosition(match.start())
                cursor.setPosition(match.end(), cursor.KeepAnchor)
                cursor.setCharFormat(fmt)

    def _save_data(self):
        """Save IRVB data to NPZ file including EFIT data"""
        if self.irvb_data is None or self.region_prad is None:
            QMessageBox.warning(self.frame, "Warning", "No data to save. Plot data first.")
            return

        if self.efit_2d is None:
            QMessageBox.warning(self.frame, "Warning", "No EFIT data available.")
            return

        # Default filename
        t_min = self.irvb_data.time[0]
        t_max = self.irvb_data.time[-1]
        default_name = f"irvb_{self.shot_number}_{t_min:.1f}-{t_max:.1f}s.npz"

        # File dialog
        filepath, _ = QFileDialog.getSaveFileName(
            self.frame,
            'Save IRVB Data',
            os.path.join(os.path.expanduser("~"), default_name),
            'NumPy NPZ (*.npz)'
        )

        if not filepath:
            return

        try:
            # Build region labels
            n_regions = len(self.psi_boundaries) + 1
            boundaries = [0] + self.psi_boundaries + [np.inf]
            region_labels = []
            for i in range(n_regions):
                psi_min = boundaries[i]
                psi_max = boundaries[i + 1]
                if psi_max == np.inf:
                    region_labels.append(f'psi_N > {psi_min}')
                else:
                    region_labels.append(f'{psi_min} < psi_N < {psi_max}')

            # Prepare metadata
            metadata = {
                'shot': self.shot_number,
                'psi_boundaries': self.psi_boundaries,
                'efit_tree': self.efit_tree_combo.currentText(),
                'time_range': [t_min, t_max],
                'region_labels': region_labels,
            }

            # Prepare EFIT data for each IRVB time point
            n_frames = len(self.irvb_data.time)
            max_bdry = self.efit_2d.bdry_r.shape[1]

            efit_bdry_r = np.zeros((n_frames, max_bdry))
            efit_bdry_z = np.zeros((n_frames, max_bdry))
            efit_nbdry = np.zeros(n_frames, dtype=int)
            efit_rmaxis = np.zeros(n_frames)
            efit_zmaxis = np.zeros(n_frames)
            efit_psi_n = np.zeros((n_frames, len(self.efit_2d.z_grid), len(self.efit_2d.r_grid)))

            for i, t in enumerate(self.irvb_data.time):
                efit_idx = self.efit_2d.find_time_index(t)
                bdry_r, bdry_z = self.efit_2d.get_boundary(efit_idx)
                nbdry = len(bdry_r)

                efit_bdry_r[i, :nbdry] = bdry_r
                efit_bdry_z[i, :nbdry] = bdry_z
                efit_nbdry[i] = nbdry
                efit_rmaxis[i], efit_zmaxis[i] = self.efit_2d.get_magnetic_axis(efit_idx)
                efit_psi_n[i] = self.efit_2d.get_psi_normalized(efit_idx)

            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=self.irvb_data.time,
                     R=self.irvb_data.x_grid,
                     Z=self.irvb_data.y_grid,
                     prad_2d=self.irvb_data.recon,
                     region_prad=self.region_prad,
                     ptot=self.irvb_data.ptot[:len(self.irvb_data.time)],
                     # EFIT data
                     efit_r_grid=self.efit_2d.r_grid,
                     efit_z_grid=self.efit_2d.z_grid,
                     efit_psi_n=efit_psi_n,
                     efit_bdry_r=efit_bdry_r,
                     efit_bdry_z=efit_bdry_z,
                     efit_nbdry=efit_nbdry,
                     efit_rmaxis=efit_rmaxis,
                     efit_zmaxis=efit_zmaxis,
                     efit_limiter_r=self.efit_2d.limiter_r,
                     efit_limiter_z=self.efit_2d.limiter_z
            )

            print(f"[IRVB] Data saved to: {filepath}")
            QMessageBox.information(self.frame, "Saved", f"Data saved to:\n{filepath}")

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to save data: {str(e)}")

    def _parse_psi_boundaries(self):
        """Parse psi boundary input string"""
        try:
            text = self.psi_entry.text().strip()
            if not text:
                return []

            values = [float(v.strip()) for v in text.split(',')]

            # Validate count
            if len(values) > self.MAX_BOUNDARIES:
                QMessageBox.warning(self.frame, "Warning",
                    f"Maximum {self.MAX_BOUNDARIES} boundaries allowed. Using first {self.MAX_BOUNDARIES}.")
                values = values[:self.MAX_BOUNDARIES]

            # Sort and validate range
            values = sorted(values)
            if any(v <= 0 or v >= 2.0 for v in values):
                QMessageBox.critical(self.frame, "Error", "psi values should be between 0 and 2.0")
                return None

            return values

        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Invalid psi boundary format. Use comma-separated numbers.")
            return None

    def _update_status(self, message, color=None, success=None):
        """Update status label with message and color

        Args:
            message: Status text
            color: Direct color ('blue', 'green', 'red', 'gray')
            success: If provided, True='green', False='red' (legacy support)
        """
        if color is None:
            if success is not None:
                color = 'green' if success else 'red'
            else:
                color = 'blue'
        from ui.ui_constants import show_status
        show_status(self.frame, "IRVB", message, color)

    def _load_data(self):
        """Load IRVB data from server"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            self._update_status('Invalid shot number', success=False)
            return

        # Load IRVB data
        try:
            from data_loaders.irvb_loader import IRVBLoader
            loader = IRVBLoader(self.app_config)
            self.irvb_data = loader.load_data(shot_number)
            self.shot_number = shot_number

            # Update frame controls
            self.total_frames = len(self.irvb_data.time)
            self.current_frame = 0
            self.frame_slider.setMaximum(self.total_frames - 1)
            self.frame_slider.setValue(0)
            self.frame_total_entry.setReadOnly(False)
            self.frame_total_entry.setText(str(self.total_frames))
            self.frame_total_entry.setReadOnly(True)
            self.frame_entry.setText('1')
            self.time_entry.setText(f'{self.irvb_data.time[0]:.3f}')

            self._update_status(f'Shot #{shot_number}: {self.total_frames} frames loaded', success=True)
            print(f"[IRVB] Data loaded for shot #{shot_number}")

        except Exception as e:
            self._update_status(f'Failed: {str(e)[:30]}...', success=False)
            QMessageBox.critical(self.frame, "Error", f"Failed to load IRVB data: {str(e)}")
            return

    def _get_ip_fault_time(self):
        """Get IP fault time from MDS+"""
        try:
            from MDSplus import Connection
            mds = Connection(self.app_config.MDS_IP)
            mds.openTree('kstar', self.shot_number)
            self.ip_fault_time = mds.get('\\t_ip_fault').data()
            mds.closeTree('kstar', self.shot_number)
            print(f"[IRVB] IP fault time = {self.ip_fault_time:.3f} s")
        except Exception:
            self.ip_fault_time = None
            print("[IRVB] IP fault time not available")

    def _slice_by_ip_fault(self):
        """Slice IRVB data by IP fault time"""
        if self.ip_fault_time is None or self.irvb_data is None:
            return

        valid_mask = self.irvb_data.time < self.ip_fault_time
        if not np.any(valid_mask):
            return

        self.irvb_data.time = self.irvb_data.time[valid_mask]
        self.irvb_data.recon = self.irvb_data.recon[valid_mask]
        self.irvb_data.ptot = self.irvb_data.ptot[valid_mask]
        print(f"[IRVB] Data sliced to {len(self.irvb_data.time)} frames (before IP fault)")

    def _slice_by_efit_time(self):
        """Slice IRVB data by valid time range (0 to EFIT last time)"""
        if self.efit_2d is None or self.irvb_data is None:
            return

        efit_end = self.efit_2d.time[-1]
        valid_mask = (self.irvb_data.time >= 0) & (self.irvb_data.time <= efit_end)
        if not np.any(valid_mask):
            return

        n_before = len(self.irvb_data.time)
        self.irvb_data.time = self.irvb_data.time[valid_mask]
        self.irvb_data.recon = self.irvb_data.recon[valid_mask]
        self.irvb_data.ptot = self.irvb_data.ptot[valid_mask]

        if len(self.irvb_data.time) < n_before:
            print(f"[IRVB] Data sliced to {len(self.irvb_data.time)} frames (0 to EFIT end)")

    def _load_efit_data(self):
        """Load EFIT 2D equilibrium data using efit_loader"""
        if self.irvb_data is None:
            return False

        efit_display = self.efit_tree_combo.currentText()
        efit_tree = self.app_config.EFIT_TREES[efit_display]

        try:
            print(f"[IRVB] Loading 2D EFIT data from {efit_tree}...")
            self.efit_2d = self.efit_loader.load_efit_2d(self.shot_number, efit_tree)
            print(f"[IRVB] EFIT data loaded ({len(self.efit_2d.time)} timepoints)")
            return True

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load EFIT data: {str(e)}")
            return False

    def _find_xpoint(self, bdry_r, bdry_z):
        """Find X-point from boundary (lower divertor)"""
        lower_mask = bdry_z < 0
        if not np.any(lower_mask):
            return None, None

        lower_r = bdry_r[lower_mask]
        lower_z = bdry_z[lower_mask]

        min_idx = np.argmin(lower_z)
        return lower_r[min_idx], lower_z[min_idx]

    def _compute_regional_prad(self):
        """Compute Prad for each region defined by psi boundaries"""
        if self.irvb_data is None or self.efit_2d is None:
            return

        n_times = len(self.irvb_data.time)
        n_regions = len(self.psi_boundaries) + 1

        self.region_prad = np.zeros((n_regions, n_times))

        x_grid = self.irvb_data.x_grid
        y_grid = self.irvb_data.y_grid
        X, Y = np.meshgrid(x_grid, y_grid)

        # Volume element: 2*pi*R * dR * dZ
        dx = x_grid[1] - x_grid[0]
        dy = y_grid[1] - y_grid[0]
        vol_factor = 2 * np.pi * X * dx * dy

        print("[IRVB] Computing regional Prad...")

        efit_idx_prev = -1
        psi_on_grid = None

        for i, t in enumerate(self.irvb_data.time):
            # Find closest EFIT time
            efit_idx = self.efit_2d.find_time_index(t)

            # Update psi map if EFIT time changed
            if efit_idx != efit_idx_prev:
                psi_n_2d = self.efit_2d.get_psi_normalized(efit_idx)

                # Use RectBivariateSpline for fast vectorized interpolation
                spline = RectBivariateSpline(
                    self.efit_2d.z_grid, self.efit_2d.r_grid, psi_n_2d
                )
                psi_on_grid = spline(y_grid, x_grid)

                efit_idx_prev = efit_idx

            recon = self.irvb_data.recon[i]

            # Compute Prad for each region
            boundaries = [0] + self.psi_boundaries + [np.inf]

            for r in range(n_regions):
                psi_min = boundaries[r]
                psi_max = boundaries[r + 1]

                mask = (psi_on_grid >= psi_min) & (psi_on_grid < psi_max)
                self.region_prad[r, i] = np.sum(recon * mask * vol_factor)

        print("[IRVB] Regional Prad computation complete")

    def _plot_data(self):
        """Main plot function"""
        t0 = tclock.time()

        if self.irvb_data is None:
            QMessageBox.warning(self.frame, "Warning", "Please load IRVB data first")
            return

        # Parse boundaries
        self._update_status("Parsing psi boundaries...", color='blue')
        self.psi_boundaries = self._parse_psi_boundaries()
        if self.psi_boundaries is None:
            self._update_status("Invalid psi boundaries", color='red')
            return

        # Load EFIT
        self._update_status("Loading EFIT data...", color='blue')
        if not self._load_efit_data():
            self._update_status("Failed to load EFIT", color='red')
            return

        # Slice data to EFIT time range and update frame controls
        self._slice_by_efit_time()
        self._update_frame_controls()

        # Compute regional Prad
        self._update_status("Computing regional Prad...", color='blue')
        self._compute_regional_prad()

        # Setup figure layout
        self._update_status("Plotting...", color='blue')
        self._setup_figure()

        # Initial plot
        self._update_plot(self.current_frame)

        # Enable save button
        self.save_button.setEnabled(True)

        elapsed = tclock.time() - t0
        self._update_status(f"Done ({elapsed:.1f}s)", color='green')
        print(f"[IRVB] Plot completed in {elapsed:.2f}s")
        print("[IRVB] " + "=" * 60 + "\n")

    def _update_frame_controls(self):
        """Update frame controls after data slicing"""
        self.total_frames = len(self.irvb_data.time)
        self.current_frame = 0
        self.frame_slider.setMaximum(self.total_frames - 1)
        self.frame_slider.setValue(0)
        self.frame_total_entry.setReadOnly(False)
        self.frame_total_entry.setText(str(self.total_frames))
        self.frame_total_entry.setReadOnly(True)
        self.frame_entry.setText('1')
        self.time_entry.setText(f'{self.irvb_data.time[0]:.3f}')

    def _setup_figure(self):
        """Setup figure with dynamic number of subplots"""
        self.figure.clear()

        n_regions = len(self.psi_boundaries) + 1

        # GridSpec layout
        gs = GridSpec(n_regions, 2, figure=self.figure,
                     width_ratios=[1.5, 1], wspace=0.10, hspace=0.0,
                     left=0.10, right=0.93, top=0.92, bottom=0.10)

        # Create time trace axes (left column)
        self.ax_traces = []
        for i in range(n_regions):
            if i == 0:
                ax = self.figure.add_subplot(gs[i, 0])
            else:
                ax = self.figure.add_subplot(gs[i, 0], sharex=self.ax_traces[0])
            self.ax_traces.append(ax)

        # Create 2D profile axis (right column)
        self.ax_2d = self.figure.add_subplot(gs[:, 1])
        self.ax_2d.set_label('2D Profile')

        # Set ax references for each time trace (for toolbar Axes access)
        n_regions = len(self.psi_boundaries) + 1
        boundaries = [0] + self.psi_boundaries + [np.inf]
        for i, ax in enumerate(self.ax_traces):
            setattr(self, f'ax_trace_{i}', ax)
            psi_min = boundaries[i]
            psi_max = boundaries[i + 1]
            if psi_max == np.inf:
                ax.set_label(f'P_rad (psi>{psi_min:.2f})')
            else:
                ax.set_label(f'P_rad (psi={psi_min:.2f}-{psi_max:.2f})')

        # Configure toolbar with all time traces (exclude 2D plot)
        if self.toolbar:
            figures_config = []
            for i in range(n_regions):
                psi_min = boundaries[i]
                psi_max = boundaries[i + 1]
                if psi_max == np.inf:
                    label = f'psi>{psi_min:.2f}'
                else:
                    label = f'psi={psi_min:.2f}-{psi_max:.2f}'
                figures_config.append((label, f'ax_trace_{i}'))
            self.toolbar.figures_config = figures_config
            self.toolbar.share_x = True

        # Plot time traces
        self._plot_time_traces()

        # Reset plot references
        self.im_2d = None
        self.contour_psi_bg = None
        self.contour_psi_bounds = []
        self.contour_bdry = None
        self.limiter_line = None
        self.maxis_marker = None
        self.colorbar = None

    def _plot_time_traces(self):
        """Plot regional Prad time traces"""
        time = self.irvb_data.time
        n_regions = len(self.psi_boundaries) + 1
        boundaries = [0] + self.psi_boundaries + [np.inf]

        # Get colors from Plot Options
        self._trace_colors = self._get_trace_colors(n_regions)

        self.time_vlines = []

        for i, ax in enumerate(self.ax_traces):
            ax.clear()

            # Title on first axis with shot number and EFIT tree
            if i == 0:
                efit_display = self.efit_tree_combo.currentText()
                ax.set_title(f'#{self.shot_number} ({efit_display})',
                            fontsize=self.trace_label_fontsize)

            # Region label
            psi_min = boundaries[i]
            psi_max = boundaries[i + 1]

            if psi_max == np.inf:
                label = f'$\\psi_N$ > {psi_min}'
            else:
                label = f'{psi_min} < $\\psi_N$ < {psi_max}'

            color = self._trace_colors[i]

            ax.plot(time, self.region_prad[i], color=color, linewidth=2)

            # Text annotation instead of legend (no line, text only)
            import matplotlib as _mpl
            _text_color = _mpl.rcParams.get('text.color', 'black')
            _bg_color = _mpl.rcParams.get('axes.facecolor', 'white')
            ax.text(0.02, 0.95, label, transform=ax.transAxes,
                   fontsize=self.trace_legend_fontsize, verticalalignment='top',
                   color=_text_color,
                   bbox=dict(boxstyle='round', facecolor=_bg_color, alpha=0.7,
                             edgecolor='gray', linewidth=0.5))

            # ylabel without psi range, fixed position
            ax.set_ylabel(r'$P_{rad}$ [MW]', fontsize=self.trace_label_fontsize)
            ax.yaxis.set_label_coords(-0.08, 0.5)

            # Set ylim with bottom=0
            y_max = np.nanmax(self.region_prad[i]) * 1.1
            if y_max <= 0 or np.isnan(y_max):
                y_max = 1.0
            ax.set_ylim(0, y_max)

            # Only show x-label on bottom plot
            if i == n_regions - 1:
                ax.set_xlabel('Time [s]', fontsize=self.trace_label_fontsize)
            else:
                plt.setp(ax.get_xticklabels(), visible=False)

            ax.set_xlim(0, self.efit_2d.time[-1])
            import matplotlib as mpl
            ax.grid(ls='--', lw=0.3, c=mpl.rcParams.get('grid.color', '#444444'))
            ax.tick_params(labelsize=self.trace_tick_fontsize)

            # Vertical line for current time
            vline = ax.axvline(time[self.current_frame], color='k',
                              linewidth=1.5, linestyle='--')
            self.time_vlines.append(vline)

    def _update_plot(self, frame_idx):
        """Update 2D plot and time markers for given frame"""
        if self.irvb_data is None:
            return

        time = self.irvb_data.time[frame_idx]
        recon = self.irvb_data.recon[frame_idx]
        x_grid = self.irvb_data.x_grid
        y_grid = self.irvb_data.y_grid

        # Get EFIT data at this time
        efit_idx = self.efit_2d.find_time_index(time)
        psi_n = self.efit_2d.get_psi_normalized(efit_idx)
        bdry_r, bdry_z = self.efit_2d.get_boundary(efit_idx)
        maxis_r, maxis_z = self.efit_2d.get_magnetic_axis(efit_idx)

        # Update 2D plot
        if self.im_2d is None:
            self.ax_2d.clear()

            # Prad contour fill
            levels = np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, self.PRAD_LEVELS)
            self.im_2d = self.ax_2d.contourf(x_grid, y_grid, recon,
                                             levels=levels,
                                             vmin=self.PRAD_VMIN,
                                             vmax=self.PRAD_VMAX,
                                             cmap=self.plot2d_colormap,
                                             extend='max')

            # Colorbar
            self.colorbar = self.figure.colorbar(
                self.im_2d, ax=self.ax_2d,
                ticks=np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, 5),
                shrink=0.8
            )
            self.colorbar.set_label(r'$P_{rad}$ [MW/m$^{3}$]',
                                    fontsize=self.plot2d_label_fontsize)

            self.ax_2d.set_xlabel('R [m]', fontsize=self.plot2d_label_fontsize)
            self.ax_2d.set_ylabel('Z [m]', fontsize=self.plot2d_label_fontsize)
            self.ax_2d.set_xlim(x_grid[0], x_grid[-1])
            self.ax_2d.set_ylim(y_grid[0], y_grid[-1])
            self.ax_2d.set_aspect('equal')
            self.ax_2d.tick_params(labelsize=self.plot2d_tick_fontsize)
        else:
            # Update contourf
            for coll in self.im_2d.collections:
                coll.remove()
            levels = np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, self.PRAD_LEVELS)
            self.im_2d = self.ax_2d.contourf(x_grid, y_grid, recon,
                                             levels=levels,
                                             vmin=self.PRAD_VMIN,
                                             vmax=self.PRAD_VMAX,
                                             cmap=self.plot2d_colormap,
                                             extend='max')

        # Remove old contours
        if self.contour_psi_bg is not None:
            for coll in self.contour_psi_bg.collections:
                coll.remove()
        for contour in self.contour_psi_bounds:
            for coll in contour.collections:
                coll.remove()
        self.contour_psi_bounds = []

        # Background psi contours (9 levels)
        self.contour_psi_bg = self.ax_2d.contour(
            self.efit_2d.r_grid, self.efit_2d.z_grid, psi_n,
            levels=self.PSI_BG_LEVELS, colors=self.plot2d_flux_color,
            linewidths=0.5, linestyles='-', alpha=0.5, zorder=2
        )

        # Update plasma boundary (draw before psi bounds so bounds appear on top)
        if self.contour_bdry is not None:
            self.contour_bdry[0].remove()
        self.contour_bdry = self.ax_2d.plot(bdry_r, bdry_z, '-',
                                            color=self.plot2d_lcfs_color,
                                            linewidth=2, zorder=3)

        # Psi boundary contours with matching colors (drawn on top of bdry)
        trace_colors = getattr(self, '_trace_colors', None)
        for idx, psi_level in enumerate(self.psi_boundaries):
            if trace_colors and idx < len(trace_colors):
                color = trace_colors[idx]
            else:
                color = self.REGION_COLORS[idx % len(self.REGION_COLORS)]
            contour = self.ax_2d.contour(
                self.efit_2d.r_grid, self.efit_2d.z_grid, psi_n,
                levels=[psi_level], colors=[color],
                linewidths=2, linestyles='--', zorder=4
            )
            self.contour_psi_bounds.append(contour)

        # Update limiter
        if self.limiter_line is not None:
            self.limiter_line[0].remove()
        if self.efit_2d.limiter_r is not None:
            self.limiter_line = self.ax_2d.plot(
                self.efit_2d.limiter_r, self.efit_2d.limiter_z,
                '-', color=self.plot2d_limiter_color, linewidth=1.5
            )

        # Update magnetic axis marker
        if self.maxis_marker is not None:
            self.maxis_marker[0].remove()
        if maxis_r is not None:
            self.maxis_marker = self.ax_2d.plot(maxis_r, maxis_z, 'x',
                                                 color=self.plot2d_maxis_color,
                                                 markersize=10, markeredgewidth=2)

        # Update title
        self.ax_2d.set_title(f'#{self.shot_number} t = {time:.3f} s',
                            fontsize=self.plot2d_label_fontsize)

        # Update time markers on traces
        for vline in self.time_vlines:
            vline.set_xdata([time, time])

        # Update UI elements
        self.current_frame = frame_idx
        self.frame_entry.setText(str(frame_idx + 1))
        self.time_entry.setText(f'{time:.3f}')

        self.canvas.draw_idle()

    # ===== Frame navigation methods =====

    def _on_slider_change(self, value):
        """Handle slider change"""
        frame_idx = value
        if frame_idx != self.current_frame and self.region_prad is not None:
            self.current_frame = frame_idx
            self._update_plot(frame_idx)

    def _on_frame_entry(self):
        """Handle frame entry"""
        self._goto_frame()

    def _on_mouse_wheel(self, event):
        """Handle mouse wheel for frame navigation"""
        if self.total_frames == 0 or self.region_prad is None:
            return

        current = self.frame_slider.value()

        if event.button == 'up':
            new_frame = max(0, current - 1)
        elif event.button == 'down':
            new_frame = min(self.total_frames - 1, current + 1)
        else:
            return

        if new_frame != current:
            self.current_frame = new_frame
            self.frame_slider.setValue(new_frame)
            self._update_plot(new_frame)

    def _goto_frame(self):
        """Go to specified frame"""
        if self.region_prad is None:
            return

        try:
            frame_num = int(self.frame_entry.text())
            frame_idx = frame_num - 1

            if 0 <= frame_idx < self.total_frames:
                self.current_frame = frame_idx
                self.frame_slider.setValue(frame_idx)
                self._update_plot(frame_idx)
            else:
                QMessageBox.warning(self.frame, "Warning",
                    f"Frame number must be between 1 and {self.total_frames}")
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid frame number")

    def _on_time_entry(self):
        """Handle time entry"""
        self._goto_time()

    def _goto_time(self):
        """Go to frame at specified time"""
        if self.region_prad is None or self.irvb_data is None:
            return

        try:
            time_sec = float(self.time_entry.text())
            # Find closest frame index
            frame_idx = np.argmin(np.abs(self.irvb_data.time - time_sec))

            self.current_frame = frame_idx
            self.frame_slider.setValue(frame_idx)
            self._update_plot(frame_idx)
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid time in seconds")

    def _step_frame(self, delta):
        """Step by delta frames"""
        if self.total_frames == 0 or self.region_prad is None:
            return

        current = self.frame_slider.value()
        new_frame = current + delta
        new_frame = max(0, min(new_frame, self.total_frames - 1))

        self.current_frame = new_frame
        self.frame_slider.setValue(new_frame)
        self._update_plot(new_frame)

    def _goto_first(self):
        """Go to first frame"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        self.current_frame = 0
        self.frame_slider.setValue(0)
        self._update_plot(0)

    def _goto_last(self):
        """Go to last frame"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        last_idx = self.total_frames - 1
        self.current_frame = last_idx
        self.frame_slider.setValue(last_idx)
        self._update_plot(last_idx)

    # ===== Playback methods =====

    def _toggle_play(self):
        """Toggle play/pause"""
        if self.total_frames == 0 or self.region_prad is None:
            return

        if self.is_playing:
            self._pause_play()
        else:
            self._start_play()

    def _start_play(self):
        """Start playback"""
        self.is_playing = True
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPause))
        self._last_frame_time = tclock.time()
        self._fps_history = []

        # Calculate initial delay
        speed_str = self.speed_combo.currentText()
        if speed_str == 'Max':
            target_delay = 1
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 300
            target_delay = max(1, int(base_delay / speed_mult))

        self._play_timer.start(target_delay)

    def _pause_play(self):
        """Pause playback"""
        self.is_playing = False
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPlay))
        self._play_timer.stop()

        self.actual_fps_label.setText('Actual: -- FPS')

    def _stop_play(self):
        """Stop playback and reset to first frame"""
        self._pause_play()
        if self.total_frames > 0 and self.region_prad is not None:
            self.current_frame = 0
            self.frame_slider.setValue(0)
            self._update_plot(0)

    def _play_next_frame(self):
        """Play next frame with speed control"""
        if not self.is_playing:
            return

        current = self.frame_slider.value()
        next_frame = current + 1

        # Handle end of sequence
        if next_frame >= self.total_frames:
            if self.loop_checkbox.isChecked():
                next_frame = 0
            else:
                self._pause_play()
                return

        # Measure actual frame time
        now = tclock.time()
        actual_frame_time = now - self._last_frame_time
        self._last_frame_time = now

        # Calculate and display actual FPS
        if actual_frame_time > 0:
            instant_fps = 1.0 / actual_frame_time
            self._fps_history.append(instant_fps)
            if len(self._fps_history) > 10:
                self._fps_history.pop(0)
            avg_fps = sum(self._fps_history) / len(self._fps_history)
            self.actual_fps_label.setText(f'Actual: {avg_fps:.1f} FPS')

        # Update state and display
        self.current_frame = next_frame
        self.frame_slider.setValue(next_frame)
        self._update_plot(next_frame)

        # Update timer interval for speed changes during playback
        speed_str = self.speed_combo.currentText()
        if speed_str == 'Max':
            target_delay = 1
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 300
            target_delay = max(1, int(base_delay / speed_mult))

        if self._play_timer.interval() != target_delay:
            self._play_timer.setInterval(target_delay)

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.text(),
            "efit_tree": self.efit_tree_combo.currentText(),
            "psi_bounds": self.psi_entry.text(),
            # Plot Options - Time Trace
            "trace_color_mode": self.trace_color_mode,
            "trace_label_fontsize": self.trace_label_fontsize,
            "trace_legend_fontsize": self.trace_legend_fontsize,
            "trace_tick_fontsize": self.trace_tick_fontsize,
            # Plot Options - 2D Plot
            "plot2d_colormap": self.plot2d_colormap,
            "plot2d_lcfs_color": self.plot2d_lcfs_color,
            "plot2d_maxis_color": self.plot2d_maxis_color,
            "plot2d_limiter_color": self.plot2d_limiter_color,
            "plot2d_flux_color": self.plot2d_flux_color,
            "plot2d_label_fontsize": self.plot2d_label_fontsize,
            "plot2d_tick_fontsize": self.plot2d_tick_fontsize,
        }
        set_tab_settings("irvb", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("irvb")

        if settings.get("shot"):
            self.shot_entry.setText(settings["shot"])

        if settings.get("efit_tree"):
            self.efit_tree_combo.setCurrentText(settings["efit_tree"])

        if settings.get("psi_bounds"):
            self.psi_entry.setText(settings["psi_bounds"])

        # Plot Options - Time Trace
        self.trace_color_mode = settings.get("trace_color_mode", "Fixed(tab10)")
        self.trace_label_fontsize = settings.get("trace_label_fontsize", 12)
        self.trace_legend_fontsize = settings.get("trace_legend_fontsize", 8)
        self.trace_tick_fontsize = settings.get("trace_tick_fontsize", 10)

        # Plot Options - 2D Plot
        self.plot2d_colormap = settings.get("plot2d_colormap", "viridis")
        self.plot2d_lcfs_color = settings.get("plot2d_lcfs_color", "black")
        self.plot2d_maxis_color = settings.get("plot2d_maxis_color", "black")
        self.plot2d_limiter_color = settings.get("plot2d_limiter_color", "white")
        self.plot2d_flux_color = settings.get("plot2d_flux_color", "gray")
        self.plot2d_label_fontsize = settings.get("plot2d_label_fontsize", 12)
        self.plot2d_tick_fontsize = settings.get("plot2d_tick_fontsize", 10)
