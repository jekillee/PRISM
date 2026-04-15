"""
TV (Visible Camera) tab for viewing image sequences from ZIP files
With line drawing feature for paper figures and TV1/TV2 compare mode
"""

import zipfile
import io
import re
import os
import time
import threading
import numpy as np
from PIL import Image
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QGroupBox,
    QCheckBox, QSlider, QSplitter, QMessageBox, QFileDialog,
    QApplication, QStyle, QScrollArea,
)
from PySide6.QtCore import Qt, QTimer

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.tv_utils import (
    TV_FPS, TV_OFFSET, get_year_from_shot, get_campaign_from_year,
    get_tv_zip_path, find_available_tvs, frame_to_time, time_to_frame
)
from config.user_settings import get_tab_settings, set_tab_settings


class TVTab:
    """TV image sequence viewer tab with line drawing and compare mode"""

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = QWidget()
        self.toolbar = None

        # Single TV mode data storage
        self.zip_file = None
        self.image_files = []
        self.current_frame = 0
        self.total_frames = 0

        # Compare mode data storage
        self.compare_mode = False
        self.tv1_zip = None
        self.tv2_zip = None
        self.tv1_images = []
        self.tv2_images = []
        self.tv1_cache = {}
        self.tv2_cache = {}

        # Image cache (store nearby frames for smooth playback)
        self.cache = {}
        self.cache_size = 100

        # Playback state
        self.is_playing = False
        self._play_timer = QTimer()
        self._play_timer.timeout.connect(self._play_next_frame)
        self._last_frame_time = 0

        # Prefetch thread control
        self._prefetch_lock = threading.Lock()
        self._prefetch_thread = None
        self._prefetch_stop = False

        # Plot references
        self.im = None
        self.im1 = None  # For compare mode (TV01)
        self.im2 = None  # For compare mode (TV02)
        self.ax = None
        self.ax1 = None  # For compare mode
        self.ax2 = None  # For compare mode

        # Line drawing
        self.line_points = []
        self.drawn_lines = []  # List of (line_obj, ax) tuples for finalized lines
        self.current_line_obj = None  # Line object currently being drawn (not finalized)
        self.preview_line = None
        self.draw_mode = False
        self.click_cid = None
        self.motion_cid = None
        self._draw_background = None  # For blitting optimization
        self.current_draw_ax = None  # Track which axis current line is being drawn on

        # Slider debounce timer
        self._slider_timer = QTimer()
        self._slider_timer.setSingleShot(True)
        self._slider_timer.timeout.connect(self._do_slider_update)

        # Slider value tracking
        self._slider_value = 0

    def create_widgets(self):
        """Create TV tab widgets"""
        main_layout = QHBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)

        # Left: Image display area
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)

        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=False)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("No image loaded")
        self.ax.set_label('TV Image')
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()
        canvas_layout.addWidget(self.canvas)

        # Bind mouse wheel event for frame navigation
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

        self._create_file_controls(control_layout)
        self._create_frame_controls(control_layout)
        self._create_playback_controls(control_layout)
        self._create_draw_line_controls(control_layout)
        control_layout.addStretch()

        scroll_area.setWidget(control_widget)
        scroll_area.viewport().setAutoFillBackground(False)
        control_widget.setAutoFillBackground(False)
        splitter.addWidget(scroll_area)

        # Load saved settings
        self.load_settings()

    def _create_file_controls(self, parent_layout):
        """Create file loading section"""
        group = QGroupBox("1. Load TV Data")
        grid = QGridLayout(group)

        # Row 0: Shot label, entry with up/down, Search button, ... button
        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_layout = QHBoxLayout()
        self.shot_entry = QLineEdit()
        self.shot_entry.setMinimumWidth(80)
        self.shot_entry.returnPressed.connect(self._search_available_tvs)
        shot_layout.addWidget(self.shot_entry, 1)

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

        grid.addLayout(shot_layout, 0, 1)

        fetch_btn = QPushButton('Fetch')
        fetch_btn.setFixedWidth(70)
        fetch_btn.clicked.connect(self._search_available_tvs)
        grid.addWidget(fetch_btn, 0, 2)

        # Row 1: TV dropdown and Load button
        grid.addWidget(QLabel('TV'), 1, 0)

        self.tv_dropdown = QComboBox()
        self.tv_dropdown.addItem('-- Select --')
        grid.addWidget(self.tv_dropdown, 1, 1)

        load_btn = QPushButton('Load')
        load_btn.setFixedWidth(70)
        load_btn.clicked.connect(self._load_selected_tv)
        grid.addWidget(load_btn, 1, 2)

        # Open TV File button
        browse_btn = QPushButton("Open TV File...")
        browse_btn.clicked.connect(self._load_zip_file)
        grid.addWidget(browse_btn, 2, 0, 1, 3)

        # File label
        self.file_label = QLabel("No file loaded")
        self.file_label.setWordWrap(True)
        grid.addWidget(self.file_label, 3, 0, 1, 3)

        # Loading status label
        self.status_label = QLabel("")
        self.status_label.setStyleSheet("color: blue; font-weight: bold; font-size: 9pt;")
        grid.addWidget(self.status_label, 4, 0, 1, 3)

        parent_layout.addWidget(group)

    def _search_available_tvs(self):
        """Search for available TVs and update dropdown"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        available_tvs = find_available_tvs(shot_number)

        if not available_tvs:
            QMessageBox.information(self.frame, "Not Found",
                f"No TV data found for shot #{shot_number}")
            self.tv_dropdown.clear()
            self.tv_dropdown.addItem('-- Select --')
            self.file_label.setText("No TV files found")
            return

        # Build dropdown values based on availability
        dropdown_values = []
        has_both = 'TV01' in available_tvs and 'TV02' in available_tvs

        if has_both:
            dropdown_values.append('TV01 + TV02')
        if 'TV01' in available_tvs:
            dropdown_values.append('TV01')
        if 'TV02' in available_tvs:
            dropdown_values.append('TV02')

        self.tv_dropdown.clear()
        self.tv_dropdown.addItems(dropdown_values)

        # Set default: TV01 + TV02 if both available, otherwise first available
        if has_both:
            self.tv_dropdown.setCurrentText('TV01 + TV02')
        else:
            self.tv_dropdown.setCurrentIndex(0)

        self.file_label.setText(f"Found: {', '.join(available_tvs)} for #{shot_number}")
        print(f"[TV] Found {available_tvs} for shot #{shot_number}")

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
        except ValueError:
            pass

    def _create_frame_controls(self, parent_layout):
        """Create frame navigation controls"""
        group = QGroupBox("2. Frame Control")
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

        self.frame_entry = QLineEdit('0')
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
        self.time_entry.setFixedWidth(60)
        self.time_entry.returnPressed.connect(self._goto_time)
        time_input_layout.addWidget(self.time_entry)

        go_time_btn = QPushButton()
        go_time_btn.setIcon(self.frame.style().standardIcon(QStyle.SP_DialogOkButton))
        go_time_btn.setFixedSize(24, 24)
        go_time_btn.setToolTip("Go to time")
        go_time_btn.clicked.connect(self._goto_time)
        time_input_layout.addStretch()

        layout.addLayout(time_input_layout)

        # Current filename display
        self.filename_label = QLabel("")
        self.filename_label.setWordWrap(True)
        layout.addWidget(self.filename_label)

        # Mouse wheel hint
        hint_label = QLabel("(Mouse wheel: navigate frames)")
        hint_label.setStyleSheet("color: gray; font-size: 8pt;")
        hint_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(hint_label)

        parent_layout.addWidget(group)

    def _create_playback_controls(self, parent_layout):
        """Create playback control section"""
        group = QGroupBox("3. Playback Control")
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

    def _create_draw_line_controls(self, parent_layout):
        """Create line drawing controls"""
        group = QGroupBox("4. Draw Line")
        layout = QVBoxLayout(group)

        # Draw mode button and Clear button
        btn_layout = QHBoxLayout()

        self.draw_btn = QPushButton('Draw Mode: OFF')
        self.draw_btn.clicked.connect(self._toggle_draw_mode)
        btn_layout.addWidget(self.draw_btn)

        clear_btn = QPushButton('Clear')
        clear_btn.clicked.connect(self._clear_line)
        btn_layout.addWidget(clear_btn)

        layout.addLayout(btn_layout)

        # Show line checkbox and smooth curve checkbox
        check_layout = QHBoxLayout()

        self.show_line_checkbox = QCheckBox('Show Line')
        self.show_line_checkbox.setChecked(True)
        self.show_line_checkbox.stateChanged.connect(lambda: self._update_line_display())
        check_layout.addWidget(self.show_line_checkbox)

        self.smooth_checkbox = QCheckBox('Smooth')
        self.smooth_checkbox.setChecked(True)
        self.smooth_checkbox.stateChanged.connect(lambda: self._update_line_display())
        check_layout.addWidget(self.smooth_checkbox)

        check_layout.addStretch()

        # Points count label
        self.points_label = QLabel('Points: 0')
        check_layout.addWidget(self.points_label)

        layout.addLayout(check_layout)

        # Line style options (2 rows for compact layout)
        style_row1 = QHBoxLayout()
        style_row1.addWidget(QLabel('Color:'))
        self.line_color_combo = QComboBox()
        self.line_color_combo.addItems(['white', 'black', 'red', 'blue', 'yellow', 'green'])
        self.line_color_combo.setCurrentText('white')
        self.line_color_combo.currentIndexChanged.connect(lambda: self._update_line_display())
        style_row1.addWidget(self.line_color_combo, 1)

        style_row1.addWidget(QLabel('Width:'))
        self.line_width_combo = QComboBox()
        self.line_width_combo.addItems(['1', '2', '3', '4', '5'])
        self.line_width_combo.setCurrentText('2')
        self.line_width_combo.setFixedWidth(50)
        self.line_width_combo.currentIndexChanged.connect(lambda: self._update_line_display())
        style_row1.addWidget(self.line_width_combo)

        layout.addLayout(style_row1)

        style_row2 = QHBoxLayout()
        style_row2.addWidget(QLabel('Style:'))
        self.line_style_combo = QComboBox()
        self.line_style_combo.addItems(['dashed', 'solid', 'dotted'])
        self.line_style_combo.setCurrentText('dashed')
        self.line_style_combo.currentIndexChanged.connect(lambda: self._update_line_display())
        style_row2.addWidget(self.line_style_combo, 1)
        style_row2.addStretch()

        layout.addLayout(style_row2)

        # Hint label
        hint_label = QLabel("(Left-click: add point, Right-click: finish)")
        hint_label.setStyleSheet("color: gray; font-size: 8pt;")
        hint_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(hint_label)

        parent_layout.addWidget(group)

    def _toggle_draw_mode(self):
        """Toggle line drawing mode"""
        self.draw_mode = not self.draw_mode

        if self.draw_mode:
            self.draw_btn.setText('Draw Mode: ON')
            self.draw_btn.setStyleSheet('background-color: #90EE90; color: black;')
            self.click_cid = self.canvas.mpl_connect('button_press_event',
                                                      self._on_line_click)
            self.motion_cid = self.canvas.mpl_connect('motion_notify_event',
                                                       self._on_mouse_motion)
            self.canvas.setCursor(Qt.CrossCursor)
            # Capture background for blitting optimization
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)
        else:
            self.draw_btn.setText('Draw Mode: OFF')
            self.draw_btn.setStyleSheet('')
            if self.click_cid:
                self.canvas.mpl_disconnect(self.click_cid)
                self.click_cid = None
            if self.motion_cid:
                self.canvas.mpl_disconnect(self.motion_cid)
                self.motion_cid = None
            self.canvas.setCursor(Qt.ArrowCursor)
            # Clear blitting background
            self._draw_background = None
            # Remove preview line
            if self.preview_line:
                self.preview_line.remove()
                self.preview_line = None
                self.canvas.draw_idle()

    def _on_line_click(self, event):
        """Handle mouse click for line drawing"""
        if event.inaxes is None:
            return

        # Determine target axis
        if self.compare_mode:
            if event.inaxes == self.ax1:
                target_ax = self.ax1
            elif event.inaxes == self.ax2:
                target_ax = self.ax2
            else:
                return
        else:
            if event.inaxes != self.ax:
                return
            target_ax = self.ax

        if event.button == 1:  # Left click - add point
            # If switching to a different axis, finalize current line first
            if self.current_draw_ax is not None and self.current_draw_ax != target_ax:
                self._finalize_current_line()

            self.current_draw_ax = target_ax
            self.line_points.append((event.xdata, event.ydata))
            self.points_label.setText(f'Points: {len(self.line_points)}')
            self._update_line_display()
        elif event.button == 3:  # Right click - finish current line
            self._finalize_current_line()

    def _on_mouse_motion(self, event):
        """Handle mouse motion for preview line"""
        if not self.draw_mode:
            return

        # Must have started drawing (current_draw_ax set)
        if self.current_draw_ax is None or len(self.line_points) == 0:
            return

        # Only show preview on the axis where we're drawing
        if event.inaxes != self.current_draw_ax:
            return

        # Update preview line from last point to cursor
        last_point = self.line_points[-1]

        if self.preview_line:
            self.preview_line.set_data([last_point[0], event.xdata],
                                        [last_point[1], event.ydata])
        else:
            self.preview_line, = self.current_draw_ax.plot([last_point[0], event.xdata],
                                                            [last_point[1], event.ydata],
                                                            'r--', linewidth=1, alpha=0.5)

        # Use blitting for smooth animation without flickering
        if self._draw_background is not None:
            self.canvas.restore_region(self._draw_background)
            self.current_draw_ax.draw_artist(self.preview_line)
            self.canvas.blit(self.figure.bbox)
        else:
            self.canvas.draw_idle()

    def _finalize_current_line(self):
        """Finalize the current line and prepare for a new one"""
        # Move current line to finalized lines list
        if self.current_line_obj is not None and self.current_draw_ax is not None:
            self.drawn_lines.append((self.current_line_obj, self.current_draw_ax))
            self.current_line_obj = None
        # Clear points for new line
        self.line_points = []
        self.current_draw_ax = None
        self.points_label.setText('Points: 0')
        # Remove preview line
        if self.preview_line:
            self.preview_line.remove()
            self.preview_line = None
        # Update background for blitting
        self.canvas.draw()
        self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

    def _clear_line(self):
        """Clear all drawn lines"""
        self.line_points = []
        self.current_draw_ax = None
        self.points_label.setText('Points: 0')

        # Remove current line being drawn
        if self.current_line_obj is not None:
            try:
                self.current_line_obj.remove()
            except ValueError:
                pass
            self.current_line_obj = None

        # Remove all finalized lines
        for line_obj, _ in self.drawn_lines:
            try:
                line_obj.remove()
            except ValueError:
                pass  # Already removed
        self.drawn_lines = []

        if self.preview_line:
            self.preview_line.remove()
            self.preview_line = None

        self.canvas.draw_idle()

        # Update background if in draw mode
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

    def _update_line_display(self):
        """Update line display based on current settings"""
        if not self.show_line_checkbox.isChecked() or len(self.line_points) < 2:
            self.canvas.draw_idle()
            return

        if self.current_draw_ax is None:
            return

        points = np.array(self.line_points)
        x, y = points[:, 0], points[:, 1]

        # Apply smoothing if enabled
        if self.smooth_checkbox.isChecked() and len(points) >= 4:
            try:
                from scipy.interpolate import splprep, splev
                tck, u = splprep([x, y], s=0)
                u_new = np.linspace(0, 1, 100)
                x, y = splev(u_new, tck)
            except Exception:
                pass  # Fall back to non-smooth

        color = self.line_color_combo.currentText()
        width = int(self.line_width_combo.currentText())
        style_map = {'dashed': '--', 'solid': '-', 'dotted': ':'}
        linestyle = style_map.get(self.line_style_combo.currentText(), '--')

        # Remove previous version of current line if exists
        if self.current_line_obj is not None:
            try:
                self.current_line_obj.remove()
            except ValueError:
                pass
            self.current_line_obj = None

        # Draw new line
        self.current_line_obj, = self.current_draw_ax.plot(x, y, color=color, linewidth=width, linestyle=linestyle)
        self.canvas.draw_idle()

        # Update background for blitting (include newly drawn line)
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

    # =========================================================================
    # Time/Frame conversion utilities
    # =========================================================================

    def _find_matching_frame(self, master_frame, master_total, slave_total):
        """Find slave frame that matches master frame's time. Returns None if out of range."""
        time_sec = frame_to_time(master_frame)

        # Calculate the time range of slave TV
        slave_min_time = frame_to_time(0)
        slave_max_time = frame_to_time(slave_total - 1)

        # Check if master time is within slave's time range
        if time_sec < slave_min_time or time_sec > slave_max_time:
            return None

        return time_to_frame(time_sec, slave_total)

    def _load_selected_tv(self):
        """Load selected TV ZIP file(s)"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        tv_selection = self.tv_dropdown.currentText()

        if not tv_selection or tv_selection == '-- Select --':
            QMessageBox.critical(self.frame, "Error", "Please search and select a TV first")
            return

        # Handle TV01 + TV02 compare mode
        if tv_selection == 'TV01 + TV02':
            tv1_path = get_tv_zip_path(shot_number, 'TV01')
            tv2_path = get_tv_zip_path(shot_number, 'TV02')

            if tv1_path is None or tv2_path is None:
                QMessageBox.critical(self.frame, "Error", "Failed to get TV file paths")
                return

            self._load_compare_mode(tv1_path, tv2_path)
        else:
            # Single TV mode
            zip_path = get_tv_zip_path(shot_number, tv_selection)

            if zip_path is None:
                QMessageBox.critical(self.frame, "Error", "Campaign folder not found")
                return

            self.compare_mode = False
            self._setup_single_mode()
            self._load_zip_from_path(zip_path)

    def _setup_single_mode(self):
        """Setup figure for single TV display"""
        self._cleanup_compare_mode()

        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("No image loaded")
        self.ax.set_label('TV Image')

        self.im = None
        self.ax1 = None
        self.ax2 = None
        self.im1 = None
        self.im2 = None

        self.canvas.draw()

    def _setup_compare_mode(self):
        """Setup figure for TV1/TV2 comparison (side by side)"""
        self.figure.clear()

        # Left: TV02, Right: TV01
        self.ax2 = self.figure.add_subplot(121)
        self.ax1 = self.figure.add_subplot(122)

        self.ax1.set_xticks([])
        self.ax1.set_yticks([])
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])

        self.ax1.set_title("TV01")
        self.ax1.set_label('TV01')
        self.ax2.set_title("TV02")
        self.ax2.set_label('TV02')

        self.im = None
        self.ax = None
        self.im1 = None
        self.im2 = None

        self.canvas.draw()

    def _cleanup_compare_mode(self):
        """Cleanup compare mode resources"""
        if self.tv1_zip:
            try:
                self.tv1_zip.close()
            except Exception:
                pass
            self.tv1_zip = None

        if self.tv2_zip:
            try:
                self.tv2_zip.close()
            except Exception:
                pass
            self.tv2_zip = None

        self.tv1_images = []
        self.tv2_images = []
        self.tv1_cache.clear()
        self.tv2_cache.clear()

    def _load_compare_mode(self, tv1_path, tv2_path):
        """Load both TV01 and TV02 for comparison"""
        self._stop_play()
        self._stop_prefetch()
        self._cleanup_compare_mode()

        # Close single mode zip if open
        if self.zip_file:
            try:
                self.zip_file.close()
            except Exception:
                pass
            self.zip_file = None
        self.cache.clear()

        self._set_status("Loading TV01...", 'blue')
        print(f"[TV] Loading TV01 from {tv1_path}")

        try:
            self.tv1_zip = zipfile.ZipFile(tv1_path, 'r')
            self.tv1_images = self._get_sorted_images(self.tv1_zip)
            print(f"[TV] TV01 has {len(self.tv1_images)} frames")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load TV01:\n{str(e)}")
            self._set_status("Failed", 'red')
            return

        self._set_status("Loading TV02...", 'blue')
        print(f"[TV] Loading TV02 from {tv2_path}")

        try:
            self.tv2_zip = zipfile.ZipFile(tv2_path, 'r')
            self.tv2_images = self._get_sorted_images(self.tv2_zip)
            print(f"[TV] TV02 has {len(self.tv2_images)} frames")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load TV02:\n{str(e)}")
            self._set_status("Failed", 'red')
            self._cleanup_compare_mode()
            return

        if not self.tv1_images or not self.tv2_images:
            QMessageBox.critical(self.frame, "Error", "No images found in one or both ZIP files")
            self._cleanup_compare_mode()
            return

        # Enable compare mode
        self.compare_mode = True
        self._setup_compare_mode()

        # Use TV01 as master for frame count
        self.total_frames = len(self.tv1_images)
        self.image_files = self.tv1_images  # For compatibility

        # Update UI
        self.file_label.setText(f"TV01: {len(self.tv1_images)} frames, TV02: {len(self.tv2_images)} frames")
        self.frame_slider.setMaximum(self.total_frames - 1)
        self.frame_total_entry.setReadOnly(False)
        self.frame_total_entry.setText(str(self.total_frames))
        self.frame_total_entry.setReadOnly(True)

        # Display first frame
        self.current_frame = 0
        self.frame_slider.setValue(0)
        self._display_compare_frame(0)

        self._start_prefetch(0)
        self._set_status("Ready", 'green')
        print(f"[TV] Compare mode ready")

    def _get_sorted_images(self, zip_file):
        """Get sorted list of image files from ZIP"""
        all_files = zip_file.namelist()

        image_extensions = ('.png', '.jpg', '.jpeg', '.tif', '.tiff', '.bmp')
        image_files = [f for f in all_files
                      if f.lower().endswith(image_extensions) and not f.startswith('__MACOSX')]

        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower()
                   for text in re.split('([0-9]+)', s)]

        image_files.sort(key=natural_sort_key)
        return image_files

    def _load_zip_file(self):
        """Open file dialog and load selected ZIP file"""
        initial_dir = '/Diag_TV' if os.path.exists('/Diag_TV') else ''
        file_path, _ = QFileDialog.getOpenFileName(
            self.frame,
            "Select TV ZIP file",
            initial_dir,
            "ZIP files (*.zip);;All files (*.*)"
        )

        if file_path:
            self.compare_mode = False
            self._setup_single_mode()
            self._load_zip_from_path(file_path)

    def _set_status(self, text, color='blue'):
        """Update status label with color"""
        self.status_label.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 9pt;")
        self.status_label.setText(text)
        QApplication.processEvents()

    def _load_zip_from_path(self, file_path):
        """Load images from ZIP file (single TV mode)"""
        self._stop_play()
        self._stop_prefetch()
        self._cleanup_compare_mode()

        if self.zip_file:
            try:
                self.zip_file.close()
            except Exception:
                pass

        self.cache.clear()
        self.im = None

        self._set_status("Opening ZIP file...", 'blue')
        print(f"[TV] Opening ZIP file: {file_path}")

        try:
            self.zip_file = zipfile.ZipFile(file_path, 'r')
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to open ZIP file:\n{str(e)}")
            self._set_status("Failed to open ZIP", 'red')
            print(f"[TV] Failed to open ZIP: {str(e)}")
            return

        self._set_status("Reading file list...", 'blue')
        self.image_files = self._get_sorted_images(self.zip_file)

        if not self.image_files:
            QMessageBox.critical(self.frame, "Error", "No image files found in ZIP")
            self._set_status("No images found", 'red')
            return

        self.total_frames = len(self.image_files)
        print(f"[TV] Found {self.total_frames} image files")

        self._set_status("Updating UI...", 'blue')
        self.file_label.setText(file_path.split('/')[-1])
        self.frame_slider.setMaximum(self.total_frames - 1)
        self.frame_total_entry.setReadOnly(False)
        self.frame_total_entry.setText(str(self.total_frames))
        self.frame_total_entry.setReadOnly(True)

        self._set_status("Loading first frame...", 'blue')
        self.current_frame = 0
        self.frame_slider.setValue(0)

        if not self._display_frame(0):
            self._set_status("Failed to load first frame", 'red')
            QMessageBox.critical(self.frame, "Error", "Failed to load first frame from ZIP")
            return

        self._start_prefetch(0)
        self._set_status("Ready", 'green')
        print(f"[TV] Successfully loaded {self.total_frames} images")

    # =========================================================================
    # Image loading and display
    # =========================================================================

    def _get_image(self, frame_idx):
        """Get image from cache or load from ZIP (single mode)"""
        if frame_idx < 0 or frame_idx >= self.total_frames:
            return None

        if self.zip_file is None:
            return None

        with self._prefetch_lock:
            if frame_idx in self.cache:
                return self.cache[frame_idx]

        try:
            filename = self.image_files[frame_idx]
            with self.zip_file.open(filename) as f:
                img_data = f.read()
                img = Image.open(io.BytesIO(img_data))
                img_array = np.array(img)

            with self._prefetch_lock:
                self.cache[frame_idx] = img_array
                self._cleanup_cache(frame_idx)

            return img_array

        except Exception as e:
            print(f"[TV] Error loading frame {frame_idx}: {str(e)}")
            return None

    def _get_image_from_tv(self, tv_num, frame_idx):
        """Get image from specific TV cache or load from ZIP"""
        if tv_num == 1:
            zip_file = self.tv1_zip
            images = self.tv1_images
            cache = self.tv1_cache
        else:
            zip_file = self.tv2_zip
            images = self.tv2_images
            cache = self.tv2_cache

        if frame_idx < 0 or frame_idx >= len(images):
            return None

        if zip_file is None:
            return None

        with self._prefetch_lock:
            if frame_idx in cache:
                return cache[frame_idx]

        try:
            filename = images[frame_idx]
            with zip_file.open(filename) as f:
                img_data = f.read()
                img = Image.open(io.BytesIO(img_data))
                img_array = np.array(img)

            with self._prefetch_lock:
                cache[frame_idx] = img_array
                # Limit cache size
                if len(cache) > self.cache_size // 2:
                    keys = list(cache.keys())
                    for k in keys[:len(cache) - self.cache_size // 2]:
                        if k != frame_idx:
                            del cache[k]

            return img_array

        except Exception as e:
            print(f"[TV] Error loading TV{tv_num} frame {frame_idx}: {str(e)}")
            return None

    def _cleanup_cache(self, center_frame):
        """Remove frames far from current position"""
        if len(self.cache) > self.cache_size:
            keys_to_remove = []
            for key in self.cache.keys():
                if abs(key - center_frame) > self.cache_size // 2:
                    keys_to_remove.append(key)
            for key in keys_to_remove[:len(self.cache) - self.cache_size]:
                del self.cache[key]

    def _start_prefetch(self, center_frame, direction=1):
        """Start background thread to prefetch frames"""
        self._stop_prefetch()

        self._prefetch_stop = False
        self._prefetch_thread = threading.Thread(
            target=self._prefetch_worker,
            args=(center_frame, direction),
            daemon=True
        )
        self._prefetch_thread.start()

    def _stop_prefetch(self):
        """Stop prefetch thread"""
        self._prefetch_stop = True
        if self._prefetch_thread and self._prefetch_thread.is_alive():
            self._prefetch_thread.join(timeout=0.1)

    def _prefetch_worker(self, center_frame, direction):
        """Background worker to prefetch frames"""
        prefetch_count = min(30, self.total_frames)

        for i in range(1, prefetch_count + 1):
            if self._prefetch_stop:
                break

            frame_idx = center_frame + (i * direction)
            if 0 <= frame_idx < self.total_frames:
                if self.compare_mode:
                    # Prefetch both TVs
                    self._get_image_from_tv(1, frame_idx)
                    tv2_frame = self._find_matching_frame(
                        frame_idx, len(self.tv1_images), len(self.tv2_images))
                    self._get_image_from_tv(2, tv2_frame)
                else:
                    with self._prefetch_lock:
                        if frame_idx not in self.cache:
                            try:
                                filename = self.image_files[frame_idx]
                                with self.zip_file.open(filename) as f:
                                    img_data = f.read()
                                    img = Image.open(io.BytesIO(img_data))
                                    img_array = np.array(img)
                                self.cache[frame_idx] = img_array
                            except Exception:
                                pass

    def _display_frame(self, frame_idx, update_ui=True):
        """Display specified frame (single TV mode)"""
        if self.compare_mode:
            return self._display_compare_frame(frame_idx, update_ui)

        img_array = self._get_image(frame_idx)

        if img_array is None:
            return False

        filename = self.image_files[frame_idx]

        if self.im is None:
            self.ax.clear()
            self.im = self.ax.imshow(img_array, cmap='gray' if len(img_array.shape) == 2 else None)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            # Clear lines associated with this axis
            self.drawn_lines = [(l, ax) for l, ax in self.drawn_lines if ax != self.ax]
            if self.current_draw_ax == self.ax:
                self.current_line_obj = None
        else:
            self.im.set_data(img_array)

        time_sec = frame_to_time(frame_idx)
        self.ax.set_title(f"Frame {frame_idx + 1}/{self.total_frames} (t = {time_sec:.3f} s)")

        self.canvas.draw_idle()
        self.canvas.flush_events()

        # Update blitting background if draw mode is active
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

        self.current_frame = frame_idx

        if update_ui:
            self.frame_entry.setText(str(frame_idx + 1))
            self.time_entry.setText(f'{time_sec:.3f}')
            self.filename_label.setText(filename)

        return True

    def _display_compare_frame(self, frame_idx, update_ui=True):
        """Display frames in compare mode (TV01 and TV02 side by side)"""
        # TV01 frame (master)
        tv1_img = self._get_image_from_tv(1, frame_idx)
        tv1_time = frame_to_time(frame_idx)

        # Find matching TV02 frame by time (returns None if out of range)
        tv2_frame = self._find_matching_frame(
            frame_idx, len(self.tv1_images), len(self.tv2_images))

        # Get TV02 image only if frame is valid
        tv2_img = None
        tv2_time = tv1_time  # Use same time for display
        if tv2_frame is not None:
            tv2_img = self._get_image_from_tv(2, tv2_frame)
            tv2_time = frame_to_time(tv2_frame)

        # Check if at least one has data
        if tv1_img is None and tv2_img is None:
            return False

        # Display TV01 (right)
        if tv1_img is not None:
            if self.im1 is None:
                self.ax1.clear()
                self.im1 = self.ax1.imshow(tv1_img, cmap='gray' if len(tv1_img.shape) == 2 else None)
                self.ax1.set_xticks([])
                self.ax1.set_yticks([])
            else:
                self.im1.set_data(tv1_img)
            self.ax1.set_title(f"TV01 Frame {frame_idx + 1}/{len(self.tv1_images)} (t = {tv1_time:.3f} s)")
        else:
            # No data for TV01 at this time
            self.ax1.clear()
            self.ax1.set_facecolor('black')
            self.ax1.text(0.5, 0.5, 'No Data', ha='center', va='center',
                         fontsize=16, color='white', transform=self.ax1.transAxes)
            self.ax1.set_xticks([])
            self.ax1.set_yticks([])
            self.ax1.set_title(f"TV01 (t = {tv1_time:.3f} s) - No Data")
            self.im1 = None

        # Display TV02 (left)
        if tv2_img is not None and tv2_frame is not None:
            if self.im2 is None:
                self.ax2.clear()
                self.im2 = self.ax2.imshow(tv2_img, cmap='gray' if len(tv2_img.shape) == 2 else None)
                self.ax2.set_xticks([])
                self.ax2.set_yticks([])
            else:
                self.im2.set_data(tv2_img)
            self.ax2.set_title(f"TV02 Frame {tv2_frame + 1}/{len(self.tv2_images)} (t = {tv2_time:.3f} s)")
        else:
            # No data for TV02 at this time
            self.ax2.clear()
            self.ax2.set_facecolor('black')
            self.ax2.text(0.5, 0.5, 'No Data', ha='center', va='center',
                         fontsize=16, color='white', transform=self.ax2.transAxes)
            self.ax2.set_xticks([])
            self.ax2.set_yticks([])
            self.ax2.set_title(f"TV02 (t = {tv1_time:.3f} s) - No Data")
            self.im2 = None

        self.canvas.draw_idle()
        self.canvas.flush_events()

        # Update blitting background if draw mode is active
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

        self.current_frame = frame_idx

        # Always update filename display
        tv1_name = self.tv1_images[frame_idx] if frame_idx < len(self.tv1_images) else "N/A"
        tv2_name = self.tv2_images[tv2_frame] if tv2_frame is not None and tv2_frame < len(self.tv2_images) else "N/A"
        self.filename_label.setText(f"TV01: {tv1_name}\nTV02: {tv2_name}")

        if update_ui:
            self.frame_entry.setText(str(frame_idx + 1))
            self.time_entry.setText(f'{tv1_time:.3f}')

        return True

    # =========================================================================
    # Frame navigation
    # =========================================================================

    def _on_slider_change(self, value):
        """Handle slider value change with debounce"""
        self._slider_value = value
        self._slider_timer.start(10)

    def _do_slider_update(self):
        """Actually perform the slider update"""
        frame_idx = self._slider_value
        self._display_frame(frame_idx, update_ui=True)

    def _on_frame_entry(self):
        """Handle frame entry input"""
        self._goto_frame()

    def _on_mouse_wheel(self, event):
        """Handle mouse wheel event for frame navigation"""
        if self.total_frames == 0:
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
            self._display_frame(new_frame)

    def _goto_frame(self):
        """Go to specified frame number"""
        try:
            frame_num = int(self.frame_entry.text())
            frame_idx = frame_num - 1

            if 0 <= frame_idx < self.total_frames:
                self.current_frame = frame_idx
                self.frame_slider.setValue(frame_idx)
                self._display_frame(frame_idx)
                self._start_prefetch(frame_idx)
            else:
                QMessageBox.warning(self.frame, "Warning",
                    f"Frame number must be between 1 and {self.total_frames}")
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid frame number")

    def _on_time_entry(self):
        """Handle time entry input"""
        self._goto_time()

    def _goto_time(self):
        """Go to frame at specified time"""
        if self.total_frames == 0:
            return

        try:
            time_sec = float(self.time_entry.text())
            frame_idx = time_to_frame(time_sec, self.total_frames)

            self.current_frame = frame_idx
            self.frame_slider.setValue(frame_idx)
            self._display_frame(frame_idx)
            self._start_prefetch(frame_idx)
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid time in seconds")

    def _step_frame(self, delta):
        """Step forward/backward by delta frames"""
        if self.total_frames == 0:
            return
        current = self.frame_slider.value()
        new_frame = current + delta
        new_frame = max(0, min(new_frame, self.total_frames - 1))

        self.current_frame = new_frame
        self.frame_slider.setValue(new_frame)
        self._display_frame(new_frame)

    def _goto_first(self):
        """Go to first frame"""
        if self.total_frames == 0:
            return
        self.current_frame = 0
        self.frame_slider.setValue(0)
        self._display_frame(0)
        self._start_prefetch(0)

    def _goto_last(self):
        """Go to last frame"""
        if self.total_frames == 0:
            return
        last_idx = self.total_frames - 1
        self.current_frame = last_idx
        self.frame_slider.setValue(last_idx)
        self._display_frame(last_idx)
        self._start_prefetch(last_idx, direction=-1)

    # =========================================================================
    # Playback controls
    # =========================================================================

    def _toggle_play(self):
        """Toggle play/pause"""
        if self.total_frames == 0:
            return

        if self.is_playing:
            self._pause_play()
        else:
            self._start_play()

    def _start_play(self):
        """Start playback"""
        self.is_playing = True
        self.play_btn.setIcon(get_icon(QStyle.SP_MediaPause))
        self._last_frame_time = time.time()
        self._fps_history = []

        self._start_prefetch(self.current_frame, direction=1)

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

        if self.total_frames > 0:
            self.frame_entry.setText(str(self.current_frame + 1))
            if self.current_frame < len(self.image_files):
                self.filename_label.setText(self.image_files[self.current_frame])

    def _stop_play(self):
        """Stop playback and reset to first frame"""
        self._pause_play()
        if self.total_frames > 0:
            self.current_frame = 0
            self.frame_slider.setValue(0)
            self._display_frame(0)

    def _play_next_frame(self):
        """Play next frame with speed multiplier"""
        if not self.is_playing:
            return

        current = self.frame_slider.value()
        next_frame = current + 1

        if next_frame >= self.total_frames:
            if self.loop_checkbox.isChecked():
                next_frame = 0
            else:
                self._pause_play()
                return

        now = time.time()
        actual_frame_time = now - self._last_frame_time
        self._last_frame_time = now

        if actual_frame_time > 0:
            instant_fps = 1.0 / actual_frame_time
            if not hasattr(self, '_fps_history'):
                self._fps_history = []
            self._fps_history.append(instant_fps)
            if len(self._fps_history) > 10:
                self._fps_history.pop(0)
            avg_fps = sum(self._fps_history) / len(self._fps_history)
            self.actual_fps_label.setText(f'Actual: {avg_fps:.1f} FPS')

        self.current_frame = next_frame
        self.frame_slider.setValue(next_frame)
        self._display_frame(next_frame, update_ui=False)

        # Update timer interval for speed changes during playback
        speed_str = self.speed_combo.currentText()
        if speed_str == 'Max':
            target_delay = 1
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 300
            target_delay = max(1, int(base_delay / speed_mult))

        # Update the timer interval if changed
        if self._play_timer.interval() != target_delay:
            self._play_timer.setInterval(target_delay)

    # =========================================================================
    # Cleanup and settings
    # =========================================================================

    def cleanup(self):
        """Cleanup resources when tab is closed"""
        self._stop_play()
        self._stop_prefetch()

        if self.zip_file:
            try:
                self.zip_file.close()
            except Exception:
                pass
            self.zip_file = None

        self._cleanup_compare_mode()
        self.cache.clear()

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.text()
        }
        set_tab_settings("tv", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("tv")

        if settings.get("shot"):
            self.shot_entry.setText(settings["shot"])
