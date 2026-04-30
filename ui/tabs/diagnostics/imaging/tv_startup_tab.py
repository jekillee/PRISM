"""
TV Startup Comparison tab for comparing startup sequences across multiple shots
Supports both TV01 and TV02 camera selection
Based on standalone TV01 startup montage viewer
"""

import zipfile
import os
import io
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QLabel, QLineEdit, QPushButton, QComboBox, QGroupBox,
    QRadioButton, QButtonGroup, QListWidget, QSplitter,
    QMessageBox, QApplication, QStyle, QScrollArea,
)
from PySide6.QtCore import Qt

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from data_loaders.tv_loader import (
    TV_FPS, TV_OFFSET, get_tv_startup_zip_path, find_available_tvs
)
from ui.tabs.diagnostics.imaging.tv_utils import frame_to_time
from config.user_settings import get_tab_settings, set_tab_settings


class TVStartupTab:
    """TV startup comparison viewer tab - compare startup sequences across shots"""

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90

    # Fixed image dimensions
    FRAME_HEIGHT = 320
    FRAME_WIDTH = 135
    LABEL_WIDTH = 110

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = QWidget()
        self.toolbar = None

        # Shot list storage: list of (shot_number, tv_name) tuples
        self.shot_list = []
        self.custom_zip_paths = {}  # (shot, tv_name) -> custom file path
        self._last_open_dir = None
        self.last_canvas = None  # Combined stacked canvas
        self.available_tvs = []  # Available TVs for current shot

        # Plot references
        self.ax = None
        self.im = None

        # Track last fetched shot for validation
        self.last_fetched_shot = None

    def create_widgets(self):
        """Create TV Startup tab widgets"""
        main_layout = QHBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)

        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)

        # Left: Image display area
        self._create_display(splitter)

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

        self._create_shot_controls(control_layout)
        self._create_display_controls(control_layout)
        self._create_plot_controls(control_layout)
        control_layout.addStretch()

        scroll_area.setWidget(control_widget)
        scroll_area.viewport().setAutoFillBackground(False)
        control_widget.setAutoFillBackground(False)
        splitter.addWidget(scroll_area)

        # Load saved settings
        self.load_settings()

    def _create_shot_controls(self, parent_layout):
        """Create shot input section"""
        group = QGroupBox("1. Shot Selection")
        grid = QGridLayout(group)

        # Shot entry with up/down buttons and Fetch button
        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_layout = QHBoxLayout()
        self.shot_entry = QLineEdit()
        self.shot_entry.setMinimumWidth(80)
        self.shot_entry.returnPressed.connect(self._fetch_available_tvs)
        self.shot_entry.textChanged.connect(self._on_shot_changed)
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

        fetch_btn = QPushButton('Fetch')
        fetch_btn.setFixedWidth(70)
        fetch_btn.clicked.connect(self._fetch_available_tvs)
        shot_layout.addWidget(fetch_btn)

        grid.addLayout(shot_layout, 0, 1)

        # TV Selection (TV01 / TV02) with Add button
        grid.addWidget(QLabel('TV'), 1, 0)

        tv_layout = QHBoxLayout()

        self.tv_button_group = QButtonGroup()
        self.tv01_radio = QRadioButton('TV01')
        self.tv01_radio.setChecked(True)
        self.tv01_radio.setEnabled(False)
        self.tv_button_group.addButton(self.tv01_radio)
        tv_layout.addWidget(self.tv01_radio)

        self.tv02_radio = QRadioButton('TV02')
        self.tv02_radio.setEnabled(False)
        self.tv_button_group.addButton(self.tv02_radio)
        tv_layout.addWidget(self.tv02_radio)

        self.add_btn = QPushButton('Add')
        self.add_btn.setFixedWidth(70)
        self.add_btn.setEnabled(False)
        self.add_btn.clicked.connect(self._add_shot)
        tv_layout.addWidget(self.add_btn)

        grid.addLayout(tv_layout, 1, 1)

        # "or" divider between Fetch+Add and Open File
        from ui.ui_constants import make_or_divider
        grid.addWidget(make_or_divider(), 2, 0, 1, 2)

        # Open Startup TV File button
        open_btn = QPushButton("Open Startup TV File...")
        open_btn.clicked.connect(self._open_startup_files)
        grid.addWidget(open_btn, 3, 0, 1, 2)

        parent_layout.addWidget(group)

    def _create_display_controls(self, parent_layout):
        """Create display options section with shot list and action buttons"""
        group = QGroupBox("2. Shot List")
        h_layout = QHBoxLayout(group)

        # Left: Shots list (QListWidget) showing shot + TV type
        self.shots_listbox = QListWidget()
        self.shots_listbox.setMaximumHeight(150)
        h_layout.addWidget(self.shots_listbox, 1)

        # Right: action buttons stacked vertically (compact)
        btn_widget = QWidget()
        btn_layout = QVBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        remove_btn = QPushButton('Remove')
        remove_btn.clicked.connect(self._remove_selected)
        btn_layout.addWidget(remove_btn)

        clear_btn = QPushButton('Clear All')
        clear_btn.clicked.connect(self._clear_all)
        btn_layout.addWidget(clear_btn)

        btn_layout.addStretch()

        hint_label = QLabel("Del key\nto remove")
        hint_label.setStyleSheet("color: gray; font-size: 8pt;")
        hint_label.setAlignment(Qt.AlignCenter)
        btn_layout.addWidget(hint_label)

        h_layout.addWidget(btn_widget, 1)

        parent_layout.addWidget(group)

    def _create_plot_controls(self, parent_layout):
        """Create plot section"""
        group = QGroupBox("3. Plot")
        grid = QGridLayout(group)

        # Frame range (start/end)
        grid.addWidget(QLabel('Frame Range'), 0, 0)

        range_layout = QHBoxLayout()

        self.frame_start_entry = QLineEdit('1')
        self.frame_start_entry.setFixedWidth(60)
        self.frame_start_entry.editingFinished.connect(self._sync_time_from_frame)
        range_layout.addWidget(self.frame_start_entry)

        range_layout.addWidget(QLabel('~'))

        self.frame_end_entry = QLineEdit('24')
        self.frame_end_entry.setFixedWidth(60)
        self.frame_end_entry.editingFinished.connect(self._sync_time_from_frame)
        range_layout.addWidget(self.frame_end_entry)
        range_layout.addStretch()

        grid.addLayout(range_layout, 0, 1)

        # Time range (synced with frame range)
        grid.addWidget(QLabel('Time Range [ms]'), 1, 0)

        time_range_layout = QHBoxLayout()

        self.time_start_entry = QLineEdit()
        self.time_start_entry.setFixedWidth(60)
        self.time_start_entry.editingFinished.connect(self._sync_frame_from_time)
        time_range_layout.addWidget(self.time_start_entry)

        time_range_layout.addWidget(QLabel('~'))

        self.time_end_entry = QLineEdit()
        self.time_end_entry.setFixedWidth(60)
        self.time_end_entry.editingFinished.connect(self._sync_frame_from_time)
        time_range_layout.addWidget(self.time_end_entry)

        time_range_layout.addStretch()

        grid.addLayout(time_range_layout, 1, 1)

        # Initialize time values from default frame values
        self._sync_time_from_frame()

        # Plot button
        plot_btn = QPushButton('Plot')
        plot_btn.clicked.connect(self._plot)
        grid.addWidget(plot_btn, 2, 0, 1, 2)

        parent_layout.addWidget(group)

    def _sync_time_from_frame(self):
        """Update time range entries from current frame range values"""
        try:
            frame_start = int(self.frame_start_entry.text().strip())
            frame_idx_start = frame_start - 1  # 0-based
            time_start_ms = frame_to_time(frame_idx_start) * 1000
            self.time_start_entry.setText(f'{time_start_ms:.1f}')
        except (ValueError, AttributeError):
            pass

        try:
            frame_end = int(self.frame_end_entry.text().strip())
            frame_idx_end = frame_end - 1  # 0-based
            time_end_ms = frame_to_time(frame_idx_end) * 1000
            self.time_end_entry.setText(f'{time_end_ms:.1f}')
        except (ValueError, AttributeError):
            pass

    def _sync_frame_from_time(self):
        """Update frame range entries from current time range values"""
        try:
            time_start_ms = float(self.time_start_entry.text().strip())
            time_start_sec = time_start_ms / 1000
            frame_num_start = max(1, round((time_start_sec + TV_OFFSET) * TV_FPS))
            self.frame_start_entry.setText(str(frame_num_start))
        except (ValueError, AttributeError):
            pass

        try:
            time_end_ms = float(self.time_end_entry.text().strip())
            time_end_sec = time_end_ms / 1000
            frame_num_end = max(1, round((time_end_sec + TV_OFFSET) * TV_FPS))
            self.frame_end_entry.setText(str(frame_num_end))
        except (ValueError, AttributeError):
            pass

    def _create_display(self, splitter):
        """Create matplotlib display area (no scroll)"""
        # Container widget
        container = QWidget()
        container_layout = QVBoxLayout(container)
        container_layout.setContentsMargins(0, 0, 0, 0)

        # Matplotlib figure
        self.figure = Figure(figsize=self.app_config.FIGURE_SIZE, tight_layout=False)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_axis_off()
        self.ax.set_label('TV Startup')
        self.ax.text(0.5, 0.5, 'Add shots and click Plot to compare startup sequences',
                     ha='center', va='center', fontsize=12,
                     transform=self.ax.transAxes, color='gray')
        apply_dark_figure_style(self.figure)

        # Matplotlib widget
        self.canvas = FigureCanvasQTAgg(self.figure)
        container_layout.addWidget(self.canvas)
        self.canvas.draw()

        splitter.addWidget(container)

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.text())
            new_shot = max(1, current + delta)
            self.shot_entry.setText(str(new_shot))
            # Disable Add and TV buttons when shot changes
            self._disable_add_controls()
        except ValueError:
            pass

    def _on_shot_changed(self, text=None):
        """Handle shot entry value change - disable Add if shot differs from last fetched"""
        try:
            current_shot = int(self.shot_entry.text().strip())
        except ValueError:
            current_shot = None

        # If shot changed from last fetched, disable Add button and TV radios
        if current_shot != self.last_fetched_shot:
            self._disable_add_controls()

    def _disable_add_controls(self):
        """Disable Add button and TV radio buttons"""
        self.add_btn.setEnabled(False)
        self.tv01_radio.setEnabled(False)
        self.tv02_radio.setEnabled(False)

    def _fetch_available_tvs(self):
        """Fetch available TV startup files for current shot"""
        try:
            shot = int(self.shot_entry.text().strip())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        # Store the fetched shot number
        self.last_fetched_shot = shot

        # Use tv_utils to find available TVs
        self.available_tvs = find_available_tvs(shot, startup=True)

        # Enable radio buttons based on availability (no TV type lock)
        self.tv01_radio.setEnabled('TV01' in self.available_tvs)
        self.tv02_radio.setEnabled('TV02' in self.available_tvs)

        # Auto-select first available if current selection is disabled
        selected_tv = 'TV01' if self.tv01_radio.isChecked() else 'TV02'
        if self.available_tvs and selected_tv not in self.available_tvs:
            if self.available_tvs[0] == 'TV01':
                self.tv01_radio.setChecked(True)
            else:
                self.tv02_radio.setChecked(True)

        # Enable Add button if selected TV is available
        selected_tv = 'TV01' if self.tv01_radio.isChecked() else 'TV02'
        if self.available_tvs and selected_tv in self.available_tvs:
            self.add_btn.setEnabled(True)
            self._set_status(f"Found: {', '.join(self.available_tvs)} for #{shot}", 'green')
            print(f"[TV Startup] Found {self.available_tvs} for shot #{shot}")
        else:
            self.add_btn.setEnabled(False)
            self._set_status(f"No TV startup data for #{shot}", 'red')
            print(f"[TV Startup] No TV startup data for shot #{shot}")

    def _add_time_label(self, canvas, frame_start, frame_end):
        """Add time labels to the top of canvas

        Uses frame_to_time (same as TV tab) for consistent time display.
        Frame 1 = first frame in ZIP (frame_idx=0).
        """
        h, w, _ = canvas.shape
        ncol = frame_end - frame_start + 1

        try:
            font = ImageFont.truetype("DejaVuSans.ttf", size=25)
        except Exception:
            font = ImageFont.load_default()

        img = Image.fromarray((canvas * 255).astype(np.uint8), mode='RGB')
        draw = ImageDraw.Draw(img)

        # Add label every 4 frames
        for i in range(0, ncol, 4):
            frame_num = frame_start + i
            # frame_num is 1-based UI frame, frame_idx is 0-based
            frame_idx = frame_num - 1
            time_sec = frame_to_time(frame_idx)
            time_ms = time_sec * 1000
            text = f'{time_ms:3.0f} ms'
            draw.text((i * self.FRAME_WIDTH + 10, 30), text, fill=(255, 255, 255), font=font)

        return np.asarray(img, dtype=np.float32) / 255.0

    def _add_shot_label(self, canvas, shot, tv_name):
        """Add shot number and TV name label to left side of canvas"""
        h, w, _ = canvas.shape

        img = Image.new('RGB', (self.LABEL_WIDTH, h), (0, 0, 0))
        draw = ImageDraw.Draw(img)

        text_shot = str(shot)
        text_tv = tv_name

        try:
            font_shot = ImageFont.truetype("DejaVuSans.ttf", size=max(16, h // 12))
            font_tv = ImageFont.truetype("DejaVuSans.ttf", size=max(12, h // 16))
        except Exception:
            font_shot = ImageFont.load_default()
            font_tv = font_shot

        # Measure text sizes
        bbox_shot = draw.textbbox((0, 0), text_shot, font=font_shot)
        tw_shot = bbox_shot[2] - bbox_shot[0]
        th_shot = bbox_shot[3] - bbox_shot[1]

        bbox_tv = draw.textbbox((0, 0), text_tv, font=font_tv)
        tw_tv = bbox_tv[2] - bbox_tv[0]
        th_tv = bbox_tv[3] - bbox_tv[1]

        spacing = 4
        total_h = th_shot + spacing + th_tv

        # Center both lines vertically
        y_start = (h - total_h) // 2
        x_shot = (self.LABEL_WIDTH - tw_shot) // 2
        x_tv = (self.LABEL_WIDTH - tw_tv) // 2

        draw.text((x_shot, y_start), text_shot, fill=(255, 255, 255), font=font_shot)
        draw.text((x_tv, y_start + th_shot + spacing), text_tv, fill=(180, 180, 180), font=font_tv)

        label_arr = np.asarray(img, dtype=np.float32) / 255.0
        return np.concatenate([label_arr, canvas], axis=1)

    def _load_startup_images(self, shot, tv_name, frame_start, frame_end):
        """Load startup images from ZIP file for specified frame range

        Frame numbering matches TV tab: frame 1 = first BMP in ZIP (index 0).
        """
        ncol = frame_end - frame_start + 1
        canvas = np.zeros((self.FRAME_HEIGHT, self.FRAME_WIDTH * ncol, 3), dtype=float)

        zip_name = self.custom_zip_paths.get((shot, tv_name))
        if zip_name is None:
            zip_name = get_tv_startup_zip_path(shot, tv_name)
        if zip_name is None or not os.path.exists(zip_name):
            raise FileNotFoundError(f"ZIP not found: {zip_name}")

        print(f"[TV Startup] Loading {zip_name} (frames {frame_start}-{frame_end})")

        with zipfile.ZipFile(zip_name, 'r') as zf:
            for i, frame_num in enumerate(range(frame_start, frame_end + 1)):
                # Frame 1 = BMP index 0 (first frame in ZIP)
                bmp_index = frame_num - 1
                bmp_name = f"{shot:06d}-{bmp_index:05d}.bmp"

                try:
                    data = zf.read(bmp_name)
                except KeyError:
                    continue

                # Read image from memory
                im = Image.open(io.BytesIO(data))
                out = im.resize((320, 320))
                rsize_arr = np.asarray(out, dtype=np.float32) / 255.0

                # Crop center portion
                crop = rsize_arr[:, 75:210, :]

                c = i * self.FRAME_WIDTH
                canvas[:, c:c + self.FRAME_WIDTH, :] = crop

        return canvas

    def _open_startup_files(self):
        """Open one or more startup ZIP files and add them to the shot list."""
        from PySide6.QtWidgets import QFileDialog
        import re

        start_dir = self._last_open_dir or os.path.expanduser("~")
        paths, _ = QFileDialog.getOpenFileNames(
            self.frame, "Open Startup TV ZIP", start_dir,
            "Startup ZIP files (*_startup.zip);;ZIP files (*.zip);;All files (*)"
        )
        if not paths:
            return

        self._last_open_dir = os.path.dirname(paths[0])

        added = 0
        for path in paths:
            basename = os.path.basename(path)
            # Expected pattern: {shot:06d}_tv{N}_startup.zip
            m = re.match(r'^(\d+)_tv(\d+)_startup\.zip$', basename, re.IGNORECASE)
            if not m:
                QMessageBox.warning(self.frame, "Invalid Filename",
                    f"Could not parse shot and TV from '{basename}'.\n"
                    f"Expected format: SHOT_tvN_startup.zip")
                continue
            shot = int(m.group(1))
            tv_name = f"TV{int(m.group(2)):02d}"

            if (shot, tv_name) in self.shot_list:
                continue

            self.shot_list.append((shot, tv_name))
            self.custom_zip_paths[(shot, tv_name)] = path
            self.shots_listbox.addItem(f"#{shot} ({tv_name})")
            added += 1

        if added:
            self._set_status(
                f"Opened {added} file(s). Total: {len(self.shot_list)} shots", 'green')
            print(f"[TV Startup] Opened {added} startup ZIP file(s)")

    def _add_shot(self):
        """Add shot to list (allows mixing TV01 and TV02)"""
        try:
            shot = int(self.shot_entry.text().strip())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        tv_name = 'TV01' if self.tv01_radio.isChecked() else 'TV02'

        # Check if same shot+TV combination already added
        if (shot, tv_name) in self.shot_list:
            QMessageBox.warning(self.frame, "Warning",
                                f"Shot #{shot} ({tv_name}) is already added")
            return

        # Add to list
        self.shot_list.append((shot, tv_name))
        self.shots_listbox.addItem(f"#{shot} ({tv_name})")

        # Disable Add and TV buttons after adding
        self._disable_add_controls()

        self._set_status(f"Added #{shot} ({tv_name}). Total: {len(self.shot_list)} shots", 'green')
        print(f"[TV Startup] Added shot #{shot} ({tv_name})")

    def _remove_selected(self):
        """Remove selected shot from list"""
        current_row = self.shots_listbox.currentRow()
        if current_row < 0:
            return

        shot, tv_name = self.shot_list[current_row]

        # Remove from storage
        del self.shot_list[current_row]
        self.custom_zip_paths.pop((shot, tv_name), None)
        self.shots_listbox.takeItem(current_row)

        if len(self.shot_list) == 0:
            self._set_status("Ready", 'gray')
        else:
            self._set_status(f"Removed #{shot}. Total: {len(self.shot_list)} shots", 'green')

        print(f"[TV Startup] Removed shot #{shot} ({tv_name})")

    def _clear_all(self):
        """Clear all shots"""
        self.shot_list = []
        self.custom_zip_paths = {}
        self.last_canvas = None
        self.shots_listbox.clear()

        # Reset display
        self.ax.clear()
        self.ax.set_axis_off()
        self.ax.text(0.5, 0.5, 'Add shots and click Plot to compare startup sequences',
                     ha='center', va='center', fontsize=12,
                     transform=self.ax.transAxes, color='gray')

        self.canvas.draw()
        self._set_status("Cleared", 'gray')

        print("[TV Startup] Cleared all shots")

    def _plot(self):
        """Generate and display the comparison plot"""
        if not self.shot_list:
            QMessageBox.warning(self.frame, "Warning", "Please add at least one shot")
            return

        # Validate frame range
        try:
            frame_start = int(self.frame_start_entry.text().strip())
            frame_end = int(self.frame_end_entry.text().strip())
            if frame_start < 1:
                frame_start = 1
            if frame_end < frame_start:
                raise ValueError("End frame must be >= start frame")
            # Limit to max 50 frames
            if frame_end - frame_start + 1 > 50:
                frame_end = frame_start + 49
                self.frame_end_entry.setText(str(frame_end))
                QMessageBox.information(self.frame, "Info", "Frame range limited to 50 frames maximum")
        except ValueError as e:
            QMessageBox.critical(self.frame, "Error", f"Invalid frame range: {e}")
            return

        ncol = frame_end - frame_start + 1

        self._set_status("Loading images...", 'blue')
        QApplication.processEvents()

        images = []
        try:
            for i, (shot, tv_name) in enumerate(self.shot_list):
                # Load images for this shot with its own TV type
                canvas = self._load_startup_images(shot, tv_name, frame_start, frame_end)

                # Add time labels only for first shot
                if i == 0:
                    canvas = self._add_time_label(canvas, frame_start, frame_end)

                # Add shot + TV name label
                canvas = self._add_shot_label(canvas, shot, tv_name)

                images.append(canvas)

            # Stack all canvases vertically
            stacked = np.vstack(images)
            self.last_canvas = stacked

            # Update plot
            self.ax.clear()
            self.ax.imshow(stacked)
            self.ax.set_axis_off()

            self.figure.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)
            self.canvas.draw()

            self._set_status(f"Plotted {len(self.shot_list)} shots", 'green')
            print(f"[TV Startup] Plotted {len(self.shot_list)} shots")

        except FileNotFoundError as e:
            self._set_status("File not found", 'red')
            QMessageBox.critical(self.frame, "Error", str(e))
        except Exception as e:
            self._set_status("Error", 'red')
            QMessageBox.critical(self.frame, "Error", f"Failed to plot:\n{str(e)}")

    def _set_status(self, text, color='gray', timeout=0):
        from ui.ui_constants import show_status
        show_status(self.frame, "TV Startup", text, color, timeout)

    def cleanup(self):
        """Cleanup resources when tab is closed"""
        self.shot_list = []
        self.last_canvas = None

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.text(),
            "frame_start": self.frame_start_entry.text(),
            "frame_end": self.frame_end_entry.text(),
            "tv": 'TV01' if self.tv01_radio.isChecked() else 'TV02'
        }
        set_tab_settings("tv_startup", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("tv_startup")

        if settings.get("shot"):
            self.shot_entry.setText(settings["shot"])

        if settings.get("frame_start"):
            self.frame_start_entry.setText(settings["frame_start"])

        if settings.get("frame_end"):
            self.frame_end_entry.setText(settings["frame_end"])

        if settings.get("tv"):
            if settings["tv"] == 'TV02':
                self.tv02_radio.setChecked(True)
            else:
                self.tv01_radio.setChecked(True)

        # Sync time range after restoring frame range
        self._sync_time_from_frame()
