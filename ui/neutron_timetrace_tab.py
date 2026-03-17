"""
Neutron Time Trace tab
Fusion Neutron Diagnostics (Near J-port): Fission Chamber, He3 Counter, Diamond Detector
"""

import os
import tempfile
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QComboBox,
    QListWidget, QAbstractItemView, QMessageBox, QFileDialog,
    QSplitter, QScrollArea, QFrame, QStyle, QSpinBox,
    QDialog, QDialogButtonBox, QTableWidget, QTableWidgetItem,
    QHeaderView, QApplication,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QKeySequence, QShortcut

from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style, get_icon
from ui.widgets.custom_toolbar import QuietNavigationToolbar
from ui.theme import ThemeManager
from config.user_settings import get_tab_settings, set_tab_settings
from plotting.plot_manager import ColorManager


class NeutronTimeTraceTab:
    """Neutron Time Trace tab with 3x1 subplot layout and shot overplot"""

    DETECTORS = ['fission', 'he3', 'diamond']
    DETECTOR_NAMES = {
        'fission': 'Fission Chamber',
        'he3': 'He3 Counter',
        'diamond': 'Diamond Detector',
    }
    DETECTOR_YLABELS = {
        'fission': 'Fission [a.u.]',
        'he3': 'He3 [a.u.]',
        'diamond': 'Diamond [a.u.]',
    }

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = QWidget()
        self.toolbar = None

        # Data cache: {shot_number: DiagnosticData}
        self.data = {}

        # Color manager
        self._color_manager = ColorManager()

        # Plot style defaults
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

    def create_widgets(self):
        """Create tab UI — canvas left, controls right"""
        # Canvas + figure (3x1 subplots)
        self.figure = Figure(figsize=(8, 10))

        self.ax_fission = self.figure.add_subplot(311)
        self.ax_he3 = self.figure.add_subplot(312, sharex=self.ax_fission)
        self.ax_diamond = self.figure.add_subplot(313, sharex=self.ax_fission)

        self.axes = {
            'fission': self.ax_fission,
            'he3': self.ax_he3,
            'diamond': self.ax_diamond,
        }

        # Set initial labels
        self.figure.suptitle('Fusion Neutron Diagnostics (Near J-port)',
                             fontsize=self.label_fontsize)
        for key, ax in self.axes.items():
            ax.set_ylabel(self.DETECTOR_YLABELS[key], fontsize=self.label_fontsize)
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.grid(True, linestyle='--', linewidth=0.3)
            ax.set_label(self.DETECTOR_NAMES[key])

        self.ax_diamond.set_xlabel('Time [s]', fontsize=self.label_fontsize)
        self.figure.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.07,
                                    hspace=0.20)
        apply_dark_figure_style(self.figure)

        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.draw()

        # Left side: canvas + toolbar
        canvas_widget = QWidget()
        canvas_layout = QVBoxLayout(canvas_widget)
        canvas_layout.setContentsMargins(0, 0, 0, 0)
        canvas_layout.addWidget(self.canvas)

        # Right side: scrollable control panel
        scroll_area = QScrollArea()
        scroll_area.setFixedWidth(CONTROL_PANEL_WIDTH)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)

        control_frame = QWidget()
        control_layout = QVBoxLayout(control_frame)
        control_layout.setContentsMargins(9, 9, 9, 9)
        control_layout.setSizeConstraint(QVBoxLayout.SetMinimumSize)

        scroll_area.setWidget(control_frame)
        scroll_area.viewport().setAutoFillBackground(False)
        control_frame.setAutoFillBackground(False)

        # Splitter for left/right layout
        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(canvas_widget)
        splitter.addWidget(scroll_area)

        main_layout = QVBoxLayout(self.frame)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.addWidget(splitter)

        # Create control panels
        self._create_shot_input(control_frame)
        self._create_shot_list(control_frame)
        self._create_plot_controls(control_frame)
        self._create_save_controls(control_frame)
        control_layout.addStretch()

        # Load settings
        self.load_settings()

    # ------------------------------------------------------------------
    # 1. Load Neutron Data
    # ------------------------------------------------------------------

    def _create_shot_input(self, parent):
        """Create shot input section"""
        group = QGroupBox("1. Load Neutron Data")
        grid = QGridLayout(group)
        grid.setColumnStretch(1, 1)

        grid.addWidget(QLabel('Shot'), 0, 0)

        shot_frame = QWidget()
        shot_layout = QHBoxLayout(shot_frame)
        shot_layout.setContentsMargins(0, 0, 0, 0)

        self.shot_entry = QLineEdit()
        shot_layout.addWidget(self.shot_entry, 1)
        self.shot_entry.returnPressed.connect(self.load_shot_data)

        btn_widget = QWidget()
        btn_layout = QVBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.setSpacing(0)

        mini_btn_style = "padding: 0px; border-radius: 2px;"
        up_btn = QPushButton()
        up_btn.setIcon(get_icon(QStyle.SP_ArrowUp))
        up_btn.setFixedSize(24, 15)
        up_btn.setStyleSheet(mini_btn_style)
        up_btn.clicked.connect(lambda: self._adjust_shot(1))
        btn_layout.addWidget(up_btn)

        down_btn = QPushButton()
        down_btn.setIcon(get_icon(QStyle.SP_ArrowDown))
        down_btn.setFixedSize(24, 15)
        down_btn.setStyleSheet(mini_btn_style)
        down_btn.clicked.connect(lambda: self._adjust_shot(-1))
        btn_layout.addWidget(down_btn)

        shot_layout.addWidget(btn_widget)
        grid.addWidget(shot_frame, 0, 1)

        self.fetch_button = QPushButton('Fetch')
        self.fetch_button.setFixedWidth(70)
        self.fetch_button.clicked.connect(self.load_shot_data)
        grid.addWidget(self.fetch_button, 0, 2)

        parent.layout().addWidget(group)

    def _adjust_shot(self, delta):
        """Adjust shot number"""
        try:
            current = int(self.shot_entry.text())
            self.shot_entry.setText(str(max(1, current + delta)))
        except ValueError:
            pass

    # ------------------------------------------------------------------
    # 2. Shot List (like TV Startup)
    # ------------------------------------------------------------------

    def _create_shot_list(self, parent):
        """Create shot list section — Fetch adds directly, Remove/Clear to manage"""
        group = QGroupBox("2. Shot List")
        h_layout = QHBoxLayout(group)

        # Left: shot list
        self.shots_listbox = QListWidget()
        self.shots_listbox.setMaximumHeight(150)
        h_layout.addWidget(self.shots_listbox, 1)

        # Right: action buttons
        btn_widget = QWidget()
        btn_layout = QVBoxLayout(btn_widget)
        btn_layout.setContentsMargins(0, 0, 0, 0)

        remove_btn = QPushButton('Remove')
        remove_btn.clicked.connect(self._remove_shot)
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

        # Keyboard shortcut
        del_shortcut = QShortcut(QKeySequence(Qt.Key_Delete), self.shots_listbox)
        del_shortcut.activated.connect(self._remove_shot)
        bs_shortcut = QShortcut(QKeySequence(Qt.Key_Backspace), self.shots_listbox)
        bs_shortcut.activated.connect(self._remove_shot)

        parent.layout().addWidget(group)

    def _remove_shot(self):
        """Remove selected shot from list"""
        current_row = self.shots_listbox.currentRow()
        if current_row < 0:
            return
        self.shots_listbox.takeItem(current_row)

    def _clear_all(self):
        """Clear all shots from list"""
        self.shots_listbox.clear()

    # ------------------------------------------------------------------
    # 3. Plot
    # ------------------------------------------------------------------

    def _create_plot_controls(self, parent):
        """Create plot control buttons — Plot + Option in same row"""
        group = QGroupBox("3. Plot")
        group_layout = QVBoxLayout(group)

        plot_row = QHBoxLayout()

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        plot_row.addWidget(plot_button, 3)

        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_style_dialog)
        plot_row.addWidget(style_btn, 1)

        group_layout.addLayout(plot_row)

        # Hidden combo for color mode (for settings save/restore)
        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        self.color_mode_combo.setCurrentText("Gradient(viridis)")
        self.color_mode_combo.hide()

        parent.layout().addWidget(group)

    def _show_style_dialog(self):
        """Show Plot Options settings dialog"""
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
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        color_combo.setCurrentText(self.color_mode_combo.currentText())
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
        btn_box = QDialogButtonBox(
            QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)

        def reset_defaults():
            color_combo.setCurrentText("Gradient(viridis)")
            label_spin.setValue(12)
            legend_spin.setValue(8)
            tick_spin.setValue(10)
        btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(reset_defaults)
        dlg_layout.addWidget(btn_box)

        if dialog.exec() == QDialog.Accepted:
            self.color_mode_combo.setCurrentText(color_combo.currentText())
            self.label_fontsize = label_spin.value()
            self.legend_fontsize = legend_spin.value()
            self.tick_fontsize = tick_spin.value()
            # Auto-apply if plot exists
            if any(ax.lines for ax in self.axes.values()):
                self.plot_data()

    # ------------------------------------------------------------------
    # 4. Save Data
    # ------------------------------------------------------------------

    def _create_save_controls(self, parent):
        """Create save controls — Preview & Save button"""
        frame = QGroupBox("4. Save Data")
        layout = QVBoxLayout(frame)

        preview_button = QPushButton('Preview && Save')
        preview_button.clicked.connect(self.preview_data)
        layout.addWidget(preview_button)

        parent.layout().addWidget(frame)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _get_data_loader(self):
        """Lazy-load neutron loader"""
        if not hasattr(self, '_loader'):
            from data_loaders.neutron_loader import NeutronLoader
            self._loader = NeutronLoader(self.app_config, {
                'name': 'Neutron',
                'mds_tree': 'kstar',
            })
        return self._loader

    def load_shot_data(self):
        """Load neutron data — add to shot list on success"""
        try:
            shot_number = int(self.shot_entry.text())
        except ValueError:
            QMessageBox.critical(self.frame, "Error", "Please enter a valid shot number")
            return

        # Check duplicate
        shot_str = f'#{shot_number}'
        for i in range(self.shots_listbox.count()):
            if self.shots_listbox.item(i).text() == shot_str:
                QMessageBox.warning(self.frame, "Warning",
                                    f"Shot {shot_number} is already in the list")
                return

        try:
            self.fetch_button.setText('Loading...')
            self.fetch_button.setEnabled(False)
            QApplication.setOverrideCursor(Qt.WaitCursor)
            QApplication.processEvents()

            loader = self._get_data_loader()
            data = loader.load_data(shot_number)
            self.data[shot_number] = data

            # Add to shot list directly
            self.shots_listbox.addItem(shot_str)

            det_names = [self.DETECTOR_NAMES[k]
                         for k in self.DETECTORS if k in data.measurements]
            print(f"[Neutron] Shot {shot_number}: {', '.join(det_names)}")

        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to load data: {str(e)}")
        finally:
            self.fetch_button.setText('Fetch')
            self.fetch_button.setEnabled(True)
            QApplication.restoreOverrideCursor()

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_color_mode(text):
        """Extract colormap name from dropdown text like 'Fixed(tab10)' -> 'tab10'"""
        start = text.find('(')
        end = text.find(')')
        if start != -1 and end != -1:
            return text[start + 1:end]
        return 'tab10'

    def _get_plot_colors(self, entries):
        """Get plot colors for entries using selected color mode"""
        mode = self.color_mode_combo.currentText()
        cmap_name = self._parse_color_mode(mode)
        return self._color_manager.get_colors_for_entries(entries, colormap=cmap_name)

    def plot_data(self):
        """Plot neutron time traces"""
        n_shots = self.shots_listbox.count()
        if n_shots == 0:
            return

        entries = [self.shots_listbox.item(i).text() for i in range(n_shots)]

        for ax in self.axes.values():
            ax.clear()

        colors = self._get_plot_colors(entries)

        for i, entry in enumerate(entries):
            shot_number = int(entry.replace('#', ''))
            color = colors[i]

            if shot_number not in self.data:
                try:
                    loader = self._get_data_loader()
                    self.data[shot_number] = loader.load_data(shot_number)
                except Exception as e:
                    print(f"[Neutron] Error loading shot {shot_number}: {e}")
                    continue

            data = self.data[shot_number]

            for key in self.DETECTORS:
                ax = self.axes[key]
                if key in data.measurements:
                    meas = data.measurements[key]
                    ax.plot(meas['time'], meas['data'], '-',
                            color=color, linewidth=0.8,
                            label=f'#{shot_number}')

        # Finalize axes
        for key, ax in self.axes.items():
            ax.set_ylabel(self.DETECTOR_YLABELS[key], fontsize=self.label_fontsize)
            ax.set_ylim(bottom=0)
            ax.tick_params(labelsize=self.tick_fontsize)
            ax.grid(True, linestyle='--', linewidth=0.3)
            ax.set_label(self.DETECTOR_NAMES[key])

            if ax.get_legend_handles_labels()[1]:
                ax.legend(fontsize=self.legend_fontsize, loc='best', frameon=False)

        self.ax_diamond.set_xlabel('Time [s]', fontsize=self.label_fontsize)

        self.figure.suptitle('Fusion Neutron Diagnostics (Near J-port)',
                             fontsize=self.label_fontsize)
        self.figure.subplots_adjust(left=0.10, right=0.98, top=0.93, bottom=0.07,
                                    hspace=0.20)

        ThemeManager.apply_theme_to_figure(self.figure)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    # ------------------------------------------------------------------
    # Data export (Preview & Save)
    # ------------------------------------------------------------------

    def _get_selected_shots(self):
        """Get shot numbers from shot list"""
        shots = []
        for i in range(self.shots_listbox.count()):
            text = self.shots_listbox.item(i).text()
            shots.append(int(text.replace('#', '')))
        return shots

    def _write_data_to_file(self, file_path, shots):
        """Write neutron data to text file"""
        with open(file_path, 'w') as f:
            f.write("# Fusion Neutron Diagnostics (Near J-port)\n")
            f.write("#%10s,%15s,%15s,%15s\n" % (
                "Shot", "Detector", "Time[s]", "Intensity[a.u.]"))

            for shot_number in shots:
                if shot_number not in self.data:
                    continue
                data = self.data[shot_number]

                for key in self.DETECTORS:
                    if key not in data.measurements:
                        continue
                    meas = data.measurements[key]
                    for t, d in zip(meas['time'], meas['data']):
                        f.write(" %10d,%15s,%15.6f,%15.6e\n" % (
                            shot_number, key, t, d))

    def preview_data(self):
        """Preview data in a spreadsheet dialog"""
        shots = self._get_selected_shots()
        if not shots:
            QMessageBox.warning(self.frame, "Warning", "No shots in the list")
            return

        try:
            fd, tmp_path = tempfile.mkstemp(suffix='.txt')
            os.close(fd)
            self._write_data_to_file(tmp_path, shots)

            with open(tmp_path, 'r') as f:
                lines = f.readlines()
            os.remove(tmp_path)

            self._show_data_preview_dialog(lines)
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to preview data: {str(e)}")

    def _show_data_preview_dialog(self, lines):
        """Show data in a spreadsheet-style dialog"""
        headers = []
        data_rows = []
        title_line = ""

        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                content = stripped.lstrip('#').strip()
                if ',' in content:
                    headers = [h.strip() for h in content.split(',')]
                else:
                    title_line = content
            else:
                cells = [c.strip() for c in stripped.split(',')]
                data_rows.append(cells)

        if not data_rows:
            QMessageBox.information(self.frame, "Preview", "No data to display")
            return

        n_cols = len(headers) if headers else len(data_rows[0])
        n_rows = len(data_rows)

        dialog = QDialog(self.frame)
        dialog.setWindowTitle(f"Data Preview — {title_line}" if title_line else "Data Preview")
        dialog.resize(min(120 * n_cols, 1200), min(30 * n_rows + 100, 700))

        layout = QVBoxLayout(dialog)

        info = QLabel(f"{n_rows} rows \u00d7 {n_cols} columns")
        info.setStyleSheet("color: gray;")
        layout.addWidget(info)

        table = QTableWidget(n_rows, n_cols)
        table.setFont(QFont('Courier', 9))
        table.setEditTriggers(QTableWidget.NoEditTriggers)

        if headers:
            table.setHorizontalHeaderLabels(headers)
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        table.verticalHeader().setDefaultSectionSize(24)

        for r, row in enumerate(data_rows):
            for c, val in enumerate(row):
                table.setItem(r, c, QTableWidgetItem(val))

        table.resizeColumnsToContents()
        for c in range(n_cols):
            if table.columnWidth(c) > 150:
                table.setColumnWidth(c, 150)

        layout.addWidget(table)

        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        save_btn = QPushButton("Save as .csv")
        save_btn.clicked.connect(lambda: (dialog.close(), self._save_data_dialog()))
        btn_layout.addWidget(save_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)

        layout.addLayout(btn_layout)
        dialog.exec()

    def _save_data_dialog(self):
        """Save data via file dialog"""
        shots = self._get_selected_shots()
        if not shots:
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self.frame, "Save Data", os.path.expanduser("~"),
            "CSV files (*.csv);;Text files (*.txt);;All files (*.*)")
        if not file_path:
            return

        try:
            self._write_data_to_file(file_path, shots)
            print(f"[Neutron] Data saved to {file_path}")
            QMessageBox.information(self.frame, "Success", f"Data saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to save data: {str(e)}")

    # ------------------------------------------------------------------
    # Settings
    # ------------------------------------------------------------------

    def save_settings(self):
        """Save tab settings"""
        settings = {
            "shot": self.shot_entry.text(),
            "color_mode": self.color_mode_combo.currentText(),
            "label_fontsize": self.label_fontsize,
            "legend_fontsize": self.legend_fontsize,
            "tick_fontsize": self.tick_fontsize,
        }
        set_tab_settings("neutron", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("neutron")

        if settings.get("shot"):
            self.shot_entry.setText(settings["shot"])
        if settings.get("color_mode"):
            self.color_mode_combo.setCurrentText(settings["color_mode"])
        if settings.get("label_fontsize"):
            self.label_fontsize = settings["label_fontsize"]
        if settings.get("legend_fontsize"):
            self.legend_fontsize = settings["legend_fontsize"]
        if settings.get("tick_fontsize"):
            self.tick_fontsize = settings["tick_fontsize"]
