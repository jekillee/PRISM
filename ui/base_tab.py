"""
Abstract base class for all diagnostic tabs
Contains only common functionality shared across all diagnostics
"""

from abc import ABC, abstractmethod
import os
import tempfile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.text import Annotation

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QGroupBox, QLabel, QLineEdit, QPushButton, QComboBox,
    QListWidget, QAbstractItemView, QCheckBox, QRadioButton,
    QButtonGroup, QMessageBox, QFileDialog, QApplication, QStyle,
    QDialog, QTableWidget, QTableWidgetItem, QHeaderView
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont

from ui.ui_constants import CONTROL_PANEL_WIDTH, get_icon


class BaseTab(ABC):
    """Abstract base class for all diagnostic tabs"""

    # Height of Select Data listboxes (available / selected)
    LISTBOX_HEIGHT = 180

    # Mapping from diagnostic_name to settings key prefix
    _SETTINGS_KEY_MAP = {
        'CES': 'tivt',
        'Thomson': 'nete',
        'MSE': 'mse',
    }

    def __init__(self, parent, app_config, diagnostic_name, tab_type,
                 data_loader, efit_loader, plot_manager):
        self.parent = parent
        self.app_config = app_config
        self.diagnostic_name = diagnostic_name
        self.tab_type = tab_type
        self.data_loader = data_loader
        self.efit_loader = efit_loader
        self.plot_manager = plot_manager

        self.frame = QWidget(parent)
        self.data = {}
        self.efit_data = {}
        self.computed_efit_tree = None
        self.toolbar = None

        from config.diagnostic_config import DIAGNOSTICS
        self.diag_config = DIAGNOSTICS[diagnostic_name]

        self.param1 = self.diag_config['parameters'][0]
        if len(self.diag_config['parameters']) > 1:
            self.param2 = self.diag_config['parameters'][1]
        else:
            self.param2 = None

        # Determine settings key for this tab (e.g. "mse_profile", "tivt_timetrace")
        prefix = self._SETTINGS_KEY_MAP.get(diagnostic_name, diagnostic_name.lower())
        self._settings_key = f"{prefix}_{tab_type}"

    def _restore_shot_from_settings(self):
        """Restore shot number and UI state from saved settings (called after create_widgets)"""
        from config.user_settings import get_tab_settings
        tab_settings = get_tab_settings(self._settings_key)
        saved_shot = tab_settings.get("shot", "")
        if saved_shot and hasattr(self, 'shot_entry'):
            self.shot_entry.setText(str(saved_shot))
        # Restore color mode
        if hasattr(self, 'color_mode_combo'):
            saved_color = tab_settings.get("color_mode", "Fixed(tab10)")
            idx = self.color_mode_combo.findText(saved_color)
            if idx >= 0:
                self.color_mode_combo.setCurrentIndex(idx)
        # Restore show nodes checkbox
        if hasattr(self, 'show_channel_checkbox'):
            saved_show = tab_settings.get("show_nodes", False)
            self.show_channel_checkbox.setChecked(saved_show)
        # Restore font sizes
        if hasattr(self, 'label_fontsize'):
            self.label_fontsize = tab_settings.get("label_fontsize", 12)
        if hasattr(self, 'legend_fontsize'):
            self.legend_fontsize = tab_settings.get("legend_fontsize", 8)
        if hasattr(self, 'tick_fontsize'):
            self.tick_fontsize = tab_settings.get("tick_fontsize", 10)

    def save_settings(self):
        """Save current tab state to settings (called on app close)"""
        from config.user_settings import get_tab_settings, set_tab_settings
        tab_settings = get_tab_settings(self._settings_key)
        if hasattr(self, 'shot_entry'):
            tab_settings["shot"] = self.shot_entry.text()
        if hasattr(self, 'color_mode_combo'):
            tab_settings["color_mode"] = self.color_mode_combo.currentText()
        if hasattr(self, 'show_channel_checkbox'):
            tab_settings["show_nodes"] = self.show_channel_checkbox.isChecked()
        if hasattr(self, 'label_fontsize'):
            tab_settings["label_fontsize"] = self.label_fontsize
        if hasattr(self, 'legend_fontsize'):
            tab_settings["legend_fontsize"] = self.legend_fontsize
        if hasattr(self, 'tick_fontsize'):
            tab_settings["tick_fontsize"] = self.tick_fontsize
        set_tab_settings(self._settings_key, tab_settings)

    @abstractmethod
    def create_widgets(self):
        """Create tab widgets - must be implemented by each diagnostic tab"""
        pass

    @abstractmethod
    def load_shot_data(self):
        """Load shot data - must be implemented by each diagnostic tab"""
        pass

    @abstractmethod
    def plot_data(self):
        """Plot data - must be implemented by each diagnostic tab"""
        pass

    # ===== Common widget creation methods =====

    def _create_selection_listboxes(self, parent):
        """Create selection listboxes (common for all tabs)"""
        outer_frame = QGroupBox("2. Select Data", parent)
        outer_layout = QVBoxLayout(outer_frame)

        # Main content layout: available list | buttons | selected list
        content_layout = QHBoxLayout()
        outer_layout.addLayout(content_layout)

        # Available list column
        available_column = QVBoxLayout()
        if self.tab_type == 'profile':
            label_text = 'Time [ms] (Diag)'
        else:
            label_text = 'R [mm] (Ch)'
        available_column.addWidget(QLabel(label_text))

        self.available_listbox = QListWidget()
        self.available_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.available_listbox.setFixedHeight(self.LISTBOX_HEIGHT)
        available_column.addWidget(self.available_listbox)
        content_layout.addLayout(available_column, stretch=1)

        # Button column
        button_layout = QVBoxLayout()
        button_layout.addStretch()
        add_button = QPushButton()
        add_button.setIcon(get_icon(QStyle.SP_ArrowForward))
        add_button.setFixedWidth(30)
        add_button.clicked.connect(self.add_selected_items)
        button_layout.addWidget(add_button)
        remove_button = QPushButton()
        remove_button.setIcon(get_icon(QStyle.SP_ArrowBack))
        remove_button.setFixedWidth(30)
        remove_button.clicked.connect(self.remove_selected_items)
        button_layout.addWidget(remove_button)
        button_layout.addStretch()
        content_layout.addLayout(button_layout)

        # Selected list column
        selected_column = QVBoxLayout()
        selected_column.addWidget(QLabel('Selected'))

        self.selected_listbox = QListWidget()
        self.selected_listbox.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.selected_listbox.setFixedHeight(self.LISTBOX_HEIGHT)
        selected_column.addWidget(self.selected_listbox)
        content_layout.addLayout(selected_column, stretch=1)

        # Keyboard shortcuts for selected listbox
        from PySide6.QtGui import QShortcut, QKeySequence
        delete_shortcut = QShortcut(QKeySequence(Qt.Key_Delete), self.selected_listbox)
        delete_shortcut.activated.connect(self.remove_selected_items)
        backspace_shortcut = QShortcut(QKeySequence(Qt.Key_Backspace), self.selected_listbox)
        backspace_shortcut.activated.connect(self.remove_selected_items)

        # Add the outer frame to the parent's layout
        parent_layout = parent.layout()
        if parent_layout is not None:
            parent_layout.addWidget(outer_frame)

        return outer_frame

    def _create_efit_controls(self, parent):
        """Create EFIT mapping controls (common for profile tabs)"""
        self.selected_x_axis_value = "R"

        frame = QGroupBox("4. EFIT Mapping (Optional)", parent)
        grid = QGridLayout(frame)

        # EFIT Tree label and dropdown
        grid.addWidget(QLabel("EFIT Tree"), 0, 0)

        efit_display_values = list(self.app_config.EFIT_TREES.keys())

        self.efit_dropdown = QComboBox()
        self.efit_dropdown.addItems(efit_display_values)
        self.efit_dropdown.setCurrentIndex(0)
        grid.addWidget(self.efit_dropdown, 0, 1, 1, 2)

        mapping_button = QPushButton('Mapping')
        mapping_button.clicked.connect(self.compute_efit)
        grid.addWidget(mapping_button, 0, 3)

        # Radio buttons for x-axis selection
        self.x_axis_button_group = QButtonGroup()

        radio_psi_n = QRadioButton("\u03C8\u2099")     # ψₙ (psi_N)
        radio_rho_pol = QRadioButton("\u03C1\u209A\u2092\u2097")   # ρₚₒₗ (rho_pol)
        radio_rho_tor = QRadioButton("\u03C1\u209C\u2092\u1D63")   # ρₜₒᵣ (rho_tor)

        self.x_axis_button_group.addButton(radio_psi_n)
        self.x_axis_button_group.addButton(radio_rho_pol)
        self.x_axis_button_group.addButton(radio_rho_tor)

        # Store mapping from button to value string
        radio_psi_n.setProperty("axis_value", "psi_N")
        radio_rho_pol.setProperty("axis_value", "rho_pol")
        radio_rho_tor.setProperty("axis_value", "rho_tor")

        # Connect to update the stored value
        self.x_axis_button_group.buttonClicked.connect(self._on_x_axis_changed)

        grid.addWidget(radio_psi_n, 1, 0)
        grid.addWidget(radio_rho_pol, 1, 1)
        grid.addWidget(radio_rho_tor, 1, 2)

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_efit_profiles)
        grid.addWidget(plot_button, 1, 3)

        # Set equal column stretch
        for i in range(4):
            grid.setColumnStretch(i, 1)

        # Add the frame to the parent's layout
        parent_layout = parent.layout()
        if parent_layout is not None:
            parent_layout.addWidget(frame)

        return frame

    def _on_x_axis_changed(self, button):
        """Handle x-axis radio button change"""
        self.selected_x_axis_value = button.property("axis_value")

    def _get_selected_x_axis(self):
        """Get the currently selected x-axis value from the radio buttons"""
        checked = self.x_axis_button_group.checkedButton()
        if checked is not None:
            return checked.property("axis_value")
        return "R"

    def _create_save_controls(self, parent, section_num=None):
        """Create save controls (common for all tabs)"""
        if section_num:
            frame = QGroupBox(f"{section_num}. Save Data", parent)
        else:
            frame = QGroupBox("Save Data", parent)
        layout = QVBoxLayout(frame)

        preview_button = QPushButton('Preview && Save')
        preview_button.clicked.connect(self.preview_data)
        layout.addWidget(preview_button)

        # Add the frame to the parent's layout
        parent_layout = parent.layout()
        if parent_layout is not None:
            parent_layout.addWidget(frame)

        return frame

    # ===== Loading state helper methods =====

    def _set_loading_state(self, is_loading, button=None, original_text='Fetch'):
        """Set loading state with visual feedback"""
        if button is None and hasattr(self, 'fetch_button'):
            button = self.fetch_button

        if button is None:
            return

        if is_loading:
            button.setText('Loading...')
            button.setEnabled(False)
            QApplication.setOverrideCursor(Qt.WaitCursor)
            QApplication.processEvents()
        else:
            button.setText(original_text)
            button.setEnabled(True)
            QApplication.restoreOverrideCursor()

    # ===== Channel label helper methods =====

    def _add_channel_labels(self, ax, x_data, y_data, node_prefix, channels):
        """Add channel labels above data points if show_channel is enabled

        node_prefix: str like 'CES_TI', 'CESNN_TI', 'TS_CORE', 'TS_EDGE', 'ECE',
                     'TGAMMA', 'pmse_qv', 'pmse_jv', 'TXCS_TI'
        channels: list of channel numbers (int)
        """
        if not hasattr(self, 'show_channel_checkbox') or not self.show_channel_checkbox.isChecked():
            return

        for x, y, ch in zip(x_data, y_data, channels):
            label = f'{node_prefix}{ch:02d}'
            ax.annotate(label, (x, y), textcoords='offset points',
                       xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                       clip_on=True, annotation_clip=True)

    def _clear_channel_labels(self, ax):
        """Clear existing channel label annotations"""
        for child in list(ax.get_children()):
            if isinstance(child, Annotation):
                child.remove()

    # ===== Common action methods =====

    def add_selected_items(self):
        """Add selected items to the selected list"""
        existing_items = [self.selected_listbox.item(i).text()
                         for i in range(self.selected_listbox.count())]
        for item in self.available_listbox.selectedItems():
            if item.text() not in existing_items:
                self.selected_listbox.addItem(item.text())

    def remove_selected_items(self):
        """Remove selected items from the selected list"""
        for item in reversed(self.selected_listbox.selectedItems()):
            self.selected_listbox.takeItem(self.selected_listbox.row(item))

    def compute_efit(self):
        """Compute EFIT equilibrium data (common logic)"""
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            QMessageBox.warning(self.frame, "Warning", "No timepoints selected")
            return

        display_name = self.efit_dropdown.currentText()
        efit_tree = self.app_config.EFIT_TREES[display_name]

        self.efit_data.clear()

        print("\n" + "="*50)
        print(f"[EFIT] Mapping ({efit_tree})")
        print("-" * 50)

        # Group by shot
        shot_groups = {}
        for entry in selected_entries:
            shot_number, param_value = self._parse_entry_for_efit(entry)
            if shot_number not in shot_groups:
                shot_groups[shot_number] = []
            shot_groups[shot_number].append((entry, param_value))

        valid_entries = []

        for shot_number, entry_list in shot_groups.items():
            try:
                print(f"[EFIT] Loading data for shot #{shot_number}...")
                efit_data = self.efit_loader.load_efit_data(shot_number, efit_tree)

                for entry, time_point in entry_list:
                    if not (np.min(efit_data.time) <= time_point <= np.max(efit_data.time)):
                        print(f"[EFIT]   Skipping {entry}: out of range")
                        continue

                    time_idx = np.argmin(np.abs(efit_data.time - time_point))
                    closest_time = efit_data.time[time_idx]

                    if abs(closest_time - time_point) > 0.05:
                        print(f"[EFIT]   Skipping {entry}: No data within ±50ms")
                        continue

                    valid_entries.append(entry)

                    self.efit_data[entry] = {
                        'R': efit_data.radius,
                        'psi_N': efit_data.psi_n[time_idx],
                        'rho_pol': efit_data.rho_pol[time_idx],
                        'rho_tor': efit_data.rho_tor[time_idx]
                    }

                    print(f"[EFIT]   {entry} processed")

            except Exception as e:
                print(f"[EFIT] Error: {str(e)}")

        if len(valid_entries) != len(selected_entries):
            missing = set(selected_entries) - set(valid_entries)
            QMessageBox.critical(self.frame, "Incomplete EFIT Data",
                               f"Missing EFIT data for:\n" + '\n'.join(['- ' + e for e in missing]))
            self.efit_data.clear()
            self.computed_efit_tree = None
            return

        self.computed_efit_tree = efit_tree
        print("      EFIT Mapping Completed!")
        print("=" * 50 + "\n")

    @abstractmethod
    def _parse_entry_for_efit(self, entry):
        """Parse entry to extract shot and time for EFIT - implemented by each tab"""
        pass

    @abstractmethod
    def plot_efit_profiles(self):
        """Plot profiles with EFIT mapping - implemented by each tab"""
        pass

    def save_data(self, default_ext='csv'):
        """Save selected data to file"""
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            QMessageBox.warning(self.frame, "Warning", "No data selected to save")
            return

        # Get save file path
        file_path, _ = QFileDialog.getSaveFileName(
            self.frame,
            "Save Data",
            os.path.expanduser("~"),
            "CSV files (*.csv);;Text files (*.txt);;All files (*.*)"
        )

        if not file_path:
            return

        try:
            self._write_data_to_file(file_path, selected_entries)
            print(f"[{self.TAB_NAME}] Data saved to {file_path}")
            QMessageBox.information(self.frame, "Success", f"Data saved to {file_path}")
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to save data: {str(e)}")

    def preview_data(self):
        """Preview selected data in a spreadsheet dialog"""
        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            QMessageBox.warning(self.frame, "Warning", "No data selected to preview")
            return

        try:
            # Write to temp file using existing _write_data_to_file
            fd, tmp_path = tempfile.mkstemp(suffix='.txt')
            os.close(fd)
            self._write_data_to_file(tmp_path, selected_entries)

            with open(tmp_path, 'r') as f:
                lines = f.readlines()
            os.remove(tmp_path)

            self._show_data_preview_dialog(lines)
        except Exception as e:
            QMessageBox.critical(self.frame, "Error", f"Failed to preview data: {str(e)}")

    def _show_data_preview_dialog(self, lines):
        """Show data in a spreadsheet-style dialog"""
        # Parse header and data rows
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

        # Create dialog
        dialog = QDialog(self.frame)
        dialog.setWindowTitle(f"Data Preview — {title_line}" if title_line else "Data Preview")
        dialog.resize(min(120 * n_cols, 1200), min(30 * n_rows + 100, 700))

        layout = QVBoxLayout(dialog)

        # Info label
        info = QLabel(f"{n_rows} rows \u00d7 {n_cols} columns")
        info.setStyleSheet("color: gray;")
        layout.addWidget(info)

        # Table
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

        # Auto-resize columns to content (capped)
        table.resizeColumnsToContents()
        for c in range(n_cols):
            if table.columnWidth(c) > 150:
                table.setColumnWidth(c, 150)

        layout.addWidget(table)

        # Bottom buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        save_btn = QPushButton("Save as .csv")
        save_btn.clicked.connect(lambda: (dialog.close(), self.save_data()))
        btn_layout.addWidget(save_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)

        layout.addLayout(btn_layout)
        dialog.exec()

    @abstractmethod
    def _write_data_to_file(self, file_path, selected_entries):
        """Write data to file - implemented by each tab"""
        pass

    def _get_efit_values_at_R(self, entry, R_positions):
        """Get EFIT coordinate values at given R positions

        Returns tuple of (psi_N, rho_pol, rho_tor) arrays.
        Returns (None, None, None) arrays if EFIT not computed.
        """
        from scipy.interpolate import interp1d

        n_points = len(R_positions)

        # Check if EFIT is computed for this entry
        if not hasattr(self, 'efit_data') or not self.efit_data:
            return [None]*n_points, [None]*n_points, [None]*n_points

        if entry not in self.efit_data:
            return [None]*n_points, [None]*n_points, [None]*n_points

        efit_entry = self.efit_data[entry]
        efit_R = efit_entry['R']

        # Interpolate EFIT values at measurement R positions
        psi_n_interp = interp1d(efit_R, efit_entry['psi_N'],
                                fill_value='extrapolate', bounds_error=False)
        rho_pol_interp = interp1d(efit_R, efit_entry['rho_pol'],
                                  fill_value='extrapolate', bounds_error=False)
        rho_tor_interp = interp1d(efit_R, efit_entry['rho_tor'],
                                  fill_value='extrapolate', bounds_error=False)

        psi_n_vals = psi_n_interp(R_positions)
        rho_pol_vals = rho_pol_interp(R_positions)
        rho_tor_vals = rho_tor_interp(R_positions)

        return psi_n_vals, rho_pol_vals, rho_tor_vals

    def _format_value(self, val, fmt="%10.3f"):
        """Format value, handling None"""
        if val is None:
            return "%10s" % "None"
        return fmt % val
