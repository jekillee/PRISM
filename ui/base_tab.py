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
        # Restore show nodes toggle
        if hasattr(self, 'show_channel_checkbox'):
            saved_show = tab_settings.get("show_nodes", False)
            self.show_channel_checkbox.setChecked(saved_show, animate=False)
        # Restore font sizes
        if hasattr(self, 'label_fontsize'):
            self.label_fontsize = tab_settings.get("label_fontsize", 12)
        if hasattr(self, 'legend_fontsize'):
            self.legend_fontsize = tab_settings.get("legend_fontsize", 8)
        if hasattr(self, 'tick_fontsize'):
            self.tick_fontsize = tab_settings.get("tick_fontsize", 10)
        # Restore fit functions (per param type)
        if hasattr(self, 'fit_func_combos') and self.fit_func_combos:
            saved_funcs = tab_settings.get("fit_funcs", {})
            for ptype, combo in self.fit_func_combos.items():
                func_name = saved_funcs.get(ptype, "mtanh")
                combo.setCurrentText(func_name)
                self.fit_func_names[ptype] = func_name
        if hasattr(self, 'fit_nonparam_options'):
            self.fit_nonparam_options = tab_settings.get("fit_nonparam_options", {})
        if hasattr(self, 'fit_xmin_entry'):
            self.fit_xmin_entry.setText(tab_settings.get("fit_xmin", "0.0"))
            self.fit_xmax_entry.setText(tab_settings.get("fit_xmax", "1.05"))
        # Restore R-shift entries
        if hasattr(self, 'r_shift_entries'):
            saved_rshifts = tab_settings.get("r_shifts", {})
            for diag_name, entry in self.r_shift_entries.items():
                entry.setText(str(saved_rshifts.get(diag_name, "0")))

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
        if hasattr(self, 'fit_func_combos') and self.fit_func_combos:
            tab_settings["fit_funcs"] = {
                ptype: combo.currentText()
                for ptype, combo in self.fit_func_combos.items()
            }
        if hasattr(self, 'fit_nonparam_options') and self.fit_nonparam_options:
            tab_settings["fit_nonparam_options"] = self.fit_nonparam_options
        if hasattr(self, 'fit_xmin_entry'):
            tab_settings["fit_xmin"] = self.fit_xmin_entry.text()
            tab_settings["fit_xmax"] = self.fit_xmax_entry.text()
        # Save R-shift entries
        if hasattr(self, 'r_shift_entries') and self.r_shift_entries:
            tab_settings["r_shifts"] = {
                diag_name: entry.text()
                for diag_name, entry in self.r_shift_entries.items()
            }
        if hasattr(self, 'tci_validate_toggle'):
            tab_settings["tci_validation"] = self.tci_validate_toggle.isChecked()
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

    def _create_efit_controls(self, parent, section_num=3):
        """Create EFIT mapping controls (dropdown + Mapping button only)"""
        frame = QGroupBox(f"{section_num}. EFIT Mapping (Optional)", parent)
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

        # Set equal column stretch
        for i in range(4):
            grid.setColumnStretch(i, 1)

        # Add the frame to the parent's layout
        parent_layout = parent.layout()
        if parent_layout is not None:
            parent_layout.addWidget(frame)

        return frame

    def _create_x_axis_radios(self, parent_layout):
        """Create R/ψₙ/ρₚₒₗ/ρₜₒᵣ radio buttons for x-axis selection.
        Flux radios are disabled until EFIT mapping is computed."""
        self.x_axis_button_group = QButtonGroup()

        radio_R = QRadioButton("R")
        radio_psi_n = QRadioButton("\u03C8\u2099")       # ψₙ
        radio_rho_pol = QRadioButton("\u03C1\u209A\u2092\u2097")  # ρₚₒₗ
        radio_rho_tor = QRadioButton("\u03C1\u209C\u2092\u1D63")  # ρₜₒᵣ

        radio_R.setProperty("axis_value", "R")
        radio_psi_n.setProperty("axis_value", "psi_N")
        radio_rho_pol.setProperty("axis_value", "rho_pol")
        radio_rho_tor.setProperty("axis_value", "rho_tor")

        self.x_axis_button_group.addButton(radio_R)
        self.x_axis_button_group.addButton(radio_psi_n)
        self.x_axis_button_group.addButton(radio_rho_pol)
        self.x_axis_button_group.addButton(radio_rho_tor)

        radio_R.setChecked(True)

        # Disable flux radios until EFIT mapped
        self._flux_radios = [radio_psi_n, radio_rho_pol, radio_rho_tor]
        for r in self._flux_radios:
            r.setEnabled(False)

        row = QHBoxLayout()
        row.addWidget(QLabel("X-axis"))
        row.addStretch(2)
        for r in [radio_R, radio_psi_n, radio_rho_pol, radio_rho_tor]:
            row.addWidget(r, 1)
        parent_layout.addLayout(row)

    def _enable_flux_radios(self, enabled: bool = True):
        """Enable or disable flux coordinate radio buttons"""
        if hasattr(self, '_flux_radios'):
            for r in self._flux_radios:
                r.setEnabled(enabled)

    def _get_selected_x_axis(self):
        """Get the currently selected x-axis value from the radio buttons"""
        if not hasattr(self, 'x_axis_button_group'):
            return "R"
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
        added = False
        for item in self.available_listbox.selectedItems():
            if item.text() not in existing_items:
                self.selected_listbox.addItem(item.text())
                added = True
        if added:
            self._invalidate_efit_and_fit()

    def remove_selected_items(self):
        """Remove selected items from the selected list"""
        items_to_remove = self.selected_listbox.selectedItems()
        if not items_to_remove:
            return
        for item in reversed(items_to_remove):
            self.selected_listbox.takeItem(self.selected_listbox.row(item))
        self._invalidate_efit_and_fit()

    def _invalidate_efit_and_fit(self):
        """Reset EFIT mapping, fit results, and x-axis when selected data changes"""
        self.efit_data.clear()
        self.computed_efit_tree = None
        self._enable_flux_radios(False)
        self._enable_fitting_group(False)
        # Reset x-axis radio to R
        if hasattr(self, 'x_axis_button_group'):
            for btn in self.x_axis_button_group.buttons():
                if btn.property("axis_value") == "R":
                    btn.setChecked(True)
                    break
        # Clear fit results if present
        if hasattr(self, 'fit_results'):
            self.fit_results.clear()

    def _enable_fitting_group(self, enabled: bool = True):
        """Enable or disable the fitting controls group"""
        if hasattr(self, '_fitting_group'):
            self._fitting_group.setEnabled(enabled)

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
        self._enable_flux_radios(True)
        self._enable_fitting_group(True)
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

        # Auto-append .csv if no extension
        if not os.path.splitext(file_path)[1]:
            file_path += '.csv'

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

        # Check if fit results are available
        has_fit = hasattr(self, 'fit_results') and bool(self.fit_results)

        if has_fit:
            self._show_preview_with_modes(selected_entries)
        else:
            try:
                fd, tmp_path = tempfile.mkstemp(suffix='.txt')
                os.close(fd)
                self._write_data_to_file(tmp_path, selected_entries)

                with open(tmp_path, 'r') as f:
                    lines = f.readlines()
                os.remove(tmp_path)

                self._show_data_preview_dialog(lines)
            except Exception as e:
                QMessageBox.critical(self.frame, "Error", f"Failed to preview data: {str(e)}")

    def _show_preview_with_modes(self, selected_entries):
        """Show preview dialog with mode selector (Raw Data / Fitted Profile / Fit Parameters)"""
        dialog = QDialog(self.frame)
        dialog.setWindowTitle("Data Preview")
        dialog.resize(800, 600)
        dlg_layout = QVBoxLayout(dialog)

        # Mode selector
        mode_row = QHBoxLayout()
        mode_row.addWidget(QLabel("View:"))
        mode_combo = QComboBox()
        mode_combo.addItems(["Raw Data", "Fitted Profile", "Fit Parameters"])
        mode_row.addWidget(mode_combo)
        mode_row.addStretch()
        dlg_layout.addLayout(mode_row)

        # Content area (stacked)
        content_stack = QWidget()
        content_layout = QVBoxLayout(content_stack)
        content_layout.setContentsMargins(0, 0, 0, 0)
        dlg_layout.addWidget(content_stack, stretch=1)

        # Table widget that gets replaced on mode change
        table_holder = [None]  # mutable reference

        def update_view():
            mode = mode_combo.currentText()
            # Clear old content
            for i in reversed(range(content_layout.count())):
                w = content_layout.itemAt(i).widget()
                if w:
                    w.deleteLater()

            if mode == "Raw Data":
                lines = self._get_raw_data_lines(selected_entries)
                widget = self._create_preview_table(lines)
            elif mode == "Fitted Profile":
                lines = self._get_fit_profile_lines(selected_entries)
                widget = self._create_preview_table(lines)
            else:
                lines = self._get_fit_params_lines(selected_entries)
                widget = self._create_text_preview(lines)

            content_layout.addWidget(widget)
            table_holder[0] = widget

        mode_combo.currentTextChanged.connect(lambda: update_view())
        update_view()

        # Bottom buttons
        btn_layout = QHBoxLayout()
        btn_layout.addStretch()

        save_btn = QPushButton("Save as .csv")

        def _update_save_label(mode_text):
            if mode_text == "Fit Parameters":
                save_btn.setText("Save as .txt")
            else:
                save_btn.setText("Save as .csv")

        mode_combo.currentTextChanged.connect(_update_save_label)

        def _save_current_mode():
            mode = mode_combo.currentText()
            if mode == "Fit Parameters":
                ext_filter = "Text files (*.txt);;All files (*.*)"
                default_ext = ".txt"
            else:
                ext_filter = "CSV files (*.csv);;Text files (*.txt);;All files (*.*)"
                default_ext = ".csv"

            file_path, _ = QFileDialog.getSaveFileName(
                dialog, "Save Data", os.path.expanduser("~"), ext_filter)
            if not file_path:
                return

            # Auto-append extension if missing
            if not os.path.splitext(file_path)[1]:
                file_path += default_ext

            try:
                if mode == "Raw Data":
                    self._write_data_to_file(file_path, selected_entries)
                elif mode == "Fitted Profile":
                    lines = self._get_fit_profile_lines(selected_entries)
                    with open(file_path, 'w') as f:
                        f.writelines(lines)
                else:  # Fit Parameters
                    lines = self._get_fit_params_lines(selected_entries)
                    with open(file_path, 'w') as f:
                        f.writelines(lines)

                print(f"[{self.TAB_NAME}] Data saved to {file_path}")
                QMessageBox.information(dialog, "Success", f"Data saved to {file_path}")
            except Exception as e:
                QMessageBox.critical(dialog, "Error", f"Failed to save data: {str(e)}")

        save_btn.clicked.connect(_save_current_mode)
        btn_layout.addWidget(save_btn)

        copy_btn = QPushButton("Copy to Clipboard")
        def _copy_current_mode():
            mode = mode_combo.currentText()
            if mode == "Raw Data":
                lines = self._get_raw_data_lines(selected_entries)
            elif mode == "Fitted Profile":
                lines = self._get_fit_profile_lines(selected_entries)
            else:
                lines = self._get_fit_params_lines(selected_entries)
            QApplication.clipboard().setText("".join(lines))
            copy_btn.setText("Copied!")
            from PySide6.QtCore import QTimer
            QTimer.singleShot(1500, lambda: copy_btn.setText("Copy to Clipboard"))
        copy_btn.clicked.connect(_copy_current_mode)
        btn_layout.addWidget(copy_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)

        dlg_layout.addLayout(btn_layout)
        dialog.exec()

    def _get_raw_data_lines(self, selected_entries):
        """Get raw data as lines for preview"""
        try:
            fd, tmp_path = tempfile.mkstemp(suffix='.txt')
            os.close(fd)
            self._write_data_to_file(tmp_path, selected_entries)
            with open(tmp_path, 'r') as f:
                lines = f.readlines()
            os.remove(tmp_path)
            return lines
        except Exception:
            return ["# No data available"]

    def _get_fit_profile_lines(self, selected_entries):
        """Get fitted profile data as lines for preview"""
        if not hasattr(self, 'fit_results') or not self.fit_results:
            return ["# No fit results available"]

        lines = ["# Fitted Profile Data\n"]
        param_types = self._get_fit_param_types() if hasattr(self, '_get_fit_param_types') else []

        # Build header
        header_parts = ["x"]
        for ptype in param_types:
            header_parts.append(f"{ptype}_fit")
        lines.append("#" + ",".join(f"{h:>12s}" for h in header_parts) + "\n")

        for entry in selected_entries:
            if entry not in self.fit_results:
                continue
            lines.append(f"# {entry}\n")
            entry_results = self.fit_results[entry]

            # Find the first successful result to get x_fit
            first_result = None
            for ptype in param_types:
                if ptype in entry_results and entry_results[ptype].success:
                    first_result = entry_results[ptype]
                    break

            if first_result is None:
                lines.append("# Fit failed\n")
                continue

            x_fit = first_result.x_fit
            for j in range(len(x_fit)):
                parts = [f"{x_fit[j]:12.6f}"]
                for ptype in param_types:
                    if ptype in entry_results and entry_results[ptype].success:
                        parts.append(f"{entry_results[ptype].y_fit[j]:12.6f}")
                    else:
                        parts.append(f"{'NaN':>12s}")
                lines.append(",".join(parts) + "\n")

        return lines

    def _get_fit_params_lines(self, selected_entries):
        """Get fit parameters as lines for preview"""
        if not hasattr(self, 'fit_results') or not self.fit_results:
            return ["# No fit results available"]

        lines = ["# Fit Parameters\n"]
        param_types = self._get_fit_param_types() if hasattr(self, '_get_fit_param_types') else []

        for entry in selected_entries:
            if entry not in self.fit_results:
                continue
            lines.append(f"# {entry}\n")
            entry_results = self.fit_results[entry]

            for ptype in param_types:
                if ptype not in entry_results:
                    continue
                result = entry_results[ptype]
                lines.append(f"# {ptype} ({result.func_name}): chi2={result.chi_squared:.4f}, {'OK' if result.success else 'FAILED'}\n")
                if result.success:
                    lines.append(f"#{'Param':>10s},{'Value':>12s},{'Error':>12s}\n")
                    for pname, val in result.params.items():
                        err = result.param_errors.get(pname, 0.0)
                        lines.append(f" {pname:>10s},{val:12.6f},{err:12.6f}\n")

        return lines

    def _create_text_preview(self, lines):
        """Create a read-only text widget for structured text preview (e.g. Fit Parameters)"""
        from PySide6.QtWidgets import QTextEdit
        text_widget = QTextEdit()
        text_widget.setReadOnly(True)
        text_widget.setFont(QFont('Courier', 10))
        text_widget.setPlainText("".join(lines))
        return text_widget

    def _create_preview_table(self, lines):
        """Create a QTableWidget from text lines"""
        headers = []
        data_rows = []
        title_line = ""

        last_comment_with_comma = None
        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                content = stripped.lstrip('#').strip()
                if ',' in content:
                    last_comment_with_comma = content
                else:
                    title_line = content
            else:
                # Use the last comma-containing comment as header
                if not headers and last_comment_with_comma:
                    headers = [h.strip() for h in last_comment_with_comma.split(',')]
                cells = [c.strip() for c in stripped.split(',')]
                data_rows.append(cells)

        if not data_rows:
            table = QTableWidget(1, 1)
            table.setItem(0, 0, QTableWidgetItem("No data available"))
            return table

        n_cols = len(headers) if headers else len(data_rows[0])
        n_rows = len(data_rows)

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
                if c < n_cols:
                    table.setItem(r, c, QTableWidgetItem(val))

        table.resizeColumnsToContents()
        for c in range(n_cols):
            if table.columnWidth(c) > 150:
                table.setColumnWidth(c, 150)

        # Enable Ctrl+C to copy selected cells
        from PySide6.QtGui import QShortcut, QKeySequence
        def _copy_selection():
            sel = table.selectedRanges()
            if not sel:
                return
            rows = range(sel[0].topRow(), sel[0].bottomRow() + 1)
            cols = range(sel[0].leftColumn(), sel[0].rightColumn() + 1)
            text = '\n'.join(
                '\t'.join((table.item(r, c).text() if table.item(r, c) else '')
                          for c in cols)
                for r in rows)
            QApplication.clipboard().setText(text)
        shortcut = QShortcut(QKeySequence.Copy, table)
        shortcut.activated.connect(_copy_selection)

        return table

    def _show_data_preview_dialog(self, lines):
        """Show data in a spreadsheet-style dialog"""
        # Parse header and data rows
        headers = []
        data_rows = []
        title_line = ""
        last_comment_with_comma = None

        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith('#'):
                content = stripped.lstrip('#').strip()
                if ',' in content:
                    last_comment_with_comma = content
                else:
                    title_line = content
            else:
                if not headers and last_comment_with_comma:
                    headers = [h.strip() for h in last_comment_with_comma.split(',')]
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

        copy_btn = QPushButton("Copy to Clipboard")
        def _copy_all():
            text = '\t'.join(headers) + '\n' if headers else ''
            text += '\n'.join('\t'.join(row) for row in data_rows)
            QApplication.clipboard().setText(text)
            copy_btn.setText("Copied!")
            from PySide6.QtCore import QTimer
            QTimer.singleShot(1500, lambda: copy_btn.setText("Copy to Clipboard"))
        copy_btn.clicked.connect(_copy_all)
        btn_layout.addWidget(copy_btn)

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
