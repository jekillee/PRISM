"""
Base class for profile tabs (Ion, Electron, MSE)
Extracts common functionality from concrete profile tabs
"""

from abc import abstractmethod
from typing import Dict, List, Tuple, Any, Optional
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QLayout,
    QGroupBox, QCheckBox, QPushButton, QMessageBox, QLabel,
    QComboBox, QSpinBox, QDialog, QScrollArea, QDialogButtonBox, QFrame,
    QTableWidget, QTableWidgetItem, QHeaderView, QTabWidget,
)
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QColor

from ui.tabs.base_tab import BaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style


class ProfileBaseTab(BaseTab):
    """Base class for all profile visualization tabs

    Provides common functionality for profile-type tabs:
    - Dual-axis plot setup (param1 on left, param2 on right)
    - R profile plotting with EFIT mapping
    - Secondary diagnostic overlay (XICS for Ion, ECE for Electron)
    """

    # Override in subclass for console logging prefix
    TAB_NAME = "Profile"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.secondary_data_cache: Dict[str, Any] = {}
        self._disabled_channels: set = set()  # Channel indices to plot in gray
        self.fit_results: Dict[str, Dict[str, Any]] = {}  # {entry: {param_name: FitResult}}
        self.fit_x_axis: Optional[str] = None  # x-axis used during last fit
        self.fit_dt_s: float = 0.0  # dt (seconds) used during last fit
        self.fit_func_names: Dict[str, str] = {}  # {param_type: func_name}, populated in _create_fitting_controls
        self.fit_func_combos: Dict[str, Any] = {}  # {param_type: QComboBox}
        self.fit_user_params: Dict[str, Dict[str, Any]] = {}  # {param_name: {a1: (val, min, max, fixed), ...}}
        self.fit_nonparam_options: Dict[str, Dict[str, Any]] = {}  # {param_type: {option: value}}
        self.fit_sol_enable: Dict[str, bool] = {}  # {param_type: SOL post-process on/off}
        self._click_points: List[Tuple] = []  # [(x, y, channel_key, ax), ...]
        self.r_shift_entries: Dict[str, Any] = {}  # {diag_name: QLineEdit}
        # Per-timepoint R-shift overrides {entry_str: mm}. Session-only (never
        # saved to settings); wins over the uniform per-diagnostic value.
        self._entry_rshifts: Dict[str, float] = {}

    def create_widgets(self) -> None:
        """Create common profile tab layout"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)

        self.ax1, self.ax2 = self.plot_manager.setup_profile_plot(
            self.figure, self.param1['label'], self.param2['label'])
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.mpl_connect('button_press_event', self._on_canvas_dblclick)
        self.canvas.draw()

        # Left side: canvas + toolbar in a vertical layout
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
        control_layout.setSizeConstraint(QLayout.SetMinimumSize)

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

        # Create control panels
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_efit_controls(control_frame, section_num=3)
        self._create_plot_controls(control_frame)
        self._create_fitting_controls(control_frame)
        self._create_save_controls(control_frame, section_num=6)
        control_layout.addStretch()

    def _create_plot_controls(self, parent: QWidget) -> None:
        """Create the '4. Plot' group.

        Layout: a column-aligned grid (X-axis radios, optional Y-axis dropdown —
        Ion: Rotation, MSE: q or j —, Time avg / R-shift, and a Channels row with
        a Show-Nodes on/off toggle + Select button), then a separator, then the
        Plot | Option action row at the bottom.
        """
        from PySide6.QtWidgets import QGridLayout
        group = QGroupBox("4. Plot")
        group_layout = QVBoxLayout(group)

        # 'Show Nodes' state holder (toggle lives inside the Select Channels
        # dialog, not in this group).
        self._ensure_show_nodes_state()

        # Hidden combo for color mode (used by _get_plot_colors, saved to settings)
        self.color_mode_combo = QComboBox()
        self.color_mode_combo.addItems([
            "Gradient(viridis)", "Gradient(hot)", "Gradient(jet)", "Gradient(coolwarm)",
            "Fixed(tab10)", "Fixed(tab20)", "Fixed(Set1)", "Fixed(Set2)", "Fixed(Set3)",
        ])
        self.color_mode_combo.setCurrentText("Gradient(viridis)")
        self.color_mode_combo.hide()

        # Default font sizes
        self.label_fontsize = 12
        self.legend_fontsize = 8
        self.tick_fontsize = 10

        # One column-aligned grid so every label (col 0) and field (col 1-3) lines
        # up. Order: X-axis, Y-axis, Plot|Option (no separator above), separator,
        # R-shift, Channels, Time avg.
        pgrid = QGridLayout()
        pgrid.setColumnStretch(0, 0)
        pgrid.setColumnStretch(1, 1)
        pgrid.setColumnStretch(2, 0)
        pgrid.setColumnStretch(3, 1)
        content_rows = []
        prow = 0
        prow = self._add_axes_row(pgrid, prow)
        content_rows += list(range(prow))

        # Plot | Option (directly under Y-axis, no separator above)
        btn_row = QHBoxLayout()
        btn_row.setContentsMargins(0, 0, 0, 0)
        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self._on_plot_clicked)
        btn_row.addWidget(plot_button, 3)
        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_style_dialog)
        btn_row.addWidget(style_btn, 1)
        pgrid.addLayout(btn_row, prow, 0, 1, 4)
        content_rows.append(prow)
        prow += 1

        # Separator (thin row — not pinned to the uniform height)
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        pgrid.addWidget(separator, prow, 0, 1, 4)
        prow += 1

        # R-shift, Channels, Time avg (below the separator)
        start = prow
        prow = self._add_rshift_channels_row(pgrid, prow)
        if self._place_time_avg_rshift_in_plot():
            prow = self._build_time_avg_row(pgrid, prow)
        content_rows += list(range(start, prow))

        # Uniform, static row height = a push button's height, so radio/label
        # rows match the (taller) button rows instead of each sizing to its own
        # content. Computed once; applied to every content row (not the separator).
        row_h = QPushButton().sizeHint().height()
        for _i in content_rows:
            pgrid.setRowMinimumHeight(_i, row_h)
        group_layout.addLayout(pgrid)

        parent.layout().addWidget(group)

    def _add_axes_row(self, grid, grid_row) -> int:
        """One row split into equal quarters: 'X-axis' + X-axis combo, and (for
        tabs with a right-axis selector) 'Y-axis (right)' + Y-axis combo. Returns
        the next free grid row."""
        from PySide6.QtWidgets import QComboBox
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(QLabel("X-axis"), 1)
        row.addWidget(self._make_x_axis_combo(), 1)

        spec = self._yaxis_param_spec()
        if spec is not None:
            label, options, default, auto_replot = spec
            row.addWidget(QLabel(label), 1)
            self.y_axis_combo = QComboBox()
            self.y_axis_combo.addItems(options)
            idx = self.y_axis_combo.findText(default)
            if idx >= 0:
                self.y_axis_combo.setCurrentIndex(idx)
            if auto_replot:
                self.y_axis_combo.currentTextChanged.connect(
                    lambda _t: self._on_plot_clicked())
            row.addWidget(self.y_axis_combo, 1)

        grid.addLayout(row, grid_row, 0, 1, 4)
        return grid_row + 1

    def _yaxis_param_spec(self):
        """Hook: (label, options, default, auto_replot) for the right-axis combo,
        or None if the tab has no right-axis selector. Overridden by Ion/MSE."""
        return None

    def _y_param_text(self) -> str:
        """Text of the selected Y-axis (right) combo item ('' if none)."""
        combo = getattr(self, 'y_axis_combo', None)
        return combo.currentText() if combo is not None else ''

    def _set_y_param_text(self, text) -> None:
        """Select the Y-axis (right) combo item whose text matches `text`."""
        combo = getattr(self, 'y_axis_combo', None)
        if combo is None:
            return
        idx = combo.findText(text)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _ensure_show_nodes_state(self) -> None:
        """Create the headless 'Show Nodes' toggle state holder. The toggle lives
        inside the Select Channels dialog; kept as a ToggleSwitch so settings
        save/restore and `_add_channel_labels` keep working unchanged."""
        if getattr(self, 'show_channel_checkbox', None) is None:
            from ui.widgets.toggle_switch import ToggleSwitch
            self.show_channel_checkbox = ToggleSwitch()
            self.show_channel_checkbox.toggled.connect(self._on_show_nodes_toggled)

    def _add_rshift_channels_row(self, grid, grid_row) -> int:
        """Two full-width rows, each with a left label + a right-aligned button:
          'R-shift'  [Adjust]   (only for tabs that support R-shift)
          'Channels' [Select]   (opens the channel enable/dim + Show-Nodes dialog)
        Returns the next free grid row."""
        def _label_button_row(text, button):
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.addWidget(QLabel(text), 1)   # left half
            row.addWidget(button, 1)          # right half
            return row

        if self._get_rshift_diagnostics():
            adjust_btn = QPushButton("Adjust")
            adjust_btn.setToolTip("Adjust R-shift per timepoint [mm]")
            adjust_btn.clicked.connect(self._show_rshift_dialog)
            grid.addLayout(_label_button_row("R-shift", adjust_btn), grid_row, 0, 1, 4)
            grid_row += 1

        channels_btn = QPushButton("Select / Deselect")
        channels_btn.setToolTip("Select channels to enable/dim, and toggle node labels")
        channels_btn.clicked.connect(self._show_channel_selector)
        grid.addLayout(_label_button_row("Channels", channels_btn), grid_row, 0, 1, 4)
        grid_row += 1
        return grid_row

    def _on_plot_clicked(self) -> None:
        """Unified plot handler: dispatches to R-space or flux-space plot"""
        x_axis = self._get_selected_x_axis()
        if x_axis == "R":
            self.plot_data()
        else:
            self.plot_efit_profiles()

    def _on_show_nodes_toggled(self, state) -> None:
        """Toggle node labels immediately when checkbox is changed"""
        if not (self.ax1.lines or self.ax1.collections):
            return
        if state:
            self._on_plot_clicked()
        else:
            from matplotlib.text import Annotation
            for ax in [self.ax1, self.ax2]:
                for child in list(ax.get_children()):
                    if isinstance(child, Annotation):
                        child.remove()
            self.canvas.draw_idle()

    def _get_channel_info(self) -> List[Tuple[str, str]]:
        """Return list of (key, label) for diagnostic channels.

        Override in subclasses for diagnostic-specific channel info.
        key: unique string identifier (e.g., "CES_0", "TS_5", "ECE_10")
        label: display label (e.g., "CES Ch1 (R=1.850m)")
        """
        return []

    def _show_channel_selector(self) -> None:
        """Show dialog with checkboxes for each diagnostic channel"""
        channel_info = self._get_channel_info()
        if not channel_info:
            QMessageBox.information(self.frame, "Channels",
                "No channel data available.\nLoad and plot data first.")
            return

        from ui.widgets.toggle_switch import ToggleSwitch
        dialog = QDialog(self.frame)
        dialog.setWindowTitle("Select / Deselect Channels")
        dialog.setMinimumWidth(300)
        dlg_layout = QVBoxLayout(dialog)

        hint = QLabel("Uncheck channels to dim them (gray) on the plot. "
                      "You can also double-click a data point on the plot to toggle. "
                      "Unchecked channels are excluded from fitting.")
        hint.setWordWrap(True)
        dlg_layout.addWidget(hint)

        # Show Nodes on/off toggle (draws channel node labels on the plot),
        # synced to the headless show_channel_checkbox state holder.
        self._ensure_show_nodes_state()
        sn_row = QHBoxLayout()
        show_nodes_toggle = ToggleSwitch()
        show_nodes_toggle.setChecked(self.show_channel_checkbox.isChecked(), animate=False)
        sn_row.addWidget(show_nodes_toggle)
        sn_row.addWidget(QLabel("Show Nodes"))
        sn_row.addStretch()
        dlg_layout.addLayout(sn_row)

        # Scrollable area for checkboxes
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMaximumHeight(400)
        cb_widget = QWidget()
        cb_layout = QVBoxLayout(cb_widget)

        checkboxes = []
        for idx, label in channel_info:
            cb = QCheckBox(label)
            cb.setChecked(idx not in self._disabled_channels)
            cb_layout.addWidget(cb)
            checkboxes.append((idx, cb))

        cb_layout.addStretch()
        scroll.setWidget(cb_widget)
        dlg_layout.addWidget(scroll)

        # Button row (3 equal parts): Select All | Deselect All | Apply
        btn_row = QHBoxLayout()
        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(
            lambda: [cb.setChecked(True) for _, cb in checkboxes])
        btn_row.addWidget(select_all_btn, 1)
        deselect_all_btn = QPushButton("Deselect All")
        deselect_all_btn.clicked.connect(
            lambda: [cb.setChecked(False) for _, cb in checkboxes])
        btn_row.addWidget(deselect_all_btn, 1)
        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(dialog.accept)
        btn_row.addWidget(apply_btn, 1)
        dlg_layout.addLayout(btn_row)

        def _apply():
            self._disabled_channels.clear()
            for idx, cb in checkboxes:
                if not cb.isChecked():
                    self._disabled_channels.add(idx)
            # Update Show-Nodes state without its own redraw; one replot below.
            self.show_channel_checkbox.blockSignals(True)
            self.show_channel_checkbox.setChecked(show_nodes_toggle.isChecked())
            self.show_channel_checkbox.blockSignals(False)
            if self.ax1.lines or self.ax1.collections:
                self._on_plot_clicked()
        dialog.accepted.connect(_apply)
        self._channel_selector_dialog = dialog
        dialog.show()

    def _get_channel_mask(self, channel_keys: List[str]) -> np.ndarray:
        """Return boolean mask: True = enabled, False = disabled (gray).

        Args:
            channel_keys: List of channel key strings (e.g., ["CES_0", "CES_1", ...])
        """
        return np.array([k not in self._disabled_channels for k in channel_keys])

    @staticmethod
    def _parse_color_mode(text):
        """Extract colormap name from dropdown text like 'Fixed(tab10)' → 'tab10'"""
        start = text.find('(')
        end = text.find(')')
        if start != -1 and end != -1:
            return text[start + 1:end]
        return 'tab10'

    def _get_plot_colors(self, entries: List[str]) -> List:
        """Get plot colors for entries using selected color mode."""
        combo = getattr(self, 'color_mode_combo', None)
        mode = combo.currentText() if combo else "Fixed(tab10)"
        cmap_name = self._parse_color_mode(mode)
        return self.plot_manager.color_manager.get_colors_for_entries(entries, colormap=cmap_name)

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
        btn_box = QDialogButtonBox(QDialogButtonBox.RestoreDefaults | QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)

        def reset_defaults():
            color_combo.setCurrentText("Gradient(viridis)")
            label_spin.setValue(12)
            legend_spin.setValue(8)
            tick_spin.setValue(10)
        btn_box.button(QDialogButtonBox.RestoreDefaults).clicked.connect(reset_defaults)
        dlg_layout.addWidget(btn_box)

        def _apply():
            self.color_mode_combo.setCurrentText(color_combo.currentText())
            self.label_fontsize = label_spin.value()
            self.legend_fontsize = legend_spin.value()
            self.tick_fontsize = tick_spin.value()
            # Auto-apply if plot exists
            if self.ax1.lines or self.ax1.collections:
                self._on_plot_clicked()
        dialog.accepted.connect(_apply)
        self._style_dialog = dialog
        dialog.show()

    def _create_fitting_controls(self, parent: QWidget) -> None:
        """Create fitting controls (section 5)"""
        from PySide6.QtWidgets import QLineEdit, QGridLayout

        group = QGroupBox("5. Fitting (EFIT Mapping Required)")
        group.setEnabled(False)
        self._fitting_group = group
        # 4 columns: label | entry1 | ~ | entry2
        grid = QGridLayout(group)
        grid.setColumnStretch(0, 0)  # label - fixed
        grid.setColumnStretch(1, 1)  # entry1 - stretch
        grid.setColumnStretch(2, 0)  # separator - fixed
        grid.setColumnStretch(3, 1)  # entry2 - stretch

        # Row 0+: Per-param-type function selectors
        param_types = self._get_fit_param_types()
        row = 0
        for ptype in param_types:
            self.fit_func_names.setdefault(ptype, 'mtanh')
            grid.addWidget(QLabel(ptype), row, 0)
            func_row = QHBoxLayout()
            combo = QComboBox()
            combo.addItems(['mtanh', 'ptanh', 'EPED', 'spline', 'RBF', 'GPR'])
            combo.setCurrentText(self.fit_func_names[ptype])
            combo.currentTextChanged.connect(
                lambda text, pt=ptype: self._on_fit_func_changed(pt, text))
            self.fit_func_combos[ptype] = combo
            self._update_fit_func_tooltip_for(ptype)
            func_row.addWidget(combo)
            params_btn = QPushButton("Option")
            params_btn.clicked.connect(
                lambda checked=False, pt=ptype: self._show_fit_options_dialog(pt))
            func_row.addWidget(params_btn)
            grid.addLayout(func_row, row, 1, 1, 3)
            row += 1

        # X Range — aligned [min] ~ [max]
        grid.addWidget(QLabel("X Range"), row, 0)
        self.fit_xmin_entry = QLineEdit("0.0")
        grid.addWidget(self.fit_xmin_entry, row, 1)
        grid.addWidget(QLabel("~"), row, 2, Qt.AlignCenter)
        self.fit_xmax_entry = QLineEdit("1.05")
        grid.addWidget(self.fit_xmax_entry, row, 3)
        row += 1

        # Time avg row — kept here unless the subclass (Ion/Electron) hosts it in
        # the '4. Plot' group. (R-shift is a button in the Plot group, not here.)
        if not self._place_time_avg_rshift_in_plot():
            row = self._build_time_avg_row(grid, row)

        # Subclass hook for extra controls (e.g. TCI validation checkbox)
        row = self._add_extra_fitting_controls(grid, row)

        # Fit button (full width, at bottom)
        fit_button = QPushButton("Fit")
        fit_button.clicked.connect(self.perform_fitting)
        grid.addWidget(fit_button, row, 0, 1, 4)

        parent.layout().addWidget(group)

    def _place_time_avg_rshift_in_plot(self) -> bool:
        """Whether Time avg + R-shift live in the '4. Plot' group instead of the
        '5. Fitting' group. Ion/Electron override to True so the controls are
        always available (Plot is enabled without EFIT mapping); other profile
        tabs (e.g. MSE) keep them in the Fitting group."""
        return False

    def _build_time_avg_row(self, grid, row: int) -> int:
        """Create the 'Time avg. [ms]' label + toggle + entry into `grid` at
        `row`. Returns the next free row."""
        from PySide6.QtWidgets import QLineEdit
        from ui.widgets.toggle_switch import ToggleSwitch

        grid.addWidget(QLabel("Time avg. [ms]"), row, 0)
        self.dt_toggle = ToggleSwitch()
        self.dt_toggle.toggled.connect(self._on_dt_toggled)
        grid.addWidget(self.dt_toggle, row, 1, 1, 2, Qt.AlignCenter)
        self.fit_dt_entry = QLineEdit("0")
        self.fit_dt_entry.setToolTip("Time averaging window: fit on mean of [t-dt, t+dt].")
        self.fit_dt_entry.setEnabled(False)
        grid.addWidget(self.fit_dt_entry, row, 3)
        self._dt_enabled = False
        return row + 1

    def _add_extra_fitting_controls(self, grid, row: int) -> int:
        """Hook for subclasses to add extra controls to the fitting section.
        Returns the next available row number."""
        return row

    def _get_rshift_diagnostics(self) -> List[str]:
        """Return list of diagnostic names that support R-shift.
        Override in subclass. E.g., ['CES'] for Ion, ['Thomson', 'ECE'] for Electron."""
        return []

    def _get_rshift(self, diag_name: str, entry: str = None) -> float:
        """R-shift in meters (dialog input is in mm) for a diagnostic.

        R-shift is set per timepoint via the R-shift dialog and is session-only.
        Returns the override for `entry` if one exists, else 0 (no shift).
        `diag_name` is accepted for call-site clarity / future per-diagnostic use."""
        if entry is not None:
            override = self._entry_rshifts.get(entry)
            if override is not None:
                return override / 1000.0  # mm -> m
        # Legacy per-diagnostic uniform field (no longer shown); 0 if absent.
        widget = self.r_shift_entries.get(diag_name)
        if widget is None:
            return 0.0
        try:
            return float(widget.text()) / 1000.0  # mm -> m
        except ValueError:
            return 0.0

    def _show_rshift_dialog(self) -> None:
        """R-shift editor: one field per selected timepoint (in mm). Blank = no
        shift. Values are session-only and applied (re-plot) on OK."""
        from PySide6.QtWidgets import QLineEdit, QGridLayout
        entries = [self.selected_listbox.item(i).text()
                   for i in range(self.selected_listbox.count())]
        if not entries:
            QMessageBox.information(self.frame, "R-shift",
                "No timepoints selected.\nSelect data first.")
            return

        dialog = QDialog(self.frame)
        dialog.setWindowTitle("R-shift [mm]")
        dialog.setMinimumWidth(360)
        dlg_layout = QVBoxLayout(dialog)
        hint = QLabel("R-shift [mm] per timepoint (blank = 0).")
        hint.setWordWrap(True)
        dlg_layout.addWidget(hint)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setMaximumHeight(400)
        form_widget = QWidget()
        form = QGridLayout(form_widget)
        form.setColumnStretch(0, 1)
        form.setColumnStretch(1, 0)
        fields = []
        for r, entry in enumerate(entries):
            form.addWidget(QLabel(entry), r, 0)
            le = QLineEdit()
            le.setFixedWidth(80)
            le.setPlaceholderText("0")
            override = self._entry_rshifts.get(entry)
            if override is not None:
                le.setText(str(override))
            form.addWidget(le, r, 1)
            fields.append((entry, le))
        scroll.setWidget(form_widget)
        dlg_layout.addWidget(scroll)

        # Button row: Clear all (left half) | Apply (right half)
        btn_row = QHBoxLayout()
        clear_btn = QPushButton("Clear all")
        clear_btn.clicked.connect(lambda: [le.clear() for _, le in fields])
        btn_row.addWidget(clear_btn, 1)
        apply_btn = QPushButton("Apply")
        apply_btn.clicked.connect(dialog.accept)
        btn_row.addWidget(apply_btn, 1)
        dlg_layout.addLayout(btn_row)

        def _apply():
            for entry, le in fields:
                txt = le.text().strip()
                if txt == '':
                    self._entry_rshifts.pop(entry, None)
                    continue
                try:
                    self._entry_rshifts[entry] = float(txt)
                except ValueError:
                    self._entry_rshifts.pop(entry, None)
            self._on_plot_clicked()
        dialog.accepted.connect(_apply)
        self._rshift_dialog = dialog
        dialog.show()

    def _update_fit_func_tooltip_for(self, ptype: str) -> None:
        """Update tooltip on function combo for a specific param type"""
        from core.fitting import FIT_FUNCTIONS
        combo = self.fit_func_combos.get(ptype)
        if combo is None:
            return
        func_name = combo.currentText()
        info = FIT_FUNCTIONS.get(func_name, {})
        desc = info.get('description', '')
        combo.setToolTip(desc)

    def _on_fit_func_changed(self, ptype: str, text: str) -> None:
        """Handle fit function change for a specific param type"""
        self.fit_func_names[ptype] = text
        # Clear user params for this param type when function changes
        self.fit_user_params.pop(ptype, None)
        self._update_fit_func_tooltip_for(ptype)

    def _get_fit_param_types(self) -> List[str]:
        """Return parameter types for fitting (e.g., ['Ti', 'vT'] or ['Te', 'ne']).
        Override in subclass."""
        return []

    def _render_formula_widget(self, dialog, formula_latex: str):
        """Render LaTeX formula as a matplotlib widget for embedding in dialogs."""
        if not formula_latex:
            return None
        try:
            from matplotlib.figure import Figure
            from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
            text_color = dialog.palette().text().color().name()
            formula_lines = formula_latex.split('\n')
            n_lines = len(formula_lines)
            fig_height = max(0.4 * n_lines, 0.6)
            fig = Figure(figsize=(5, fig_height))
            fig.patch.set_alpha(0.0)
            for i, line in enumerate(formula_lines):
                y_pos = 1.0 - (i + 0.5) / n_lines
                fig.text(0.02, y_pos, line, fontsize=11,
                         verticalalignment='center', color=text_color)
            canvas = FigureCanvasQTAgg(fig)
            canvas.setFixedHeight(int(fig_height * 72))
            canvas.setStyleSheet("background: transparent;")
            return canvas
        except Exception:
            return None

    def _show_fit_options_dialog(self, ptype: str) -> None:
        """Show Fit Options dialog with parameter table for a specific param type"""
        from core.fitting import FIT_FUNCTIONS, get_default_params

        func_name = self.fit_func_names.get(ptype, 'mtanh')

        if func_name == 'RBF':
            self._show_rbf_options_dialog(ptype)
            return

        # Non-parametric functions (spline, GPR): info-only dialog
        if func_name in ('spline', 'GPR'):
            func_info = FIT_FUNCTIONS.get(func_name, {})
            dialog = QDialog(self.frame)
            dialog.setWindowTitle(f"{ptype} — {func_name}")
            dialog.setMinimumWidth(400)
            dlg_layout = QVBoxLayout(dialog)

            title_label = QLabel(func_name)
            title_label.setStyleSheet("font-weight: bold; font-size: 13px;")
            dlg_layout.addWidget(title_label)

            desc_label = QLabel(func_info.get('description', ''))
            desc_label.setWordWrap(True)
            dlg_layout.addWidget(desc_label)

            formula_widget = self._render_formula_widget(dialog, func_info.get('formula_latex', ''))
            if formula_widget:
                dlg_layout.addWidget(formula_widget)

            dlg_layout.addSpacing(10)
            info_label = QLabel("No parameters to configure.")
            info_label.setStyleSheet("color: gray;")
            dlg_layout.addWidget(info_label)

            close_btn = QPushButton("Close")
            close_btn.clicked.connect(dialog.close)
            dlg_layout.addWidget(close_btn)
            self._fit_options_dialog = dialog
            dialog.show()
            return

        func_info = FIT_FUNCTIONS[func_name]
        param_names = list(func_info['param_names'])
        param_descs = dict(func_info.get('param_descriptions', {}))
        defaults = dict(get_default_params(func_name, ptype))
        user = self.fit_user_params.get(ptype, {})
        sol_on = self.fit_sol_enable.get(ptype, False)
        # Always include sol_b row; it's greyed when the SOL toggle is off.
        from core.fitting import _DEFAULT_SOL_B
        param_names.append('sol_b')
        defaults['sol_b'] = _DEFAULT_SOL_B
        param_descs['sol_b'] = 'SOL linear-decay slope (ψ_N > 1 region)'

        dialog = QDialog(self.frame)
        dialog.setWindowTitle(f"{ptype} — {func_name}")
        dialog.setMinimumWidth(500)
        dialog.setMinimumHeight(300)
        dlg_layout = QVBoxLayout(dialog)

        # Function info header
        title_label = QLabel(func_name)
        title_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        dlg_layout.addWidget(title_label)

        desc_label = QLabel(func_info.get('description', ''))
        desc_label.setWordWrap(True)
        dlg_layout.addWidget(desc_label)

        # Formula
        formula_widget = self._render_formula_widget(dialog, func_info.get('formula_latex', ''))
        if formula_widget:
            dlg_layout.addWidget(formula_widget)

        # SOL toggle — appears in every Option dialog. Greys out the sol_b row.
        sol_row_layout = QHBoxLayout()
        sol_cb = QCheckBox("Add SOL tail:")
        sol_cb.setChecked(sol_on)
        sol_cb.setToolTip(
            "Replaces the ψ_N > 1 region of the fitted curve with a linear\n"
            "tail  y = y(ψ=1) · max(1 − sol_b·(ψ−1), 0).\n"
            "sol_b is fitted on ψ_N > 1 data (or fixed via the table row).\n"
            "Toggling SOL never changes the ψ_N ≤ 1 portion of the fit.")
        sol_row_layout.addWidget(sol_cb)
        sol_formula = QLabel("y_SOL = y(ψ=1) · max(1 − sol_b·(ψ−1), 0)")
        sol_formula.setStyleSheet("color: #888; font-style: italic;")
        sol_row_layout.addWidget(sol_formula)
        sol_row_layout.addStretch()
        dlg_layout.addLayout(sol_row_layout)

        # Parameter table
        table = QTableWidget(len(param_names), 5)
        table.setHorizontalHeaderLabels(['Param', 'Value', 'Min', 'Max', 'Fix'])
        table.setVerticalHeaderLabels(['' for _ in param_names])
        table.verticalHeader().setVisible(False)
        for col in range(1, 5):
            table.horizontalHeader().setSectionResizeMode(col, QHeaderView.ResizeToContents)
        table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        table.horizontalHeader().setStretchLastSection(False)
        table.setFont(QFont('Courier', 9))

        for row, pname in enumerate(param_names):
            if pname in user:
                val, lb, ub, fixed = user[pname]
            elif pname in defaults:
                val, lb, ub = defaults[pname]
                fixed = False
            else:
                val, lb, ub, fixed = 1.0, -1e10, 1e10, False

            # Param name + description
            desc = param_descs.get(pname, '')
            param_item = QTableWidgetItem(f"{pname}: {desc}" if desc else pname)
            param_item.setFlags(Qt.ItemIsEnabled)  # read-only
            table.setItem(row, 0, param_item)

            table.setItem(row, 1, QTableWidgetItem(f"{val:.4g}"))

            lb_str = "-inf" if lb == -np.inf else f"{lb:.4g}"
            table.setItem(row, 2, QTableWidgetItem(lb_str))

            ub_str = "inf" if ub == np.inf else f"{ub:.4g}"
            table.setItem(row, 3, QTableWidgetItem(ub_str))

            fix_cb = QCheckBox()
            fix_cb.setChecked(fixed)
            fix_container = QWidget()
            fix_layout = QHBoxLayout(fix_container)
            fix_layout.addWidget(fix_cb)
            fix_layout.setAlignment(Qt.AlignCenter)
            fix_layout.setContentsMargins(0, 0, 0, 0)
            table.setCellWidget(row, 4, fix_container)

        dlg_layout.addWidget(table)

        # Enable/disable the sol_b row visually based on SOL toggle.
        # Pick colors that read as disabled vs normal in both themes.
        from ui.theme import ThemeManager
        if ThemeManager.current_theme == 'dark':
            _enabled_brush = QColor('#cccccc')
            _disabled_brush = QColor('#666666')
        else:
            _enabled_brush = QColor('#202020')
            _disabled_brush = QColor('#a0a0a0')
        def _set_sol_row_enabled(enabled: bool):
            try:
                row = param_names.index('sol_b')
            except ValueError:
                return
            for col in range(4):
                it = table.item(row, col)
                if it is None:
                    continue
                flags = it.flags()
                if enabled:
                    if col != 0:
                        it.setFlags(flags | Qt.ItemIsEditable)
                    it.setForeground(_enabled_brush)
                else:
                    it.setFlags(flags & ~Qt.ItemIsEditable)
                    it.setForeground(_disabled_brush)
            cell = table.cellWidget(row, 4)
            if cell is not None:
                cb = cell.findChild(QCheckBox)
                if cb is not None:
                    cb.setEnabled(enabled)
        _set_sol_row_enabled(sol_on)
        sol_cb.toggled.connect(_set_sol_row_enabled)

        # Buttons
        btn_layout = QHBoxLayout()

        restore_btn = QPushButton("Restore Defaults")
        def restore_defaults():
            for row, pname in enumerate(param_names):
                if pname in defaults:
                    val, lb, ub = defaults[pname]
                else:
                    val, lb, ub = 1.0, -1e10, 1e10
                table.item(row, 1).setText(f"{val:.4g}")
                lb_str = "-inf" if lb == -np.inf else f"{lb:.4g}"
                table.item(row, 2).setText(lb_str)
                ub_str = "inf" if ub == np.inf else f"{ub:.4g}"
                table.item(row, 3).setText(ub_str)
                table.cellWidget(row, 4).findChild(QCheckBox).setChecked(False)
        restore_btn.clicked.connect(restore_defaults)
        btn_layout.addWidget(restore_btn)

        btn_layout.addStretch()

        apply_btn = QPushButton("Apply")
        def apply_params():
            params = {}
            for row, pname in enumerate(param_names):
                try:
                    val_text = table.item(row, 1).text()
                    val = float(val_text)
                    lb_text = table.item(row, 2).text()
                    lb = -np.inf if lb_text.strip().lower() in ('-inf', '') else float(lb_text)
                    ub_text = table.item(row, 3).text()
                    ub = np.inf if ub_text.strip().lower() in ('inf', '') else float(ub_text)
                    fixed = table.cellWidget(row, 4).findChild(QCheckBox).isChecked()
                    params[pname] = (val, lb, ub, fixed)
                except ValueError:
                    QMessageBox.warning(dialog, "Invalid Value",
                        f"Invalid value for {ptype}/{pname}")
                    return
            self.fit_user_params[ptype] = params
            self.fit_sol_enable[ptype] = sol_cb.isChecked()
            dialog.accept()
        apply_btn.clicked.connect(apply_params)
        btn_layout.addWidget(apply_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)

        dlg_layout.addLayout(btn_layout)
        self._fit_options_dialog = dialog
        dialog.show()

    def _show_rbf_options_dialog(self, ptype: str) -> None:
        """Show RBF options dialog with n_bases setting"""
        from core.fitting import FIT_FUNCTIONS

        func_info = FIT_FUNCTIONS['RBF']
        opts = self.fit_nonparam_options.get(ptype, {})
        cur_n_bases = opts.get('n_bases', 5)

        dialog = QDialog(self.frame)
        dialog.setWindowTitle(f"{ptype} — RBF")
        dialog.setMinimumWidth(400)
        dlg_layout = QVBoxLayout(dialog)

        # Title (bold)
        title_label = QLabel("RBF")
        title_label.setStyleSheet("font-weight: bold; font-size: 13px;")
        dlg_layout.addWidget(title_label)

        # Description
        desc_label = QLabel(func_info.get('description', ''))
        desc_label.setWordWrap(True)
        dlg_layout.addWidget(desc_label)

        # Formula
        formula_widget = self._render_formula_widget(dialog, func_info.get('formula_latex', ''))
        if formula_widget:
            dlg_layout.addWidget(formula_widget)

        dlg_layout.addSpacing(10)

        # n_bases setting
        from PySide6.QtWidgets import QSpinBox, QGridLayout
        grid = QGridLayout()
        grid.addWidget(QLabel("Number of bases:"), 0, 0)
        n_bases_spin = QSpinBox()
        n_bases_spin.setRange(0, 200)
        n_bases_spin.setSpecialValueText("Auto (= data points)")
        n_bases_spin.setValue(cur_n_bases)
        n_bases_spin.setToolTip("0 = use all data points as centers.\n"
                                "Fewer bases → smoother fit.\n"
                                "More bases → more detailed fit.")
        grid.addWidget(n_bases_spin, 0, 1)
        dlg_layout.addLayout(grid)

        # Buttons
        btn_layout = QHBoxLayout()
        apply_btn = QPushButton("Apply")
        def apply_opts():
            val = n_bases_spin.value()
            if ptype not in self.fit_nonparam_options:
                self.fit_nonparam_options[ptype] = {}
            self.fit_nonparam_options[ptype]['n_bases'] = val
            dialog.accept()
        apply_btn.clicked.connect(apply_opts)
        btn_layout.addWidget(apply_btn)

        close_btn = QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        btn_layout.addWidget(close_btn)

        dlg_layout.addLayout(btn_layout)
        self._rbf_options_dialog = dialog
        dialog.show()

    def perform_fitting(self) -> None:
        """Perform profile fitting on current x-axis coordinate"""
        x_axis = self._get_selected_x_axis()

        if x_axis == "R":
            QMessageBox.warning(self.frame, "Warning",
                "Fitting is available only in flux coordinates (ψₙ, ρₚₒₗ, ρₜₒᵣ).\n"
                "Please compute EFIT mapping and select a flux coordinate.")
            return

        if not self.efit_data or self.computed_efit_tree is None:
            QMessageBox.warning(self.frame, "Warning",
                "Please compute EFIT mapping first before fitting.")
            return

        from core.fitting import fit_profile

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        param_types = self._get_fit_param_types()
        self.fit_results.clear()
        self.fit_x_axis = x_axis

        try:
            x_fit_min = float(self.fit_xmin_entry.text())
            x_fit_max = float(self.fit_xmax_entry.text())
        except ValueError:
            x_fit_min, x_fit_max = 0.0, 1.05
        x_fit_range = (x_fit_min, x_fit_max)

        success_count = 0
        fail_count = 0

        # dt validation
        dt_s = self._get_dt_seconds()
        dt_ms = dt_s * 1000.0
        self.fit_dt_s = dt_s

        if self._dt_enabled:
            dt_text = self.fit_dt_entry.text().strip()
            try:
                val = float(dt_text)
            except ValueError:
                QMessageBox.warning(self.frame, "Invalid dt",
                    f"'{dt_text}' is not a valid number. Please enter dt in ms (e.g., 10, 50).")
                return
            if val < 0:
                QMessageBox.warning(self.frame, "Invalid dt",
                    f"dt cannot be negative ({val} ms). Please enter a positive value.")
                return

            # Check if dt is too small to actually average anything
            if dt_s > 0:
                first_entry = selected_entries[0]
                sampling_dt = self._get_data_sampling_dt(first_entry)
                if sampling_dt is not None and dt_s < sampling_dt:
                    QMessageBox.information(self.frame, "dt too small",
                        f"dt={dt_ms:.1f}ms is smaller than the data sampling interval "
                        f"({sampling_dt*1000:.1f}ms).\n\n"
                        f"No averaging will be performed. "
                        f"Use dt \u2265 {sampling_dt*1000:.1f}ms to average at least 2 time slices.")

        func_summary = ", ".join(f"{pt}={self.fit_func_names.get(pt, 'mtanh')}" for pt in param_types)
        dt_info = f", dt={dt_ms:.0f}ms" if dt_s > 0 else ""
        print(f"\n[Fitting] {func_summary} on {len(selected_entries)} entries (x={x_axis}{dt_info})")

        for entry in selected_entries:
            if entry not in self.efit_data:
                continue

            efit_entry = self.efit_data[entry]
            efit_R = efit_entry['R']

            if x_axis == "psi_N":
                interp_func = interp1d(efit_R, efit_entry['psi_N'], fill_value='extrapolate')
            elif x_axis == "rho_pol":
                interp_func = interp1d(efit_R, efit_entry['rho_pol'], fill_value='extrapolate')
            else:
                interp_func = interp1d(efit_R, efit_entry['rho_tor'], fill_value='extrapolate')

            # Get profile data from subclass (with optional time averaging)
            profile_data = self._get_efit_profile_data(entry, interp_func, dt_s=dt_s)
            if profile_data is None:
                continue

            entry_results = {}
            for ptype in param_types:
                if ptype not in profile_data:
                    continue

                func_name = self.fit_func_names.get(ptype, 'mtanh')
                x_data = profile_data[ptype]['x']
                y_data = profile_data[ptype]['y']
                sigma = profile_data[ptype].get('err', None)

                user_params = self.fit_user_params.get(ptype, None)
                nonparam_opts = self.fit_nonparam_options.get(ptype, {})
                sol_enable = self.fit_sol_enable.get(ptype, False)

                result = fit_profile(
                    x_data, y_data,
                    func_name=func_name,
                    param_type=ptype,
                    user_params=user_params,
                    sigma=sigma,
                    x_fit_range=x_fit_range,
                    n_bases=nonparam_opts.get('n_bases', 5),
                    sol_enable=sol_enable,
                )

                entry_results[ptype] = result

                if result.success:
                    success_count += 1
                    print(f"[Fitting]   {entry} {ptype} ({func_name}): chi2={result.chi_squared:.4f}")
                    # Parameter table for parametric functions
                    if result.params and func_name in ('mtanh', 'ptanh', 'EPED'):
                        from core.fitting import FIT_FUNCTIONS, get_default_params
                        fi = FIT_FUNCTIONS.get(func_name, {})
                        pnames = list(fi.get('param_names', []))
                        defaults = dict(get_default_params(func_name, ptype))
                        up = self.fit_user_params.get(ptype, {})
                        print(f"[Fitting]     {'Param':>5}  {'Value':>10}  {'Error':>10}  {'Min':>10}  {'Max':>10}  Fixed")
                        print(f"[Fitting]     {'─'*5}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  ─────")
                        for pn in pnames:
                            val = result.params.get(pn, 0.0)
                            err = result.param_errors.get(pn, 0.0)
                            if pn in up:
                                _, lb, ub, fixed = up[pn]
                            elif pn in defaults:
                                _, lb, ub = defaults[pn]
                                fixed = False
                            else:
                                lb, ub, fixed = -np.inf, np.inf, False
                            lb_s = '-inf' if lb == -np.inf else f'{lb:.4g}'
                            ub_s = 'inf' if ub == np.inf else f'{ub:.4g}'
                            fx_s = 'Y' if fixed else ''
                            print(f"[Fitting]     {pn:>5}  {val:10.4f}  {err:10.4f}  {lb_s:>10}  {ub_s:>10}  {fx_s:>5}")
                    # SOL slope (any function), printed once after the core table
                    if sol_enable and 'sol_b' in (result.params or {}):
                        b = result.params['sol_b']
                        b_err = (result.param_errors or {}).get('sol_b', 0.0)
                        up = self.fit_user_params.get(ptype, {})
                        if 'sol_b' in up:
                            _, lb, ub, fx = up['sol_b']
                        else:
                            lb, ub, fx = 0.01, 50.0, False
                        lb_s = '-inf' if lb == -np.inf else f'{lb:.4g}'
                        ub_s = 'inf' if ub == np.inf else f'{ub:.4g}'
                        fx_s = 'Y' if fx else ''
                        # Reuse the parametric table format for sol_b
                        if not (result.params and func_name in ('mtanh', 'ptanh', 'EPED')):
                            print(f"[Fitting]     {'Param':>5}  {'Value':>10}  {'Error':>10}  {'Min':>10}  {'Max':>10}  Fixed")
                            print(f"[Fitting]     {'─'*5}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*10}  ─────")
                        print(f"[Fitting]     {'sol_b':>5}  {b:10.4f}  {b_err:10.4f}  {lb_s:>10}  {ub_s:>10}  {fx_s:>5}")
                        print(f"[Fitting]     SOL: y_SOL = y(ψ=1)·max(1 − {b:.3f}·(ψ−1), 0)")
                    # Pedestal summary (PRISM style — sentence-case, ± sign,
                    # parenthesized physical units).
                    ped_map = self._PEDESTAL_PARAMS.get(func_name, {})
                    if ped_map:
                        for line in self._format_ped_lines(
                                func_name, result, ped_map, ptype, x_axis,
                                entry, prefix="[Fitting]     "):
                            print(line)
                else:
                    fail_count += 1
                    print(f"[Fitting]   {entry} {ptype} ({func_name}): FAILED - {result.message}")

            self.fit_results[entry] = entry_results

        print(f"[Fitting] Done: {success_count} succeeded, {fail_count} failed")

        # Post-fitting hook (e.g. TCI validation in ne/Te tab)
        self._post_fitting(selected_entries, x_axis)

        # Re-plot with fit curves
        self._on_plot_clicked()

    def _post_fitting(self, selected_entries: List[str], x_axis: str) -> None:
        """Hook called after fitting completes. Override in subclass."""
        pass

    def _on_dt_toggled(self, checked):
        self._dt_enabled = checked
        self.fit_dt_entry.setEnabled(checked)

    def _get_dt_seconds(self) -> float:
        """Get effective dt in seconds (0 if toggle off or invalid)."""
        if not getattr(self, '_dt_enabled', False):
            return 0.0
        try:
            return max(0.0, float(self.fit_dt_entry.text())) / 1000.0
        except ValueError:
            return 0.0

    def _get_data_sampling_dt(self, entry: str) -> Optional[float]:
        """Get the time sampling interval (in seconds) for the data in this entry.
        Override in subclasses. Returns None if unknown."""
        return None

    def _get_efit_profile_data(self, entry: str, interp_func, dt_s: float = 0.0) -> Optional[Dict]:
        """Get profile data for fitting in EFIT coordinates.

        Must be implemented by subclasses.

        Args:
            entry: Listbox entry string
            interp_func: R -> EFIT coordinate interpolation function
            dt_s: time averaging half-window in seconds. 0 = no averaging.

        Returns:
            Dict with structure:
            {
                'Ti': {'x': ndarray, 'y': ndarray, 'err': ndarray or None},
                'vT': {'x': ndarray, 'y': ndarray, 'err': ndarray or None},
            }
            or None if data not available.
        """
        return None

    def _overlay_fit_curves(self, ax, param_type: str, selected_entries: List[str],
                           colors: List) -> None:
        """Overlay fitted curves on the current EFIT plot.

        Args:
            ax: matplotlib axes to plot on
            param_type: 'Ti', 'vT', 'Te', or 'ne'
            selected_entries: list of entry strings
            colors: list of colors matching entries
        """
        # Fit results are tied to the x-axis used during fitting
        # (psi_N / rho_pol / rho_tor). Skip overlay if user switched axes.
        current_x_axis = self._get_selected_x_axis()
        if self.fit_x_axis is not None and self.fit_x_axis != current_x_axis:
            return

        for i, entry in enumerate(selected_entries):
            if entry not in self.fit_results:
                continue
            entry_results = self.fit_results[entry]
            if param_type not in entry_results:
                continue

            result = entry_results[param_type]
            if not result.success or len(result.x_fit) == 0:
                continue

            color = colors[i]
            ax.plot(result.x_fit, result.y_fit, '-', color=color, linewidth=1.5,
                    alpha=0.8, zorder=15)

            # GPR uncertainty band
            if result.y_fit_std is not None:
                ax.fill_between(result.x_fit,
                                result.y_fit - 2 * result.y_fit_std,
                                result.y_fit + 2 * result.y_fit_std,
                                color=color, alpha=0.15, zorder=14)

    # ===== Fit preview overrides =====

    def _get_fit_profile_lines(self, selected_entries):
        """Get fitted profile data on psi_N grid with shot/time/coordinate columns"""
        from core.fitting import FIT_FUNCTIONS

        if not hasattr(self, 'fit_results') or not self.fit_results:
            return ["# No fit results available"]

        param_types = self._get_fit_param_types() if hasattr(self, '_get_fit_param_types') else []
        func_summary = " | ".join(f"{pt}={self.fit_func_names.get(pt, '?')}" for pt in param_types)
        lines = [f"# Fitted Profile Data ({func_summary})\n"]

        # Header with shot/time/source columns
        header_parts = ["shot", "time[ms]", "source", "psi_N", "rho_pol", "rho_tor", "R[m]"]
        for ptype in param_types:
            header_parts.append(f"{ptype}_fit")
        lines.append("#" + ",".join(f"{h:>12s}" for h in header_parts) + "\n")

        # Fixed psi_N grid: 0 to 1, 101 points
        psi_n_grid = np.linspace(0.0, 1.0, 101)

        for entry in selected_entries:
            if entry not in self.fit_results:
                continue

            # Parse shot, time, source
            try:
                shot, time_s = self._parse_entry_for_efit(entry)
                time_ms = time_s * 1e3
            except Exception:
                shot, time_ms = 0, 0.0
            try:
                _, _, source = self._parse_entry(entry)
            except Exception:
                source = ''

            entry_results = self.fit_results[entry]

            # Compute coordinate mappings from EFIT data
            R_vals = np.full_like(psi_n_grid, np.nan)
            rho_pol_vals = np.full_like(psi_n_grid, np.nan)
            rho_tor_vals = np.full_like(psi_n_grid, np.nan)

            if entry in self.efit_data:
                efit_entry = self.efit_data[entry]
                efit_R = efit_entry['R']
                efit_psi = efit_entry['psi_N']
                try:
                    psi_to_R = interp1d(efit_psi, efit_R, fill_value='extrapolate')
                    R_vals = psi_to_R(psi_n_grid)
                    rho_pol_interp = interp1d(efit_R, efit_entry['rho_pol'], fill_value='extrapolate')
                    rho_tor_interp = interp1d(efit_R, efit_entry['rho_tor'], fill_value='extrapolate')
                    rho_pol_vals = rho_pol_interp(R_vals)
                    rho_tor_vals = rho_tor_interp(R_vals)
                except Exception:
                    pass

            # Evaluate fit functions on psi_N grid (or current x_axis)
            x_axis = self._get_selected_x_axis()
            if x_axis == "psi_N":
                eval_x = psi_n_grid
            elif x_axis == "rho_pol":
                eval_x = rho_pol_vals
            else:
                eval_x = rho_tor_vals

            fit_vals = {}
            for ptype in param_types:
                if ptype not in entry_results or not entry_results[ptype].success:
                    fit_vals[ptype] = np.full_like(psi_n_grid, np.nan)
                    continue
                result = entry_results[ptype]
                func_name = result.func_name
                func_info = FIT_FUNCTIONS.get(func_name)
                if func_info and func_info['func'] is not None:
                    func = func_info['func']
                    param_vals = [result.params[pn] for pn in func_info['param_names']]
                    try:
                        fit_vals[ptype] = func(eval_x, *param_vals)
                    except Exception:
                        fit_vals[ptype] = np.full_like(psi_n_grid, np.nan)
                else:
                    # Spline: interpolate from stored x_fit/y_fit
                    try:
                        spl_interp = interp1d(result.x_fit, result.y_fit,
                                              fill_value='extrapolate')
                        fit_vals[ptype] = spl_interp(eval_x)
                    except Exception:
                        fit_vals[ptype] = np.full_like(psi_n_grid, np.nan)

            for j in range(len(psi_n_grid)):
                parts = [
                    f"{shot:>12d}",
                    f"{time_ms:12.1f}",
                    f"{source:>12s}",
                    f"{psi_n_grid[j]:12.6f}",
                    f"{rho_pol_vals[j]:12.6f}",
                    f"{rho_tor_vals[j]:12.6f}",
                    f"{R_vals[j]:12.6f}",
                ]
                for ptype in param_types:
                    parts.append(f"{fit_vals[ptype][j]:12.6f}")
                lines.append(",".join(parts) + "\n")

        return lines

    # Pedestal parameter mapping: {func_name: {role: param_name}}
    # For EPED, the pedestal symmetry point is implicit at psi_N = 1 - a3/2
    # (computed via _resolve_ped_position below).
    _PEDESTAL_PARAMS = {
        'mtanh': {'width': 'a3', 'height': 'a2', 'position': 'a4', 'foot': 'a1'},
        'ptanh': {'width': 'a6', 'height': 'a2', 'position': 'a7', 'foot': 'a1'},
        'EPED':  {'width': 'a3', 'height': 'a2', 'foot': 'a1'},
    }

    @staticmethod
    def _resolve_ped_position(func_name, params, param_errors):
        """Return (position, position_error) for the given fit function.

        For mtanh/ptanh this is a fitted parameter. For EPED the position is
        implicit at psi_N = 1 - a3/2 (a3 = width), with error propagated as
        |d position / d a3| * sigma(a3) = 0.5 * sigma(a3).
        """
        if func_name == 'mtanh':
            return params.get('a4'), param_errors.get('a4', 0.0)
        if func_name == 'ptanh':
            return params.get('a7'), param_errors.get('a7', 0.0)
        if func_name == 'EPED':
            w = params.get('a3')
            if w is None:
                return None, 0.0
            return 1.0 - 0.5 * w, 0.5 * param_errors.get('a3', 0.0)
        return None, 0.0

    def _format_ped_lines(self, func_name, result, ped_map, ptype, x_axis,
                          entry, prefix=""):
        """Format pedestal summary in PRISM style.

        Returns list of strings (no trailing newlines). Caller adds prefix
        like '[Fitting]     ' or '#   '. Output:

            Pedestal:
              Height      1.0663 keV
              Top         0.9912 keV  @ ψ_N = 0.9024  (R = 2.2050 m)
              Foot        0.3000 keV  @ ψ_N = 0.9720  (R = 2.2316 m)
              Width       0.0540 ± 0.0043 ψ_N    (1.14 cm)
              Location    0.9512 ± 0.0017 ψ_N    (R = 2.2278 m)
        """
        unit = self._PARAM_UNITS.get(ptype, '')
        coord_label = {'psi_N': 'ψ_N', 'rho_pol': 'ρ_pol',
                       'rho_tor': 'ρ_tor'}.get(x_axis, x_axis)

        w_pn = ped_map.get('width')
        w_val = result.params.get(w_pn) if w_pn else None
        w_err = result.param_errors.get(w_pn, 0.0) if w_pn else 0.0
        p_val, p_err = self._resolve_ped_position(
            func_name, result.params, result.param_errors)
        height_v = self._compute_ped_height(func_name, result.params, p_val)

        # Top/Foot edges (function evaluated at position ± width/2)
        top_val = foot_val = top_psi = foot_psi = None
        R_top = R_foot = None
        if p_val is not None and w_val is not None:
            top_psi = p_val - 0.5 * w_val
            foot_psi = p_val + 0.5 * w_val
            top_val = self._compute_ped_height(func_name, result.params, top_psi)
            foot_val = self._compute_ped_height(func_name, result.params, foot_psi)
            R_top, _ = self._get_ped_R_info(entry, top_psi, 0.0, x_axis)
            R_foot, _ = self._get_ped_R_info(entry, foot_psi, 0.0, x_axis)

        out = [f"{prefix}Pedestal:"]
        if height_v is not None:
            out.append(f"{prefix}  Height      {height_v:.4f} {unit}")
        if top_val is not None and top_psi is not None:
            line = (f"{prefix}  Top         {top_val:.4f} {unit}"
                    f"  @ {coord_label} = {top_psi:.4f}")
            if R_top is not None:
                line += f"  (R = {R_top:.4f} m)"
            out.append(line)
        if foot_val is not None and foot_psi is not None:
            line = (f"{prefix}  Foot        {foot_val:.4f} {unit}"
                    f"  @ {coord_label} = {foot_psi:.4f}")
            if R_foot is not None:
                line += f"  (R = {R_foot:.4f} m)"
            out.append(line)
        if w_val is not None:
            line = f"{prefix}  Width       {w_val:.4f} ± {w_err:.4f} {coord_label}"
            if p_val is not None:
                _, R_width = self._get_ped_R_info(entry, p_val, w_val, x_axis)
                if R_width is not None:
                    line += f"    ({R_width*100:.2f} cm)"
            out.append(line)
        if p_val is not None:
            line = f"{prefix}  Location    {p_val:.4f} ± {p_err:.4f} {coord_label}"
            R_pos, _ = self._get_ped_R_info(entry, p_val, w_val or 0, x_axis)
            if R_pos is not None:
                line += f"    (R = {R_pos:.4f} m)"
            out.append(line)
        return out

    @staticmethod
    def _compute_ped_height(func_name, params, position):
        """Pedestal height = function value at the pedestal symmetry point.

        Pedestal height = function value at the symmetry point — well defined
        and directly comparable across mtanh / ptanh / EPED parameterizations.
        """
        from core.fitting import FIT_FUNCTIONS
        info = FIT_FUNCTIONS.get(func_name)
        if info is None or info.get('func') is None:
            return None
        func = info['func']
        param_names = info['param_names']
        try:
            args = [params[n] for n in param_names]
        except KeyError:
            return None
        if position is None:
            return None
        try:
            import numpy as _np
            val = func(_np.asarray([position], dtype=float), *args)[0]
            return float(val)
        except Exception:
            return None

    # Unit labels per param type
    _PARAM_UNITS = {
        'Ti': 'keV', 'Te': 'keV',
        'ne': '10^19/m3', 'vT': 'km/s',
    }

    def _get_ped_R_info(self, entry: str, position: float, width: float,
                        x_axis: str = 'psi_N'):
        """Convert pedestal position/width from x_axis coordinate to R [m].

        Args:
            entry: Data entry key
            position: Pedestal position in x_axis coordinate
            width: Pedestal width in x_axis coordinate
            x_axis: Coordinate system ('psi_N', 'rho_pol', 'rho_tor')

        Returns:
            (R_pos, R_width) in meters, or (None, None) if EFIT data unavailable
        """
        if not self.efit_data or entry not in self.efit_data:
            return None, None
        efit_entry = self.efit_data[entry]
        efit_R = efit_entry['R']
        coord_key = x_axis if x_axis in efit_entry else 'psi_N'
        efit_coord = efit_entry[coord_key]
        try:
            coord_to_R = interp1d(efit_coord, efit_R, fill_value='extrapolate')
            R_pos = float(coord_to_R(position))
            R_top = float(coord_to_R(position - 0.5 * width))
            R_end = float(coord_to_R(position + 0.5 * width))
            R_width = abs(R_top - R_end)
            return R_pos, R_width
        except Exception:
            return None, None

    def _get_fit_params_lines(self, selected_entries):
        """Get fit parameters with function name per parameter type"""
        if not hasattr(self, 'fit_results') or not self.fit_results:
            return ["# No fit results available"]

        x_axis = self._get_selected_x_axis() if hasattr(self, '_get_selected_x_axis') else 'psi_N'
        lines = [f"# Fit Parameters\n"]
        param_types = self._get_fit_param_types() if hasattr(self, '_get_fit_param_types') else []

        for entry in selected_entries:
            if entry not in self.fit_results:
                continue
            # Include shot number and time (no commas to avoid header misparse)
            try:
                shot, time_s = self._parse_entry_for_efit(entry)
                lines.append(f"# {entry} | shot={shot} | time={time_s*1e3:.1f} ms\n")
            except Exception:
                lines.append(f"# {entry}\n")
            entry_results = self.fit_results[entry]

            for ptype in param_types:
                if ptype not in entry_results:
                    continue
                result = entry_results[ptype]
                ped_map = self._PEDESTAL_PARAMS.get(result.func_name, {})
                status = 'OK' if result.success else 'FAILED'
                lines.append(f"# {ptype} — function: {result.func_name} | chi2={result.chi_squared:.4f} | {status}\n")
                if result.success:
                    lines.append(f"#{'Param':>10s},{'Value':>12s},{'Error':>12s}\n")
                    for pname, val in result.params.items():
                        err = result.param_errors.get(pname, 0.0)
                        lines.append(f" {pname:>10s},{val:12.6f},{err:12.6f}\n")
                    # Pedestal summary (PRISM style).
                    if ped_map:
                        ped_lines_raw = self._format_ped_lines(
                            result.func_name, result, ped_map, ptype,
                            x_axis, entry, prefix="#   ")
                        ped_lines = [ln + "\n" for ln in ped_lines_raw]
                        for pl in ped_lines:
                            lines.append(pl)

        return lines

    # ===== p-File format export =====

    # Mapping from PRISM param types to p-file variable names and units
    _PFILE_VARNAMES = {
        'Te': ('te', 'keV', 'dte'),
        'ne': ('ne', '10^20/m^3', 'dne'),
        'Ti': ('ti', 'keV', 'dti'),
        'vT': ('vtor1', 'km/s', 'dvtor'),
    }

    # ne scale: PRISM uses 10^19/m^3 internally, p-file uses 10^20/m^3
    _PFILE_NE_SCALE = 0.1

    def _get_pfile_lines(self, selected_entries):
        """Generate p-file (PEQDSK) format text from fitted profiles"""
        from core.fitting import FIT_FUNCTIONS

        if not hasattr(self, 'fit_results') or not self.fit_results:
            return ["# No fit results available\n"]

        param_types = self._get_fit_param_types() if hasattr(self, '_get_fit_param_types') else []
        psi_n_grid = np.linspace(0.0, 1.0, 401)
        npts = len(psi_n_grid)

        lines = []

        for entry in selected_entries:
            if entry not in self.fit_results:
                continue

            try:
                shot, time_s = self._parse_entry_for_efit(entry)
                time_ms = int(round(time_s * 1e3))
            except Exception:
                shot, time_ms = 0, 0

            lines.append(f"# p-file format: shot={shot}, time={time_ms}ms\n")

            entry_results = self.fit_results[entry]

            # Evaluate each param on psi_N grid
            x_axis = self._get_selected_x_axis()
            eval_x = psi_n_grid
            if x_axis != "psi_N" and entry in self.efit_data:
                efit_entry = self.efit_data[entry]
                try:
                    psi_to_R = interp1d(efit_entry['psi_N'], efit_entry['R'], fill_value='extrapolate')
                    R_vals = psi_to_R(psi_n_grid)
                    coord_interp = interp1d(efit_entry['R'], efit_entry[x_axis], fill_value='extrapolate')
                    eval_x = coord_interp(R_vals)
                except Exception:
                    pass

            for ptype in param_types:
                if ptype not in entry_results or not entry_results[ptype].success:
                    continue

                result = entry_results[ptype]
                varname, unit, dname = self._PFILE_VARNAMES.get(ptype, (ptype.lower(), '', 'd' + ptype.lower()))

                # Evaluate fit on grid
                func_info = FIT_FUNCTIONS.get(result.func_name)
                if func_info and func_info['func'] is not None:
                    func = func_info['func']
                    param_vals = [result.params[pn] for pn in func_info['param_names']]
                    try:
                        y_vals = func(eval_x, *param_vals)
                    except Exception:
                        continue
                else:
                    try:
                        spl = interp1d(result.x_fit, result.y_fit, fill_value='extrapolate')
                        y_vals = spl(eval_x)
                    except Exception:
                        continue

                # Convert ne from 10^19/m^3 to 10^20/m^3
                if ptype == 'ne':
                    y_vals = y_vals * self._PFILE_NE_SCALE

                # Compute derivative d(var)/d(psiN)
                dy_dpsi = np.gradient(y_vals, psi_n_grid)

                # Write block header
                lines.append(f" {npts} psinorm {varname}({unit}) {dname}/dpsiN\n")

                # Write data rows: 1 space + psiN + 3 spaces + value + 3 spaces + derivative
                for j in range(npts):
                    lines.append(f" {psi_n_grid[j]:.6f}   {y_vals[j]:.6f}   {dy_dpsi[j]:.6f}\n")

            lines.append("\n")

        return lines

    # ===== Canvas click-toggle for channel enable/disable =====

    def _clear_click_points(self) -> None:
        """Clear registered click points (call at start of each plot method)"""
        self._click_points.clear()

    def _register_click_points(self, ax, x_data, y_data, channel_keys: List[str]) -> None:
        """Register data points for interactive channel toggle via double-click.

        Args:
            ax: matplotlib axes the points are plotted on
            x_data: x coordinates (array-like)
            y_data: y coordinates (array-like)
            channel_keys: list of channel key strings (e.g., ["CES_0", "CES_1", ...])
        """
        x_arr = np.atleast_1d(x_data)
        y_arr = np.atleast_1d(y_data)
        for x, y, key in zip(x_arr, y_arr, channel_keys):
            if np.isfinite(x) and np.isfinite(y):
                self._click_points.append((float(x), float(y), key, ax))

    def _on_canvas_dblclick(self, event) -> None:
        """Toggle channel enabled/disabled on double-click near a data point"""
        if not event.dblclick or event.inaxes not in [self.ax1, self.ax2]:
            return
        if not self._click_points:
            return

        ax = event.inaxes
        click_disp = ax.transData.transform((event.xdata, event.ydata))

        min_dist = float('inf')
        nearest_key = None

        for x, y, key, pt_ax in self._click_points:
            if pt_ax is not ax:
                continue
            pt_disp = ax.transData.transform((x, y))
            dist = ((click_disp[0] - pt_disp[0])**2 + (click_disp[1] - pt_disp[1])**2)**0.5
            if dist < min_dist:
                min_dist = dist
                nearest_key = key

        if min_dist > 15 or nearest_key is None:  # 15 pixel threshold
            return

        # Resolve display label from channel info
        label_map = dict(self._get_channel_info())
        display = label_map.get(nearest_key, nearest_key)

        if nearest_key in self._disabled_channels:
            self._disabled_channels.discard(nearest_key)
            print(f"[Channel] Enabled: {display}")
        else:
            self._disabled_channels.add(nearest_key)
            print(f"[Channel] Disabled: {display}")

        # Re-plot current view, preserving axis limits
        if self.ax1.lines or self.ax1.collections:
            xlim1, ylim1 = self.ax1.get_xlim(), self.ax1.get_ylim()
            xlim2, ylim2 = self.ax2.get_xlim(), self.ax2.get_ylim()
            self._on_plot_clicked()
            self.ax1.set_xlim(xlim1); self.ax1.set_ylim(ylim1)
            self.ax2.set_xlim(xlim2); self.ax2.set_ylim(ylim2)
            self.canvas.draw_idle()

    @abstractmethod
    def _create_shot_input(self, parent: QWidget) -> None:
        """Create diagnostic-specific shot input controls

        Must be implemented by subclasses to provide diagnostic-specific
        options (e.g., analysis type dropdown for CES, diagnostic selection
        for Thomson/ECE).
        """
        pass

    @abstractmethod
    def _get_secondary_loader(self) -> Any:
        """Get secondary data loader (XICS for Ion, ECE for Electron)

        Returns:
            Loader instance for secondary diagnostic, or None
        """
        pass

    @abstractmethod
    def _parse_entry(self, entry: str) -> Tuple[int, float, str]:
        """Parse listbox entry to extract shot, time, and source

        Args:
            entry: Listbox entry string (e.g., "012345_001000 (mod)")

        Returns:
            Tuple of (shot_number, time_point, source)
        """
        pass

    def plot_data(self) -> None:
        """Plot R profiles with primary and secondary data"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_xlabel('R [m]')
        self.ax2.set_xlabel('R [m]')
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        y1_max, y2_max, y2_min = 0.0, 0.0, 0.0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                y1_max, y2_max, y2_min = self._plot_single_entry(
                    entry, colors[i], y1_max, y2_max, y2_min)
            except Exception as e:
                print(f"[{self.TAB_NAME}] Error plotting {entry}: {str(e)}")

        # Apply axis limits
        self._apply_plot_limits(y1_max, y2_max, y2_min)

        self.plot_manager.apply_common_styling(
            self.ax1, self.ax2,
            legend_fontsize=self.legend_fontsize,
            label_fontsize=self.label_fontsize,
            tick_fontsize=self.tick_fontsize)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def _plot_single_entry(
        self,
        entry: str,
        color: str,
        y1_max: float,
        y2_max: float,
        y2_min: float
    ) -> Tuple[float, float, float]:
        """Plot single entry - diagnostic specific logic

        Args:
            entry: Listbox entry string
            color: Plot color for this entry
            y1_max: Current maximum value for y1 axis
            y2_max: Current maximum value for y2 axis
            y2_min: Current minimum value for y2 axis

        Returns:
            Updated tuple of (y1_max, y2_max, y2_min)

        Note:
            Override this method OR override plot_data() entirely.
        """
        raise NotImplementedError("Subclass must implement _plot_single_entry or override plot_data")

    def _apply_plot_limits(
        self,
        y1_max: float,
        y2_max: float,
        y2_min: float
    ) -> None:
        """Apply axis limits to both plots

        Can be overridden by subclasses for diagnostic-specific handling
        (e.g., velocity plots with negative values).
        """
        self.ax1.set_xlim(self.app_config.R_LIMITS)
        self.ax2.set_xlim(self.app_config.R_LIMITS)

        y1_margin = y1_max * 0.1 if y1_max > 0 else 0.1
        self.ax1.set_ylim(0, y1_max + y1_margin)

        self._apply_y2_limits(y2_max, y2_min)

    def _apply_y2_limits(self, y2_max: float, y2_min: float) -> None:
        """Apply y2 axis limits - can be overridden by subclasses

        Default implementation assumes non-negative values (like ne, Te).
        Override for parameters that can be negative (like vT).
        """
        y2_margin = y2_max * 0.1 if y2_max > 0 else 0.1
        self.ax2.set_ylim(0, y2_max + y2_margin)

    def plot_efit_profiles(self) -> None:
        """Plot profiles with EFIT mapping"""
        if not self.efit_data or self.computed_efit_tree is None:
            QMessageBox.warning(self.frame, "Warning", "Please compute EFIT first.")
            return

        efit_tree = self.computed_efit_tree

        self.ax1.clear()
        self.ax2.clear()

        selected_entries = [self.selected_listbox.item(i).text()
                           for i in range(self.selected_listbox.count())]
        if not selected_entries:
            return

        x_axis = self._get_selected_x_axis()

        if x_axis == "psi_N":
            x_label = rf"$\psi_N$ ({efit_tree})"
        elif x_axis == "rho_pol":
            x_label = rf"$\rho_{{pol}}$ ({efit_tree})"
        else:
            x_label = rf"$\rho_{{tor}}$ ({efit_tree})"

        y1_max, y2_max, y2_min = 0.0, 0.0, 0.0
        colors = self._get_plot_colors(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                if entry not in self.efit_data:
                    continue

                efit_entry = self.efit_data[entry]
                efit_R = efit_entry['R']

                if x_axis == "psi_N":
                    interp_func = interp1d(efit_R, efit_entry['psi_N'], fill_value='extrapolate')
                elif x_axis == "rho_pol":
                    interp_func = interp1d(efit_R, efit_entry['rho_pol'], fill_value='extrapolate')
                else:
                    interp_func = interp1d(efit_R, efit_entry['rho_tor'], fill_value='extrapolate')

                y1_max, y2_max, y2_min = self._plot_single_efit_entry(
                    entry, colors[i], interp_func, y1_max, y2_max, y2_min)

            except Exception as e:
                print(f"[{self.TAB_NAME}] Error plotting {entry}: {str(e)}")

        # Set axis labels and limits
        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])

        # Expand x-axis if fitting range exceeds default limits
        x_left, x_right = -0.1, 1.1
        try:
            fit_left = float(self.fit_xmin_entry.text()) - 0.1
            fit_right = float(self.fit_xmax_entry.text()) + 0.1
            if fit_left < x_left:
                x_left = fit_left
            if fit_right > x_right:
                x_right = fit_right
        except (ValueError, AttributeError):
            pass
        self.ax1.set_xlim(x_left, x_right)
        self.ax2.set_xlim(x_left, x_right)
        self.ax1.set_ylim(0, y1_max * 1.1)
        self._apply_y2_limits_efit(y2_max, y2_min)

        # Shade outside plasma region (core & SOL)
        for ax in [self.ax1, self.ax2]:
            ax.axvspan(-1, 0, color='gray', alpha=0.15, zorder=0)
            ax.axvspan(1, 2, color='gray', alpha=0.15, zorder=0)

        self.plot_manager.apply_common_styling(
            self.ax1, self.ax2,
            legend_fontsize=self.legend_fontsize,
            label_fontsize=self.label_fontsize,
            tick_fontsize=self.tick_fontsize)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def _plot_single_efit_entry(
        self,
        entry: str,
        color: str,
        interp_func: Any,
        y1_max: float,
        y2_max: float,
        y2_min: float
    ) -> Tuple[float, float, float]:
        """Plot single entry with EFIT mapping

        Args:
            entry: Listbox entry string
            color: Plot color for this entry
            interp_func: Interpolation function R -> coordinate
            y1_max: Current maximum value for y1 axis
            y2_max: Current maximum value for y2 axis
            y2_min: Current minimum value for y2 axis

        Returns:
            Updated tuple of (y1_max, y2_max, y2_min)

        Note:
            Override this method OR override plot_efit_profiles() entirely.
        """
        raise NotImplementedError("Subclass must implement _plot_single_efit_entry or override plot_efit_profiles")

    def _apply_y2_limits_efit(self, y2_max: float, y2_min: float) -> None:
        """Apply y2 axis limits for EFIT plots - can be overridden"""
        self.ax2.set_ylim(0, y2_max * 1.1)

    def _get_marker_style(self, source: str) -> Dict[str, Any]:
        """Get marker style based on data source

        Args:
            source: Data source identifier

        Returns:
            Dictionary with marker style properties
        """
        # Default styles - can be overridden by subclasses
        styles = {
            'default': {
                'marker': 'o',
                'markerfacecolor': None,
                'markeredgewidth': 1,
                'fillstyle': 'full'
            },
            'empty': {
                'marker': 'o',
                'markerfacecolor': 'none',
                'markeredgewidth': 1.5,
                'fillstyle': 'full'
            },
            'square': {
                'marker': 's',
                'markerfacecolor': None,
                'markeredgewidth': 1,
                'fillstyle': 'full'
            },
            'half': {
                'marker': 'o',
                'markerfacecolor': None,
                'markeredgewidth': 1,
                'fillstyle': 'right'
            }
        }
        return styles.get('default', styles['default'])
