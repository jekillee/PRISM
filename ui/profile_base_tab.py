"""
Base class for profile tabs (TiVT, NeTe, MSE)
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
)
from PySide6.QtCore import Qt

from ui.base_tab import BaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style


class ProfileBaseTab(BaseTab):
    """Base class for all profile visualization tabs

    Provides common functionality for profile-type tabs:
    - Dual-axis plot setup (param1 on left, param2 on right)
    - R profile plotting with EFIT mapping
    - Secondary diagnostic overlay (XICS for TiVT, ECE for NeTe)
    """

    # Override in subclass for console logging prefix
    TAB_NAME = "Profile"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.secondary_data_cache: Dict[str, Any] = {}
        self._disabled_channels: set = set()  # Channel indices to plot in gray

    def create_widgets(self) -> None:
        """Create common profile tab layout"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)

        self.ax1, self.ax2 = self.plot_manager.setup_profile_plot(
            self.figure, self.param1['label'], self.param2['label'])
        apply_dark_figure_style(self.figure)

        # Create canvas
        self.canvas = FigureCanvasQTAgg(self.figure)
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
        self._create_plot_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_save_controls(control_frame, section_num=5)
        control_layout.addStretch()

    def _create_plot_controls(self, parent: QWidget) -> None:
        """Create plot control buttons (common for all profile tabs)"""
        group = QGroupBox("3. Plot")
        group_layout = QVBoxLayout(group)

        # Plot + Option buttons in same row
        plot_row = QHBoxLayout()

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        plot_row.addWidget(plot_button, 3)

        style_btn = QPushButton("Option")
        style_btn.clicked.connect(self._show_style_dialog)
        plot_row.addWidget(style_btn, 1)

        group_layout.addLayout(plot_row)

        # Separator line
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        group_layout.addWidget(separator)

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

        # Show Nodes | Select Channels
        row2 = QHBoxLayout()

        self.show_channel_checkbox = QCheckBox("Show Nodes")
        self.show_channel_checkbox.setChecked(False)
        row2.addWidget(self.show_channel_checkbox)

        channels_btn = QPushButton("Select Channels")
        channels_btn.setToolTip("Select which channels to enable or dim")
        channels_btn.clicked.connect(self._show_channel_selector)
        row2.addWidget(channels_btn)

        group_layout.addLayout(row2)

        parent.layout().addWidget(group)

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

        dialog = QDialog(self.frame)
        dialog.setWindowTitle("Select Channels")
        dialog.setMinimumWidth(300)
        dlg_layout = QVBoxLayout(dialog)

        hint = QLabel("Uncheck channels to dim them (gray) on the plot:")
        hint.setWordWrap(True)
        dlg_layout.addWidget(hint)

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

        # Select All / Deselect All buttons
        sel_layout = QHBoxLayout()
        select_all_btn = QPushButton("Select All")
        select_all_btn.clicked.connect(
            lambda: [cb.setChecked(True) for _, cb in checkboxes])
        sel_layout.addWidget(select_all_btn)
        deselect_all_btn = QPushButton("Deselect All")
        deselect_all_btn.clicked.connect(
            lambda: [cb.setChecked(False) for _, cb in checkboxes])
        sel_layout.addWidget(deselect_all_btn)
        dlg_layout.addLayout(sel_layout)

        # OK / Cancel
        btn_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        btn_box.accepted.connect(dialog.accept)
        btn_box.rejected.connect(dialog.reject)
        dlg_layout.addWidget(btn_box)

        if dialog.exec() == QDialog.Accepted:
            self._disabled_channels.clear()
            for idx, cb in checkboxes:
                if not cb.isChecked():
                    self._disabled_channels.add(idx)

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

        if dialog.exec() == QDialog.Accepted:
            self.color_mode_combo.setCurrentText(color_combo.currentText())
            self.label_fontsize = label_spin.value()
            self.legend_fontsize = legend_spin.value()
            self.tick_fontsize = tick_spin.value()
            # Auto-apply if plot exists
            if self.ax1.lines or self.ax1.collections:
                self.plot_data()

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
        """Get secondary data loader (XICS for TiVT, ECE for NeTe)

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

        x_axis = self.selected_x_axis.get()

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

        self.ax1.set_xlim(0, 1.05)
        self.ax2.set_xlim(0, 1.05)
        self.ax1.set_ylim(0, y1_max * 1.1)
        self._apply_y2_limits_efit(y2_max, y2_min)

        # Add vertical lines at core and LCFS
        for ax in [self.ax1, self.ax2]:
            ax.axvline(x=0, c='k', ls='--')
            ax.axvline(x=1, c='k', ls='--')

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
