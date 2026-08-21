"""
Base class for time trace tabs (Ion, Electron, MSE)
Extracts common functionality from concrete time trace tabs
"""

from abc import abstractmethod
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QSplitter, QScrollArea,
    QGroupBox, QComboBox, QLabel, QPushButton, QMessageBox, QFrame,
    QSpinBox, QDialog, QDialogButtonBox,
)
from PySide6.QtCore import Qt

from ui.tabs.base_tab import BaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, apply_dark_figure_style
from plotting.plot_manager import apply_legend_with_limit, TIMETRACE_LEGEND_LIMIT


class TimeTraceBaseTab(BaseTab):
    """Base class for all time trace visualization tabs

    Provides common functionality for time trace-type tabs:
    - Dual-axis plot setup (param1 on top, param2 on bottom)
    - Time series plotting
    - Common plot controls
    """

    def create_widgets(self) -> None:
        """Create common time trace tab layout"""
        self.figure = Figure(self.app_config.FIGURE_SIZE)

        self.ax1, self.ax2 = self.plot_manager.setup_timetrace_plot(
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
        control_layout.setSizeConstraint(QVBoxLayout.SetMinimumSize)

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
        self._create_save_controls(control_frame, section_num=4)
        control_layout.addStretch()

    def _create_plot_controls(self, parent: QWidget) -> None:
        """Create the '3. Plot' group.

        An optional Y-axis (bottom) parameter combo (Ion: vT/ωT, MSE: q/j — via
        the shared `_yaxis_param_spec` hook) sits above the Plot | Style row. Tabs
        without a bottom-panel selector (Electron, Neutron) show only the buttons.
        """
        group = QGroupBox("3. Plot")
        group_layout = QVBoxLayout(group)

        # Optional Y-axis (bottom) selector row — shared mechanism with profiles
        label, combo = self._build_yaxis_combo()
        if combo is not None:
            row_y = QHBoxLayout()
            row_y.setContentsMargins(0, 0, 0, 0)
            row_y.addWidget(QLabel(label), 1)
            row_y.addWidget(combo, 1)
            group_layout.addLayout(row_y)

        # Plot + Style buttons in same row
        plot_row = QHBoxLayout()

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        plot_row.addWidget(plot_button, 3)

        style_btn = QPushButton("Style")
        style_btn.clicked.connect(self._show_style_dialog)
        plot_row.addWidget(style_btn, 1)

        group_layout.addLayout(plot_row)

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

        parent.layout().addWidget(group)

    @abstractmethod
    def _create_shot_input(self, parent: QWidget) -> None:
        """Create diagnostic-specific shot input controls

        Must be implemented by subclasses to provide diagnostic-specific
        options (e.g., analysis type dropdown for CES, diagnostic selection
        for Thomson/ECE).
        """
        pass

    @abstractmethod
    def load_shot_data(self) -> None:
        """Load shot data from MDS+ or file

        Must be implemented by subclasses with diagnostic-specific
        loading logic.
        """
        pass

    @abstractmethod
    def _parse_entry(self, entry: str):
        """Parse listbox entry to extract shot, radius, and source info

        Must be implemented by subclasses with diagnostic-specific
        parsing logic.
        """
        pass

    @abstractmethod
    def plot_data(self) -> None:
        """Plot time trace data

        Must be implemented by subclasses with diagnostic-specific
        plotting logic.
        """
        pass

    @abstractmethod
    def _write_data_to_file(self, file_path: str, selected_entries: list) -> None:
        """Write data to file

        Must be implemented by subclasses with diagnostic-specific
        export format.
        """
        pass

    def _parse_entry_for_efit(self, entry: str):
        """Not used for time trace tabs"""
        return None, None

    def plot_efit_profiles(self) -> None:
        """Not applicable for time trace tabs"""
        QMessageBox.information(self.frame, "Info", "EFIT mapping not available for time trace tabs")

    @staticmethod
    def _parse_color_mode(text):
        """Extract colormap name from dropdown text like 'Fixed(tab10)' → 'tab10'"""
        start = text.find('(')
        end = text.find(')')
        if start != -1 and end != -1:
            return text[start + 1:end]
        return 'tab10'

    def _get_plot_colors(self, entries):
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
            if self.ax1.lines:
                self.plot_data()
        dialog.accepted.connect(_apply)
        self._style_dialog = dialog
        dialog.show()

    def _finalize_plot(self) -> None:
        """Apply common styling and update canvas after plotting

        Call this at the end of plot_data() in subclasses.
        """
        legend_fs = getattr(self, 'legend_fontsize', 8)
        label_fs = getattr(self, 'label_fontsize', 12)
        tick_fs = getattr(self, 'tick_fontsize', 10)

        for ax in [self.ax1, self.ax2]:
            apply_legend_with_limit(ax, TIMETRACE_LEGEND_LIMIT,
                                    frameon=False, fontsize=legend_fs)

        self.plot_manager.apply_common_styling(
            self.ax1, self.ax2, skip_legend=True,
            label_fontsize=label_fs, tick_fontsize=tick_fs)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()
