"""
Base class for time trace tabs (TiVT, NeTe, MSE)
Extracts common functionality from concrete time trace tabs
"""

from abc import abstractmethod
from matplotlib.figure import Figure
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg

from PySide6.QtWidgets import (
    QWidget, QVBoxLayout, QSplitter, QScrollArea,
    QGroupBox, QPushButton, QMessageBox,
)
from PySide6.QtCore import Qt

from ui.base_tab import BaseTab
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

        # Create control panels
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_save_controls(control_frame, section_num=4)

    def _create_plot_controls(self, parent: QWidget) -> None:
        """Create plot control buttons (default implementation)

        Can be overridden by subclasses that need additional controls
        (e.g., MSE with q/j parameter selection).
        """
        group = QGroupBox("3. Plot")
        group_layout = QVBoxLayout(group)

        plot_button = QPushButton("Plot")
        plot_button.clicked.connect(self.plot_data)
        group_layout.addWidget(plot_button)

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

    def _finalize_plot(self) -> None:
        """Apply common styling and update canvas after plotting

        Call this at the end of plot_data() in subclasses.
        """
        for ax in [self.ax1, self.ax2]:
            apply_legend_with_limit(ax, TIMETRACE_LEGEND_LIMIT,
                                    frameon=False, fontsize=8)

        self.plot_manager.apply_common_styling(self.ax1, self.ax2, skip_legend=True)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()
