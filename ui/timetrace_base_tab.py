#!/usr/bin/python3.8

"""
Base class for time trace tabs (TiVT, NeTe, MSE)
Extracts common functionality from concrete time trace tabs
"""

from abc import abstractmethod
import tkinter as tk
from tkinter import ttk, messagebox
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ui.base_tab import BaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, PAD_X, PAD_Y
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
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)

        self.ax1, self.ax2 = self.plot_manager.setup_timetrace_plot(
            self.figure, self.param1['label'], self.param2['label'])

        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()

        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill='both', expand=True)

        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)

        # Create control panels
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_axis_controls(control_frame)
        self._create_save_controls(control_frame)

    def _create_plot_controls(self, parent: ttk.Frame) -> None:
        """Create plot control buttons (default implementation)

        Can be overridden by subclasses that need additional controls
        (e.g., MSE with q/j parameter selection).
        """
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        ttk.Button(frame, text='Plot time traces', command=self.plot_data).pack(
            fill='x', padx=PAD_X, pady=PAD_Y)

    @abstractmethod
    def _create_shot_input(self, parent: ttk.Frame) -> None:
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
        messagebox.showinfo("Info", "EFIT mapping not available for time trace tabs")

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
