#!/usr/bin/python3.8

"""
Base class for profile tabs (TiVT, NeTe, MSE)
Extracts common functionality from concrete profile tabs
"""

from abc import abstractmethod
from typing import Dict, List, Tuple, Any, Optional
import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ui.base_tab import BaseTab
from ui.ui_constants import CONTROL_PANEL_WIDTH, PAD_X, PAD_Y
from ui.widgets.custom_toolbar import AxisControlToolbar


class ProfileBaseTab(BaseTab):
    """Base class for all profile visualization tabs

    Provides common functionality for profile-type tabs:
    - Dual-axis plot setup (param1 on left, param2 on right)
    - R profile plotting with EFIT mapping
    - Secondary diagnostic overlay (XICS for TiVT, ECE for NeTe)
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.secondary_data_cache: Dict[str, Any] = {}

    def create_widgets(self) -> None:
        """Create common profile tab layout"""
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)

        self.ax1, self.ax2 = self.plot_manager.setup_profile_plot(
            self.figure, self.param1['label'], self.param2['label'])

        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()

        # Create toolbar frame to hold canvas and toolbar
        plot_frame = ttk.Frame(self.frame)
        plot_frame.pack(side=tk.LEFT, fill='both', expand=True)

        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill='both', expand=True, in_=plot_frame)

        # Add axis control toolbar
        self.toolbar = AxisControlToolbar(self.canvas, plot_frame, tab_instance=self)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill='x', in_=plot_frame)
        self.toolbar.configure_axes(
            has_y2=self.param2 is not None,
            ax1_label=self.param1['label'],
            ax2_label=self.param2['label'] if self.param2 else 'Axes 2'
        )

        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)

        # Create control panels
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_save_controls(control_frame)

    def _create_plot_controls(self, parent: ttk.Frame) -> None:
        """Create plot control buttons (common for all profile tabs)"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        row_frame = ttk.Frame(frame)
        row_frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        row_frame.grid_columnconfigure(1, weight=1)

        self.show_channel_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            row_frame, text='Show Nodes', variable=self.show_channel_var
        ).grid(row=0, column=0, sticky='w')

        ttk.Button(
            row_frame, text='Plot R profiles', command=self.plot_data
        ).grid(row=0, column=1, sticky='ew', padx=(10, 0))

    @abstractmethod
    def _create_shot_input(self, parent: ttk.Frame) -> None:
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

        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return

        y1_max, y2_max, y2_min = 0.0, 0.0, 0.0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                y1_max, y2_max, y2_min = self._plot_single_entry(
                    entry, colors[i], y1_max, y2_max, y2_min)
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")

        # Apply axis limits
        self._apply_plot_limits(y1_max, y2_max, y2_min)

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
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
            messagebox.showwarning("Warning", "Please compute EFIT first.")
            return

        efit_tree = self.computed_efit_tree

        self.ax1.clear()
        self.ax2.clear()

        selected_entries = list(self.selected_listbox.get(0, tk.END))
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
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)

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
                print(f"Error plotting {entry}: {str(e)}")

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

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
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
