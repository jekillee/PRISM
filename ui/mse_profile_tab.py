#!/usr/bin/python3.8

"""
MSE Profile tab with j/q selection
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ui.profile_base_tab import ProfileBaseTab
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y, ENTRY_WIDTH_AXIS, LABEL_WIDTH_SHORT
)
from ui.widgets.custom_toolbar import AxisControlToolbar


class MSEProfileTab(ProfileBaseTab):
    """MSE Profile tab with TGAMMA and j/q profiles"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _get_secondary_loader(self):
        """MSE has no secondary diagnostic"""
        return None

    def create_widgets(self):
        """Create MSE profile tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)

        # Setup 1x2 plot with shared x-axis: TGAMMA (left), j or q (right)
        self.ax1 = self.figure.add_subplot(121)
        self.ax2 = self.figure.add_subplot(122, sharex=self.ax1)

        self.ax1.set_xlabel('R [m]')
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('R [m]')

        # Create plot frame to hold canvas and toolbar
        plot_frame = ttk.Frame(self.frame)
        plot_frame.pack(side=tk.LEFT, fill='both', expand=True)

        # Create canvas inside plot_frame
        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.draw()

        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.TOP, fill='both', expand=True)

        # Add axis control toolbar at bottom of plot_frame
        self.toolbar = AxisControlToolbar(self.canvas, plot_frame, tab_instance=self)
        self.toolbar.update()
        self.toolbar.pack(side=tk.BOTTOM, fill='x')
        self.toolbar.configure_axes(has_y2=True, ax1_label='γ [rad]', ax2_label='q / j')

        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)

        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_save_controls(control_frame)

    def _create_shot_input(self, parent):
        """Create data loading section"""
        frame = ttk.LabelFrame(parent, text="1. Load MSE Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        frame.grid_columnconfigure(1, weight=1)

        # Row 0: Shot label, entry, Fetch
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_SHORT, anchor='e').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='e')
        self.shot_entry = ttk.Entry(frame)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self.load_shot_data())

        self.fetch_button = ttk.Button(frame, text='Fetch', command=self.load_shot_data, width=8)
        self.fetch_button.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')

    def _create_plot_controls(self, parent):
        """Create plot control buttons"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        row_frame = ttk.Frame(frame)
        row_frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        row_frame.grid_columnconfigure(2, weight=1)

        self.show_channel_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row_frame, text='Show Nodes', variable=self.show_channel_var).grid(
            row=0, column=0, sticky='w')

        self.selected_param = tk.StringVar(value='q')
        param_dropdown = ttk.Combobox(row_frame, textvariable=self.selected_param,
                                      values=['q', 'j'], state='readonly', width=3)
        param_dropdown.grid(row=0, column=1, sticky='w', padx=(10, 0))

        ttk.Button(row_frame, text='Plot R profiles', command=self.plot_data).grid(
            row=0, column=2, sticky='ew', padx=(10, 0))

    def load_shot_data(self):
        """Load MSE shot data from MDS+"""
        try:
            shot_number = int(self.shot_entry.get())

            data = self.data_loader.load_data(shot_number)

            cache_key = f'{shot_number}_MSE'
            self.data[cache_key] = data

            self.available_listbox.delete(0, tk.END)

            # Use profile time for time selection
            for tp in data.time_prof:
                item_str = f'{shot_number:06d}_{tp*1e3:06.0f} (MSE)'
                self.available_listbox.insert(tk.END, item_str)

            print(f"MSE data loaded: {len(data.time_prof)} timepoints")
            print(f"  NB source: {data.nb_source}")

        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data: {str(e)}")

    def _parse_entry(self, entry):
        """Parse entry to get shot and time"""
        # entry format: '034142_000730 (MSE)'
        # Remove ' (MSE)' suffix first
        entry_clean = entry.replace(' (MSE)', '')
        parts = entry_clean.split('_')
        shot_number = int(parts[0])
        time_point = float(parts[1]) / 1e3  # ms to s
        return shot_number, time_point, 'MSE'

    def _parse_entry_for_efit(self, entry):
        """Parse entry for EFIT (extract shot and time)"""
        shot_number, time_point, _ = self._parse_entry(entry)
        return shot_number, time_point

    def _get_r_edge_at_time(self, data, time_point):
        """Get R_edge at specific time by interpolation"""
        r_edge_info = data.measurements['r_edge']
        interp_func = interp1d(r_edge_info['time'], r_edge_info['data'],
                               fill_value='extrapolate')
        return float(interp_func(time_point))

    def plot_data(self):
        """Plot R profiles"""
        self.ax1.clear()
        self.ax2.clear()

        self.ax1.set_xlabel('R [m]')
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('R [m]')

        # Set y-label based on selected parameter
        param = self.selected_param.get()
        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')

        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return

        gamma_min, gamma_max = np.inf, -np.inf
        param_min, param_max = np.inf, -np.inf
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, _ = self._parse_entry(entry)
                color = colors[i]

                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                # Get time indices
                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))

                # Get R_edge and magx at this time
                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                # =============================================================
                # Plot TGAMMA vs R (left plot)
                # =============================================================
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                label = f'#{shot_number} {int(time_point*1e3):06d}ms (MSE)'

                self.ax1.errorbar(R_raw, tgamma, xerr=drr*0.5, yerr=sgamma*0.5,
                                 fmt='o-', capsize=3, markersize=5, color=color,
                                 linewidth=0.8, label=label)

                # Add channel labels for TGAMMA (node: TGAMMA01~25)
                raw_channels = [i+1 for i in range(len(good_mask)) if good_mask[i]]
                self._add_channel_labels(self.ax1, R_raw, tgamma, 'TGAMMA', raw_channels)

                gamma_min = min(gamma_min, np.nanpercentile(tgamma, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma, 98))

                # =============================================================
                # Plot j or q vs R (right plot)
                # - Solid line: pmse profile (20 points)
                # - Markers: interpolated values at TGAMMA R positions
                # =============================================================
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]

                # Convert r/a to R
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)

                # Sort by R for proper line plotting
                sort_idx = np.argsort(R_prof)
                R_prof_sorted = R_prof[sort_idx]
                param_data_sorted = param_data[sort_idx]
                param_err_sorted = param_err[sort_idx]

                # Plot pmse profile as solid line with error band
                self.ax2.plot(R_prof_sorted, param_data_sorted, '-', color=color,
                             linewidth=1.5, label=label)
                self.ax2.fill_between(R_prof_sorted,
                                      param_data_sorted - param_err_sorted,
                                      param_data_sorted + param_err_sorted,
                                      color=color, alpha=0.2)

                # Interpolate param values at TGAMMA R positions
                # Use R_raw from TGAMMA (good channels only)
                param_interp_func = interp1d(R_prof_sorted, param_data_sorted,
                                             fill_value='extrapolate', bounds_error=False)
                param_at_tgamma = param_interp_func(R_raw)

                # Plot markers at TGAMMA positions
                self.ax2.plot(R_raw, param_at_tgamma, 'o', color=color,
                             markersize=5, markeredgecolor='white', markeredgewidth=0.5)

                # Add channel labels at TGAMMA positions (node: TGAMMA01~25)
                self._add_channel_labels(self.ax2, R_raw, param_at_tgamma, 'TGAMMA', raw_channels)

                param_max = max(param_max, np.nanmax(param_data))
                param_min = min(param_min, np.nanmin(param_data))

            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")

        # Set limits
        self.ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
        gamma_margin = (gamma_max - gamma_min) * 0.1
        self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)

        if param == 'q':
            self.ax2.axhline(1, color='gray', linestyle='--', alpha=0.5)
        else:
            self.ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        self.ax2.set_ylim(0, 6)

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def plot_efit_profiles(self):
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

        param = self.selected_param.get()

        gamma_max, gamma_min, param_max, param_min = 0, 0, 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)

        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, _ = self._parse_entry(entry)

                if entry not in self.efit_data:
                    continue

                color = colors[i]
                efit_entry = self.efit_data[entry]
                efit_R = efit_entry['R']

                if x_axis == "psi_N":
                    interp_func = interp1d(efit_R, efit_entry['psi_N'], fill_value='extrapolate')
                elif x_axis == "rho_pol":
                    interp_func = interp1d(efit_R, efit_entry['rho_pol'], fill_value='extrapolate')
                else:
                    interp_func = interp1d(efit_R, efit_entry['rho_tor'], fill_value='extrapolate')

                cache_key = f'{shot_number}_MSE'
                data = self.data[cache_key]

                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))

                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                # TGAMMA
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                x_raw = interp_func(R_raw)
                label = f'#{shot_number} {int(time_point*1e3):06d}ms (MSE)'

                self.ax1.errorbar(x_raw, tgamma, yerr=sgamma*0.5,
                                 fmt='o-', capsize=3, markersize=5, color=color,
                                 linewidth=0.8, label=label)

                # Add channel labels for TGAMMA
                raw_channels = [j+1 for j in range(len(good_mask)) if good_mask[j]]
                self._add_channel_labels(self.ax1, x_raw, tgamma, 'TGAMMA', raw_channels)

                gamma_min = min(gamma_min, np.nanpercentile(tgamma, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma, 98))

                # j or q
                # - Solid line: pmse profile (20 points)
                # - Markers: interpolated values at TGAMMA positions
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]

                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)
                x_prof = interp_func(R_prof)

                # Sort by x coordinate for proper line plotting
                sort_idx = np.argsort(x_prof)
                x_prof_sorted = x_prof[sort_idx]
                param_data_sorted = param_data[sort_idx]
                param_err_sorted = param_err[sort_idx]

                # Plot pmse profile as solid line with error band
                self.ax2.plot(x_prof_sorted, param_data_sorted, '-', color=color,
                             linewidth=1.5, label=label)
                self.ax2.fill_between(x_prof_sorted,
                                      param_data_sorted - param_err_sorted,
                                      param_data_sorted + param_err_sorted,
                                      color=color, alpha=0.2)

                # Interpolate param values at TGAMMA positions (in EFIT coordinates)
                param_interp_func = interp1d(x_prof_sorted, param_data_sorted,
                                             fill_value='extrapolate', bounds_error=False)
                param_at_tgamma = param_interp_func(x_raw)

                # Plot markers at TGAMMA positions
                self.ax2.plot(x_raw, param_at_tgamma, 'o', color=color,
                             markersize=5, markeredgecolor='white', markeredgewidth=0.5)

                # Add channel labels at TGAMMA positions
                self._add_channel_labels(self.ax2, x_raw, param_at_tgamma, 'TGAMMA', raw_channels)

                param_max = max(param_max, np.nanmax(param_data))
                param_min = min(param_min, np.nanmin(param_data))

            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")

        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(r'$\gamma$ [rad]')

        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')

        self.ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
        gamma_margin = (gamma_max - gamma_min) * 0.1
        self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)

        if param == 'q':
            self.ax2.axhline(1, color='gray', linestyle='--', alpha=0.5)
        else:
            self.ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        self.ax2.set_ylim(0, 6)

        for ax in [self.ax1, self.ax2]:
            ax.axvline(x=0, c='k', ls='--')
            ax.axvline(x=1, c='k', ls='--')

        self.ax1.set_xlim(-0.5, 1.05)
        self.ax2.set_xlim(0, 1.05)

        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
        self.canvas.draw()

        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()

    def _write_data_to_file(self, file_path, selected_entries):
        """Write MSE profile data to text file (gamma, q, j at TGAMMA positions)"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# MSE Profile Data (q and j interpolated at TGAMMA positions)\n")
            if self.computed_efit_tree:
                f.write(f"# EFIT Source: {self.computed_efit_tree}\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "gamma", "gamma_err", "drr[m]", "q", "q_err", "j[MA/m2]", "j_err"
            ))

            for entry in selected_entries:
                shot_number, time_point, _ = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'

                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]

                idx_raw = np.argmin(np.abs(data.time - time_point))
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))
                actual_time = data.time[idx_raw]

                # TGAMMA data
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']

                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]

                # Get q and j profiles for interpolation
                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]

                q_meas = data.measurements['q']
                j_meas = data.measurements['j']
                roa = q_meas['roa'][:, idx_prof]
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)

                # Sort by R for interpolation
                sort_idx = np.argsort(R_prof)
                R_prof_sorted = R_prof[sort_idx]

                q_data_sorted = q_meas['data'][:, idx_prof][sort_idx]
                q_err_sorted = q_meas['error'][:, idx_prof][sort_idx]
                j_data_sorted = j_meas['data'][:, idx_prof][sort_idx]
                j_err_sorted = j_meas['error'][:, idx_prof][sort_idx]

                # Interpolate q and j at TGAMMA R positions
                q_at_tgamma = np.interp(R_raw, R_prof_sorted, q_data_sorted)
                q_err_at_tgamma = np.interp(R_raw, R_prof_sorted, q_err_sorted)
                j_at_tgamma = np.interp(R_raw, R_prof_sorted, j_data_sorted)
                j_err_at_tgamma = np.interp(R_raw, R_prof_sorted, j_err_sorted)

                # Get EFIT values
                psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_raw)

                # Write data rows (gamma, q, j at same R positions)
                for i in range(len(R_raw)):
                    f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, actual_time, R_raw[i],
                        self._format_value(psi_n[i]),
                        self._format_value(rho_pol[i]),
                        self._format_value(rho_tor[i]),
                        tgamma[i], sgamma[i], drr[i],
                        q_at_tgamma[i], q_err_at_tgamma[i],
                        j_at_tgamma[i], j_err_at_tgamma[i]
                    ))
