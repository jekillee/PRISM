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
from ui.base_tab import BaseTab
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y, ENTRY_WIDTH_AXIS,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT, LABEL_WIDTH_MEDIUM
)


class MSEProfileTab(BaseTab):
    """MSE Profile tab with TGAMMA and j/q profiles"""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def create_widgets(self):
        """Create MSE profile tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)
        
        # Setup 1x2 plot with shared x-axis: TGAMMA (left), j or q (right)
        self.ax1 = self.figure.add_subplot(121)
        self.ax2 = self.figure.add_subplot(122, sharex=self.ax1)
        
        self.ax1.set_xlabel('R [m]')
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('R [m]')
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill='both', expand=True)
        
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)
        
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_axis_controls(control_frame)
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
    
    def _create_axis_controls(self, parent):
        """Create custom axis control panel for MSE (gamma + q/j)"""
        frame = ttk.LabelFrame(parent, text="Axis Control Panel", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        for i in range(4):
            frame.grid_columnconfigure(i, weight=1)
        
        self.axis_entries = {}
        
        # Headers
        ttk.Label(frame, text="", anchor="center").grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky="ew")
        ttk.Label(frame, text="x", anchor="center").grid(
            row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky="ew")
        ttk.Label(frame, text="gamma [rad]", anchor="center").grid(
            row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky="ew")
        ttk.Label(frame, text="q or j [MA/m2]", anchor="center").grid(
            row=0, column=3, padx=PAD_X, pady=PAD_Y, sticky="ew")
        
        # Min entries
        ttk.Label(frame, text="min", anchor="center").grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y)
        
        x_min_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        x_min_entry.grid(row=1, column=1, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['xmin'] = x_min_entry
        
        y1_min_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        y1_min_entry.grid(row=1, column=2, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['y1min'] = y1_min_entry
        
        y2_min_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        y2_min_entry.grid(row=1, column=3, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['y2min'] = y2_min_entry
        
        # Max entries
        ttk.Label(frame, text="max", anchor="center").grid(
            row=2, column=0, padx=PAD_X, pady=PAD_Y)
        
        x_max_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        x_max_entry.grid(row=2, column=1, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['xmax'] = x_max_entry
        
        y1_max_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        y1_max_entry.grid(row=2, column=2, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['y1max'] = y1_max_entry
        
        y2_max_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
        y2_max_entry.grid(row=2, column=3, padx=PAD_X, pady=PAD_Y, sticky="ew")
        self.axis_entries['y2max'] = y2_max_entry
        
        # Apply button
        ttk.Button(frame, text="Apply", command=self.apply_axis_limits).grid(
            row=0, column=4, rowspan=3, padx=PAD_X, pady=PAD_Y, sticky="nsew")
    
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
                # =============================================================
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]
                
                # Convert r/a to R
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)
                
                self.ax2.errorbar(R_prof, param_data, yerr=param_err,
                                 fmt='o-', capsize=3, markersize=5, color=color,
                                 label=label)
                
                # Add channel labels for profile (node: pmse_qv01~20 or pmse_jv01~20)
                n_prof = len(R_prof)
                prof_channels = list(range(1, n_prof + 1))
                node_prefix = 'pmse_qv' if param == 'q' else 'pmse_jv'
                self._add_channel_labels(self.ax2, R_prof, param_data, node_prefix, prof_channels)
                
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
                
                gamma_min = min(gamma_min, np.nanpercentile(tgamma, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma, 98))
                
                # j or q
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]
                
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)
                x_prof = interp_func(R_prof)
                
                self.ax2.errorbar(x_prof, param_data, yerr=param_err,
                                 fmt='o-', capsize=3, markersize=5, color=color,
                                 label=label)
                
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
        """Write MSE profile data to text file"""
        param = self.selected_param.get()
        
        with open(file_path, 'w') as f:
            # Write header for TGAMMA data (10 char width each)
            f.write("# MSE Profile Data - Section 1: TGAMMA (raw measurements)\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "gamma", "gamma_err", "drr[m]"
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
                actual_time = data.time[idx_raw]
                
                # TGAMMA data
                tgamma_meas = data.measurements['tgamma']
                good_mask = tgamma_meas['good_mask']
                
                R_raw = tgamma_meas['R'][good_mask]
                drr = tgamma_meas['drr'][good_mask]
                tgamma = tgamma_meas['data'][good_mask, idx_raw]
                sgamma = tgamma_meas['error'][good_mask, idx_raw]
                
                # Get EFIT values
                psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_raw)
                
                # Write TGAMMA data rows (10 char width each)
                for i in range(len(R_raw)):
                    f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, actual_time, R_raw[i],
                        self._format_value(psi_n[i]),
                        self._format_value(rho_pol[i]),
                        self._format_value(rho_tor[i]),
                        tgamma[i], sgamma[i], drr[i]
                    ))
            
            # Write header for q/j profile data (10 char width each)
            f.write("#\n")
            f.write("# MSE Profile Data - Section 2: %s profile (processed)\n" % param)
            if param == 'q':
                f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                    "Shot", "Time[s]", "r/a", "R[m]", "psi_N", "rho_pol", "rho_tor", "q", "q_err"
                ))
            else:
                f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                    "Shot", "Time[s]", "r/a", "R[m]", "psi_N", "rho_pol", "rho_tor", "j[MA/m2]", "j_err"
                ))
            
            for entry in selected_entries:
                shot_number, time_point, _ = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'
                data = self.data[cache_key]
                
                idx_prof = np.argmin(np.abs(data.time_prof - time_point))
                actual_time = data.time_prof[idx_prof]
                
                r_edge = self._get_r_edge_at_time(data, time_point)
                magx = data.measurements['magx']['data'][idx_prof]
                
                param_meas = data.measurements[param]
                roa = param_meas['roa'][:, idx_prof]
                param_data = param_meas['data'][:, idx_prof]
                param_err = param_meas['error'][:, idx_prof]
                
                R_prof = self.data_loader.get_R_from_roa(roa, magx, r_edge)
                
                # Get EFIT values for profile R positions
                psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_prof)
                
                # Write profile data rows (10 char width each)
                for i in range(len(roa)):
                    f.write(" %10d,%10.3f,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f\n" % (
                        shot_number, actual_time, roa[i], R_prof[i],
                        self._format_value(psi_n[i]),
                        self._format_value(rho_pol[i]),
                        self._format_value(rho_tor[i]),
                        param_data[i], param_err[i]
                    ))