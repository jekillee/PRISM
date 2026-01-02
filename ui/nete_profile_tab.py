#!/usr/bin/python3.8

"""
ne, Te Profile tab with Thomson/ECE selection
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ui.base_tab import BaseTab
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT
)


class NeTeProfileTab(BaseTab):
    """ne, Te Profile tab supporting Thomson and/or ECE"""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ece_loader = None
        self.ece_data_cache = {}
    
    def create_widgets(self):
        """Create ne, Te profile tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)
        
        self.ax1, self.ax2 = self.plot_manager.setup_profile_plot(
            self.figure, self.param1['label'], self.param2['label'])
        
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
        """Create data loading section with diagnostic selection"""
        frame = ttk.LabelFrame(parent, text="1. Load ne, Te Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Row 0: Shot label, entry, diagnostic dropdown, Fetch
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_SHORT, anchor='e').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='e')
        self.shot_entry = ttk.Entry(frame)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self.load_shot_data())
        
        diag_options = ['TS+ECE', 'TS only', 'ECE only (100Hz)', 'ECE only (1kHz)']
        self.selected_diagnostic = tk.StringVar(value='TS+ECE')
        
        diag_dropdown = ttk.Combobox(frame, textvariable=self.selected_diagnostic,
                                    values=diag_options, state="readonly", width=15)
        diag_dropdown.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.fetch_button = ttk.Button(frame, text='Fetch', command=self.load_shot_data, width=8)
        self.fetch_button.grid(row=0, column=3, padx=PAD_X, pady=PAD_Y, sticky='e')
    
    def _create_plot_controls(self, parent):
        """Create plot control buttons"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        row_frame = ttk.Frame(frame)
        row_frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        row_frame.grid_columnconfigure(1, weight=1)
        
        self.show_channel_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(row_frame, text='Show Nodes', variable=self.show_channel_var).grid(
            row=0, column=0, sticky='w')
        
        ttk.Button(row_frame, text='Plot R profiles', command=self.plot_data).grid(
            row=0, column=1, sticky='ew', padx=(10, 0))
    
    def _get_ece_loader(self):
        """On-demand initialization of ECE loader"""
        if self.ece_loader is None:
            from data_loaders.ece_loader import ECELoader
            from config.diagnostic_config import DIAGNOSTICS
            self.ece_loader = ECELoader(self.app_config, DIAGNOSTICS['ECE'])
        return self.ece_loader
    
    def load_shot_data(self):
        """Load shot data based on diagnostic selection"""
        try:
            shot_number = int(self.shot_entry.get())
            selection = self.selected_diagnostic.get()
            
            self.available_listbox.delete(0, tk.END)
            
            # Parse selection to determine mode and ECE sampling
            if selection.startswith('ECE only'):
                load_thomson = False
                load_ece = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'
            elif selection == 'TS only':
                load_thomson = True
                load_ece = False
                sampling_rate = None
                sampling_key = None
            else:  # TS+ECE
                load_thomson = True
                load_ece = True
                sampling_rate = 0.01
                sampling_key = '100Hz'
            
            if load_thomson:
                # Load Thomson data
                data = self.data_loader.load_data(shot_number)
                cache_key = f'{shot_number}_TS'
                self.data[cache_key] = data
                
                # Label depends on mode
                label = 'TS+ECE' if selection == 'TS+ECE' else 'TS'
                
                for tp in data.time:
                    item_str = f'{shot_number:06d}_{tp*1e3:06.0f} ({label})'
                    self.available_listbox.insert(tk.END, item_str)
                
                print(f"Thomson data loaded: {len(data.radius)} channels, {len(data.time)} timepoints")
            
            if load_ece:
                loader = self._get_ece_loader()
                ece_data = loader.load_data(shot_number, sampling_rate=sampling_rate)
                
                cache_key = f'{shot_number}_ECE_{sampling_key}'
                self.ece_data_cache[cache_key] = ece_data
                
                # Add ECE timepoints to listbox only for "ECE only" mode
                if selection.startswith('ECE only'):
                    for tp in ece_data.time:
                        item_str = f'{shot_number:06d}_{tp*1e3:06.0f} (ECE)'
                        self.available_listbox.insert(tk.END, item_str)
                
                n_valid = np.sum(ece_data.measurements['Te']['valid_mask'])
                n_overlap = np.sum(ece_data.measurements['Te']['overlap_mask'])
                print(f"ECE data loaded: {len(ece_data.radius)} channels, {len(ece_data.time)} timepoints @ {sampling_key}")
                print(f"  Valid: {n_valid}, Overlap: {n_overlap}")
                
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data: {str(e)}")
    
    def _parse_entry(self, entry):
        """Parse entry to get shot, time, and source"""
        if '(' in entry and ')' in entry:
            main_part = entry.split('(')[0].strip()
            source = entry.split('(')[1].split(')')[0].strip()
            parts = main_part.split('_')
        else:
            parts = entry.split('_')
            source = 'TS'
        
        shot_number = int(parts[0])
        time_point = float(parts[1]) / 1e3  # ms to s
        
        return shot_number, time_point, source
    
    def _parse_entry_for_efit(self, entry):
        """Parse entry for EFIT (extract shot and time)"""
        shot_number, time_point, _ = self._parse_entry(entry)
        return shot_number, time_point
    
    def _get_ece_cache_key(self, shot_number):
        """Find ECE cache key for given shot"""
        for key in self.ece_data_cache.keys():
            if key.startswith(f'{shot_number}_ECE'):
                return key
        return None
    
    def plot_data(self):
        """Plot R profiles"""
        self.ax1.clear()
        self.ax2.clear()
        
        self.ax1.set_xlabel('R [m]')
        self.ax2.set_xlabel('R [m]')
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])
        
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return
        
        te_max, ne_max = 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)
        
        has_thomson = False
        
        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, source = self._parse_entry(entry)
                color = colors[i]
                
                if source in ['TS', 'TS+ECE']:
                    # Plot Thomson data
                    has_thomson = True
                    cache_key = f'{shot_number}_TS'
                    data = self.data.get(cache_key) or self.data_loader.load_data(shot_number)
                    if cache_key not in self.data:
                        self.data[cache_key] = data
                    
                    time_idx = np.argmin(np.abs(data.time - time_point))
                    R_data = data.radius
                    
                    Te_data, Te_err = data.get_parameter('Te')
                    ne_data, ne_err = data.get_parameter('ne')
                    
                    Te_profile = Te_data[:, time_idx]
                    Te_err_profile = Te_err[:, time_idx]
                    ne_profile = ne_data[:, time_idx]
                    ne_err_profile = ne_err[:, time_idx]
                    
                    valid_idx = np.argmin(np.abs(R_data - self.app_config.R_EDGE))
                    te_max = max(te_max, np.nanmax(Te_profile[:valid_idx]))
                    ne_max = max(ne_max, np.nanmax(ne_profile[:valid_idx]))
                    
                    label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (TS)'
                    
                    self.ax1.errorbar(R_data, Te_profile, Te_err_profile,
                                     fmt='o', capsize=5, label=label,
                                     color=color, markersize=5, zorder=10)
                    
                    self.ax2.errorbar(R_data, ne_profile, ne_err_profile,
                                     fmt='o', capsize=5, label=label,
                                     color=color, markersize=5, zorder=10)
                    
                    # Add Thomson channel labels (TS_CORE1-14 for core, TS_EDGE1-17 for edge)
                    n_channels = len(R_data)
                    for ch_idx in range(n_channels):
                        if self.show_channel_var.get():
                            if ch_idx < 14:
                                lbl = f'TS_CORE{ch_idx + 1}'
                            else:
                                lbl = f'TS_EDGE{ch_idx - 14 + 1}'
                            self.ax1.annotate(lbl, (R_data[ch_idx], Te_profile[ch_idx]), 
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8)
                            self.ax2.annotate(lbl, (R_data[ch_idx], ne_profile[ch_idx]), 
                                             textcoords='offset points',
                                             xytext=(0, 5), ha='center', fontsize=7, alpha=0.8)
                    
                    # Plot ECE data only for "TS+ECE" mode
                    if source == 'TS+ECE':
                        ece_cache_key = self._get_ece_cache_key(shot_number)
                        if ece_cache_key is not None:
                            ece_data = self.ece_data_cache[ece_cache_key]
                            
                            ece_time_idx = np.argmin(np.abs(ece_data.time - time_point))
                            ece_R_data = ece_data.radius
                            
                            ece_Te_data = ece_data.measurements['Te']['data']
                            ece_Te_profile = ece_Te_data[:, ece_time_idx]
                            
                            valid_mask = ece_data.measurements['Te']['valid_mask']
                            
                            valid_in_range = valid_mask & (ece_R_data <= self.app_config.R_EDGE)
                            if np.any(valid_in_range):
                                te_max = max(te_max, np.nanmax(ece_Te_profile[valid_in_range]))
                            
                            ece_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (ECE)'
                            
                            # Plot ECE valid channels only (2nd harmonic) - empty square
                            self.ax1.plot(ece_R_data[valid_mask], ece_Te_profile[valid_mask],
                                         's', color=color, markersize=5,
                                         markerfacecolor='none', markeredgewidth=1.5,
                                         label=ece_label, zorder=5)
                            
                            # Add ECE channel labels (ECE01~76)
                            ece_channels = ece_data.measurements['Te'].get('channels', 
                                           list(range(1, len(ece_R_data) + 1)))
                            valid_R = ece_R_data[valid_mask]
                            valid_Te = ece_Te_profile[valid_mask]
                            valid_ch = [ece_channels[i] for i in range(len(valid_mask)) if valid_mask[i]]
                            self._add_channel_labels(self.ax1, valid_R, valid_Te, 'ECE', valid_ch)
                
                elif source == 'ECE':
                    # Plot ECE data (valid channels only)
                    cache_key = self._get_ece_cache_key(shot_number)
                    
                    if cache_key is None:
                        # Load with default sampling rate if not found
                        loader = self._get_ece_loader()
                        ece_data = loader.load_data(shot_number, sampling_rate=0.01)  # 100Hz default
                        cache_key = f'{shot_number}_ECE_100Hz'
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]
                    
                    time_idx = np.argmin(np.abs(ece_data.time - time_point))
                    R_data = ece_data.radius
                    
                    Te_data = ece_data.measurements['Te']['data']
                    Te_profile = Te_data[:, time_idx]
                    
                    valid_mask = ece_data.measurements['Te']['valid_mask']
                    
                    valid_in_range = valid_mask & (R_data <= self.app_config.R_EDGE)
                    if np.any(valid_in_range):
                        te_max = max(te_max, np.nanmax(Te_profile[valid_in_range]))
                    
                    label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (ECE)'
                    
                    # Plot valid channels only (2nd harmonic) - empty square
                    self.ax1.plot(R_data[valid_mask], Te_profile[valid_mask],
                                 's', color=color, markersize=5,
                                 markerfacecolor='none', markeredgewidth=1.5,
                                 label=label, zorder=5)
                    
                    # Add ECE channel labels
                    ece_channels = ece_data.measurements['Te'].get('channels',
                                   list(range(1, len(R_data) + 1)))
                    valid_R = R_data[valid_mask]
                    valid_Te = Te_profile[valid_mask]
                    valid_ch = [ece_channels[i] for i in range(len(valid_mask)) if valid_mask[i]]
                    self._add_channel_labels(self.ax1, valid_R, valid_Te, 'ECE', valid_ch)
                
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")
        
        # Set limits
        self.ax1.set_xlim(self.app_config.R_LIMITS)
        te_margin = te_max * 0.1
        self.ax1.set_ylim(0, te_max + te_margin)
        
        if has_thomson:
            self.ax2.set_xlim(self.app_config.R_LIMITS)
            ne_margin = ne_max * 0.1
            self.ax2.set_ylim(0, ne_max + ne_margin)
        
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
        
        te_max, ne_max = 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)
        
        has_thomson = False
        
        for i, entry in enumerate(selected_entries):
            try:
                shot_number, time_point, source = self._parse_entry(entry)
                
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
                
                if source in ['TS', 'TS+ECE']:
                    has_thomson = True
                    cache_key = f'{shot_number}_TS'
                    data = self.data[cache_key]
                    
                    time_idx = np.argmin(np.abs(data.time - time_point))
                    x_data = interp_func(data.radius)
                    
                    Te_data, Te_err = data.get_parameter('Te')
                    ne_data, ne_err = data.get_parameter('ne')
                    
                    Te_profile = Te_data[:, time_idx]
                    Te_err_profile = Te_err[:, time_idx]
                    ne_profile = ne_data[:, time_idx]
                    ne_err_profile = ne_err[:, time_idx]
                    
                    lcfs_idx = np.argmin(np.abs(x_data - 1))
                    te_max = max(te_max, np.nanmax(Te_profile[:lcfs_idx]))
                    ne_max = max(ne_max, np.nanmax(ne_profile[:lcfs_idx]))
                    
                    label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (TS)'
                    
                    self.ax1.errorbar(x_data, Te_profile, Te_err_profile,
                                     fmt='o', capsize=5, label=label,
                                     color=color, markersize=5, zorder=10)
                    
                    self.ax2.errorbar(x_data, ne_profile, ne_err_profile,
                                     fmt='o', capsize=5, label=label,
                                     color=color, markersize=5, zorder=10)
                    
                    # Plot ECE data only for "TS+ECE" mode
                    if source == 'TS+ECE':
                        ece_cache_key = self._get_ece_cache_key(shot_number)
                        if ece_cache_key is not None:
                            ece_data = self.ece_data_cache[ece_cache_key]
                            ece_time_idx = np.argmin(np.abs(ece_data.time - time_point))
                            ece_x_data = interp_func(ece_data.radius)
                            
                            ece_Te_data = ece_data.measurements['Te']['data']
                            ece_Te_profile = ece_Te_data[:, ece_time_idx]
                            
                            valid_mask = ece_data.measurements['Te']['valid_mask']
                            
                            in_range = (ece_x_data <= 1.0)
                            if np.any(valid_mask & in_range):
                                te_max = max(te_max, np.nanmax(ece_Te_profile[valid_mask & in_range]))
                            
                            ece_label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (ECE)'
                            
                            # Plot ECE valid channels only (2nd harmonic) - empty square
                            self.ax1.plot(ece_x_data[valid_mask], ece_Te_profile[valid_mask],
                                         's', color=color, markersize=5,
                                         markerfacecolor='none', markeredgewidth=1.5,
                                         label=ece_label, zorder=5)
                
                elif source == 'ECE':
                    cache_key = self._get_ece_cache_key(shot_number)
                    
                    if cache_key:
                        ece_data = self.ece_data_cache[cache_key]
                        time_idx = np.argmin(np.abs(ece_data.time - time_point))
                        x_data = interp_func(ece_data.radius)
                        
                        Te_data = ece_data.measurements['Te']['data']
                        Te_profile = Te_data[:, time_idx]
                        
                        valid_mask = ece_data.measurements['Te']['valid_mask']
                        
                        in_range = (x_data <= 1.0)
                        if np.any(valid_mask & in_range):
                            te_max = max(te_max, np.nanmax(Te_profile[valid_mask & in_range]))
                        
                        label = f'#{shot_number} {entry.split("_")[1].split()[0]}ms (ECE)'
                        
                        # Plot valid channels only (2nd harmonic) - empty square
                        self.ax1.plot(x_data[valid_mask], Te_profile[valid_mask],
                                     's', color=color, markersize=5,
                                     markerfacecolor='none', markeredgewidth=1.5,
                                     label=label, zorder=5)
                
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")
        
        self.ax1.set_xlabel(x_label)
        self.ax2.set_xlabel(x_label)
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_ylabel(self.param2['label'])
        
        self.ax1.set_xlim(0, 1.05)
        self.ax2.set_xlim(0, 1.05)
        self.ax1.set_ylim(0, te_max * 1.1)
        
        if has_thomson:
            self.ax2.set_ylim(0, ne_max * 1.1)
        
        for ax in [self.ax1, self.ax2]:
            ax.axvline(x=0, c='k', ls='--')
            ax.axvline(x=1, c='k', ls='--')
        
        self.plot_manager.apply_common_styling(self.ax1, self.ax2)
        self.canvas.draw()
        
        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()
    
    def _write_data_to_file(self, file_path, selected_entries):
        """Write ne, Te profile data to text file"""
        with open(file_path, 'w') as f:
            # Write header with aligned columns (10 char width each)
            f.write("# ne, Te Profile Data\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "psi_N", "rho_pol", "rho_tor",
                "Te[keV]", "Te_err", "ne[e19]", "ne_err", "Source"
            ))
            
            for entry in selected_entries:
                shot_number, time_point, source = self._parse_entry(entry)
                
                if source in ['TS', 'TS+ECE']:
                    cache_key = f'{shot_number}_TS'
                    data = self.data.get(cache_key)
                    if data is None:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data
                    
                    time_idx = np.argmin(np.abs(data.time - time_point))
                    actual_time = data.time[time_idx]
                    R_data = data.radius
                    
                    Te_data, Te_err = data.get_parameter('Te')
                    ne_data, ne_err = data.get_parameter('ne')
                    
                    Te_profile = Te_data[:, time_idx]
                    Te_err_profile = Te_err[:, time_idx]
                    ne_profile = ne_data[:, time_idx]
                    ne_err_profile = ne_err[:, time_idx]
                    
                    # Get EFIT values
                    psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_data)
                    
                    # Write Thomson data rows (10 char width each)
                    for i in range(len(R_data)):
                        f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, actual_time, R_data[i],
                            self._format_value(psi_n[i]),
                            self._format_value(rho_pol[i]),
                            self._format_value(rho_tor[i]),
                            Te_profile[i], Te_err_profile[i],
                            ne_profile[i], ne_err_profile[i], 'TS'
                        ))
                    
                    # Check if ECE data exists for "TS+ECE" mode
                    ece_cache_key = self._get_ece_cache_key(shot_number)
                    if ece_cache_key is not None:
                        ece_data = self.ece_data_cache[ece_cache_key]
                        ece_time_idx = np.argmin(np.abs(ece_data.time - time_point))
                        ece_R_data = ece_data.radius
                        
                        ece_Te_data = ece_data.measurements['Te']['data']
                        ece_Te_profile = ece_Te_data[:, ece_time_idx]
                        valid_mask = ece_data.measurements['Te']['valid_mask']
                        
                        # Get EFIT values for ECE
                        ece_psi_n, ece_rho_pol, ece_rho_tor = self._get_efit_values_at_R(entry, ece_R_data)
                        
                        # Write ECE data rows (valid channels only, 10 char width each)
                        for i in range(len(ece_R_data)):
                            if valid_mask[i]:
                                f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10s,%10s,%10s,%10s\n" % (
                                    shot_number, actual_time, ece_R_data[i],
                                    self._format_value(ece_psi_n[i]),
                                    self._format_value(ece_rho_pol[i]),
                                    self._format_value(ece_rho_tor[i]),
                                    ece_Te_profile[i], 'NaN', 'NaN', 'NaN', 'ECE'
                                ))
                
                elif source == 'ECE':
                    cache_key = self._get_ece_cache_key(shot_number)
                    if cache_key is None:
                        loader = self._get_ece_loader()
                        ece_data = loader.load_data(shot_number, sampling_rate=0.01)
                        cache_key = f'{shot_number}_ECE_100Hz'
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]
                    
                    time_idx = np.argmin(np.abs(ece_data.time - time_point))
                    actual_time = ece_data.time[time_idx]
                    R_data = ece_data.radius
                    
                    Te_data = ece_data.measurements['Te']['data']
                    Te_profile = Te_data[:, time_idx]
                    valid_mask = ece_data.measurements['Te']['valid_mask']
                    
                    # Get EFIT values
                    psi_n, rho_pol, rho_tor = self._get_efit_values_at_R(entry, R_data)
                    
                    # Write ECE data rows (valid channels only, 10 char width each)
                    for i in range(len(R_data)):
                        if valid_mask[i]:
                            f.write(" %10d,%10.3f,%10.3f,%s,%s,%s,%10.3f,%10s,%10s,%10s,%10s\n" % (
                                shot_number, actual_time, R_data[i],
                                self._format_value(psi_n[i]),
                                self._format_value(rho_pol[i]),
                                self._format_value(rho_tor[i]),
                                Te_profile[i], 'NaN', 'NaN', 'NaN', 'ECE'
                            ))