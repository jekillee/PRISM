#!/usr/bin/python3.8

"""
MSE Time Trace tab with j/q selection
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
    ENTRY_WIDTH_SHOT, LABEL_WIDTH_SHORT, LABEL_WIDTH_MEDIUM
)
from plotting.plot_manager import apply_legend_with_limit, TIMETRACE_LEGEND_LIMIT


class MSETimeTraceTab(BaseTab):
    """MSE Time Trace tab with TGAMMA and j/q time traces"""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    
    def create_widgets(self):
        """Create MSE time trace tab widgets"""
        self.figure = Figure(self.app_config.FIGURE_SIZE, tight_layout=True)
        
        # Setup 2x1 plot: TGAMMA (top), j or q (bottom)
        self.ax1, self.ax2 = self.plot_manager.setup_timetrace_plot(
            self.figure, r'$\gamma$ [rad]', 'q')
        
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
        row_frame.grid_columnconfigure(1, weight=1)
        
        self.selected_param = tk.StringVar(value='q')
        param_dropdown = ttk.Combobox(row_frame, textvariable=self.selected_param,
                                      values=['q', 'j'], state='readonly', width=3)
        param_dropdown.grid(row=0, column=0, sticky='w')
        
        ttk.Button(row_frame, text='Plot time traces', command=self.plot_data).grid(
            row=0, column=1, sticky='ew', padx=(10, 0))
    
    def _get_r_edge_at_time(self, data, time_point):
        """Get R_edge at specific time by interpolation"""
        r_edge_info = data.measurements['r_edge']
        interp_func = interp1d(r_edge_info['time'], r_edge_info['data'], 
                               fill_value='extrapolate')
        return float(interp_func(time_point))
    
    def load_shot_data(self):
        """Load MSE shot data from MDS+"""
        try:
            shot_number = int(self.shot_entry.get())
            
            data = self.data_loader.load_data(shot_number)
            
            cache_key = f'{shot_number}_MSE'
            self.data[cache_key] = data
            
            self.available_listbox.delete(0, tk.END)
            
            # Get good channels for TGAMMA
            tgamma_meas = data.measurements['tgamma']
            good_mask = tgamma_meas['good_mask']
            R_raw = tgamma_meas['R']
            raw_channels = tgamma_meas.get('raw_channels', list(range(1, 26)))
            
            # Add TGAMMA channels (good channels only)
            for ch_idx in np.where(good_mask)[0]:
                R_val = R_raw[ch_idx]
                ch_num = raw_channels[ch_idx] if ch_idx < len(raw_channels) else ch_idx + 1
                item_str = f'{shot_number:06d}_{R_val*1e3:.0f} (MSE{ch_num:02d})'
                self.available_listbox.insert(tk.END, item_str)
            
            print(f"MSE data loaded: {np.sum(good_mask)} good channels")
            print(f"  NB source: {data.nb_source}")
            
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data: {str(e)}")
    
    def _parse_entry(self, entry):
        """Parse time trace entry: XXXXXX_YYYY (MSE##)"""
        # Format: 036053_1800 (MSE01)
        main_part = entry.split('(')[0].strip()
        ch_label = entry.split('(')[1].split(')')[0].strip()
        parts = main_part.split('_')
        
        shot_number = int(parts[0])
        radius = float(parts[1]) / 1e3  # mm to m
        
        return shot_number, radius, ch_label
    
    def _parse_entry_for_efit(self, entry):
        """Not used for time trace tabs"""
        return None, None
    
    def plot_data(self):
        """Plot time traces"""
        self.ax1.clear()
        self.ax2.clear()
        
        self.ax1.set_ylabel(r'$\gamma$ [rad]')
        self.ax2.set_xlabel('Time [s]')
        
        param = self.selected_param.get()
        if param == 'q':
            self.ax2.set_ylabel('q')
        else:
            self.ax2.set_ylabel('j [MA/m$^2$]')
        
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return
        
        gamma_max, param_max, param_min = 0, 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)
        
        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius, ch_label = self._parse_entry(entry)
                color = colors[i]
                
                cache_key = f'{shot_number}_MSE'
                
                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]
                
                # Find channel index by R position
                tgamma_meas = data.measurements['tgamma']
                R_raw = tgamma_meas['R']
                ch_idx = np.argmin(np.abs(R_raw - radius))
                actual_R = R_raw[ch_idx]
                
                # Plot TGAMMA time trace with fill_between
                tgamma_trace = tgamma_meas['data'][ch_idx, :]
                sgamma_trace = tgamma_meas['error'][ch_idx, :]
                
                label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'
                
                self.ax1.plot(data.time, tgamma_trace, '-', color=color,
                             linewidth=0.8, label=label)
                self.ax1.fill_between(data.time, 
                                      tgamma_trace - sgamma_trace*0.5,
                                      tgamma_trace + sgamma_trace*0.5,
                                      color=color, alpha=0.3)
                
                gamma_min = min(gamma_min, np.nanpercentile(tgamma_trace, 2))
                gamma_max = max(gamma_max, np.nanpercentile(tgamma_trace, 98))
                
                # Plot j or q time trace interpolated at gamma R position
                param_meas = data.measurements[param]
                roa_data = param_meas['roa']
                param_data_all = param_meas['data']
                param_err_all = param_meas['error']
                
                # Interpolate at each time step
                param_trace = np.zeros(len(data.time_prof))
                param_err_trace = np.zeros(len(data.time_prof))
                
                for t_idx in range(len(data.time_prof)):
                    magx = data.measurements['magx']['data'][t_idx]
                    r_edge = self._get_r_edge_at_time(data, data.time_prof[t_idx])
                    
                    # Convert r/a to R for this time step
                    roa_profile = roa_data[:, t_idx]
                    R_profile = magx + roa_profile * (r_edge - magx)
                    
                    # Interpolate param at actual_R
                    param_profile = param_data_all[:, t_idx]
                    param_err_profile = param_err_all[:, t_idx]
                    
                    # Sort by R for interpolation
                    sort_idx = np.argsort(R_profile)
                    R_sorted = R_profile[sort_idx]
                    param_sorted = param_profile[sort_idx]
                    param_err_sorted = param_err_profile[sort_idx]
                    
                    param_trace[t_idx] = np.interp(actual_R, R_sorted, param_sorted)
                    param_err_trace[t_idx] = np.interp(actual_R, R_sorted, param_err_sorted)
                
                label2 = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'
                
                self.ax2.errorbar(data.time_prof, param_trace, yerr=param_err_trace,
                                 fmt='o-', capsize=3, markersize=3, color=color,
                                 linewidth=0.8, label=label2)
                
                param_max = max(param_max, np.nanmax(param_trace))
                param_min = min(param_min, np.nanmin(param_trace))
                
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")
        
        # Set limits and styling
        self.ax1.axhline(0, color='gray', linestyle='--', alpha=0.5)
        if gamma_max > 0:
            gamma_margin = (gamma_max - gamma_min) * 0.1
            self.ax1.set_ylim(gamma_min - gamma_margin, gamma_max + gamma_margin)
        
        if param == 'q':
            self.ax2.axhline(1, color='gray', linestyle='--', alpha=0.5)
        else:
            self.ax2.axhline(0, color='gray', linestyle='--', alpha=0.5)
        self.ax2.set_ylim(0, 6)
        
        # Apply styling
        for ax in [self.ax1, self.ax2]:
            apply_legend_with_limit(ax, TIMETRACE_LEGEND_LIMIT,
                                    frameon=False, fontsize=8)
        
        self.plot_manager.apply_common_styling(self.ax1, self.ax2, skip_legend=True)
        self.canvas.draw()
        
        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()
    
    def _create_axis_controls(self, parent):
        """Create axis control panel with dynamic y2 label for MSE"""
        frame = ttk.LabelFrame(parent, text="Axis Control Panel", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        for i in range(4):
            frame.grid_columnconfigure(i, weight=1)
        
        self.axis_entries = {}
        
        # Headers
        ttk.Label(frame, text="", anchor="center").grid(row=0, column=0, padx=5, pady=5, sticky="ew")
        ttk.Label(frame, text='Time [s]', anchor="center").grid(row=0, column=1, padx=5, pady=5, sticky="ew")
        ttk.Label(frame, text=r'gamma[rad]', anchor="center").grid(row=0, column=2, padx=5, pady=5, sticky="ew")
        ttk.Label(frame, text='q or j [MA/m2]', anchor="center").grid(row=0, column=3, padx=5, pady=5, sticky="ew")
        
        # Min entries
        ttk.Label(frame, text="min", anchor="center").grid(row=1, column=0, padx=5, pady=5)
        
        x_min_entry = ttk.Entry(frame, width=6)
        x_min_entry.grid(row=1, column=1, padx=5, pady=5, sticky="ew")
        self.axis_entries['xmin'] = x_min_entry
        
        y1_min_entry = ttk.Entry(frame, width=6)
        y1_min_entry.grid(row=1, column=2, padx=5, pady=5, sticky="ew")
        self.axis_entries['y1min'] = y1_min_entry
        
        y2_min_entry = ttk.Entry(frame, width=6)
        y2_min_entry.grid(row=1, column=3, padx=5, pady=5, sticky="ew")
        self.axis_entries['y2min'] = y2_min_entry
        
        # Max entries
        ttk.Label(frame, text="max", anchor="center").grid(row=2, column=0, padx=5, pady=5)
        
        x_max_entry = ttk.Entry(frame, width=6)
        x_max_entry.grid(row=2, column=1, padx=5, pady=5, sticky="ew")
        self.axis_entries['xmax'] = x_max_entry
        
        y1_max_entry = ttk.Entry(frame, width=6)
        y1_max_entry.grid(row=2, column=2, padx=5, pady=5, sticky="ew")
        self.axis_entries['y1max'] = y1_max_entry
        
        y2_max_entry = ttk.Entry(frame, width=6)
        y2_max_entry.grid(row=2, column=3, padx=5, pady=5, sticky="ew")
        self.axis_entries['y2max'] = y2_max_entry
        
        # Apply button
        ttk.Button(frame, text="Apply", command=self.apply_axis_limits).grid(
            row=0, column=4, rowspan=3, padx=5, pady=5, sticky="nsew")
    
    def plot_efit_profiles(self):
        """Not applicable for time trace tabs"""
        messagebox.showinfo("Info", "EFIT mapping not available for time trace tabs")
    
    def _write_data_to_file(self, file_path, selected_entries):
        """Write MSE time trace data to text file"""
        param = self.selected_param.get()
        
        with open(file_path, 'w') as f:
            # Write header for TGAMMA time trace (10 char width each)
            f.write("# MSE Time Trace Data - Section 1: TGAMMA time trace\n")
            f.write("#%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "gamma", "gamma_err"
            ))
            
            for entry in selected_entries:
                shot_number, radius, ch_label = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'
                
                if cache_key not in self.data:
                    data = self.data_loader.load_data(shot_number)
                    self.data[cache_key] = data
                else:
                    data = self.data[cache_key]
                
                # Find channel index by R position
                tgamma_meas = data.measurements['tgamma']
                R_raw = tgamma_meas['R']
                ch_idx = np.argmin(np.abs(R_raw - radius))
                actual_R = R_raw[ch_idx]
                
                tgamma_trace = tgamma_meas['data'][ch_idx, :]
                sgamma_trace = tgamma_meas['error'][ch_idx, :]
                
                # Write TGAMMA time trace data (10 char width each)
                for i in range(len(data.time)):
                    f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, data.time[i], actual_R,
                        tgamma_trace[i], sgamma_trace[i]
                    ))
            
            # Write header for q/j time trace (10 char width each)
            f.write("#\n")
            f.write("# MSE Time Trace Data - Section 2: %s time trace (interpolated at gamma R)\n" % param)
            if param == 'q':
                f.write("#%10s,%10s,%10s,%10s,%10s\n" % (
                    "Shot", "Time[s]", "R[m]", "q", "q_err"
                ))
            else:
                f.write("#%10s,%10s,%10s,%10s,%10s\n" % (
                    "Shot", "Time[s]", "R[m]", "j[MA/m2]", "j_err"
                ))
            
            for entry in selected_entries:
                shot_number, radius, ch_label = self._parse_entry(entry)
                cache_key = f'{shot_number}_MSE'
                data = self.data[cache_key]
                
                # Find actual gamma R position
                tgamma_meas = data.measurements['tgamma']
                R_raw = tgamma_meas['R']
                ch_idx = np.argmin(np.abs(R_raw - radius))
                actual_R = R_raw[ch_idx]
                
                param_meas = data.measurements[param]
                roa_data = param_meas['roa']
                param_data_all = param_meas['data']
                param_err_all = param_meas['error']
                
                # Interpolate at each time step
                for i in range(len(data.time_prof)):
                    magx = data.measurements['magx']['data'][i]
                    r_edge = self._get_r_edge_at_time(data, data.time_prof[i])
                    
                    # Convert r/a to R for this time step
                    roa_profile = roa_data[:, i]
                    R_profile = magx + roa_profile * (r_edge - magx)
                    
                    # Sort by R for interpolation
                    param_profile = param_data_all[:, i]
                    param_err_profile = param_err_all[:, i]
                    sort_idx = np.argsort(R_profile)
                    R_sorted = R_profile[sort_idx]
                    param_sorted = param_profile[sort_idx]
                    param_err_sorted = param_err_profile[sort_idx]
                    
                    param_val = np.interp(actual_R, R_sorted, param_sorted)
                    param_err_val = np.interp(actual_R, R_sorted, param_err_sorted)
                    
                    f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f\n" % (
                        shot_number, data.time_prof[i], actual_R,
                        param_val, param_err_val
                    ))