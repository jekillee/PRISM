#!/usr/bin/python3.8

"""
ne, Te Time Trace tab with unified Thomson/ECE/TCI loading
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ui.base_tab import BaseTab
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT
)
from plotting.plot_manager import apply_legend_with_limit, TIMETRACE_LEGEND_LIMIT


class NeTeTimeTraceTab(BaseTab):
    """ne, Te Time Trace tab with unified Thomson/ECE/TCI loading"""
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.ece_loader = None
        self.ece_data_cache = {}
        self.tci_loader = None
        self.tci_data_cache = {}
    
    def create_widgets(self):
        """Create ne, Te time trace tab widgets"""
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
        
        self._create_shot_input(control_frame)
        self._create_selection_listboxes(control_frame)
        self._create_plot_controls(control_frame)
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
        
        diag_options = ['TS', 'ECE (100Hz)', 'ECE (1kHz)', 'TCI (100Hz)', 'TCI (1kHz)']
        self.selected_diagnostic = tk.StringVar(value='TS')
        
        diag_dropdown = ttk.Combobox(frame, textvariable=self.selected_diagnostic,
                                    values=diag_options, state="readonly", width=15)
        diag_dropdown.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.fetch_button = ttk.Button(frame, text='Fetch', command=self.load_shot_data, width=8)
        self.fetch_button.grid(row=0, column=3, padx=PAD_X, pady=PAD_Y, sticky='e')
    
    def _create_plot_controls(self, parent):
        """Create plot control buttons"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        ttk.Button(frame, text='Plot time traces', command=self.plot_data).pack(
            fill='x', padx=PAD_X, pady=PAD_Y)
    
    def _get_ece_loader(self):
        """On-demand initialization of ECE loader"""
        if self.ece_loader is None:
            from data_loaders.ece_loader import ECELoader
            from config.diagnostic_config import DIAGNOSTICS
            self.ece_loader = ECELoader(self.app_config, DIAGNOSTICS['ECE'])
        return self.ece_loader
    
    def _get_tci_loader(self):
        """On-demand initialization of TCI loader"""
        if self.tci_loader is None:
            from data_loaders.tci_loader import TCILoader
            from config.diagnostic_config import DIAGNOSTICS
            self.tci_loader = TCILoader(self.app_config, DIAGNOSTICS['TCI'])
        return self.tci_loader
    
    def load_shot_data(self):
        """Load Thomson/ECE/TCI data based on selection"""
        try:
            shot_number = int(self.shot_entry.get())
            selection = self.selected_diagnostic.get()
            
            all_channels = []
            
            # Parse selection to determine mode
            load_thomson = False
            load_ece = False
            load_tci = False
            sampling_rate = None
            sampling_key = None
            
            if selection.startswith('ECE'):
                load_ece = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'
            elif selection == 'TS':
                load_thomson = True
            elif selection.startswith('TCI'):
                load_tci = True
                if '100Hz' in selection:
                    sampling_rate = 0.01
                    sampling_key = '100Hz'
                else:  # 1kHz
                    sampling_rate = 0.001
                    sampling_key = '1kHz'
            
            if load_thomson:
                # Load Thomson data
                try:
                    ts_data = self.data_loader.load_data(shot_number)
                    cache_key = f'{shot_number}_TS'
                    self.data[cache_key] = ts_data
                    
                    # Thomson Core: TSC01-14, Edge: TSE01-17
                    for i, radius in enumerate(ts_data.radius):
                        if i < 14:
                            ch_label = f'TSC{i+1:02d}'
                        else:
                            ch_label = f'TSE{i-14+1:02d}'
                        
                        all_channels.append({
                            'R': radius,
                            'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})',
                            'source': 'TS',
                            'channel_idx': i,
                            'ch_label': ch_label
                        })
                    
                    print(f"Thomson data loaded: {len(ts_data.radius)} channels")
                except Exception as e:
                    print(f"Thomson not available: {str(e)}")
            
            if load_ece:
                # Load ECE data
                try:
                    loader = self._get_ece_loader()
                    ece_data = loader.load_data(shot_number, sampling_rate=sampling_rate)
                    
                    cache_key = f'{shot_number}_ECE_{sampling_key}'
                    self.ece_data_cache[cache_key] = ece_data
                    
                    channels = ece_data.measurements['Te'].get('channels', 
                               list(range(1, len(ece_data.radius) + 1)))
                    
                    for i, radius in enumerate(ece_data.radius):
                        ch_num = channels[i] if i < len(channels) else i + 1
                        ch_label = f'ECE{ch_num:02d}'
                        all_channels.append({
                            'R': radius,
                            'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label}) [{sampling_key}]',
                            'source': 'ECE',
                            'channel_idx': i,
                            'ch_label': ch_label,
                            'sampling_key': sampling_key
                        })
                    
                    n_valid = np.sum(ece_data.measurements['Te']['valid_mask'])
                    n_overlap = np.sum(ece_data.measurements['Te']['overlap_mask'])
                    print(f"ECE data loaded: {len(ece_data.radius)} channels @ {sampling_key}")
                    print(f"  Valid: {n_valid}, Overlap: {n_overlap}")
                except Exception as e:
                    print(f"ECE not available: {str(e)}")
            
            if load_tci:
                # Load TCI data with resampling
                try:
                    loader = self._get_tci_loader()
                    tci_data = loader.load_data(shot_number, sampling_rate=sampling_rate)
                    cache_key = f'{shot_number}_TCI_{sampling_key}'
                    self.tci_data_cache[cache_key] = tci_data
                    
                    channels = tci_data.measurements['ne']['channels']
                    for i, ch in enumerate(channels):
                        ch_label = f'TCI{ch:02d}'
                        all_channels.append({
                            'R': ch,  # Use channel number for sorting
                            'label': f'{shot_number:06d} ({ch_label}) [{sampling_key}]',
                            'source': 'TCI',
                            'channel_idx': i,
                            'ch_label': ch_label,
                            'sampling_key': sampling_key
                        })
                    
                    print(f"TCI data loaded: {len(channels)} channels @ {sampling_key}")
                except Exception as e:
                    print(f"TCI not available: {str(e)}")
            
            if not all_channels:
                messagebox.showerror("Error", "No data available for this shot")
                return
            
            # Sort by R position (or channel number for TCI)
            all_channels.sort(key=lambda x: x['R'])
            
            # Update available listbox
            self.available_listbox.delete(0, tk.END)
            for channel in all_channels:
                self.available_listbox.insert(tk.END, channel['label'])
            
            print(f"Total channels loaded: {len(all_channels)}")
            
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load data: {str(e)}")
    
    def _parse_entry(self, entry):
        """Parse time trace entry with sampling rate info
        
        Format: XXXXXX_YYYY (DIAG##) [rate] or XXXXXX (TCI##) [rate]
        Returns: shot_number, radius_or_channel, source, channel_label, sampling_key
        """
        # Extract sampling key if present
        sampling_key = None
        if '[' in entry and ']' in entry:
            sampling_key = entry.split('[')[1].split(']')[0].strip()
            entry_base = entry.split('[')[0].strip()
        else:
            entry_base = entry
        
        # Format: 012345_1800 (ECE01) or 012345 (TCI01)
        main_part = entry_base.split('(')[0].strip()
        diag_label = entry_base.split('(')[1].split(')')[0].strip()
        
        parts = main_part.split('_')
        shot_number = int(parts[0])
        
        # Determine source from diag_label
        if diag_label.startswith('TSC') or diag_label.startswith('TSE'):
            source = 'TS'
            radius = float(parts[1]) / 1e3  # mm to m
            return shot_number, radius, source, diag_label, sampling_key
        elif diag_label.startswith('ECE'):
            source = 'ECE'
            radius = float(parts[1]) / 1e3  # mm to m
            return shot_number, radius, source, diag_label, sampling_key
        elif diag_label.startswith('TCI'):
            source = 'TCI'
            channel = int(diag_label.replace('TCI', ''))
            return shot_number, channel, source, diag_label, sampling_key
        else:
            # Fallback
            radius = float(parts[1]) / 1e3 if len(parts) > 1 and parts[1] else 0
            return shot_number, radius, 'TS', diag_label, sampling_key
    
    def _parse_entry_for_efit(self, entry):
        """Not used for time trace tabs"""
        return None, None
    
    def plot_data(self):
        """Plot time traces with markers only (no lines)"""
        self.ax1.clear()
        self.ax2.clear()
        
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_xlabel('Time [s]')
        self.ax2.set_ylabel(self.param2['label'])
        
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return
        
        te_max, ne_max = 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)
        
        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius_or_ch, source, ch_label, sampling_key = self._parse_entry(entry)
                color = colors[i]
                
                if source == 'TS':
                    # Plot Thomson time trace
                    cache_key = f'{shot_number}_TS'
                    
                    if cache_key not in self.data:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]
                    
                    radius_idx = np.argmin(np.abs(data.radius - radius_or_ch))
                    actual_R = data.radius[radius_idx]
                    
                    x_data = data.time
                    
                    Te_data, Te_err = data.get_parameter('Te')
                    ne_data, ne_err = data.get_parameter('ne')
                    
                    Te_trace = Te_data[radius_idx, :]
                    Te_err_trace = Te_err[radius_idx, :]
                    ne_trace = ne_data[radius_idx, :]
                    ne_err_trace = ne_err[radius_idx, :]
                    
                    te_max = max(te_max, np.nanpercentile(Te_trace, 98))
                    ne_max = max(ne_max, np.nanpercentile(ne_trace, 98))
                    
                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'
                    
                    # TS: filled circle marker with error bars
                    self.ax1.errorbar(x_data, Te_trace, Te_err_trace,
                                     fmt='o', capsize=3, color=color, 
                                     markersize=3, label=label, zorder=10)
                    self.ax2.errorbar(x_data, ne_trace, ne_err_trace,
                                     fmt='o', capsize=3, color=color, 
                                     markersize=3, label=label, zorder=10)
                
                elif source == 'ECE':
                    # Plot ECE time trace with specific sampling rate
                    cache_key = f'{shot_number}_ECE_{sampling_key}'
                    
                    if cache_key not in self.ece_data_cache:
                        # Load with specified sampling rate
                        loader = self._get_ece_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        ece_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]
                    
                    radius_idx = np.argmin(np.abs(ece_data.radius - radius_or_ch))
                    actual_R = ece_data.radius[radius_idx]
                    
                    x_data = ece_data.time
                    Te_data = ece_data.measurements['Te']['data']
                    Te_trace = Te_data[radius_idx, :]
                    
                    valid_mask = ece_data.measurements['Te']['valid_mask']
                    is_valid = valid_mask[radius_idx]
                    
                    te_max = max(te_max, np.nanmax(Te_trace))
                    
                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label}) [{sampling_key}]'
                    if not is_valid:
                        label += ' ovlp'
                    
                    # ECE valid: empty square, ECE overlap: x marker
                    if is_valid:
                        self.ax1.plot(x_data, Te_trace, 's', color=color, 
                                     markersize=3, linestyle='',
                                     markerfacecolor='none', markeredgewidth=1.5,
                                     label=label, zorder=5)
                    else:
                        self.ax1.plot(x_data, Te_trace, 'x', color=color, 
                                     markersize=3, linestyle='',
                                     alpha=0.5, markeredgewidth=1.5,
                                     label=label, zorder=5)
                
                elif source == 'TCI':
                    # Plot TCI time trace with specific sampling rate
                    cache_key = f'{shot_number}_TCI_{sampling_key}'
                    
                    if cache_key not in self.tci_data_cache:
                        # Load with specified sampling rate
                        loader = self._get_tci_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        tci_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.tci_data_cache[cache_key] = tci_data
                    else:
                        tci_data = self.tci_data_cache[cache_key]
                    
                    channel = int(radius_or_ch)  # For TCI, this is channel number
                    channels = tci_data.measurements['ne']['channels']
                    
                    if channel in channels:
                        ch_idx = channels.index(channel)
                        x_data = tci_data.time
                        ne_data = tci_data.measurements['ne']['data']
                        ne_trace = ne_data[ch_idx, :]
                        
                        ne_max = max(ne_max, np.nanmax(ne_trace))
                        
                        # TCI has no R position
                        label = f'#{shot_number} ({ch_label}) [{sampling_key}]'
                        
                        # TCI: solid line (no markers)
                        self.ax2.plot(x_data, ne_trace, '-', color=color, 
                                     linewidth=1.5, label=label, zorder=5)
                
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")
        
        # Set limits
        if te_max > 0:
            te_margin = te_max * 0.1
            self.ax1.set_ylim(0, te_max + te_margin)
        
        if ne_max > 0:
            ne_margin = ne_max * 0.1
            self.ax2.set_ylim(0, ne_max + ne_margin)
        
        # Apply styling
        for ax in [self.ax1, self.ax2]:
            apply_legend_with_limit(ax, TIMETRACE_LEGEND_LIMIT,
                                    frameon=False, fontsize=8)
        
        self.plot_manager.apply_common_styling(self.ax1, self.ax2, skip_legend=True)
        self.canvas.draw()
        
        if self.toolbar:
            self.toolbar.update()
            self.toolbar.push_current()
    
    def plot_efit_profiles(self):
        """Not applicable for time trace tabs"""
        messagebox.showinfo("Info", "EFIT mapping not available for time trace tabs")
    
    def _write_data_to_file(self, file_path, selected_entries):
        """Write ne, Te time trace data to text file"""
        with open(file_path, 'w') as f:
            # Write header with aligned columns (10 char width each)
            f.write("# ne, Te Time Trace Data\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "Te[keV]", "Te_err", "ne[e19]", "ne_err", "Source"
            ))
            
            for entry in selected_entries:
                shot_number, radius, source, ch_label, sampling_key = self._parse_entry(entry)
                
                if source == 'TS':
                    cache_key = f'{shot_number}_TS'
                    if cache_key not in self.data:
                        data = self.data_loader.load_data(shot_number)
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]
                    
                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]
                    
                    Te_data, Te_err = data.get_parameter('Te')
                    ne_data, ne_err = data.get_parameter('ne')
                    
                    Te_trace = Te_data[radius_idx, :]
                    Te_err_trace = Te_err[radius_idx, :]
                    ne_trace = ne_data[radius_idx, :]
                    ne_err_trace = ne_err[radius_idx, :]
                    
                    # Write Thomson data rows (10 char width each)
                    for i in range(len(data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, data.time[i], actual_R,
                            Te_trace[i], Te_err_trace[i],
                            ne_trace[i], ne_err_trace[i], 'TS'
                        ))
                
                elif source == 'ECE':
                    cache_key = f'{shot_number}_ECE_{sampling_key}'
                    if cache_key not in self.ece_data_cache:
                        loader = self._get_ece_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        ece_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.ece_data_cache[cache_key] = ece_data
                    else:
                        ece_data = self.ece_data_cache[cache_key]
                    
                    radius_idx = np.argmin(np.abs(ece_data.radius - radius))
                    actual_R = ece_data.radius[radius_idx]
                    
                    Te_data = ece_data.measurements['Te']['data']
                    Te_trace = Te_data[radius_idx, :]
                    
                    # Write ECE data rows (no ne for ECE, 10 char width each)
                    for i in range(len(ece_data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10s,%10s,%10s,%10s\n" % (
                            shot_number, ece_data.time[i], actual_R,
                            Te_trace[i], 'NaN', 'NaN', 'NaN', 'ECE'
                        ))
                
                elif source == 'TCI':
                    cache_key = f'{shot_number}_TCI_{sampling_key}'
                    if cache_key not in self.tci_data_cache:
                        loader = self._get_tci_loader()
                        sr = 0.01 if sampling_key == '100Hz' else 0.001
                        tci_data = loader.load_data(shot_number, sampling_rate=sr)
                        self.tci_data_cache[cache_key] = tci_data
                    else:
                        tci_data = self.tci_data_cache[cache_key]
                    
                    channel = int(radius)  # For TCI, radius is channel number
                    channels = tci_data.measurements['ne']['channels']
                    
                    if channel in channels:
                        ch_idx = channels.index(channel)
                        ne_data = tci_data.measurements['ne']['data']
                        ne_trace = ne_data[ch_idx, :]
                        
                        # Write TCI data rows (no Te, no R for TCI)
                        for i in range(len(tci_data.time)):
                            f.write(" %10d,%10.3f,%10s,%10s,%10s,%10.3f,%10s,%10s\n" % (
                                shot_number, tci_data.time[i], 'NaN',
                                'NaN', 'NaN', ne_trace[i], 'NaN', 'TCI'
                            ))