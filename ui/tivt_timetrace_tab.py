#!/usr/bin/python3.8

"""
Ti/vT Time Trace tab with CES and XICS integration
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ui.base_tab import BaseTab
from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT
)
from plotting.plot_manager import apply_legend_with_limit, TIMETRACE_LEGEND_LIMIT


class TiVTTimeTraceTab(BaseTab):
    """Ti/vT Time Trace tab with CES and XICS support"""
    
    def __init__(self, parent, app_config, diagnostic_name, tab_type,
                 data_loader, efit_loader, plot_manager, file_parser):
        super().__init__(parent, app_config, diagnostic_name, tab_type,
                        data_loader, efit_loader, plot_manager)
        self.file_parser = file_parser
        self.xics_loader = None
        self.xics_data_cache = {}
    
    def create_widgets(self):
        """Create Ti/vT time trace tab widgets"""
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
        """Create data loading section with analysis type selection"""
        frame = ttk.LabelFrame(parent, text="1. Load Ti/vT Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Row 0: Shot label, entry, dropdown, Fetch + File
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_SHORT, anchor='e').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='e')
        self.shot_entry = ttk.Entry(frame)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self.load_shot_data())
        
        analysis_types = list(self.diag_config['analysis_types'].keys())
        self.selected_analysis_type = tk.StringVar(value='mod')
        
        analysis_dropdown = ttk.Combobox(frame, textvariable=self.selected_analysis_type,
                                        values=analysis_types, state="readonly", width=15)
        analysis_dropdown.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=0, column=3, padx=PAD_X, pady=PAD_Y, sticky='e')
        
        self.fetch_button = ttk.Button(btn_frame, text='Fetch', command=self.load_shot_data, width=5)
        self.fetch_button.pack(side=tk.LEFT)
        
        ttk.Button(btn_frame, text='...', command=self.load_file_data, width=3).pack(
            side=tk.LEFT, padx=(2, 0))
    
    def _create_plot_controls(self, parent):
        """Create plot control buttons"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(frame, text='Plot time traces', command=self.plot_data).pack(
            fill='x', padx=5, pady=5)
    
    def _get_xics_loader(self):
        """Lazy initialization of XICS loader"""
        if self.xics_loader is None:
            from data_loaders.xics_loader import XICSLoader
            from config.diagnostic_config import DIAGNOSTICS
            self.xics_loader = XICSLoader(self.app_config, DIAGNOSTICS.get('XICS', {}))
        return self.xics_loader
    
    def load_shot_data(self):
        """Load CES and XICS shot data from MDS+, sorted by R"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return
        
        analysis_type = self.selected_analysis_type.get()
        ces_loaded = False
        xics_loaded = False
        data = None
        xics_data = None
        
        # Try to load CES data
        try:
            data = self.data_loader.load_data(shot_number, analysis_type)
            cache_key = f'{shot_number}_{analysis_type}'
            self.data[cache_key] = data
            ces_loaded = True
            print(f"CES {analysis_type} data loaded: {len(data.radius)} channels, {len(data.time)} timepoints")
        except Exception as e:
            print(f"CES {analysis_type} not available: {str(e)}")
        
        # Try to load XICS data
        try:
            loader = self._get_xics_loader()
            xics_data = loader.load_data(shot_number)
            if xics_data is not None:
                xics_cache_key = f'{shot_number}_XICS'
                self.xics_data_cache[xics_cache_key] = xics_data
                xics_loaded = True
                print(f"XICS data loaded: {len(xics_data.time)} timepoints")
        except Exception as e:
            print(f"XICS not available: {str(e)}")
        
        # Check if any data loaded
        if not ces_loaded and not xics_loaded:
            messagebox.showerror("Error", f"No CES or XICS data available for shot #{shot_number}")
            return
        
        # Build list of all channels with R positions
        all_channels = []
        
        # Add CES channels: mod01-32 or nn01-32
        if ces_loaded:
            for i, radius in enumerate(data.radius):
                ch_label = f'{analysis_type}{i+1:02d}'
                all_channels.append({
                    'R': radius,
                    'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})'
                })
        
        # Add XICS channel
        if xics_loaded:
            all_channels.append({
                'R': 1.8,
                'label': f'{shot_number:06d}_1800 (XICS)'
            })
        
        # Sort by R position
        all_channels.sort(key=lambda x: x['R'])
        
        # Update listbox
        self.available_listbox.delete(0, tk.END)
        for ch in all_channels:
            self.available_listbox.insert(tk.END, ch['label'])
    
    def load_file_data(self):
        """Load CES data from result file, sorted by R"""
        file_path = filedialog.askopenfilename(
            title="Select CES Result File",
            initialdir=self.app_config.CES_RESULT_PATH,
            filetypes=[("Text files", "*.txt"), ("All Files", "*.*")]
        )
        
        if not file_path:
            return
        
        try:
            data, shot_number = self.file_parser.parse_file(file_path)
            
            file_key = f'file_{shot_number}_{data.source}'
            self.data[file_key] = data
            
            print(f"CES result file loaded: {data.source}")
            
            # Try to load XICS data for this shot
            loader = self._get_xics_loader()
            xics_data = loader.load_data(shot_number)
            
            if xics_data is not None:
                xics_cache_key = f'{shot_number}_XICS'
                self.xics_data_cache[xics_cache_key] = xics_data
                print(f"XICS data loaded: {len(xics_data.time)} timepoints")
            
            # Build list of all channels with R positions
            all_channels = []
            
            # Add CES file channels
            for i, radius in enumerate(data.radius):
                ch_label = f'{data.source}{i+1:02d}'
                all_channels.append({
                    'R': radius,
                    'label': f'{shot_number:06d}_{radius*1e3:.0f} ({ch_label})'
                })
            
            # Add XICS channel
            if xics_data is not None:
                all_channels.append({
                    'R': 1.8,
                    'label': f'{shot_number:06d}_1800 (XICS)'
                })
            
            # Sort by R position
            all_channels.sort(key=lambda x: x['R'])
            
            # Update listbox
            self.available_listbox.delete(0, tk.END)
            for ch in all_channels:
                self.available_listbox.insert(tk.END, ch['label'])
            
        except ValueError as e:
            messagebox.showerror("Invalid File Format", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {str(e)}")
    
    def _parse_entry(self, entry):
        """Parse time trace entry: XXXXXX_YYYY (DIAG##)
        
        Returns: shot_number, radius, source, channel_label
        """
        # Format: 012345_1800 (mod01) or 012345_1800 (XICS)
        main_part = entry.split('(')[0].strip()
        diag_label = entry.split('(')[1].split(')')[0].strip()
        
        parts = main_part.split('_')
        shot_number = int(parts[0])
        radius = float(parts[1]) / 1e3  # mm to m
        
        # Determine source from diag_label
        if diag_label.startswith('mod'):
            source = 'mod'
        elif diag_label.startswith('nn'):
            source = 'nn'
        elif diag_label == 'XICS':
            source = 'XICS'
        else:
            # File data
            source = diag_label
        
        return shot_number, radius, source, diag_label
    
    def _parse_entry_for_efit(self, entry):
        """Not used for time trace tabs"""
        return None, None
    
    def _get_marker_style(self, source):
        """Get marker style based on data source"""
        if source == 'mod':
            # Filled circle
            return {'marker': 'o', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'full'}
        elif source == 'nn':
            # Empty circle
            return {'marker': 'o', 'markerfacecolor': 'none', 'markeredgewidth': 1.5, 'fillstyle': 'full'}
        elif source == 'XICS':
            # Filled square
            return {'marker': 's', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'full'}
        else:
            # File data: half-filled circle
            return {'marker': 'o', 'markerfacecolor': None, 'markeredgewidth': 1, 'fillstyle': 'right'}
    
    def plot_data(self):
        """Plot time traces with markers"""
        self.ax1.clear()
        self.ax2.clear()
        
        self.ax1.set_ylabel(self.param1['label'])
        self.ax2.set_xlabel('Time [s]')
        self.ax2.set_ylabel(self.param2['label'])
        
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            return
        
        ti_max, vt_max, vt_min = 0, 0, 0
        colors = self.plot_manager.color_manager.get_colors_for_entries(selected_entries)
        
        for i, entry in enumerate(selected_entries):
            try:
                shot_number, radius, source, ch_label = self._parse_entry(entry)
                color = colors[i]
                
                if source == 'XICS':
                    # Plot XICS time trace
                    xics_cache_key = f'{shot_number}_XICS'
                    
                    if xics_cache_key not in self.xics_data_cache:
                        continue
                    
                    xics_data = self.xics_data_cache[xics_cache_key]
                    
                    x_data = xics_data.time
                    Ti_data, Ti_err = xics_data.get_parameter('Ti')
                    vT_data, vT_err = xics_data.get_parameter('vT')
                    
                    Ti_trace = Ti_data[0, :]
                    vT_trace = vT_data[0, :]
                    
                    ti_max = max(ti_max, np.nanpercentile(Ti_trace, 98))
                    vt_min = min(vt_min, np.nanpercentile(vT_trace, 2))
                    vt_max = max(vt_max, np.nanpercentile(vT_trace, 98))
                    
                    label = f'#{shot_number} 1800mm ({ch_label})'
                    
                    style = self._get_marker_style('XICS')
                    
                    # XICS: closed square without error bars
                    self.ax1.plot(x_data, Ti_trace,
                                 marker=style['marker'], color=color,
                                 markersize=3, label=label,
                                 linestyle='none')
                    self.ax2.plot(x_data, vT_trace,
                                 marker=style['marker'], color=color,
                                 markersize=3, label=label,
                                 linestyle='none')
                
                else:
                    # Plot CES time trace
                    if source in ['mod', 'nn']:
                        cache_key = f'{shot_number}_{source}'
                    else:
                        cache_key = f'file_{shot_number}_{source}'
                    
                    # Load data if not cached
                    if cache_key not in self.data:
                        if source in ['mod', 'nn']:
                            data = self.data_loader.load_data(shot_number, source)
                        else:
                            continue
                        self.data[cache_key] = data
                    else:
                        data = self.data[cache_key]
                    
                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]
                    
                    x_data = data.time
                    
                    Ti_data, Ti_err = data.get_parameter('Ti')
                    vT_data, vT_err = data.get_parameter('vT')
                    
                    Ti_trace = Ti_data[radius_idx, :]
                    Ti_err_trace = Ti_err[radius_idx, :]
                    vT_trace = vT_data[radius_idx, :]
                    vT_err_trace = vT_err[radius_idx, :]
                    
                    ti_max = max(ti_max, np.nanpercentile(Ti_trace, 98))
                    vt_min = min(vt_min, np.nanpercentile(vT_trace, 2))
                    vt_max = max(vt_max, np.nanpercentile(vT_trace, 98))
                    
                    label = f'#{shot_number} {actual_R*1e3:.0f}mm ({ch_label})'
                    
                    style = self._get_marker_style(source)
                    
                    plot_kwargs = {
                        'fmt': style['marker'],
                        'capsize': 3,
                        'label': label,
                        'color': color,
                        'markersize': 3,
                        'fillstyle': style['fillstyle']
                    }
                    if style['markerfacecolor'] == 'none':
                        plot_kwargs['markerfacecolor'] = 'none'
                        plot_kwargs['markeredgewidth'] = style['markeredgewidth']
                    
                    self.ax1.errorbar(x_data, Ti_trace, Ti_err_trace, **plot_kwargs)
                    self.ax2.errorbar(x_data, vT_trace, vT_err_trace, **plot_kwargs)
                
            except Exception as e:
                print(f"Error plotting {entry}: {str(e)}")
        
        # Set limits
        ti_margin = ti_max * 0.1
        self.ax1.set_ylim(0, ti_max + ti_margin)
        
        vt_margin = (vt_max - vt_min) * 0.1
        self.ax2.set_ylim(vt_min - vt_margin, vt_max + vt_margin)
        self.ax2.axhline(y=0, c='silver', ls='--')
        
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
        """Write Ti/vT time trace data to text file"""
        with open(file_path, 'w') as f:
            # Write header
            f.write("# Ti/vT Time Trace Data\n")
            f.write("#%10s,%10s,%10s,%10s,%10s,%10s,%10s,%10s\n" % (
                "Shot", "Time[s]", "R[m]", "Ti[keV]", "Ti_err", "vT[km/s]", "vT_err", "Source"
            ))
            
            for entry in selected_entries:
                shot_number, radius, source, ch_label = self._parse_entry(entry)
                
                if source == 'XICS':
                    # Write XICS data
                    xics_cache_key = f'{shot_number}_XICS'
                    if xics_cache_key not in self.xics_data_cache:
                        continue
                    
                    xics_data = self.xics_data_cache[xics_cache_key]
                    
                    Ti_data, Ti_err = xics_data.get_parameter('Ti')
                    vT_data, vT_err = xics_data.get_parameter('vT')
                    
                    Ti_trace = Ti_data[0, :]
                    Ti_err_trace = Ti_err[0, :]
                    vT_trace = vT_data[0, :]
                    vT_err_trace = vT_err[0, :]
                    
                    for i in range(len(xics_data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, xics_data.time[i], 1.8,
                            Ti_trace[i], Ti_err_trace[i],
                            vT_trace[i], vT_err_trace[i], 'XICS'
                        ))
                
                else:
                    # Write CES data
                    if source in ['mod', 'nn']:
                        cache_key = f'{shot_number}_{source}'
                    else:
                        cache_key = f'file_{shot_number}_{source}'
                    
                    if cache_key not in self.data:
                        if source in ['mod', 'nn']:
                            data = self.data_loader.load_data(shot_number, source)
                            self.data[cache_key] = data
                        else:
                            continue
                    else:
                        data = self.data[cache_key]
                    
                    radius_idx = np.argmin(np.abs(data.radius - radius))
                    actual_R = data.radius[radius_idx]
                    
                    Ti_data, Ti_err = data.get_parameter('Ti')
                    vT_data, vT_err = data.get_parameter('vT')
                    
                    Ti_trace = Ti_data[radius_idx, :]
                    Ti_err_trace = Ti_err[radius_idx, :]
                    vT_trace = vT_data[radius_idx, :]
                    vT_err_trace = vT_err[radius_idx, :]
                    
                    for i in range(len(data.time)):
                        f.write(" %10d,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10.3f,%10s\n" % (
                            shot_number, data.time[i], actual_R,
                            Ti_trace[i], Ti_err_trace[i],
                            vT_trace[i], vT_err_trace[i], source
                        ))