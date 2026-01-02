#!/usr/bin/python3.8

"""
Abstract base class for all diagnostic tabs
Contains only common functionality shared across all diagnostics
"""

from abc import ABC, abstractmethod
import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.text import Annotation

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y, PAD_SECTION,
    ENTRY_WIDTH_SHOT, ENTRY_WIDTH_AXIS, LISTBOX_HEIGHT, LISTBOX_WIDTH,
    BUTTON_WIDTH_SMALL, BUTTON_WIDTH_MEDIUM, BUTTON_WIDTH_ARROW,
    LABEL_WIDTH_SHORT, LABEL_WIDTH_MEDIUM, LABEL_WIDTH_LONG,
    SCROLLBAR_STYLE, SCROLLBAR_WIDTH
)


class BaseTab(ABC):
    """Abstract base class for all diagnostic tabs"""
    
    def __init__(self, parent, app_config, diagnostic_name, tab_type, 
                 data_loader, efit_loader, plot_manager):
        self.parent = parent
        self.app_config = app_config
        self.diagnostic_name = diagnostic_name
        self.tab_type = tab_type
        self.data_loader = data_loader
        self.efit_loader = efit_loader
        self.plot_manager = plot_manager
        
        self.frame = ttk.Frame(parent)
        self.data = {}
        self.efit_data = {}
        self.computed_efit_tree = None
        self.toolbar = None
        
        from config.diagnostic_config import DIAGNOSTICS
        self.diag_config = DIAGNOSTICS[diagnostic_name]
        
        self.param1 = self.diag_config['parameters'][0]
        if len(self.diag_config['parameters']) > 1:
            self.param2 = self.diag_config['parameters'][1]
        else:
            self.param2 = None
    
    @abstractmethod
    def create_widgets(self):
        """Create tab widgets - must be implemented by each diagnostic tab"""
        pass
    
    @abstractmethod
    def load_shot_data(self):
        """Load shot data - must be implemented by each diagnostic tab"""
        pass
    
    @abstractmethod
    def plot_data(self):
        """Plot data - must be implemented by each diagnostic tab"""
        pass
    
    # ===== Common widget creation methods =====
    
    def _create_selection_listboxes(self, parent):
        """Create selection listboxes (common for all tabs)"""
        outer_frame = ttk.LabelFrame(parent, text="2. Select Data", labelanchor="n")
        outer_frame.pack(fill='both', expand=True, padx=PAD_X, pady=PAD_Y)
        
        style = ttk.Style()
        style.configure(SCROLLBAR_STYLE, width=SCROLLBAR_WIDTH)
        
        frame = ttk.Frame(outer_frame)
        frame.pack(fill='both', expand=True, padx=PAD_SECTION, pady=PAD_SECTION)
        
        frame.grid_rowconfigure(1, weight=1)
        frame.grid_columnconfigure(0, weight=1)
        frame.grid_columnconfigure(1, weight=0)
        frame.grid_columnconfigure(2, weight=1)
        
        # Label
        if self.tab_type == 'profile':
            label_text = 'Time [ms] (Diag)'
        else:
            label_text = 'R [mm] (Ch)'
        ttk.Label(frame, text=label_text).grid(row=0, column=0, padx=2, pady=2)
        
        # Available listbox with vertical and horizontal scrollbars
        available_frame = ttk.Frame(frame)
        available_frame.grid(row=1, column=0, padx=2, pady=2, sticky='nsew')
        
        # Horizontal scrollbar at bottom (pack first)
        available_xscroll = ttk.Scrollbar(available_frame, orient='horizontal')
        available_xscroll.pack(side='bottom', fill='x')
        
        # Vertical scrollbar on right
        available_yscroll = ttk.Scrollbar(available_frame, orient='vertical',
                                          style=SCROLLBAR_STYLE)
        available_yscroll.pack(side='right', fill='y')
        
        # Listbox
        self.available_listbox = tk.Listbox(available_frame, height=LISTBOX_HEIGHT, 
                                           width=LISTBOX_WIDTH, selectmode='extended',
                                           xscrollcommand=available_xscroll.set,
                                           yscrollcommand=available_yscroll.set)
        self.available_listbox.pack(side='left', fill='both', expand=True)
        
        available_xscroll.config(command=self.available_listbox.xview)
        available_yscroll.config(command=self.available_listbox.yview)
        
        # Buttons
        button_frame = ttk.Frame(frame)
        button_frame.grid(row=1, column=1, padx=2, pady=2, sticky='ns')
        
        ttk.Button(button_frame, text='>', width=BUTTON_WIDTH_ARROW, 
                  command=self.add_selected_items).pack(fill='both', expand=True, pady=(0, 2))
        ttk.Button(button_frame, text='<', width=BUTTON_WIDTH_ARROW, 
                  command=self.remove_selected_items).pack(fill='both', expand=True)
        
        # Selected listbox with vertical and horizontal scrollbars
        ttk.Label(frame, text='Selected').grid(row=0, column=2, padx=2, pady=2)
        
        selected_frame = ttk.Frame(frame)
        selected_frame.grid(row=1, column=2, padx=2, pady=2, sticky='nsew')
        
        # Horizontal scrollbar at bottom (pack first)
        selected_xscroll = ttk.Scrollbar(selected_frame, orient='horizontal')
        selected_xscroll.pack(side='bottom', fill='x')
        
        # Vertical scrollbar on right
        selected_yscroll = ttk.Scrollbar(selected_frame, orient='vertical',
                                         style=SCROLLBAR_STYLE)
        selected_yscroll.pack(side='right', fill='y')
        
        # Listbox
        self.selected_listbox = tk.Listbox(selected_frame, height=LISTBOX_HEIGHT, 
                                          width=LISTBOX_WIDTH, selectmode='extended',
                                          xscrollcommand=selected_xscroll.set,
                                          yscrollcommand=selected_yscroll.set)
        self.selected_listbox.pack(side='left', fill='both', expand=True)
        
        selected_xscroll.config(command=self.selected_listbox.xview)
        selected_yscroll.config(command=self.selected_listbox.yview)
        
        # Keyboard shortcuts for listboxes
        self.selected_listbox.bind('<Delete>', lambda e: self.remove_selected_items())
        self.selected_listbox.bind('<BackSpace>', lambda e: self.remove_selected_items())
    
    def _create_efit_controls(self, parent):
        """Create EFIT mapping controls (common for profile tabs)"""
        self.selected_x_axis = tk.StringVar(value="R")
        
        frame = ttk.LabelFrame(parent, text="4. EFIT Mapping (Optional)", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        for i in range(4):
            frame.grid_columnconfigure(i, weight=1)
        
        ttk.Label(frame, text="EFIT Tree:").grid(row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky="e")
        
        efit_display_values = list(self.app_config.EFIT_TREES.keys())
        self.selected_efit_tree_display = tk.StringVar(value=efit_display_values[0])
        
        efit_dropdown = ttk.Combobox(frame, textvariable=self.selected_efit_tree_display,
                                    values=efit_display_values,
                                    state="readonly", width=20)
        efit_dropdown.grid(row=0, column=1, columnspan=2, padx=PAD_X, pady=PAD_Y, sticky="ew")
        
        ttk.Button(frame, text='Mapping', command=self.compute_efit).grid(
            row=0, column=3, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        ttk.Radiobutton(frame, text="psi_N", variable=self.selected_x_axis,
                       value="psi_N").grid(row=1, column=0, padx=PAD_X, pady=PAD_Y)
        ttk.Radiobutton(frame, text="rho_pol", variable=self.selected_x_axis,
                       value="rho_pol").grid(row=1, column=1, padx=PAD_X, pady=PAD_Y)
        ttk.Radiobutton(frame, text="rho_tor", variable=self.selected_x_axis,
                       value="rho_tor").grid(row=1, column=2, padx=PAD_X, pady=PAD_Y)
        ttk.Button(frame, text="Plot", command=self.plot_efit_profiles).grid(
            row=1, column=3, padx=PAD_X, pady=PAD_Y, sticky='ew')
    
    def _create_axis_controls(self, parent):
        """Create axis control panel (common for all tabs)"""
        frame = ttk.LabelFrame(parent, text="Axis Control Panel", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        for i in range(4):
            frame.grid_columnconfigure(i, weight=1)
        
        self.axis_entries = {}
        
        # Headers
        ttk.Label(frame, text="", anchor="center").grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky="ew")
        
        if self.tab_type == 'profile':
            x_label = 'x'
        else:
            x_label = 'Time [s]'
        ttk.Label(frame, text=x_label, anchor="center").grid(
            row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky="ew")
        
        ttk.Label(frame, text=self.param1['short_label'], anchor="center").grid(
            row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky="ew")
        if self.param2:
            ttk.Label(frame, text=self.param2['short_label'], anchor="center").grid(
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
        
        if self.param2:
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
        
        if self.param2:
            y2_max_entry = ttk.Entry(frame, width=ENTRY_WIDTH_AXIS)
            y2_max_entry.grid(row=2, column=3, padx=PAD_X, pady=PAD_Y, sticky="ew")
            self.axis_entries['y2max'] = y2_max_entry
        
        # Apply button
        ttk.Button(frame, text="Apply", command=self.apply_axis_limits).grid(
            row=0, column=4, rowspan=3, padx=PAD_X, pady=PAD_Y, sticky="nsew")
    
    def _create_save_controls(self, parent):
        """Create save controls (common for all tabs)"""
        ttk.Button(parent, text='Save as .txt', command=self.save_data).pack(
            fill='x', padx=PAD_X, pady=2)
    
    # ===== Loading state helper methods =====
    
    def _set_loading_state(self, is_loading, button=None, original_text='Fetch'):
        """Set loading state with visual feedback"""
        if button is None and hasattr(self, 'fetch_button'):
            button = self.fetch_button
        
        if button is None:
            return
        
        if is_loading:
            button.config(text='Loading...', state='disabled')
            self.frame.config(cursor='watch')
            self.frame.update_idletasks()
        else:
            button.config(text=original_text, state='normal')
            self.frame.config(cursor='')
    
    # ===== Channel label helper methods =====
    
    def _add_channel_labels(self, ax, x_data, y_data, node_prefix, channels):
        """Add channel labels above data points if show_channel is enabled
        
        node_prefix: str like 'CES_TI', 'CESNN_TI', 'TS_CORE', 'TS_EDGE', 'ECE', 
                     'TGAMMA', 'pmse_qv', 'pmse_jv', 'TXCS_TI'
        channels: list of channel numbers (int)
        """
        if not hasattr(self, 'show_channel_var') or not self.show_channel_var.get():
            return
        
        for x, y, ch in zip(x_data, y_data, channels):
            label = f'{node_prefix}{ch:02d}'
            ax.annotate(label, (x, y), textcoords='offset points', 
                       xytext=(0, 5), ha='center', fontsize=7, alpha=0.8,
                       clip_on=True)
    
    def _clear_channel_labels(self, ax):
        """Clear existing channel label annotations"""
        for child in list(ax.get_children()):
            if isinstance(child, Annotation):
                child.remove()
    
    # ===== Common action methods =====
    
    def add_selected_items(self):
        """Add selected items to the selected list"""
        for index in self.available_listbox.curselection():
            item = self.available_listbox.get(index)
            if item not in self.selected_listbox.get(0, tk.END):
                self.selected_listbox.insert(tk.END, item)
    
    def remove_selected_items(self):
        """Remove selected items from the selected list"""
        for index in reversed(self.selected_listbox.curselection()):
            self.selected_listbox.delete(index)
    
    def apply_axis_limits(self):
        """Apply custom axis limits"""
        try:
            limits = {}
            for name, entry in self.axis_entries.items():
                value = entry.get()
                limits[name] = float(value) if value else None
            
            # Validate min < max for x/t axis
            xmin_key = 'xmin' if 'xmin' in limits else 'tmin'
            xmax_key = 'xmax' if 'xmax' in limits else 'tmax'
            
            if limits.get(xmin_key) is not None and limits.get(xmax_key) is not None:
                if limits[xmin_key] >= limits[xmax_key]:
                    messagebox.showerror("Error", f"min must be less than max")
                    return
            
            if limits['y1min'] is not None and limits['y1max'] is not None:
                if limits['y1min'] >= limits['y1max']:
                    messagebox.showerror("Error", "y1_min must be less than y1_max")
                    return
            
            if self.param2:
                if limits.get('y2min') is not None and limits.get('y2max') is not None:
                    if limits['y2min'] >= limits['y2max']:
                        messagebox.showerror("Error", "y2_min must be less than y2_max")
                        return
            
            # Apply x/t axis limits
            if limits.get(xmin_key) is not None or limits.get(xmax_key) is not None:
                self.ax1.set_xlim(left=limits.get(xmin_key), right=limits.get(xmax_key))
                if self.param2:
                    self.ax2.set_xlim(left=limits.get(xmin_key), right=limits.get(xmax_key))
            
            # Apply y1 axis limits
            if limits.get('y1min') is not None or limits.get('y1max') is not None:
                self.ax1.set_ylim(bottom=limits.get('y1min'), top=limits.get('y1max'))
            
            # Apply y2 axis limits
            if self.param2:
                if limits.get('y2min') is not None or limits.get('y2max') is not None:
                    self.ax2.set_ylim(bottom=limits.get('y2min'), top=limits.get('y2max'))
            
            self.canvas.draw_idle()
            
        except ValueError:
            messagebox.showerror("Error", "Please enter valid numeric values")
    
    def compute_efit(self):
        """Compute EFIT equilibrium data (common logic)"""
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            messagebox.showwarning("Warning", "No timepoints selected")
            return
        
        display_name = self.selected_efit_tree_display.get()
        efit_tree = self.app_config.EFIT_TREES[display_name]
        
        self.efit_data.clear()
        
        print("\n" + "="*50)
        print(f"      EFIT Mapping ({efit_tree})")
        print("-" * 50)
        
        # Group by shot
        shot_groups = {}
        for entry in selected_entries:
            shot_number, param_value = self._parse_entry_for_efit(entry)
            if shot_number not in shot_groups:
                shot_groups[shot_number] = []
            shot_groups[shot_number].append((entry, param_value))
        
        valid_entries = []
        
        for shot_number, entry_list in shot_groups.items():
            try:
                print(f"-> Loading EFIT data for shot #{shot_number}...")
                efit_data = self.efit_loader.load_efit_data(shot_number, efit_tree)
                
                for entry, time_point in entry_list:
                    if not (np.min(efit_data.time) <= time_point <= np.max(efit_data.time)):
                        print(f"   -> Skipping {entry}: EFIT out of range")
                        continue
                    
                    time_idx = np.argmin(np.abs(efit_data.time - time_point))
                    closest_time = efit_data.time[time_idx]
                    
                    if abs(closest_time - time_point) > 0.05:
                        print(f"   -> Skipping {entry}: No EFIT within ÃƒÆ’Ã¢â‚¬Å¡Ãƒâ€šÃ‚Â±50ms")
                        continue
                    
                    valid_entries.append(entry)
                    
                    self.efit_data[entry] = {
                        'R': efit_data.radius,
                        'psi_N': efit_data.psi_n[time_idx],
                        'rho_pol': efit_data.rho_pol[time_idx],
                        'rho_tor': efit_data.rho_tor[time_idx]
                    }
                    
                    print(f"   -> {entry} processed...")
                
            except Exception as e:
                print(f"-> Error: {str(e)}")
        
        if len(valid_entries) != len(selected_entries):
            missing = set(selected_entries) - set(valid_entries)
            messagebox.showerror("Incomplete EFIT Data",
                               f"Missing EFIT data for:\n" + '\n'.join(['- ' + e for e in missing]))
            self.efit_data.clear()
            self.computed_efit_tree = None
            return
        
        self.computed_efit_tree = efit_tree
        print("      EFIT Mapping Completed!")
        print("=" * 50 + "\n")
    
    @abstractmethod
    def _parse_entry_for_efit(self, entry):
        """Parse entry to extract shot and time for EFIT - implemented by each tab"""
        pass
    
    @abstractmethod
    def plot_efit_profiles(self):
        """Plot profiles with EFIT mapping - implemented by each tab"""
        pass
    
    def save_data(self):
        """Save selected data to text file"""
        selected_entries = list(self.selected_listbox.get(0, tk.END))
        if not selected_entries:
            messagebox.showwarning("Warning", "No data selected to save")
            return
        
        # Get save file path
        file_path = filedialog.asksaveasfilename(
            initialdir=os.path.expanduser("~"),
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            title="Save Data"
        )
        
        if not file_path:
            return
        
        try:
            self._write_data_to_file(file_path, selected_entries)
            print(f"Data saved to {file_path}")
            messagebox.showinfo("Success", f"Data saved to {file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {str(e)}")
    
    @abstractmethod
    def _write_data_to_file(self, file_path, selected_entries):
        """Write data to file - implemented by each tab"""
        pass
    
    def _get_efit_values_at_R(self, entry, R_positions):
        """Get EFIT coordinate values at given R positions
        
        Returns tuple of (psi_N, rho_pol, rho_tor) arrays.
        Returns (None, None, None) arrays if EFIT not computed.
        """
        from scipy.interpolate import interp1d
        
        n_points = len(R_positions)
        
        # Check if EFIT is computed for this entry
        if not hasattr(self, 'efit_data') or not self.efit_data:
            return [None]*n_points, [None]*n_points, [None]*n_points
        
        if entry not in self.efit_data:
            return [None]*n_points, [None]*n_points, [None]*n_points
        
        efit_entry = self.efit_data[entry]
        efit_R = efit_entry['R']
        
        # Interpolate EFIT values at measurement R positions
        psi_n_interp = interp1d(efit_R, efit_entry['psi_N'], 
                                fill_value='extrapolate', bounds_error=False)
        rho_pol_interp = interp1d(efit_R, efit_entry['rho_pol'], 
                                  fill_value='extrapolate', bounds_error=False)
        rho_tor_interp = interp1d(efit_R, efit_entry['rho_tor'], 
                                  fill_value='extrapolate', bounds_error=False)
        
        psi_n_vals = psi_n_interp(R_positions)
        rho_pol_vals = rho_pol_interp(R_positions)
        rho_tor_vals = rho_tor_interp(R_positions)
        
        return psi_n_vals, rho_pol_vals, rho_tor_vals
    
    def _format_value(self, val, fmt="%10.3f"):
        """Format value, handling None"""
        if val is None:
            return "%10s" % "None"
        return fmt % val