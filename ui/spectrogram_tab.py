#!/usr/bin/python3.8

"""
Spectrogram tab for high-frequency signal analysis
Supports ECE, Mirnov, and BES diagnostics
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
from scipy.signal import spectrogram
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from MDSplus import Connection

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT
)
from data_loaders.ecei_loader import ECEILoader
from data_loaders.ece_loader import ECELoader
from data_loaders.bes_loader import BESLoader


class SpectrogramTab:
    """Spectrogram visualization tab"""
    
    # NFFT options
    NFFT_OPTIONS = ['256', '512', '1024', '2048', '4096']
    DEFAULT_NFFT = '1024'
    
    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        
        self.frame = ttk.Frame(parent)
        self.toolbar = None
        
        # ECEI loader
        self.ecei_loader = ECEILoader(app_config.MDS_IP)
        
        # ECE loader
        self.ece_loader = ECELoader(app_config, None)
        
        # BES loader
        self.bes_loader = BESLoader(app_config, None)
        
        # Data cache
        self.ece_info = {}      # {shot: {'channels': [], 'R': [], 'Z': [], 'I_TF': float}}
        self.bes_info = {}      # {shot: {'channels': [], 'R': [], 'Z': []}}
        self.ecei_info = {}     # {(shot, device): {'channels': [], 'R': [], 'Z': []}}
        self.signal_data = None
        self.spectrogram_data = None
        
        # Current loaded shot
        self.current_shot = None
        
        # Plot references
        self.im = None
        self.colorbar = None
        
        # Widget references for enable/disable
        self.signal_widgets = []
    
    def create_widgets(self):
        """Create spectrogram tab widgets"""
        # Single canvas for spectrogram
        self.figure = Figure((10, 6), tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill='both', expand=True)
        
        # Control panel
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)
        
        self._create_shot_input(control_frame)
        self._create_signal_selection(control_frame)
        self._create_spectrogram_params(control_frame)
        self._create_color_controls(control_frame)
        self._create_save_controls(control_frame)
        
        # Initially disable signal selection
        self._set_signal_selection_state('disabled')
    
    def _create_shot_input(self, parent):
        """Create shot input section"""
        frame = ttk.LabelFrame(parent, text="1. Load Shot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_SHORT, anchor='center').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry = ttk.Entry(frame, width=ENTRY_WIDTH_SHOT)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self._load_shot_info())
        
        self.fetch_button = ttk.Button(frame, text='Fetch', command=self._load_shot_info, 
                                       width=BUTTON_WIDTH_MEDIUM)
        self.fetch_button.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')
    
    def _create_signal_selection(self, parent):
        """Create signal selection section"""
        frame = ttk.LabelFrame(parent, text="2. Select Signal", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        self.signal_frame = frame
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Shot number display (row 0)
        ttk.Label(frame, text='Shot:', anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.shot_display_label = ttk.Label(frame, text='Not loaded', anchor='w', 
                                            foreground='gray')
        self.shot_display_label.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        # Diagnostic type dropdown (row 1)
        ttk.Label(frame, text='Diagnostic:', anchor='w').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.diag_options = ['Mirnov (Toroidal)', 'Mirnov (Poloidal)', 'ECE', 'BES', 'TCI', 
                             'ECEI-GT', 'ECEI-GR', 'ECEI-HT']
        self.selected_diag = tk.StringVar(value='Mirnov (Toroidal)')
        
        self.diag_dropdown = ttk.Combobox(frame, textvariable=self.selected_diag,
                                          values=self.diag_options, state="disabled", width=18)
        self.diag_dropdown.grid(row=1, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.diag_dropdown.bind('<<ComboboxSelected>>', self._on_diag_changed)
        self.signal_widgets.append(self.diag_dropdown)
        
        # Channel dropdown (row 2)
        ttk.Label(frame, text='Channel:', anchor='w').grid(
            row=2, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.selected_channel = tk.StringVar()
        self.channel_dropdown = ttk.Combobox(frame, textvariable=self.selected_channel,
                                              state="disabled", width=18)
        self.channel_dropdown.grid(row=2, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.signal_widgets.append(self.channel_dropdown)
    
    def _create_spectrogram_params(self, parent):
        """Create spectrogram parameter section"""
        frame = ttk.LabelFrame(parent, text="3. Spectrogram Parameters", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        frame.grid_columnconfigure(3, weight=1)
        
        # Time range
        ttk.Label(frame, text='Time [s]:', anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        time_frame = ttk.Frame(frame)
        time_frame.grid(row=0, column=1, columnspan=3, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        self.time_min_entry = ttk.Entry(time_frame, width=8)
        self.time_min_entry.pack(side=tk.LEFT, padx=(0, 2))
        self.time_min_entry.insert(0, '0')
        
        ttk.Label(time_frame, text='-').pack(side=tk.LEFT)
        
        self.time_max_entry = ttk.Entry(time_frame, width=8)
        self.time_max_entry.pack(side=tk.LEFT, padx=(2, 0))
        self.time_max_entry.insert(0, '10')
        
        # Frequency range
        ttk.Label(frame, text='Freq [kHz]:', anchor='w').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        freq_frame = ttk.Frame(frame)
        freq_frame.grid(row=1, column=1, columnspan=3, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        self.freq_min_entry = ttk.Entry(freq_frame, width=8)
        self.freq_min_entry.pack(side=tk.LEFT, padx=(0, 2))
        self.freq_min_entry.insert(0, '0')
        
        ttk.Label(freq_frame, text='-').pack(side=tk.LEFT)
        
        self.freq_max_entry = ttk.Entry(freq_frame, width=8)
        self.freq_max_entry.pack(side=tk.LEFT, padx=(2, 0))
        self.freq_max_entry.insert(0, '100')
        
        # NFFT
        ttk.Label(frame, text='NFFT:', anchor='w').grid(
            row=2, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.selected_nfft = tk.StringVar(value=self.DEFAULT_NFFT)
        nfft_dropdown = ttk.Combobox(frame, textvariable=self.selected_nfft,
                                      values=self.NFFT_OPTIONS, state="readonly", width=10)
        nfft_dropdown.grid(row=2, column=1, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        # Plot button
        ttk.Button(frame, text='Plot Spectrogram', command=self._plot_spectrogram).grid(
            row=3, column=0, columnspan=4, padx=PAD_X, pady=10, sticky='ew')
    
    def _create_color_controls(self, parent):
        """Create color range control section"""
        frame = ttk.LabelFrame(parent, text="4. Color Range", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        # Dynamic range slider (1 to 11 decades)
        ttk.Label(frame, text='Dynamic Range (decades):', anchor='w').pack(
            padx=5, pady=(5, 0), anchor='w')
        
        slider_frame = ttk.Frame(frame)
        slider_frame.pack(fill='x', padx=5, pady=5)
        
        self.dyn_range_var = tk.DoubleVar(value=6.0)
        self.dyn_range_slider = ttk.Scale(slider_frame, from_=1, to=11,
                                           variable=self.dyn_range_var, orient='horizontal',
                                           command=self._on_dyn_range_changed)
        self.dyn_range_slider.pack(side=tk.LEFT, fill='x', expand=True)
        
        self.dyn_range_label = ttk.Label(slider_frame, text='6.0', width=5)
        self.dyn_range_label.pack(side=tk.RIGHT, padx=(5, 0))
        
        # Labels for scale reference
        ref_frame = ttk.Frame(frame)
        ref_frame.pack(fill='x', padx=5, pady=(0, 5))
        ttk.Label(ref_frame, text='1', font=('TkDefaultFont', 8)).pack(side=tk.LEFT)
        ttk.Label(ref_frame, text='11', font=('TkDefaultFont', 8)).pack(side=tk.RIGHT)
    
    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = ttk.LabelFrame(parent, text="5. Save Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        self.save_button = ttk.Button(frame, text='Save as NPZ', 
                                       command=self._save_data, state='disabled')
        self.save_button.pack(fill='x', padx=PAD_X, pady=PAD_Y)
    
    def _save_data(self):
        """Save spectrogram data to NPZ file (only selected time/freq range)"""
        if self.spectrogram_data is None or self.signal_data is None:
            messagebox.showwarning("Warning", "No spectrogram data to save")
            return
        
        # Get parameters for filename
        shot = self.signal_data['shot']
        channel = self.signal_data['channel'].split('(')[0].strip()  # Remove position info
        t_min = float(self.time_min_entry.get())
        t_max = float(self.time_max_entry.get())
        f_min = float(self.freq_min_entry.get())  # kHz
        f_max = float(self.freq_max_entry.get())  # kHz
        
        # Default filename
        default_name = f"spectrogram_{shot}_{channel}_{t_min:.1f}-{t_max:.1f}s.npz"
        
        # File dialog
        filepath = filedialog.asksaveasfilename(
            initialdir=os.path.expanduser("~"),
            defaultextension='.npz',
            filetypes=[('NumPy NPZ', '*.npz')],
            initialfile=default_name,
            title='Save Spectrogram Data'
        )
        
        if not filepath:
            return
        
        try:
            # Get full data
            f_full = self.spectrogram_data['f']  # Hz
            t_full = self.spectrogram_data['t']  # s
            Sxx_log_full = self.spectrogram_data['Sxx_log']  # (freq, time)
            
            # Apply frequency mask (convert kHz to Hz for comparison)
            f_mask = (f_full >= f_min * 1e3) & (f_full <= f_max * 1e3)
            
            # Slice data to selected frequency range
            f_sliced = f_full[f_mask] / 1e3  # Convert to kHz
            Sxx_log_sliced = Sxx_log_full[f_mask, :]  # (freq_sliced, time)
            
            # Prepare metadata
            metadata = {
                'shot': shot,
                'channel': self.signal_data['channel'],
                'diagnostic': self.selected_diag.get(),
                'time_range': [t_min, t_max],
                'freq_range': [f_min, f_max],
                'nfft': int(self.selected_nfft.get()),
            }
            
            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=t_full,
                     frequency=f_sliced,  # Already in kHz, sliced
                     power=Sxx_log_sliced.T  # Transpose: (time, freq)
            )
            
            print(f"Spectrogram data saved to: {filepath}")
            print(f"  Time: {len(t_full)} points, Freq: {len(f_sliced)} points")
            messagebox.showinfo("Saved", f"Data saved to:\n{filepath}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {str(e)}")
    
    def _set_signal_selection_state(self, state):
        """Enable or disable signal selection widgets"""
        # state: 'normal' or 'disabled'
        combo_state = 'readonly' if state == 'normal' else 'disabled'
        
        for widget in self.signal_widgets:
            widget.config(state=combo_state)
    
    def _load_shot_info(self):
        """Load shot info and populate channel list"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return
        
        self.current_shot = shot_number
        
        # Update shot display label
        self.shot_display_label.config(text=f'#{shot_number}', foreground='black')
        
        # Enable signal selection widgets
        self._set_signal_selection_state('normal')
        
        # Load channels for current diagnostic
        diag_type = self.selected_diag.get()
        self._update_channel_list(diag_type, shot_number)
    
    def _update_channel_list(self, diag_type, shot_number):
        """Update channel list based on diagnostic type"""
        if diag_type == 'ECE':
            self._load_ece_channels(shot_number)
        elif diag_type == 'BES':
            self._load_bes_channels(shot_number)
        elif diag_type == 'TCI':
            self._load_tci_channels()
        elif diag_type.startswith('ECEI-'):
            device = diag_type.split('-')[1]  # GT, GR, or HT
            self._load_ecei_channels(shot_number, device)
        else:
            self._load_mirnov_channels(diag_type)
    
    def _load_ece_channels(self, shot_number):
        """Load ECE channel info for given shot using ECELoader"""
        if shot_number in self.ece_info:
            info = self.ece_info[shot_number]
            self._update_channel_dropdown_ece(info)
            return
        
        try:
            print(f"Loading ECE channel info for #{shot_number}...")
            
            # Use ECELoader to get channel positions
            positions = self.ece_loader.get_channel_positions(shot_number)
            
            # Cache the info
            self.ece_info[shot_number] = positions
            
            print(f"  I_TF = {positions['I_TF']:.2f} kA, {len(positions['channels'])} channels in range")
            
            self._update_channel_dropdown_ece(positions)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load ECE info: {str(e)}")
    
    def _update_channel_dropdown_ece(self, info):
        """Update channel dropdown with ECE channels (R, Z positions)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"ECE{ch:02d} (R={R:.3f}m, Z={Z:.3f}m)")
        
        self.channel_dropdown['values'] = channel_list
        if channel_list:
            self.channel_dropdown.current(0)
    
    def _load_bes_channels(self, shot_number):
        """Load BES channel info for given shot using BESLoader"""
        if shot_number in self.bes_info:
            info = self.bes_info[shot_number]
            self._update_channel_dropdown_bes(info)
            return
        
        try:
            print(f"Loading BES channel info for #{shot_number}...")
            
            # Use BESLoader to get channel positions
            positions = self.bes_loader.get_channel_positions(shot_number)
            
            if not positions['channels']:
                messagebox.showwarning("Warning", "No BES channels available for this shot")
                return
            
            # Cache the info
            self.bes_info[shot_number] = positions
            
            print(f"  {len(positions['channels'])} BES channels loaded")
            
            self._update_channel_dropdown_bes(positions)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load BES info: {str(e)}")
    
    def _update_channel_dropdown_bes(self, info):
        """Update channel dropdown with BES channels (R, Z positions)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"BES_{ch} (R={R:.3f}m, Z={Z:.3f}m)")
        
        self.channel_dropdown['values'] = channel_list
        if channel_list:
            self.channel_dropdown.current(0)
    
    def _load_tci_channels(self):
        """Load TCI channel list (processed ne and raw signals)"""
        # Processed ne channels
        channels = [f'TCI{i:02d} (ne)' for i in range(1, 6)]
        # Raw signal channels
        channels += [f'TCI{i:02d} (raw)' for i in range(1, 6)]
        
        self.channel_dropdown['values'] = channels
        if channels:
            self.channel_dropdown.current(0)
    
    def _load_mirnov_channels(self, diag_type):
        """Load Mirnov channel list"""
        from config.diagnostic_config import DIAGNOSTICS
        mirnov_config = DIAGNOSTICS.get('Mirnov', {})
        
        if 'Toroidal' in diag_type:
            channels = mirnov_config.get('toroidal_channels', [])
        else:
            channels = mirnov_config.get('poloidal_channels', [])
        
        self.channel_dropdown['values'] = channels
        if channels:
            self.channel_dropdown.current(0)
    
    def _load_ecei_channels(self, shot_number, device):
        """Load ECEI channel info for given shot and device using ECEILoader"""
        cache_key = (shot_number, device)
        if cache_key in self.ecei_info:
            info = self.ecei_info[cache_key]
            self._update_channel_dropdown_ecei(info, device)
            return
        
        try:
            print(f"Loading ECEI-{device} channel info for #{shot_number}...")
            
            # Use ECEILoader to get channel positions
            positions = self.ecei_loader.get_channel_positions(shot_number, device)
            
            # GR excluded vertical channels
            excluded_v = ECEILoader.GR_EXCLUDED_VERTICAL if device == 'GR' else []
            
            # Build channel list with positions (flattened from 2D arrays)
            channels = []
            R_values = []
            Z_values = []
            
            for v in range(1, 25):
                if v in excluded_v:
                    continue
                
                for r in range(1, 9):
                    ch_name = f'{device}{v:02d}{r:02d}'
                    R_pos = positions['R'][v - 1, r - 1]
                    Z_pos = positions['Z'][v - 1, r - 1]
                    
                    channels.append(ch_name)
                    R_values.append(R_pos)
                    Z_values.append(Z_pos)
            
            # Cache the info
            self.ecei_info[cache_key] = {
                'channels': channels,
                'R': np.array(R_values),
                'Z': np.array(Z_values),
                'Bt': positions['Bt'],
                'mode': positions['mode'],
                'LO': positions['LO']
            }
            
            print(f"  Bt = {positions['Bt']:.1f} T, mode = {positions['mode']}, LO = {positions['LO']} GHz")
            print(f"  R range: {min(R_values):.3f} - {max(R_values):.3f} m")
            print(f"  {len(channels)} channels loaded")
            
            self._update_channel_dropdown_ecei(self.ecei_info[cache_key], device)
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load ECEI-{device} info: {str(e)}")
    
    def _update_channel_dropdown_ecei(self, info, device):
        """Update channel dropdown with ECEI channels (R, Z)"""
        channel_list = []
        for ch, R, Z in zip(info['channels'], info['R'], info['Z']):
            channel_list.append(f"{ch} (R={R:.3f}m, Z={Z:.3f}m)")
        
        self.channel_dropdown['values'] = channel_list
        if channel_list:
            self.channel_dropdown.current(0)
    
    def _on_diag_changed(self, event=None):
        """Handle diagnostic type change"""
        if self.current_shot is None:
            return
        
        diag_type = self.selected_diag.get()
        self._update_channel_list(diag_type, self.current_shot)
    
    def _format_frequency(self, freq_hz):
        """Format frequency with appropriate unit (Hz, kHz, MHz, GHz)"""
        if freq_hz >= 1e9:
            return f'{int(freq_hz / 1e9)}GHz'
        elif freq_hz >= 1e6:
            return f'{int(freq_hz / 1e6)}MHz'
        elif freq_hz >= 1e3:
            return f'{int(freq_hz / 1e3)}kHz'
        else:
            return f'{int(freq_hz)}Hz'
    
    def _on_dyn_range_changed(self, value):
        """Handle dynamic range slider change"""
        val = float(value)
        self.dyn_range_label.config(text=f'{val:.1f}')
        
        # Update colorbar if spectrogram exists
        if self.im is not None and self.spectrogram_data is not None:
            vmax = self.spectrogram_data['vmax']
            new_vmin = vmax - val
            self.im.set_clim(vmin=new_vmin, vmax=vmax)
            self.canvas.draw_idle()
    
    def _load_signal_data(self):
        """Load selected signal data"""
        if self.current_shot is None:
            messagebox.showerror("Error", "Please fetch a shot first")
            return None
        
        shot_number = self.current_shot
        channel_str = self.selected_channel.get()
        
        if not channel_str:
            messagebox.showerror("Error", "Please select a channel")
            return None
        
        diag_type = self.selected_diag.get()
        
        try:
            mds = Connection(self.app_config.MDS_IP)
            mds.openTree('kstar', shot_number)
            
            if diag_type == 'ECE':
                # Parse ECE channel number: "ECE05 (R=2.105m, Z=0.000m)"
                ch_num = int(channel_str.split('(')[0].replace('ECE', '').strip())
                node_name = f'\\ECE{ch_num:02d}'
            elif diag_type == 'BES':
                # Parse BES channel: "BES_0108 (R=...)" or "BES_0108"
                if '(' in channel_str:
                    ch_name = channel_str.split('(')[0].replace('BES_', '').strip()
                else:
                    ch_name = channel_str.replace('BES_', '').strip()
                node_name = f'\\BES_{ch_name}:FOO'
            elif diag_type == 'TCI':
                # Parse TCI channel: "TCI01 (ne)" or "TCI01 (raw)"
                ch_num = int(channel_str[3:5])
                if '(ne)' in channel_str:
                    node_name = f'\\ne_tci{ch_num:02d}'
                else:  # raw
                    node_name = f'\\tci{ch_num:02d}_ls:foo'
            elif diag_type.startswith('ECEI-'):
                # Parse ECEI channel: "GT0101 (R=1.850m, Z=-0.150m)"
                device = diag_type.split('-')[1]
                ch_name = channel_str.split('(')[0].strip()
                node_name = f'\\ECEI_{ch_name}:FOO'
            else:
                # Mirnov channel name directly
                node_name = f'\\{channel_str}'
            
            print(f"Loading {node_name}...")
            data = mds.get(node_name).data()
            time = mds.get(f'dim_of({node_name})').data()
            
            mds.closeTree('kstar', shot_number)
            
            # TCI raw signal has one extra time frame at the beginning
            if diag_type == 'TCI' and '(raw)' in channel_str:
                time = time[1:]
            
            # Check for bad channel (mean == 0)
            if np.mean(data) == 0:
                messagebox.showwarning("Bad Channel", 
                    f"{channel_str} appears to be a bad channel (mean=0).\n"
                    "Please select another channel.")
                return None
            
            # Build result dict
            result = {
                'time': np.array(time),
                'data': np.array(data),
                'shot': shot_number,
                'channel': channel_str,
                'diag_type': diag_type
            }
            
            # Add position info for ECEI
            if diag_type.startswith('ECEI-'):
                device = diag_type.split('-')[1]
                cache_key = (shot_number, device)
                if cache_key in self.ecei_info:
                    info = self.ecei_info[cache_key]
                    ch_name = channel_str.split('(')[0].strip()
                    try:
                        idx = info['channels'].index(ch_name)
                        result['R'] = info['R'][idx]
                        result['Z'] = info['Z'][idx]
                    except ValueError:
                        pass
            
            return result
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load signal: {str(e)}")
            return None
    
    def _plot_spectrogram(self):
        """Calculate and plot spectrogram"""
        # Load signal data
        signal = self._load_signal_data()
        if signal is None:
            return
        
        self.signal_data = signal
        
        # Get parameters
        try:
            t_min = float(self.time_min_entry.get())
            t_max = float(self.time_max_entry.get())
            f_min = float(self.freq_min_entry.get()) * 1e3  # kHz to Hz
            f_max = float(self.freq_max_entry.get()) * 1e3
            nfft = int(self.selected_nfft.get())
        except ValueError:
            messagebox.showerror("Error", "Invalid parameter values")
            return
        
        # Time slice
        time = signal['time']
        data = signal['data']
        
        mask = (time >= t_min) & (time <= t_max)
        if np.sum(mask) < nfft:
            messagebox.showerror("Error", "Time range too short for selected NFFT")
            return
        
        time_slice = time[mask]
        data_slice = data[mask]
        
        # Calculate sampling frequency
        fs = 1.0 / (time[1] - time[0])
        
        # Calculate spectrogram
        print(f"Calculating spectrogram (NFFT={nfft}, fs={fs/1e6:.2f}MHz)...")
        f, t_spec, Sxx = spectrogram(data_slice, fs, nperseg=nfft)
        
        # Convert to log scale (a.u.)
        Sxx_log = np.log10(Sxx + 1e-20)  # Avoid log(0)
        vmax = np.max(Sxx_log)
        
        # Store spectrogram data
        self.spectrogram_data = {
            'f': f,
            't': t_spec + t_min,  # Offset to actual time
            'Sxx_log': Sxx_log,
            'vmax': vmax
        }
        
        # Clear and plot
        self.ax.clear()
        
        # Frequency mask
        f_mask = (f >= f_min) & (f <= f_max)
        
        # Calculate vmin from dynamic range slider
        dyn_range = self.dyn_range_var.get()
        vmin = vmax - dyn_range
        
        self.im = self.ax.imshow(
            Sxx_log[f_mask, :],
            aspect='auto',
            origin='lower',
            extent=[t_spec[0] + t_min, t_spec[-1] + t_min, 
                    f[f_mask][0]/1e3, f[f_mask][-1]/1e3],
            vmin=vmin,
            vmax=vmax,
            cmap='viridis'
        )
        
        self.ax.set_xlabel('Time [s]')
        self.ax.set_ylabel('Frequency [kHz]')
        
        # Build title with position info for ECEI and sampling rate
        if signal.get('diag_type', '').startswith('ECEI-') and 'R' in signal:
            R = signal['R']
            Z = signal['Z']
            ch_name = signal['channel'].split('(')[0].strip()
            title = f"#{signal['shot']} ECEI_{ch_name} (R={R:.3f}m, Z={Z:.3f}m)"
        else:
            title = f"#{signal['shot']} {signal['channel']}"
        
        # Add sampling rate to title
        fs_str = self._format_frequency(fs)
        title += f" - {fs_str}"
        self.ax.set_title(title)
        
        # Update or create colorbar
        if self.colorbar is not None:
            self.colorbar.remove()
        self.colorbar = self.figure.colorbar(self.im, ax=self.ax, label='log$_{10}$(Power) [a.u.]')
        
        self.canvas.draw()
        print("Spectrogram plotted.")
        
        # Enable save button
        self.save_button.config(state='normal')