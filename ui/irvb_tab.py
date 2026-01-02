#!/usr/bin/python3.8

"""
IRVB (Infra-Red Video Bolometer) tab
2D radiation profile with EFIT overlay and regional Prad time traces
"""

import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.gridspec import GridSpec
from scipy.interpolate import RectBivariateSpline

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_LONG
)


class IRVBTab:
    """IRVB visualization tab"""
    
    # Default psi boundaries for region separation
    DEFAULT_PSI_BOUNDARIES = "0.7, 1.0"
    MAX_BOUNDARIES = 5
    
    # 2D plot color settings
    PRAD_VMIN = 0.0
    PRAD_VMAX = 1.5
    PRAD_LEVELS = 101
    
    # Background psi contour levels (9 levels from 0.1 to 0.9)
    PSI_BG_LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    
    # Region colors (tab10 colormap)
    REGION_COLORS = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5']
    
    def __init__(self, parent, app_config, diagnostic_config, efit_loader):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        self.efit_loader = efit_loader
        
        self.frame = ttk.Frame(parent)
        self.toolbar = None
        
        # Data storage
        self.irvb_data = None
        self.efit_2d = None
        self.region_prad = None
        self.psi_boundaries = []
        self.shot_number = None
        self.ip_fault_time = None
        
        # Frame control
        self.current_frame = 0
        self.total_frames = 0
        
        # Playback state
        self.is_playing = False
        self.play_job = None
        self._last_frame_time = 0
        self._fps_history = []
        
        # Plot references for fast update
        self.im_2d = None
        self.contour_psi_bg = None
        self.contour_psi_bounds = []
        self.contour_bdry = None
        self.limiter_line = None
        self.maxis_marker = None
        self.time_vlines = []
        self.colorbar = None
    
    def create_widgets(self):
        """Create IRVB tab widgets"""
        # Main figure
        self.figure = Figure((8, 6), tight_layout=False)
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill='both', expand=True)
        
        # Bind mouse wheel for frame navigation
        self.canvas.mpl_connect('scroll_event', self._on_mouse_wheel)
        
        # Control panel
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)
        
        self._create_load_controls(control_frame)
        self._create_efit_controls(control_frame)
        self._create_plot_controls(control_frame)
        self._create_frame_controls(control_frame)
        self._create_playback_controls(control_frame)
        self._create_save_controls(control_frame)
    
    def _create_load_controls(self, parent):
        """Create data loading section"""
        frame = ttk.LabelFrame(parent, text="1. Load IRVB Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Shot number
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_LONG, anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        self.shot_entry = ttk.Entry(frame, width=ENTRY_WIDTH_SHOT)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self._load_data())
        
        self.fetch_button = ttk.Button(frame, text='Fetch', command=self._load_data, 
                                       width=BUTTON_WIDTH_MEDIUM)
        self.fetch_button.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')
        
        # Status label (left-aligned)
        self.status_label = ttk.Label(frame, text='', foreground='gray', anchor='w')
        self.status_label.grid(row=1, column=0, columnspan=3, padx=PAD_X, pady=(0, 5), sticky='w')
    
    def _create_efit_controls(self, parent):
        """Create EFIT settings section with psi boundaries"""
        frame = ttk.LabelFrame(parent, text="2. EFIT Settings", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # EFIT Tree selection
        ttk.Label(frame, text="EFIT Tree:", width=LABEL_WIDTH_LONG, anchor='e').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='e')
        
        efit_options = list(self.app_config.EFIT_TREES.keys())
        self.efit_tree_var = tk.StringVar(value=efit_options[0])
        
        efit_dropdown = ttk.Combobox(frame, textvariable=self.efit_tree_var,
                                     values=efit_options, state="readonly", width=18)
        efit_dropdown.grid(row=0, column=1, columnspan=2, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        # Psi boundaries
        ttk.Label(frame, text='psi bounds:', width=LABEL_WIDTH_LONG, anchor='e').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='e')
        self.psi_entry = ttk.Entry(frame, width=20)
        self.psi_entry.insert(0, self.DEFAULT_PSI_BOUNDARIES)
        self.psi_entry.grid(row=1, column=1, columnspan=2, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        # Hint label
        hint = ttk.Label(frame, text="(comma-separated, max 5 values)", 
                        font=('TkDefaultFont', 8), foreground='gray')
        hint.grid(row=2, column=0, columnspan=3, padx=PAD_X, pady=(0, 5))
    
    def _create_plot_controls(self, parent):
        """Create plot button section"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        ttk.Button(frame, text='Plot', command=self._plot_data).pack(
            fill='x', padx=PAD_X, pady=PAD_Y)
    
    def _create_frame_controls(self, parent):
        """Create frame navigation controls"""
        frame = ttk.LabelFrame(parent, text="4. Frame Control", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        # Frame slider
        slider_frame = ttk.Frame(frame)
        slider_frame.pack(fill='x', padx=5, pady=5)
        
        self.frame_var = tk.IntVar(value=0)
        self.frame_slider = ttk.Scale(
            slider_frame, from_=0, to=0, orient='horizontal',
            variable=self.frame_var, command=self._on_slider_change
        )
        self.frame_slider.pack(fill='x', expand=True)
        
        # Frame number and time display
        num_frame = ttk.Frame(frame)
        num_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(num_frame, text='Frame:').pack(side=tk.LEFT)
        
        self.frame_entry = ttk.Entry(num_frame, width=8)
        self.frame_entry.pack(side=tk.LEFT, padx=5)
        self.frame_entry.insert(0, '1')
        self.frame_entry.bind('<Return>', self._on_frame_entry)
        
        self.frame_total_label = ttk.Label(num_frame, text='/ 0')
        self.frame_total_label.pack(side=tk.LEFT)
        
        ttk.Button(num_frame, text='Go', width=5, command=self._goto_frame).pack(
            side=tk.LEFT, padx=5)
        
        # Time display
        time_frame = ttk.Frame(frame)
        time_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(time_frame, text='Time:').pack(side=tk.LEFT)
        self.time_label = ttk.Label(time_frame, text='0.000 s', width=12)
        self.time_label.pack(side=tk.LEFT, padx=5)
        
        # Navigation buttons
        nav_frame = ttk.Frame(frame)
        nav_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Button(nav_frame, text='|<', width=4, command=self._goto_first).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        ttk.Button(nav_frame, text='<10', width=4, command=lambda: self._step_frame(-10)).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        ttk.Button(nav_frame, text='<', width=4, command=lambda: self._step_frame(-1)).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        ttk.Button(nav_frame, text='>', width=4, command=lambda: self._step_frame(1)).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        ttk.Button(nav_frame, text='10>', width=4, command=lambda: self._step_frame(10)).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        ttk.Button(nav_frame, text='>|', width=4, command=self._goto_last).pack(
            side=tk.LEFT, expand=True, fill='x', padx=1)
        
        # Mouse wheel hint
        hint = ttk.Label(frame, text="(Mouse wheel: navigate frames)",
                        font=('TkDefaultFont', 8), foreground='gray')
        hint.pack(pady=(0, 5))
    
    def _create_playback_controls(self, parent):
        """Create playback control section"""
        frame = ttk.LabelFrame(parent, text="5. Playback Control", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        # Play/Pause/Stop buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        self.play_btn = ttk.Button(btn_frame, text='Play', command=self._toggle_play)
        self.play_btn.pack(side=tk.LEFT, expand=True, fill='x', padx=2)
        
        ttk.Button(btn_frame, text='Stop', command=self._stop_play).pack(
            side=tk.LEFT, expand=True, fill='x', padx=2)
        
        # Speed control
        speed_frame = ttk.Frame(frame)
        speed_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(speed_frame, text='Speed:').pack(side=tk.LEFT)
        
        self.speed_var = tk.StringVar(value='1x')
        speed_combo = ttk.Combobox(speed_frame, textvariable=self.speed_var,
                                    values=['0.5x', '1x', '2x', 'Max'],
                                    state='readonly', width=6)
        speed_combo.pack(side=tk.LEFT, padx=5)
        
        # Loop checkbox
        self.loop_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(speed_frame, text='Loop', variable=self.loop_var).pack(
            side=tk.LEFT, padx=10)
        
        # Actual FPS display
        fps_frame = ttk.Frame(frame)
        fps_frame.pack(fill='x', padx=5, pady=2)
        
        self.actual_fps_label = ttk.Label(fps_frame, text='Actual: -- FPS',
                                           foreground='gray')
        self.actual_fps_label.pack(side=tk.LEFT)
    
    def _create_save_controls(self, parent):
        """Create save data section"""
        frame = ttk.LabelFrame(parent, text="6. Save Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        self.save_button = ttk.Button(frame, text='Save as NPZ', 
                                       command=self._save_data, state='disabled')
        self.save_button.pack(fill='x', padx=PAD_X, pady=PAD_Y)
    
    def _save_data(self):
        """Save IRVB data to NPZ file including EFIT data"""
        if self.irvb_data is None or self.region_prad is None:
            messagebox.showwarning("Warning", "No data to save. Plot data first.")
            return
        
        if self.efit_2d is None:
            messagebox.showwarning("Warning", "No EFIT data available.")
            return
        
        # Default filename
        t_min = self.irvb_data.time[0]
        t_max = self.irvb_data.time[-1]
        default_name = f"irvb_{self.shot_number}_{t_min:.1f}-{t_max:.1f}s.npz"
        
        # File dialog
        filepath = filedialog.asksaveasfilename(
            initialdir=os.path.expanduser("~"),
            defaultextension='.npz',
            filetypes=[('NumPy NPZ', '*.npz')],
            initialfile=default_name,
            title='Save IRVB Data'
        )
        
        if not filepath:
            return
        
        try:
            # Build region labels
            n_regions = len(self.psi_boundaries) + 1
            boundaries = [0] + self.psi_boundaries + [np.inf]
            region_labels = []
            for i in range(n_regions):
                psi_min = boundaries[i]
                psi_max = boundaries[i + 1]
                if psi_max == np.inf:
                    region_labels.append(f'psi_N > {psi_min}')
                else:
                    region_labels.append(f'{psi_min} < psi_N < {psi_max}')
            
            # Prepare metadata
            metadata = {
                'shot': self.shot_number,
                'psi_boundaries': self.psi_boundaries,
                'efit_tree': self.efit_tree_var.get(),
                'time_range': [t_min, t_max],
                'region_labels': region_labels,
            }
            
            # Prepare EFIT data for each IRVB time point
            n_frames = len(self.irvb_data.time)
            max_bdry = self.efit_2d.bdry_r.shape[1]
            
            efit_bdry_r = np.zeros((n_frames, max_bdry))
            efit_bdry_z = np.zeros((n_frames, max_bdry))
            efit_nbdry = np.zeros(n_frames, dtype=int)
            efit_rmaxis = np.zeros(n_frames)
            efit_zmaxis = np.zeros(n_frames)
            efit_psi_n = np.zeros((n_frames, len(self.efit_2d.z_grid), len(self.efit_2d.r_grid)))
            
            for i, t in enumerate(self.irvb_data.time):
                efit_idx = self.efit_2d.find_time_index(t)
                bdry_r, bdry_z = self.efit_2d.get_boundary(efit_idx)
                nbdry = len(bdry_r)
                
                efit_bdry_r[i, :nbdry] = bdry_r
                efit_bdry_z[i, :nbdry] = bdry_z
                efit_nbdry[i] = nbdry
                efit_rmaxis[i], efit_zmaxis[i] = self.efit_2d.get_magnetic_axis(efit_idx)
                efit_psi_n[i] = self.efit_2d.get_psi_normalized(efit_idx)
            
            # Save to NPZ
            np.savez(filepath,
                     metadata=np.array(metadata),
                     time=self.irvb_data.time,
                     R=self.irvb_data.x_grid,
                     Z=self.irvb_data.y_grid,
                     prad_2d=self.irvb_data.recon,
                     region_prad=self.region_prad,
                     ptot=self.irvb_data.ptot[:len(self.irvb_data.time)],
                     # EFIT data
                     efit_r_grid=self.efit_2d.r_grid,
                     efit_z_grid=self.efit_2d.z_grid,
                     efit_psi_n=efit_psi_n,
                     efit_bdry_r=efit_bdry_r,
                     efit_bdry_z=efit_bdry_z,
                     efit_nbdry=efit_nbdry,
                     efit_rmaxis=efit_rmaxis,
                     efit_zmaxis=efit_zmaxis,
                     efit_limiter_r=self.efit_2d.limiter_r,
                     efit_limiter_z=self.efit_2d.limiter_z
            )
            
            print(f"IRVB data saved to: {filepath}")
            messagebox.showinfo("Saved", f"Data saved to:\n{filepath}")
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save data: {str(e)}")
    
    def _parse_psi_boundaries(self):
        """Parse psi boundary input string"""
        try:
            text = self.psi_entry.get().strip()
            if not text:
                return []
            
            values = [float(v.strip()) for v in text.split(',')]
            
            # Validate count
            if len(values) > self.MAX_BOUNDARIES:
                messagebox.showwarning("Warning", 
                    f"Maximum {self.MAX_BOUNDARIES} boundaries allowed. Using first {self.MAX_BOUNDARIES}.")
                values = values[:self.MAX_BOUNDARIES]
            
            # Sort and validate range
            values = sorted(values)
            if any(v <= 0 or v >= 2.0 for v in values):
                messagebox.showerror("Error", "psi values should be between 0 and 2.0")
                return None
            
            return values
            
        except ValueError:
            messagebox.showerror("Error", "Invalid psi boundary format. Use comma-separated numbers.")
            return None
    
    def _update_status(self, message, success=True):
        """Update status label with message"""
        color = 'green' if success else 'red'
        self.status_label.config(text=message, foreground=color)
    
    def _load_data(self):
        """Load IRVB data from server"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            self._update_status('Invalid shot number', success=False)
            return
        
        # Load IRVB data
        try:
            from data_loaders.irvb_loader import IRVBLoader
            loader = IRVBLoader(self.app_config)
            self.irvb_data = loader.load_data(shot_number)
            self.shot_number = shot_number
            
            # Get IP fault time and slice data
            self._get_ip_fault_time()
            self._slice_by_ip_fault()
            
            # Update frame controls
            self.total_frames = len(self.irvb_data.time)
            self.current_frame = 0
            self.frame_slider.config(to=self.total_frames - 1)
            self.frame_var.set(0)
            self.frame_total_label.config(text=f'/ {self.total_frames}')
            self.frame_entry.delete(0, tk.END)
            self.frame_entry.insert(0, '1')
            self.time_label.config(text=f'{self.irvb_data.time[0]:.3f} s')
            
            self._update_status(f'Shot #{shot_number}: {self.total_frames} frames loaded', success=True)
            print(f"IRVB: Data loaded for shot #{shot_number}")
            
        except Exception as e:
            self._update_status(f'Failed: {str(e)[:30]}...', success=False)
            messagebox.showerror("Error", f"Failed to load IRVB data: {str(e)}")
            return
    
    def _get_ip_fault_time(self):
        """Get IP fault time from MDS+"""
        try:
            from MDSplus import Connection
            mds = Connection(self.app_config.MDS_IP)
            mds.openTree('kstar', self.shot_number)
            self.ip_fault_time = mds.get('\\t_ip_fault').data()
            mds.closeTree('kstar', self.shot_number)
            print(f"IRVB: IP fault time = {self.ip_fault_time:.3f} s")
        except:
            self.ip_fault_time = None
            print("IRVB: IP fault time not available")
    
    def _slice_by_ip_fault(self):
        """Slice IRVB data by IP fault time"""
        if self.ip_fault_time is None or self.irvb_data is None:
            return
        
        valid_mask = self.irvb_data.time < self.ip_fault_time
        if not np.any(valid_mask):
            return
        
        self.irvb_data.time = self.irvb_data.time[valid_mask]
        self.irvb_data.recon = self.irvb_data.recon[valid_mask]
        print(f"IRVB: Data sliced to {len(self.irvb_data.time)} frames (before IP fault)")
    
    def _slice_by_efit_time(self):
        """Slice IRVB data by valid time range (0 to EFIT last time)"""
        if self.efit_2d is None or self.irvb_data is None:
            return
        
        efit_end = self.efit_2d.time[-1]
        valid_mask = (self.irvb_data.time >= 0) & (self.irvb_data.time <= efit_end)
        if not np.any(valid_mask):
            return
        
        n_before = len(self.irvb_data.time)
        self.irvb_data.time = self.irvb_data.time[valid_mask]
        self.irvb_data.recon = self.irvb_data.recon[valid_mask]
        
        if len(self.irvb_data.time) < n_before:
            print(f"IRVB: Data sliced to {len(self.irvb_data.time)} frames (0 to EFIT end)")
    
    def _load_efit_data(self):
        """Load EFIT 2D equilibrium data using efit_loader"""
        if self.irvb_data is None:
            return False
        
        efit_display = self.efit_tree_var.get()
        efit_tree = self.app_config.EFIT_TREES[efit_display]
        
        try:
            print(f"IRVB: Loading 2D EFIT data from {efit_tree}...")
            self.efit_2d = self.efit_loader.load_efit_2d(self.shot_number, efit_tree)
            print(f"IRVB: EFIT data loaded ({len(self.efit_2d.time)} timepoints)")
            return True
            
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load EFIT data: {str(e)}")
            return False
    
    def _find_xpoint(self, bdry_r, bdry_z):
        """Find X-point from boundary (lower divertor)"""
        lower_mask = bdry_z < 0
        if not np.any(lower_mask):
            return None, None
        
        lower_r = bdry_r[lower_mask]
        lower_z = bdry_z[lower_mask]
        
        min_idx = np.argmin(lower_z)
        return lower_r[min_idx], lower_z[min_idx]
    
    def _compute_regional_prad(self):
        """Compute Prad for each region defined by psi boundaries"""
        if self.irvb_data is None or self.efit_2d is None:
            return
        
        n_times = len(self.irvb_data.time)
        n_regions = len(self.psi_boundaries) + 1
        
        self.region_prad = np.zeros((n_regions, n_times))
        
        x_grid = self.irvb_data.x_grid
        y_grid = self.irvb_data.y_grid
        X, Y = np.meshgrid(x_grid, y_grid)
        
        # Volume element: 2*pi*R * dR * dZ
        dx = x_grid[1] - x_grid[0]
        dy = y_grid[1] - y_grid[0]
        vol_factor = 2 * np.pi * X * dx * dy
        
        print("IRVB: Computing regional Prad...")
        
        efit_idx_prev = -1
        psi_on_grid = None
        
        for i, t in enumerate(self.irvb_data.time):
            # Find closest EFIT time
            efit_idx = self.efit_2d.find_time_index(t)
            
            # Update psi map if EFIT time changed
            if efit_idx != efit_idx_prev:
                psi_n_2d = self.efit_2d.get_psi_normalized(efit_idx)
                
                # Use RectBivariateSpline for fast vectorized interpolation
                spline = RectBivariateSpline(
                    self.efit_2d.z_grid, self.efit_2d.r_grid, psi_n_2d
                )
                psi_on_grid = spline(y_grid, x_grid)
                
                efit_idx_prev = efit_idx
            
            recon = self.irvb_data.recon[i]
            
            # Compute Prad for each region
            boundaries = [0] + self.psi_boundaries + [np.inf]
            
            for r in range(n_regions):
                psi_min = boundaries[r]
                psi_max = boundaries[r + 1]
                
                mask = (psi_on_grid >= psi_min) & (psi_on_grid < psi_max)
                self.region_prad[r, i] = np.sum(recon * mask * vol_factor)
        
        print("IRVB: Regional Prad computation complete")
    
    def _plot_data(self):
        """Main plot function"""
        if self.irvb_data is None:
            messagebox.showwarning("Warning", "Please load IRVB data first")
            return
        
        # Parse boundaries
        self.psi_boundaries = self._parse_psi_boundaries()
        if self.psi_boundaries is None:
            return
        
        # Load EFIT
        if not self._load_efit_data():
            return
        
        # Slice data by EFIT time range and update frame controls
        self._slice_by_efit_time()
        self._update_frame_controls()
        
        # Compute regional Prad
        self._compute_regional_prad()
        
        # Setup figure layout
        self._setup_figure()
        
        # Initial plot
        self._update_plot(self.current_frame)
        
        # Enable save button
        self.save_button.config(state='normal')
    
    def _update_frame_controls(self):
        """Update frame controls after data slicing"""
        self.total_frames = len(self.irvb_data.time)
        self.current_frame = 0
        self.frame_slider.config(to=self.total_frames - 1)
        self.frame_var.set(0)
        self.frame_total_label.config(text=f'/ {self.total_frames}')
        self.frame_entry.delete(0, tk.END)
        self.frame_entry.insert(0, '1')
        self.time_label.config(text=f'{self.irvb_data.time[0]:.3f} s')
    
    def _setup_figure(self):
        """Setup figure with dynamic number of subplots"""
        self.figure.clear()
        
        n_regions = len(self.psi_boundaries) + 1
        
        # GridSpec layout
        gs = GridSpec(n_regions, 2, figure=self.figure, 
                     width_ratios=[1.5, 1], wspace=0.0, hspace=0.0,
                     left=0.10, right=0.95, top=0.95, bottom=0.10)
        
        # Create time trace axes (left column)
        self.ax_traces = []
        for i in range(n_regions):
            if i == 0:
                ax = self.figure.add_subplot(gs[i, 0])
            else:
                ax = self.figure.add_subplot(gs[i, 0], sharex=self.ax_traces[0])
            self.ax_traces.append(ax)
        
        # Create 2D profile axis (right column)
        self.ax_2d = self.figure.add_subplot(gs[:, 1])
        
        # Plot time traces
        self._plot_time_traces()
        
        # Reset plot references
        self.im_2d = None
        self.contour_psi_bg = None
        self.contour_psi_bounds = []
        self.contour_bdry = None
        self.limiter_line = None
        self.maxis_marker = None
        self.colorbar = None
    
    def _plot_time_traces(self):
        """Plot regional Prad time traces"""
        time = self.irvb_data.time
        n_regions = len(self.psi_boundaries) + 1
        boundaries = [0] + self.psi_boundaries + [np.inf]
        
        self.time_vlines = []
        
        for i, ax in enumerate(self.ax_traces):
            ax.clear()
            
            # Title on first axis with shot number and EFIT tree
            if i == 0:
                efit_display = self.efit_tree_var.get()
                ax.set_title(f'#{self.shot_number} ({efit_display})')
            
            # Region label
            psi_min = boundaries[i]
            psi_max = boundaries[i + 1]
            
            if psi_max == np.inf:
                label = f'$\\psi_N$ > {psi_min}'
            else:
                label = f'{psi_min} < $\\psi_N$ < {psi_max}'
            
            # Last region (SOL) uses gray color, others use REGION_COLORS
            if i == n_regions - 1:
                color = 'gray'
            else:
                color = self.REGION_COLORS[i % len(self.REGION_COLORS)]
            
            ax.plot(time, self.region_prad[i], color=color, linewidth=2)
            
            # Text annotation instead of legend (no line, text only)
            ax.text(0.02, 0.95, label, transform=ax.transAxes,
                   fontsize=9, verticalalignment='top',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.5,
                             edgecolor='gray', linewidth=0.5))
            
            # ylabel without psi range, fixed position
            ax.set_ylabel(r'$P_{rad}$ [MW]')
            ax.yaxis.set_label_coords(-0.08, 0.5)
            
            # Set ylim with bottom=0
            y_max = np.nanmax(self.region_prad[i]) * 1.1
            if y_max <= 0 or np.isnan(y_max):
                y_max = 1.0
            ax.set_ylim(0, y_max)
            
            # Only show x-label on bottom plot
            if i == n_regions - 1:
                ax.set_xlabel('Time [s]')
            else:
                plt.setp(ax.get_xticklabels(), visible=False)
            
            ax.set_xlim(0, self.efit_2d.time[-1])
            ax.grid(ls='--', lw=0.3, c='lightgray')
            
            # Vertical line for current time
            vline = ax.axvline(time[self.current_frame], color='k', 
                              linewidth=1.5, linestyle='--')
            self.time_vlines.append(vline)
    
    def _update_plot(self, frame_idx):
        """Update 2D plot and time markers for given frame"""
        if self.irvb_data is None:
            return
        
        time = self.irvb_data.time[frame_idx]
        recon = self.irvb_data.recon[frame_idx]
        x_grid = self.irvb_data.x_grid
        y_grid = self.irvb_data.y_grid
        
        # Get EFIT data at this time
        efit_idx = self.efit_2d.find_time_index(time)
        psi_n = self.efit_2d.get_psi_normalized(efit_idx)
        bdry_r, bdry_z = self.efit_2d.get_boundary(efit_idx)
        maxis_r, maxis_z = self.efit_2d.get_magnetic_axis(efit_idx)
        
        # Update 2D plot
        if self.im_2d is None:
            self.ax_2d.clear()
            
            # Prad contour fill
            levels = np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, self.PRAD_LEVELS)
            self.im_2d = self.ax_2d.contourf(x_grid, y_grid, recon, 
                                             levels=levels, 
                                             vmin=self.PRAD_VMIN,
                                             vmax=self.PRAD_VMAX,
                                             extend='max')
            
            # Colorbar
            self.colorbar = self.figure.colorbar(
                self.im_2d, ax=self.ax_2d,
                ticks=np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, 5),
                shrink=0.8
            )
            self.colorbar.set_label(r'$P_{rad}$ [MW/m$^{3}$]')
            
            self.ax_2d.set_xlabel('R [m]')
            self.ax_2d.set_ylabel('Z [m]')
            self.ax_2d.set_xlim(x_grid[0], x_grid[-1])
            self.ax_2d.set_ylim(y_grid[0], y_grid[-1])
            self.ax_2d.set_aspect('equal')
        else:
            # Update contourf
            for coll in self.im_2d.collections:
                coll.remove()
            levels = np.linspace(self.PRAD_VMIN, self.PRAD_VMAX, self.PRAD_LEVELS)
            self.im_2d = self.ax_2d.contourf(x_grid, y_grid, recon,
                                             levels=levels,
                                             vmin=self.PRAD_VMIN,
                                             vmax=self.PRAD_VMAX,
                                             extend='max')
        
        # Remove old contours
        if self.contour_psi_bg is not None:
            for coll in self.contour_psi_bg.collections:
                coll.remove()
        for contour in self.contour_psi_bounds:
            for coll in contour.collections:
                coll.remove()
        self.contour_psi_bounds = []
        
        # Background psi contours (9 levels, light gray)
        self.contour_psi_bg = self.ax_2d.contour(
            self.efit_2d.r_grid, self.efit_2d.z_grid, psi_n,
            levels=self.PSI_BG_LEVELS, colors='gray',
            linewidths=0.5, linestyles='-', alpha=0.5, zorder=2
        )
        
        # Update plasma boundary (draw before psi bounds so bounds appear on top)
        if self.contour_bdry is not None:
            self.contour_bdry[0].remove()
        self.contour_bdry = self.ax_2d.plot(bdry_r, bdry_z, 'k-', linewidth=2, zorder=3)
        
        # Psi boundary contours with matching colors (drawn on top of bdry)
        for idx, psi_level in enumerate(self.psi_boundaries):
            color = self.REGION_COLORS[idx % len(self.REGION_COLORS)]
            contour = self.ax_2d.contour(
                self.efit_2d.r_grid, self.efit_2d.z_grid, psi_n,
                levels=[psi_level], colors=[color],
                linewidths=2, linestyles='--', zorder=4
            )
            self.contour_psi_bounds.append(contour)
        
        # Update limiter
        if self.limiter_line is not None:
            self.limiter_line[0].remove()
        if self.efit_2d.limiter_r is not None:
            self.limiter_line = self.ax_2d.plot(
                self.efit_2d.limiter_r, self.efit_2d.limiter_z, 
                'k-', linewidth=1.5
            )
        
        # Update magnetic axis marker
        if self.maxis_marker is not None:
            self.maxis_marker[0].remove()
        if maxis_r is not None:
            self.maxis_marker = self.ax_2d.plot(maxis_r, maxis_z, 'kx', 
                                                 markersize=10, markeredgewidth=2)
        
        # Update title
        self.ax_2d.set_title(f'#{self.shot_number} t = {time:.3f} s')
        
        # Update time markers on traces
        for vline in self.time_vlines:
            vline.set_xdata([time, time])
        
        # Update UI elements
        self.current_frame = frame_idx
        self.frame_entry.delete(0, tk.END)
        self.frame_entry.insert(0, str(frame_idx + 1))
        self.time_label.config(text=f'{time:.3f} s')
        
        self.canvas.draw_idle()
    
    # ===== Frame navigation methods =====
    
    def _on_slider_change(self, value):
        """Handle slider change"""
        frame_idx = int(float(value))
        if frame_idx != self.current_frame and self.region_prad is not None:
            self.current_frame = frame_idx
            self._update_plot(frame_idx)
    
    def _on_frame_entry(self, event):
        """Handle frame entry"""
        self._goto_frame()
    
    def _on_mouse_wheel(self, event):
        """Handle mouse wheel for frame navigation"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        
        current = self.frame_var.get()
        
        if event.button == 'up':
            new_frame = max(0, current - 1)
        elif event.button == 'down':
            new_frame = min(self.total_frames - 1, current + 1)
        else:
            return
        
        if new_frame != current:
            self.current_frame = new_frame
            self.frame_var.set(new_frame)
            self._update_plot(new_frame)
    
    def _goto_frame(self):
        """Go to specified frame"""
        if self.region_prad is None:
            return
        
        try:
            frame_num = int(self.frame_entry.get())
            frame_idx = frame_num - 1
            
            if 0 <= frame_idx < self.total_frames:
                self.current_frame = frame_idx
                self.frame_var.set(frame_idx)
                self._update_plot(frame_idx)
            else:
                messagebox.showwarning("Warning",
                    f"Frame number must be between 1 and {self.total_frames}")
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid frame number")
    
    def _step_frame(self, delta):
        """Step by delta frames"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        
        current = self.frame_var.get()
        new_frame = current + delta
        new_frame = max(0, min(new_frame, self.total_frames - 1))
        
        self.current_frame = new_frame
        self.frame_var.set(new_frame)
        self._update_plot(new_frame)
    
    def _goto_first(self):
        """Go to first frame"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        self.current_frame = 0
        self.frame_var.set(0)
        self._update_plot(0)
    
    def _goto_last(self):
        """Go to last frame"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        last_idx = self.total_frames - 1
        self.current_frame = last_idx
        self.frame_var.set(last_idx)
        self._update_plot(last_idx)
    
    # ===== Playback methods =====
    
    def _toggle_play(self):
        """Toggle play/pause"""
        if self.total_frames == 0 or self.region_prad is None:
            return
        
        if self.is_playing:
            self._pause_play()
        else:
            self._start_play()
    
    def _start_play(self):
        """Start playback"""
        self.is_playing = True
        self.play_btn.config(text='Pause')
        self._last_frame_time = time.time()
        self._fps_history = []
        self._play_next_frame()
    
    def _pause_play(self):
        """Pause playback"""
        self.is_playing = False
        self.play_btn.config(text='Play')
        if self.play_job:
            self.frame.after_cancel(self.play_job)
            self.play_job = None
        
        self.actual_fps_label.config(text='Actual: -- FPS')
    
    def _stop_play(self):
        """Stop playback and reset to first frame"""
        self._pause_play()
        if self.total_frames > 0 and self.region_prad is not None:
            self.current_frame = 0
            self.frame_var.set(0)
            self._update_plot(0)
    
    def _play_next_frame(self):
        """Play next frame with speed control"""
        if not self.is_playing:
            return
        
        current = self.frame_var.get()
        next_frame = current + 1
        
        # Handle end of sequence
        if next_frame >= self.total_frames:
            if self.loop_var.get():
                next_frame = 0
            else:
                self._pause_play()
                return
        
        # Measure actual frame time
        now = time.time()
        actual_frame_time = now - self._last_frame_time
        self._last_frame_time = now
        
        # Calculate and display actual FPS
        if actual_frame_time > 0:
            instant_fps = 1.0 / actual_frame_time
            self._fps_history.append(instant_fps)
            if len(self._fps_history) > 10:
                self._fps_history.pop(0)
            avg_fps = sum(self._fps_history) / len(self._fps_history)
            self.actual_fps_label.config(text=f'Actual: {avg_fps:.1f} FPS')
        
        # Update state and display
        self.current_frame = next_frame
        self.frame_var.set(next_frame)
        self._update_plot(next_frame)
        
        # Calculate delay from speed multiplier (1x = 10 FPS = 100ms)
        speed_str = self.speed_var.get()
        if speed_str == 'Max':
            target_delay = 1
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 100
            target_delay = max(1, int(base_delay / speed_mult))
        
        self.play_job = self.frame.after(target_delay, self._play_next_frame)