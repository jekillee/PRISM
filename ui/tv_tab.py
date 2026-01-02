#!/usr/bin/python3.8

"""
TV (Visible Camera) tab for viewing image sequences from ZIP files
With line drawing feature for paper figures
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import zipfile
import io
import re
import os
import time
import threading
import numpy as np
from PIL import Image
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM, LABEL_WIDTH_SHORT
)


class TVTab:
    """TV image sequence viewer tab with line drawing"""
    
    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        
        self.frame = ttk.Frame(parent)
        self.toolbar = None
        
        # Data storage
        self.zip_file = None
        self.image_files = []
        self.current_frame = 0
        self.total_frames = 0
        
        # Image cache (store nearby frames for smooth playback)
        self.cache = {}
        self.cache_size = 100  # Increased for better playback
        
        # Playback state
        self.is_playing = False
        self.play_job = None
        self._last_frame_time = 0  # For adaptive frame timing
        
        # Prefetch thread control
        self._prefetch_lock = threading.Lock()
        self._prefetch_thread = None
        self._prefetch_stop = False
        
        # Plot reference
        self.im = None
        
        # Line drawing
        self.line_points = []
        self.drawn_line = None
        self.preview_line = None
        self.draw_mode = False
        self.click_cid = None
        self.motion_cid = None
    
    def create_widgets(self):
        """Create TV tab widgets"""
        # Left: Image display area
        self.figure = Figure((8, 6), tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("No image loaded")
        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        
        canvas_widget = self.canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill='both', expand=True)
        
        # Bind mouse wheel event for frame navigation
        self.canvas.mpl_connect('scroll_event', self._on_mouse_wheel)
        
        # Right: Control panel
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)
        
        self._create_file_controls(control_frame)
        self._create_frame_controls(control_frame)
        self._create_playback_controls(control_frame)
        self._create_draw_line_controls(control_frame)
    
    def _create_file_controls(self, parent):
        """Create file loading section"""
        frame = ttk.LabelFrame(parent, text="1. Load TV Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)
        
        frame.grid_columnconfigure(1, weight=1)
        
        # Row 0: Shot input and Search button
        ttk.Label(frame, text='Shot', width=LABEL_WIDTH_SHORT, anchor='center').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry = ttk.Entry(frame, width=20)
        self.shot_entry.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        self.shot_entry.bind('<Return>', lambda e: self._search_available_tvs())
        
        ttk.Button(frame, text='Search', command=self._search_available_tvs, 
                  width=BUTTON_WIDTH_MEDIUM).grid(
            row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')
        
        # Row 1: TV dropdown and Load button
        ttk.Label(frame, text='TV', width=LABEL_WIDTH_SHORT, anchor='center').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        self.tv_selection_var = tk.StringVar(value='-- Select TV --')
        self.tv_dropdown = ttk.Combobox(frame, textvariable=self.tv_selection_var,
                                         state='readonly', width=17)
        self.tv_dropdown.grid(row=1, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        ttk.Button(frame, text='Load', command=self._load_selected_tv, 
                  width=BUTTON_WIDTH_MEDIUM).grid(
            row=1, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')
        
        # Row 2: Or label and Load from File button
        ttk.Label(frame, text='Or', width=LABEL_WIDTH_SHORT, anchor='center').grid(
            row=2, column=0, padx=PAD_X, pady=PAD_Y, sticky='ew')
        ttk.Button(frame, text='Load from File',
                  command=self._load_zip_file).grid(
            row=2, column=1, columnspan=2, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        # File label
        self.file_label = ttk.Label(frame, text="No file loaded", wraplength=360)
        self.file_label.grid(row=3, column=0, columnspan=3, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        # Loading status label
        self.status_label = ttk.Label(frame, text="", foreground='blue')
        self.status_label.grid(row=4, column=0, columnspan=3, padx=PAD_X, pady=2, sticky='w')
    
    def _create_frame_controls(self, parent):
        """Create frame navigation controls"""
        frame = ttk.LabelFrame(parent, text="2. Frame Control", labelanchor="n")
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
        
        # Frame number display
        num_frame = ttk.Frame(frame)
        num_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(num_frame, text='Frame:').pack(side=tk.LEFT)
        
        self.frame_entry = ttk.Entry(num_frame, width=8)
        self.frame_entry.pack(side=tk.LEFT, padx=5)
        self.frame_entry.insert(0, '0')
        self.frame_entry.bind('<Return>', self._on_frame_entry)
        
        self.frame_total_label = ttk.Label(num_frame, text='/ 0')
        self.frame_total_label.pack(side=tk.LEFT)
        
        ttk.Button(num_frame, text='Go', width=5, command=self._goto_frame).pack(
            side=tk.LEFT, padx=5)
        
        # Current filename display
        self.filename_label = ttk.Label(frame, text="", wraplength=320)
        self.filename_label.pack(fill='x', padx=5, pady=5)
        
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
        hint_label = ttk.Label(frame, text="(Mouse wheel: navigate frames)", 
                               font=('TkDefaultFont', 8), foreground='gray')
        hint_label.pack(pady=(0, 5))
    
    def _create_playback_controls(self, parent):
        """Create playback control section"""
        frame = ttk.LabelFrame(parent, text="3. Playback Control", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        # Play/Pause/Stop buttons
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        self.play_btn = ttk.Button(btn_frame, text='Play', command=self._toggle_play)
        self.play_btn.pack(side=tk.LEFT, expand=True, fill='x', padx=2)
        
        ttk.Button(btn_frame, text='Stop', command=self._stop_play).pack(
            side=tk.LEFT, expand=True, fill='x', padx=2)
        
        # Speed control (1x = 10 FPS base)
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
        fps_display_frame = ttk.Frame(frame)
        fps_display_frame.pack(fill='x', padx=5, pady=2)
        
        self.actual_fps_label = ttk.Label(fps_display_frame, text='Actual: -- FPS', 
                                           foreground='gray')
        self.actual_fps_label.pack(side=tk.LEFT)
        
        # Frame skip
        skip_frame = ttk.Frame(frame)
        skip_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(skip_frame, text='Frame skip:').pack(side=tk.LEFT)
        
        self.skip_var = tk.StringVar(value='1')
        skip_combo = ttk.Combobox(skip_frame, textvariable=self.skip_var,
                                   values=['1', '2', '5', '10', '20'],
                                   state='readonly', width=6)
        skip_combo.pack(side=tk.LEFT, padx=5)
    
    def _create_draw_line_controls(self, parent):
        """Create line drawing controls"""
        frame = ttk.LabelFrame(parent, text="4. Draw Line", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        # Draw mode button and Clear button
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', padx=5, pady=5)
        
        self.draw_btn = ttk.Button(btn_frame, text='Draw Mode: OFF', 
                                   command=self._toggle_draw_mode)
        self.draw_btn.pack(side=tk.LEFT, expand=True, fill='x', padx=2)
        
        ttk.Button(btn_frame, text='Clear', command=self._clear_line).pack(
            side=tk.LEFT, expand=True, fill='x', padx=2)
        
        # Show line checkbox and smooth curve checkbox
        check_frame = ttk.Frame(frame)
        check_frame.pack(fill='x', padx=5, pady=5)
        
        self.show_line_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(check_frame, text='Show Line', 
                       variable=self.show_line_var,
                       command=self._update_line_display).pack(side=tk.LEFT)
        
        self.smooth_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(check_frame, text='Smooth', 
                       variable=self.smooth_var,
                       command=self._update_line_display).pack(side=tk.LEFT, padx=10)
        
        # Points count label
        self.points_label = ttk.Label(check_frame, text='Points: 0')
        self.points_label.pack(side=tk.RIGHT, padx=5)
        
        # Line style options
        style_frame = ttk.Frame(frame)
        style_frame.pack(fill='x', padx=5, pady=5)
        
        ttk.Label(style_frame, text='Color:').pack(side=tk.LEFT)
        self.line_color_var = tk.StringVar(value='black')
        color_combo = ttk.Combobox(style_frame, textvariable=self.line_color_var,
                                    values=['black', 'white', 'red', 'blue', 'yellow', 'green'],
                                    state='readonly', width=8)
        color_combo.pack(side=tk.LEFT, padx=5)
        color_combo.bind('<<ComboboxSelected>>', lambda e: self._update_line_display())
        
        ttk.Label(style_frame, text='Width:').pack(side=tk.LEFT, padx=(10, 0))
        self.line_width_var = tk.StringVar(value='2')
        width_combo = ttk.Combobox(style_frame, textvariable=self.line_width_var,
                                    values=['1', '2', '3', '4', '5'],
                                    state='readonly', width=4)
        width_combo.pack(side=tk.LEFT, padx=5)
        width_combo.bind('<<ComboboxSelected>>', lambda e: self._update_line_display())
        
        # Hint label
        hint_label = ttk.Label(frame, text="(Left-click: add point, Right-click: finish)", 
                               font=('TkDefaultFont', 8), foreground='gray')
        hint_label.pack(pady=(0, 5))
    
    def _toggle_draw_mode(self):
        """Toggle line drawing mode"""
        self.draw_mode = not self.draw_mode
        
        if self.draw_mode:
            self.draw_btn.config(text='Draw Mode: ON')
            self.click_cid = self.canvas.mpl_connect('button_press_event', 
                                                      self._on_line_click)
            self.motion_cid = self.canvas.mpl_connect('motion_notify_event',
                                                       self._on_mouse_motion)
            self.canvas.get_tk_widget().config(cursor='crosshair')
        else:
            self.draw_btn.config(text='Draw Mode: OFF')
            if self.click_cid:
                self.canvas.mpl_disconnect(self.click_cid)
                self.click_cid = None
            if self.motion_cid:
                self.canvas.mpl_disconnect(self.motion_cid)
                self.motion_cid = None
            self.canvas.get_tk_widget().config(cursor='')
            # Remove preview line
            if self.preview_line:
                self.preview_line.remove()
                self.preview_line = None
                self.canvas.draw_idle()
    
    def _on_line_click(self, event):
        """Handle mouse click for line drawing"""
        if event.inaxes != self.ax:
            return
        
        if event.button == 1:  # Left click - add point
            self.line_points.append((event.xdata, event.ydata))
            self.points_label.config(text=f'Points: {len(self.line_points)}')
            self._update_line_display()
        elif event.button == 3:  # Right click - finish drawing
            self._toggle_draw_mode()
    
    def _on_mouse_motion(self, event):
        """Handle mouse motion for preview line"""
        if not self.draw_mode or event.inaxes != self.ax:
            return
        
        if len(self.line_points) == 0:
            return
        
        # Update preview line from last point to cursor
        last_point = self.line_points[-1]
        
        if self.preview_line:
            self.preview_line.set_data([last_point[0], event.xdata],
                                        [last_point[1], event.ydata])
        else:
            self.preview_line, = self.ax.plot([last_point[0], event.xdata],
                                               [last_point[1], event.ydata],
                                               'r--', linewidth=1, alpha=0.5)
        self.canvas.draw_idle()
    
    def _clear_line(self):
        """Clear drawn line"""
        self.line_points = []
        self.points_label.config(text='Points: 0')
        
        if self.drawn_line:
            self.drawn_line.remove()
            self.drawn_line = None
        if self.preview_line:
            self.preview_line.remove()
            self.preview_line = None
        
        self.canvas.draw_idle()
    
    def _update_line_display(self):
        """Update line display based on current settings"""
        # Remove existing drawn line
        if self.drawn_line:
            self.drawn_line.remove()
            self.drawn_line = None
        
        if not self.show_line_var.get() or len(self.line_points) < 2:
            self.canvas.draw_idle()
            return
        
        points = np.array(self.line_points)
        x, y = points[:, 0], points[:, 1]
        
        # Apply smoothing if enabled
        if self.smooth_var.get() and len(points) >= 4:
            try:
                from scipy.interpolate import splprep, splev
                tck, u = splprep([x, y], s=0)
                u_new = np.linspace(0, 1, 100)
                x, y = splev(u_new, tck)
            except:
                pass  # Fall back to non-smooth
        
        color = self.line_color_var.get()
        width = int(self.line_width_var.get())
        
        self.drawn_line, = self.ax.plot(x, y, color=color, linewidth=width)
        self.canvas.draw_idle()
    
    def _get_year_from_shot(self, shot_number):
        """Get year from shot number (KSTAR shot ranges)"""
        if shot_number > 37741:
            return 2025
        elif shot_number > 34836:
            return 2024
        elif shot_number > 32768:
            return 2023
        elif shot_number > 30445:
            return 2022
        elif shot_number > 27400:
            return 2021
        elif shot_number > 24081:
            return 2020
        elif shot_number > 21758:
            return 2019
        elif shot_number > 19396:
            return 2018
        elif shot_number > 17376:
            return 2017
        elif shot_number > 14407:
            return 2016
        elif shot_number > 11724:
            return 2015
        elif shot_number > 9427:
            return 2014
        elif shot_number > 8354:
            return 2013
        elif shot_number > 6470:
            return 2012
        elif shot_number > 4468:
            return 2011
        elif shot_number > 2342:
            return 2010
        elif shot_number > 1283:
            return 2009
        else:
            return 2008
    
    def _get_campaign_from_year(self, year):
        """Get campaign folder name from year"""
        base_path = '/Diag_TV'
        
        if not os.path.exists(base_path):
            print(f"TV: Base path {base_path} not found")
            return None
        
        try:
            entries = os.listdir(base_path)
            campaign_dirs = [d for d in entries if d.startswith(f'{year}C')]
            
            if campaign_dirs:
                campaign_dirs.sort()
                return campaign_dirs[-1]
            return None
        except Exception as e:
            print(f"TV: Error listing campaigns: {str(e)}")
            return None
    
    def _find_available_tvs(self, shot_number):
        """Find available TV01/TV02 ZIP files for given shot"""
        year = self._get_year_from_shot(shot_number)
        campaign = self._get_campaign_from_year(year)
        
        if campaign is None:
            return []
        
        base_path = f'/Diag_TV/{campaign}'
        available_tvs = []
        shot_str = f'{shot_number:06d}'
        
        for tv_num in ['01', '02']:
            tv_path = f'{base_path}/TV{tv_num}'
            zip_filename = f'{shot_str}_tv{tv_num}.zip'
            full_path = f'{tv_path}/{zip_filename}'
            
            try:
                if os.path.exists(full_path):
                    available_tvs.append(f'TV{tv_num}')
                    print(f"TV: Found {full_path}")
            except Exception as e:
                print(f"TV: Error checking {full_path}: {str(e)}")
        
        return available_tvs
    
    def _search_available_tvs(self):
        """Search for available TV01/TV02 files for given shot"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return
        
        available_tvs = self._find_available_tvs(shot_number)
        
        if not available_tvs:
            messagebox.showinfo("Not Found", 
                f"No TV data found for shot #{shot_number}")
            self.tv_dropdown['values'] = []
            self.tv_selection_var.set('-- Select TV --')
            self.file_label.config(text="No TV files found")
            return
        
        # Update dropdown with available TVs
        self.tv_dropdown['values'] = available_tvs
        self.tv_selection_var.set(available_tvs[0])
        self.file_label.config(text=f"Found: {', '.join(available_tvs)} for #{shot_number}")
        print(f"TV: Found {available_tvs} for shot #{shot_number}")
    
    def _load_selected_tv(self):
        """Load selected TV ZIP file"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return
        
        tv_selection = self.tv_selection_var.get()
        
        if not tv_selection or tv_selection == '-- Select TV --':
            messagebox.showerror("Error", "Please search and select a TV first")
            return
        
        year = self._get_year_from_shot(shot_number)
        campaign = self._get_campaign_from_year(year)
        
        if campaign is None:
            messagebox.showerror("Error", "Campaign folder not found")
            return
        
        shot_str = f'{shot_number:06d}'
        tv_num = tv_selection.replace('TV', '')
        
        zip_path = f'/Diag_TV/{campaign}/{tv_selection}/{shot_str}_tv{tv_num}.zip'
        
        self._load_zip_from_path(zip_path)
    
    def _load_zip_file(self):
        """Open file dialog and load selected ZIP file"""
        initial_dir = '/Diag_TV' if os.path.exists('/Diag_TV') else None
        file_path = filedialog.askopenfilename(
            title="Select TV ZIP file",
            initialdir=initial_dir,
            filetypes=[("ZIP files", "*.zip"), ("All files", "*.*")]
        )
        
        if file_path:
            self._load_zip_from_path(file_path)
    
    def _set_status(self, text):
        """Update status label"""
        self.status_label.config(text=text)
        self.status_label.update()
    
    def _load_zip_from_path(self, file_path):
        """Load images from ZIP file"""
        # Stop any playback
        self._stop_play()
        self._stop_prefetch()
        
        # Step 1: Close existing ZIP file
        if self.zip_file:
            try:
                self.zip_file.close()
            except:
                pass
        
        self.cache.clear()
        self.im = None  # Reset image reference
        
        # Step 2: Open ZIP
        self._set_status("Opening ZIP file...")
        print(f"TV: Opening ZIP file {file_path}")
        
        try:
            self.zip_file = zipfile.ZipFile(file_path, 'r')
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open ZIP file:\n{str(e)}")
            self._set_status("Failed to open ZIP")
            print(f"TV: Failed to open ZIP: {str(e)}")
            return
        
        # Step 3: Get file list
        self._set_status("Reading file list...")
        print(f"TV: Reading file list...")
        
        all_files = self.zip_file.namelist()
        
        # Step 4: Filter image files
        self._set_status("Filtering image files...")
        print(f"TV: Filtering image files from {len(all_files)} entries...")
        
        image_extensions = ('.png', '.jpg', '.jpeg', '.tif', '.tiff', '.bmp')
        image_files = [f for f in all_files 
                      if f.lower().endswith(image_extensions) and not f.startswith('__MACOSX')]
        
        if not image_files:
            messagebox.showerror("Error", "No image files found in ZIP")
            self._set_status("No images found")
            print(f"TV: No image files found")
            return
        
        # Sort files naturally
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower()
                   for text in re.split('([0-9]+)', s)]
        
        image_files.sort(key=natural_sort_key)
        
        self.image_files = image_files
        self.total_frames = len(image_files)
        
        print(f"TV: Found {self.total_frames} image files")
        
        # Step 5: Update UI
        self._set_status("Updating UI...")
        self.file_label.config(text=file_path.split('/')[-1])
        self.frame_slider.config(to=self.total_frames - 1)
        self.frame_total_label.config(text=f'/ {self.total_frames}')
        
        # Step 6: Load first frame
        self._set_status("Loading first frame...")
        print(f"TV: Loading first frame...")
        
        self.current_frame = 0
        self.frame_var.set(0)
        
        if not self._display_frame(0):
            self._set_status("Failed to load first frame")
            print(f"TV: Failed to load first frame")
            messagebox.showerror("Error", "Failed to load first frame from ZIP")
            return
        
        # Start prefetching nearby frames
        self._start_prefetch(0)
        
        self._set_status("Ready")
        print(f"TV: Successfully loaded {self.total_frames} images")
    
    def _get_image(self, frame_idx):
        """Get image from cache or load from ZIP"""
        if frame_idx < 0 or frame_idx >= self.total_frames:
            return None
        
        if self.zip_file is None:
            return None
        
        # Check cache
        with self._prefetch_lock:
            if frame_idx in self.cache:
                return self.cache[frame_idx]
        
        # Load from ZIP
        try:
            filename = self.image_files[frame_idx]
            with self.zip_file.open(filename) as f:
                img_data = f.read()
                img = Image.open(io.BytesIO(img_data))
                img_array = np.array(img)
            
            # Add to cache
            with self._prefetch_lock:
                self.cache[frame_idx] = img_array
                self._cleanup_cache(frame_idx)
            
            return img_array
            
        except Exception as e:
            print(f"TV: Error loading frame {frame_idx}: {str(e)}")
            return None
    
    def _cleanup_cache(self, center_frame):
        """Remove frames far from current position"""
        if len(self.cache) > self.cache_size:
            keys_to_remove = []
            for key in self.cache.keys():
                if abs(key - center_frame) > self.cache_size // 2:
                    keys_to_remove.append(key)
            for key in keys_to_remove[:len(self.cache) - self.cache_size]:
                del self.cache[key]
    
    def _start_prefetch(self, center_frame, direction=1):
        """Start background thread to prefetch frames"""
        self._stop_prefetch()
        
        self._prefetch_stop = False
        self._prefetch_thread = threading.Thread(
            target=self._prefetch_worker,
            args=(center_frame, direction),
            daemon=True
        )
        self._prefetch_thread.start()
    
    def _stop_prefetch(self):
        """Stop prefetch thread"""
        self._prefetch_stop = True
        if self._prefetch_thread and self._prefetch_thread.is_alive():
            self._prefetch_thread.join(timeout=0.1)
    
    def _prefetch_worker(self, center_frame, direction):
        """Background worker to prefetch frames"""
        # Prefetch frames in the playback direction
        prefetch_count = min(30, self.total_frames)
        
        for i in range(1, prefetch_count + 1):
            if self._prefetch_stop:
                break
            
            # Prefetch in playback direction first
            frame_idx = center_frame + (i * direction)
            if 0 <= frame_idx < self.total_frames:
                with self._prefetch_lock:
                    if frame_idx not in self.cache:
                        try:
                            filename = self.image_files[frame_idx]
                            with self.zip_file.open(filename) as f:
                                img_data = f.read()
                                img = Image.open(io.BytesIO(img_data))
                                img_array = np.array(img)
                            self.cache[frame_idx] = img_array
                        except:
                            pass
    
    def _display_frame(self, frame_idx, update_ui=True):
        """Display specified frame using set_data for speed"""
        img_array = self._get_image(frame_idx)
        
        if img_array is None:
            return False
        
        filename = self.image_files[frame_idx]
        
        # First frame: use imshow, subsequent frames: use set_data
        if self.im is None:
            self.ax.clear()
            self.im = self.ax.imshow(img_array, cmap='gray' if len(img_array.shape) == 2 else None)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            # Re-draw line after clear
            self.drawn_line = None
            self._update_line_display()
        else:
            self.im.set_data(img_array)
        
        self.ax.set_title(f"Frame {frame_idx + 1}/{self.total_frames}: {filename}")
        
        # Use blit-like approach for faster rendering
        self.canvas.draw_idle()
        self.canvas.flush_events()
        
        # Update internal state - this is the key fix for synchronization
        self.current_frame = frame_idx
        
        # Update UI elements (skip during fast playback for performance)
        if update_ui:
            self.frame_entry.delete(0, tk.END)
            self.frame_entry.insert(0, str(frame_idx + 1))
            self.filename_label.config(text=filename)
        
        return True
    
    def _on_slider_change(self, value):
        """Handle slider value change"""
        frame_idx = int(float(value))
        # Always sync current_frame with slider value
        if frame_idx != self.current_frame:
            self.current_frame = frame_idx  # Sync first
            self._display_frame(frame_idx)
            # Restart prefetch from new position
            self._start_prefetch(frame_idx)
    
    def _on_frame_entry(self, event):
        """Handle frame entry input"""
        self._goto_frame()
    
    def _on_mouse_wheel(self, event):
        """Handle mouse wheel event for frame navigation"""
        if self.total_frames == 0:
            return
        
        # Use frame_var.get() instead of current_frame for consistent state
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
            self._display_frame(new_frame)
    
    def _goto_frame(self):
        """Go to specified frame number"""
        try:
            frame_num = int(self.frame_entry.get())
            frame_idx = frame_num - 1
            
            if 0 <= frame_idx < self.total_frames:
                self.current_frame = frame_idx
                self.frame_var.set(frame_idx)
                self._display_frame(frame_idx)
                self._start_prefetch(frame_idx)
            else:
                messagebox.showwarning("Warning", 
                    f"Frame number must be between 1 and {self.total_frames}")
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid frame number")
    
    def _step_frame(self, delta):
        """Step forward/backward by delta frames"""
        if self.total_frames == 0:
            return
        # Use frame_var for consistent state
        current = self.frame_var.get()
        new_frame = current + delta
        new_frame = max(0, min(new_frame, self.total_frames - 1))
        
        self.current_frame = new_frame
        self.frame_var.set(new_frame)
        self._display_frame(new_frame)
    
    def _goto_first(self):
        """Go to first frame"""
        if self.total_frames == 0:
            return
        self.current_frame = 0
        self.frame_var.set(0)
        self._display_frame(0)
        self._start_prefetch(0)
    
    def _goto_last(self):
        """Go to last frame"""
        if self.total_frames == 0:
            return
        last_idx = self.total_frames - 1
        self.current_frame = last_idx
        self.frame_var.set(last_idx)
        self._display_frame(last_idx)
        self._start_prefetch(last_idx, direction=-1)
    
    def _toggle_play(self):
        """Toggle play/pause"""
        if self.total_frames == 0:
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
        self._fps_history = []  # Reset FPS history
        
        # Start prefetching in playback direction
        skip = int(self.skip_var.get())
        self._start_prefetch(self.current_frame, direction=skip)
        
        self._play_next_frame()
    
    def _pause_play(self):
        """Pause playback"""
        self.is_playing = False
        self.play_btn.config(text='Play')
        if self.play_job:
            self.frame.after_cancel(self.play_job)
            self.play_job = None
        
        # Reset actual FPS display
        self.actual_fps_label.config(text='Actual: -- FPS')
        
        # Update UI with current frame info after pause
        if self.total_frames > 0:
            self.frame_entry.delete(0, tk.END)
            self.frame_entry.insert(0, str(self.current_frame + 1))
            if self.current_frame < len(self.image_files):
                self.filename_label.config(text=self.image_files[self.current_frame])
    
    def _stop_play(self):
        """Stop playback and reset to first frame"""
        self._pause_play()
        if self.total_frames > 0:
            self.current_frame = 0
            self.frame_var.set(0)
            self._display_frame(0)
    
    def _play_next_frame(self):
        """Play next frame with speed multiplier"""
        if not self.is_playing:
            return
        
        skip = int(self.skip_var.get())
        
        # Use frame_var for consistent state
        current = self.frame_var.get()
        next_frame = current + skip
        
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
        
        # Calculate and display actual FPS (rolling average)
        if actual_frame_time > 0:
            instant_fps = 1.0 / actual_frame_time
            # Simple rolling average
            if not hasattr(self, '_fps_history'):
                self._fps_history = []
            self._fps_history.append(instant_fps)
            if len(self._fps_history) > 10:
                self._fps_history.pop(0)
            avg_fps = sum(self._fps_history) / len(self._fps_history)
            self.actual_fps_label.config(text=f'Actual: {avg_fps:.1f} FPS')
        
        # Update state and display
        self.current_frame = next_frame
        self.frame_var.set(next_frame)
        self._display_frame(next_frame, update_ui=False)
        
        # Calculate delay from speed multiplier (1x = 10 FPS = 100ms)
        speed_str = self.speed_var.get()
        if speed_str == 'Max':
            target_delay = 1  # Minimum delay for maximum speed
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 100  # 1x = 100ms = 10 FPS
            target_delay = max(1, int(base_delay / speed_mult))
        
        self.play_job = self.frame.after(target_delay, self._play_next_frame)
    
    def cleanup(self):
        """Cleanup resources when tab is closed"""
        self._stop_play()
        self._stop_prefetch()
        if self.zip_file:
            try:
                self.zip_file.close()
            except:
                pass
            self.zip_file = None
        self.cache.clear()