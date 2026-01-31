#!/usr/bin/python3.8

"""
TV (Visible Camera) tab for viewing image sequences from ZIP files
With line drawing feature for paper figures and TV1/TV2 compare mode
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
from config.user_settings import get_tab_settings, set_tab_settings


class TVTab:
    """TV image sequence viewer tab with line drawing and compare mode"""

    # TV camera parameters
    TV_FPS = 210.0
    TV_OFFSET = 0.1  # 100ms trigger offset

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90
    
    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config
        
        self.frame = ttk.Frame(parent)
        self.toolbar = None
        
        # Single TV mode data storage
        self.zip_file = None
        self.image_files = []
        self.current_frame = 0
        self.total_frames = 0
        
        # Compare mode data storage
        self.compare_mode = False
        self.tv1_zip = None
        self.tv2_zip = None
        self.tv1_images = []
        self.tv2_images = []
        self.tv1_cache = {}
        self.tv2_cache = {}
        
        # Image cache (store nearby frames for smooth playback)
        self.cache = {}
        self.cache_size = 100
        
        # Playback state
        self.is_playing = False
        self.play_job = None
        self._last_frame_time = 0
        
        # Prefetch thread control
        self._prefetch_lock = threading.Lock()
        self._prefetch_thread = None
        self._prefetch_stop = False
        
        # Plot references
        self.im = None
        self.im1 = None  # For compare mode (TV01)
        self.im2 = None  # For compare mode (TV02)
        self.ax = None
        self.ax1 = None  # For compare mode
        self.ax2 = None  # For compare mode
        
        # Line drawing
        self.line_points = []
        self.drawn_lines = []  # List of (line_obj, ax) tuples for finalized lines
        self.current_line_obj = None  # Line object currently being drawn (not finalized)
        self.preview_line = None
        self.draw_mode = False
        self.click_cid = None
        self.motion_cid = None
        self._draw_background = None  # For blitting optimization
        self.current_draw_ax = None  # Track which axis current line is being drawn on
    
    def create_widgets(self):
        """Create TV tab widgets"""
        # Left: Image display area
        self.figure = Figure((8, 6), tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("No image loaded")

        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.LEFT, fill='both', expand=True)

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

        # Load saved settings
        self.load_settings()
    
    def _create_file_controls(self, parent):
        """Create file loading section"""
        frame = ttk.LabelFrame(parent, text="1. Load TV Data", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)
        
        # Row 0: Shot label, entry with up/down, Search button, ... button
        ttk.Label(frame, text='Shot', anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')

        shot_frame = ttk.Frame(frame)
        shot_frame.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')

        self.shot_entry = ttk.Entry(shot_frame, width=10)
        self.shot_entry.pack(side=tk.LEFT, fill='x', expand=True)
        self.shot_entry.bind('<Return>', lambda e: self._search_available_tvs())

        ttk.Button(shot_frame, text='\u25B2', width=2,
                   command=lambda: self._adjust_shot(1)).pack(side=tk.LEFT, padx=(2, 0))
        ttk.Button(shot_frame, text='\u25BC', width=2,
                   command=lambda: self._adjust_shot(-1)).pack(side=tk.LEFT)

        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=0, column=2, padx=PAD_X, pady=PAD_Y, sticky='e')

        ttk.Button(btn_frame, text='Fetch', command=self._search_available_tvs, width=8).pack(
            side=tk.LEFT)
        ttk.Button(btn_frame, text='...', command=self._load_zip_file, width=3).pack(
            side=tk.LEFT, padx=(2, 0))
        
        # Row 1: TV dropdown and Load button
        ttk.Label(frame, text='TV', width=LABEL_WIDTH_SHORT, anchor='w').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        self.tv_selection_var = tk.StringVar(value='-- Select --')
        self.tv_dropdown = ttk.Combobox(frame, textvariable=self.tv_selection_var,
                                         state='readonly', width=12)
        self.tv_dropdown['values'] = []
        self.tv_dropdown.grid(row=1, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')
        
        ttk.Button(frame, text='Load', command=self._load_selected_tv, width=8).grid(
            row=1, column=2, padx=PAD_X, pady=PAD_Y, sticky='w')
        
        # File label
        self.file_label = ttk.Label(frame, text="No file loaded", wraplength=360)
        self.file_label.grid(row=2, column=0, columnspan=3, padx=PAD_X, pady=2, sticky='w')
        
        # Loading status label
        self.status_label = ttk.Label(frame, text="", foreground='blue')
        self.status_label.grid(row=3, column=0, columnspan=3, padx=PAD_X, pady=2, sticky='w')
    
    def _search_available_tvs(self):
        """Search for available TVs and update dropdown"""
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
            self.tv_selection_var.set('-- Select --')
            self.file_label.config(text="No TV files found")
            return
        
        # Build dropdown values based on availability
        dropdown_values = []
        has_both = 'TV01' in available_tvs and 'TV02' in available_tvs
        
        if has_both:
            dropdown_values.append('TV01 + TV02')
        if 'TV01' in available_tvs:
            dropdown_values.append('TV01')
        if 'TV02' in available_tvs:
            dropdown_values.append('TV02')
        
        self.tv_dropdown['values'] = dropdown_values
        
        # Set default: TV01 + TV02 if both available, otherwise first available
        if has_both:
            self.tv_selection_var.set('TV01 + TV02')
        else:
            self.tv_selection_var.set(dropdown_values[0])
        
        self.file_label.config(text=f"Found: {', '.join(available_tvs)} for #{shot_number}")
        print(f"TV: Found {available_tvs} for shot #{shot_number}")

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.get())
            new_shot = max(1, current + delta)
            self.shot_entry.delete(0, tk.END)
            self.shot_entry.insert(0, str(new_shot))
        except ValueError:
            pass

    def _create_frame_controls(self, parent):
        """Create frame navigation controls"""
        frame = ttk.LabelFrame(parent, text="2. Frame Control", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)

        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)

        row = 0

        # Frame slider with < > buttons
        slider_frame = ttk.Frame(frame)
        slider_frame.grid(row=row, column=0, columnspan=3, padx=5, pady=5, sticky='ew')

        ttk.Button(slider_frame, text='<', width=3, command=lambda: self._step_frame(-1)).pack(
            side=tk.LEFT, padx=(0, 2))

        self.frame_var = tk.IntVar(value=0)
        self.frame_slider = ttk.Scale(
            slider_frame, from_=0, to=0, orient='horizontal',
            variable=self.frame_var, command=self._on_slider_change
        )
        self.frame_slider.pack(side=tk.LEFT, fill='x', expand=True)

        ttk.Button(slider_frame, text='>', width=3, command=lambda: self._step_frame(1)).pack(
            side=tk.LEFT, padx=(2, 0))
        row += 1

        # Frame number input
        ttk.Label(frame, text='Frame', anchor='w').grid(
            row=row, column=0, padx=5, pady=2, sticky='w')

        frame_input = ttk.Frame(frame)
        frame_input.grid(row=row, column=1, columnspan=2, padx=5, pady=2, sticky='w')

        self.frame_entry = ttk.Entry(frame_input, width=8)
        self.frame_entry.pack(side=tk.LEFT)
        self.frame_entry.insert(0, '0')
        self.frame_entry.bind('<Return>', self._on_frame_entry)

        ttk.Label(frame_input, text='/').pack(side=tk.LEFT, padx=2)

        self.frame_total_entry = ttk.Entry(frame_input, width=8, state='readonly')
        self.frame_total_entry.pack(side=tk.LEFT)

        ttk.Button(frame_input, text='Go', width=5, command=self._goto_frame).pack(
            side=tk.LEFT, padx=(10, 0))
        row += 1

        # Time input
        ttk.Label(frame, text='Time [s]', anchor='w').grid(
            row=row, column=0, padx=5, pady=2, sticky='w')

        time_input = ttk.Frame(frame)
        time_input.grid(row=row, column=1, columnspan=2, padx=5, pady=2, sticky='w')

        self.time_entry = ttk.Entry(time_input, width=8)
        self.time_entry.pack(side=tk.LEFT)
        self.time_entry.insert(0, '0.0')
        self.time_entry.bind('<Return>', self._on_time_entry)

        ttk.Button(time_input, text='Go', width=5, command=self._goto_time).pack(
            side=tk.LEFT, padx=(10, 0))
        row += 1

        # Current filename display
        self.filename_label = ttk.Label(frame, text="", wraplength=320)
        self.filename_label.grid(row=row, column=0, columnspan=3, padx=5, pady=5, sticky='w')
        row += 1

        # Mouse wheel hint
        hint_label = ttk.Label(frame, text="(Mouse wheel: navigate frames)",
                               font=('TkDefaultFont', 8), foreground='gray')
        hint_label.grid(row=row, column=0, columnspan=3, pady=(0, 5))
    
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
    
    def _create_draw_line_controls(self, parent):
        """Create line drawing controls"""
        frame = ttk.LabelFrame(parent, text="4. Draw Line", labelanchor="n")
        frame.pack(fill='x', padx=5, pady=5)
        
        # Draw mode button and Clear button
        btn_frame = ttk.Frame(frame)
        btn_frame.pack(fill='x', padx=5, pady=5)

        self.draw_btn = tk.Button(btn_frame, text='Draw Mode: OFF',
                                  command=self._toggle_draw_mode, width=14)
        self.draw_btn.pack(side=tk.LEFT, expand=True, fill='x', padx=2)

        tk.Button(btn_frame, text='Clear', command=self._clear_line, width=8).pack(
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
        
        # Line style options (grid layout for even spacing)
        style_frame = ttk.Frame(frame)
        style_frame.pack(fill='x', padx=5, pady=5)
        for i in range(3):
            style_frame.columnconfigure(i, weight=1)

        # Color
        color_frame = ttk.Frame(style_frame)
        color_frame.grid(row=0, column=0, sticky='w')
        ttk.Label(color_frame, text='Color:').pack(side=tk.LEFT, padx=(0, 2))
        self.line_color_var = tk.StringVar(value='white')
        color_combo = ttk.Combobox(color_frame, textvariable=self.line_color_var,
                                    values=['white', 'black', 'red', 'blue', 'yellow', 'green'],
                                    state='readonly', width=6)
        color_combo.pack(side=tk.LEFT)
        color_combo.bind('<<ComboboxSelected>>', lambda e: self._update_line_display())

        # Linestyle
        linestyle_frame = ttk.Frame(style_frame)
        linestyle_frame.grid(row=0, column=1)
        ttk.Label(linestyle_frame, text='Linestyle:').pack(side=tk.LEFT, padx=(0, 2))
        self.line_style_var = tk.StringVar(value='dashed')
        style_combo = ttk.Combobox(linestyle_frame, textvariable=self.line_style_var,
                                    values=['dashed', 'solid', 'dotted'],
                                    state='readonly', width=6)
        style_combo.pack(side=tk.LEFT)
        style_combo.bind('<<ComboboxSelected>>', lambda e: self._update_line_display())

        # Width
        width_frame = ttk.Frame(style_frame)
        width_frame.grid(row=0, column=2, sticky='e')
        ttk.Label(width_frame, text='Width:').pack(side=tk.LEFT, padx=(0, 2))
        self.line_width_var = tk.StringVar(value='2')
        width_combo = ttk.Combobox(width_frame, textvariable=self.line_width_var,
                                    values=['1', '2', '3', '4', '5'],
                                    state='readonly', width=3)
        width_combo.pack(side=tk.LEFT)
        width_combo.bind('<<ComboboxSelected>>', lambda e: self._update_line_display())
        
        # Hint label
        hint_label = ttk.Label(frame, text="(Left-click: add point, Right-click: finish)", 
                               font=('TkDefaultFont', 8), foreground='gray')
        hint_label.pack(pady=(0, 5))
    
    def _toggle_draw_mode(self):
        """Toggle line drawing mode"""
        self.draw_mode = not self.draw_mode

        if self.draw_mode:
            self.draw_btn.config(text='Draw Mode: ON', bg='#90EE90', activebackground='#7CCD7C')
            self.click_cid = self.canvas.mpl_connect('button_press_event',
                                                      self._on_line_click)
            self.motion_cid = self.canvas.mpl_connect('motion_notify_event',
                                                       self._on_mouse_motion)
            self.canvas.get_tk_widget().config(cursor='crosshair')
            # Capture background for blitting optimization
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)
        else:
            self.draw_btn.config(text='Draw Mode: OFF', bg='#d9d9d9', activebackground='#ececec')
            if self.click_cid:
                self.canvas.mpl_disconnect(self.click_cid)
                self.click_cid = None
            if self.motion_cid:
                self.canvas.mpl_disconnect(self.motion_cid)
                self.motion_cid = None
            self.canvas.get_tk_widget().config(cursor='')
            # Clear blitting background
            self._draw_background = None
            # Remove preview line
            if self.preview_line:
                self.preview_line.remove()
                self.preview_line = None
                self.canvas.draw_idle()
    
    def _on_line_click(self, event):
        """Handle mouse click for line drawing"""
        if event.inaxes is None:
            return

        # Determine target axis
        if self.compare_mode:
            if event.inaxes == self.ax1:
                target_ax = self.ax1
            elif event.inaxes == self.ax2:
                target_ax = self.ax2
            else:
                return
        else:
            if event.inaxes != self.ax:
                return
            target_ax = self.ax

        if event.button == 1:  # Left click - add point
            # If switching to a different axis, finalize current line first
            if self.current_draw_ax is not None and self.current_draw_ax != target_ax:
                self._finalize_current_line()

            self.current_draw_ax = target_ax
            self.line_points.append((event.xdata, event.ydata))
            self.points_label.config(text=f'Points: {len(self.line_points)}')
            self._update_line_display()
        elif event.button == 3:  # Right click - finish current line
            self._finalize_current_line()
    
    def _on_mouse_motion(self, event):
        """Handle mouse motion for preview line"""
        if not self.draw_mode:
            return

        # Must have started drawing (current_draw_ax set)
        if self.current_draw_ax is None or len(self.line_points) == 0:
            return

        # Only show preview on the axis where we're drawing
        if event.inaxes != self.current_draw_ax:
            return

        # Update preview line from last point to cursor
        last_point = self.line_points[-1]

        if self.preview_line:
            self.preview_line.set_data([last_point[0], event.xdata],
                                        [last_point[1], event.ydata])
        else:
            self.preview_line, = self.current_draw_ax.plot([last_point[0], event.xdata],
                                                            [last_point[1], event.ydata],
                                                            'r--', linewidth=1, alpha=0.5)

        # Use blitting for smooth animation without flickering
        if self._draw_background is not None:
            self.canvas.restore_region(self._draw_background)
            self.current_draw_ax.draw_artist(self.preview_line)
            self.canvas.blit(self.figure.bbox)
        else:
            self.canvas.draw_idle()
    
    def _finalize_current_line(self):
        """Finalize the current line and prepare for a new one"""
        # Move current line to finalized lines list
        if self.current_line_obj is not None and self.current_draw_ax is not None:
            self.drawn_lines.append((self.current_line_obj, self.current_draw_ax))
            self.current_line_obj = None
        # Clear points for new line
        self.line_points = []
        self.current_draw_ax = None
        self.points_label.config(text='Points: 0')
        # Remove preview line
        if self.preview_line:
            self.preview_line.remove()
            self.preview_line = None
        # Update background for blitting
        self.canvas.draw()
        self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

    def _clear_line(self):
        """Clear all drawn lines"""
        self.line_points = []
        self.current_draw_ax = None
        self.points_label.config(text='Points: 0')

        # Remove current line being drawn
        if self.current_line_obj is not None:
            try:
                self.current_line_obj.remove()
            except ValueError:
                pass
            self.current_line_obj = None

        # Remove all finalized lines
        for line_obj, _ in self.drawn_lines:
            try:
                line_obj.remove()
            except ValueError:
                pass  # Already removed
        self.drawn_lines = []

        if self.preview_line:
            self.preview_line.remove()
            self.preview_line = None

        self.canvas.draw_idle()

        # Update background if in draw mode
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

    def _update_line_display(self):
        """Update line display based on current settings"""
        if not self.show_line_var.get() or len(self.line_points) < 2:
            self.canvas.draw_idle()
            return

        if self.current_draw_ax is None:
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
        style_map = {'dashed': '--', 'solid': '-', 'dotted': ':'}
        linestyle = style_map.get(self.line_style_var.get(), '--')

        # Remove previous version of current line if exists
        if self.current_line_obj is not None:
            try:
                self.current_line_obj.remove()
            except ValueError:
                pass
            self.current_line_obj = None

        # Draw new line
        self.current_line_obj, = self.current_draw_ax.plot(x, y, color=color, linewidth=width, linestyle=linestyle)
        self.canvas.draw_idle()

        # Update background for blitting (include newly drawn line)
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)
    
    # =========================================================================
    # Time/Frame conversion utilities
    # =========================================================================
    
    def _frame_to_time(self, frame_idx):
        """Convert frame index to time in seconds"""
        return (frame_idx + 1) / self.TV_FPS - self.TV_OFFSET
    
    def _time_to_frame(self, time_sec, total_frames):
        """Convert time in seconds to nearest frame index"""
        frame_idx = int(round((time_sec + self.TV_OFFSET) * self.TV_FPS - 1))
        return max(0, min(frame_idx, total_frames - 1))
    
    def _find_matching_frame(self, master_frame, master_total, slave_total):
        """Find slave frame that matches master frame's time. Returns None if out of range."""
        time_sec = self._frame_to_time(master_frame)
        
        # Calculate the time range of slave TV
        slave_min_time = self._frame_to_time(0)
        slave_max_time = self._frame_to_time(slave_total - 1)
        
        # Check if master time is within slave's time range
        if time_sec < slave_min_time or time_sec > slave_max_time:
            return None
        
        return self._time_to_frame(time_sec, slave_total)
    
    # =========================================================================
    # Shot/Campaign utilities
    # =========================================================================
    
    def _get_year_from_shot(self, shot_number):
        """Get year from shot number (KSTAR shot ranges)"""
        if shot_number >= 40464:
            return 2026
        elif shot_number > 37741:
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
    
    def _get_tv_zip_path(self, shot_number, tv_name):
        """Get full path to TV ZIP file"""
        year = self._get_year_from_shot(shot_number)
        campaign = self._get_campaign_from_year(year)
        
        if campaign is None:
            return None
        
        shot_str = f'{shot_number:06d}'
        tv_num = tv_name.replace('TV', '')
        return f'/Diag_TV/{campaign}/{tv_name}/{shot_str}_tv{tv_num}.zip'
    
    def _load_selected_tv(self):
        """Load selected TV ZIP file(s)"""
        try:
            shot_number = int(self.shot_entry.get())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return
        
        tv_selection = self.tv_selection_var.get()
        
        if not tv_selection or tv_selection == '-- Select --':
            messagebox.showerror("Error", "Please search and select a TV first")
            return
        
        # Handle TV01 + TV02 compare mode
        if tv_selection == 'TV01 + TV02':
            tv1_path = self._get_tv_zip_path(shot_number, 'TV01')
            tv2_path = self._get_tv_zip_path(shot_number, 'TV02')
            
            if tv1_path is None or tv2_path is None:
                messagebox.showerror("Error", "Failed to get TV file paths")
                return
            
            self._load_compare_mode(tv1_path, tv2_path)
        else:
            # Single TV mode
            zip_path = self._get_tv_zip_path(shot_number, tv_selection)
            
            if zip_path is None:
                messagebox.showerror("Error", "Campaign folder not found")
                return
            
            self.compare_mode = False
            self._setup_single_mode()
            self._load_zip_from_path(zip_path)
    
    def _setup_single_mode(self):
        """Setup figure for single TV display"""
        self._cleanup_compare_mode()
        
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        self.ax.set_xticks([])
        self.ax.set_yticks([])
        self.ax.set_title("No image loaded")
        
        self.im = None
        self.ax1 = None
        self.ax2 = None
        self.im1 = None
        self.im2 = None
        
        self.canvas.draw()
    
    def _setup_compare_mode(self):
        """Setup figure for TV1/TV2 comparison (side by side)"""
        self.figure.clear()
        
        # Left: TV02, Right: TV01
        self.ax2 = self.figure.add_subplot(121)
        self.ax1 = self.figure.add_subplot(122)
        
        self.ax1.set_xticks([])
        self.ax1.set_yticks([])
        self.ax2.set_xticks([])
        self.ax2.set_yticks([])
        
        self.ax1.set_title("TV01")
        self.ax2.set_title("TV02")
        
        self.im = None
        self.ax = None
        self.im1 = None
        self.im2 = None
        
        self.canvas.draw()
    
    def _cleanup_compare_mode(self):
        """Cleanup compare mode resources"""
        if self.tv1_zip:
            try:
                self.tv1_zip.close()
            except:
                pass
            self.tv1_zip = None
        
        if self.tv2_zip:
            try:
                self.tv2_zip.close()
            except:
                pass
            self.tv2_zip = None
        
        self.tv1_images = []
        self.tv2_images = []
        self.tv1_cache.clear()
        self.tv2_cache.clear()
    
    def _load_compare_mode(self, tv1_path, tv2_path):
        """Load both TV01 and TV02 for comparison"""
        self._stop_play()
        self._stop_prefetch()
        self._cleanup_compare_mode()
        
        # Close single mode zip if open
        if self.zip_file:
            try:
                self.zip_file.close()
            except:
                pass
            self.zip_file = None
        self.cache.clear()
        
        self._set_status("Loading TV01...")
        print(f"TV: Loading TV01 from {tv1_path}")
        
        try:
            self.tv1_zip = zipfile.ZipFile(tv1_path, 'r')
            self.tv1_images = self._get_sorted_images(self.tv1_zip)
            print(f"TV: TV01 has {len(self.tv1_images)} frames")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load TV01:\n{str(e)}")
            self._set_status("Failed")
            return
        
        self._set_status("Loading TV02...")
        print(f"TV: Loading TV02 from {tv2_path}")
        
        try:
            self.tv2_zip = zipfile.ZipFile(tv2_path, 'r')
            self.tv2_images = self._get_sorted_images(self.tv2_zip)
            print(f"TV: TV02 has {len(self.tv2_images)} frames")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load TV02:\n{str(e)}")
            self._set_status("Failed")
            self._cleanup_compare_mode()
            return
        
        if not self.tv1_images or not self.tv2_images:
            messagebox.showerror("Error", "No images found in one or both ZIP files")
            self._cleanup_compare_mode()
            return
        
        # Enable compare mode
        self.compare_mode = True
        self._setup_compare_mode()
        
        # Use TV01 as master for frame count
        self.total_frames = len(self.tv1_images)
        self.image_files = self.tv1_images  # For compatibility
        
        # Update UI
        self.file_label.config(text=f"TV01: {len(self.tv1_images)} frames, TV02: {len(self.tv2_images)} frames")
        self.frame_slider.config(to=self.total_frames - 1)
        self.frame_total_entry.config(state='normal')
        self.frame_total_entry.delete(0, tk.END)
        self.frame_total_entry.insert(0, str(self.total_frames))
        self.frame_total_entry.config(state='readonly')
        
        # Display first frame
        self.current_frame = 0
        self.frame_var.set(0)
        self._display_compare_frame(0)
        
        self._start_prefetch(0)
        self._set_status("Ready")
        print(f"TV: Compare mode ready")
    
    def _get_sorted_images(self, zip_file):
        """Get sorted list of image files from ZIP"""
        all_files = zip_file.namelist()
        
        image_extensions = ('.png', '.jpg', '.jpeg', '.tif', '.tiff', '.bmp')
        image_files = [f for f in all_files 
                      if f.lower().endswith(image_extensions) and not f.startswith('__MACOSX')]
        
        def natural_sort_key(s):
            return [int(text) if text.isdigit() else text.lower()
                   for text in re.split('([0-9]+)', s)]
        
        image_files.sort(key=natural_sort_key)
        return image_files
    
    def _load_zip_file(self):
        """Open file dialog and load selected ZIP file"""
        initial_dir = '/Diag_TV' if os.path.exists('/Diag_TV') else None
        file_path = filedialog.askopenfilename(
            title="Select TV ZIP file",
            initialdir=initial_dir,
            filetypes=[("ZIP files", "*.zip"), ("All files", "*.*")]
        )
        
        if file_path:
            self.compare_mode = False
            self._setup_single_mode()
            self._load_zip_from_path(file_path)
    
    def _set_status(self, text):
        """Update status label"""
        self.status_label.config(text=text)
        self.status_label.update()
    
    def _load_zip_from_path(self, file_path):
        """Load images from ZIP file (single TV mode)"""
        self._stop_play()
        self._stop_prefetch()
        self._cleanup_compare_mode()
        
        if self.zip_file:
            try:
                self.zip_file.close()
            except:
                pass
        
        self.cache.clear()
        self.im = None
        
        self._set_status("Opening ZIP file...")
        print(f"TV: Opening ZIP file {file_path}")
        
        try:
            self.zip_file = zipfile.ZipFile(file_path, 'r')
        except Exception as e:
            messagebox.showerror("Error", f"Failed to open ZIP file:\n{str(e)}")
            self._set_status("Failed to open ZIP")
            print(f"TV: Failed to open ZIP: {str(e)}")
            return
        
        self._set_status("Reading file list...")
        self.image_files = self._get_sorted_images(self.zip_file)
        
        if not self.image_files:
            messagebox.showerror("Error", "No image files found in ZIP")
            self._set_status("No images found")
            return
        
        self.total_frames = len(self.image_files)
        print(f"TV: Found {self.total_frames} image files")
        
        self._set_status("Updating UI...")
        self.file_label.config(text=file_path.split('/')[-1])
        self.frame_slider.config(to=self.total_frames - 1)
        self.frame_total_entry.config(state='normal')
        self.frame_total_entry.delete(0, tk.END)
        self.frame_total_entry.insert(0, str(self.total_frames))
        self.frame_total_entry.config(state='readonly')
        
        self._set_status("Loading first frame...")
        self.current_frame = 0
        self.frame_var.set(0)
        
        if not self._display_frame(0):
            self._set_status("Failed to load first frame")
            messagebox.showerror("Error", "Failed to load first frame from ZIP")
            return
        
        self._start_prefetch(0)
        self._set_status("Ready")
        print(f"TV: Successfully loaded {self.total_frames} images")
    
    # =========================================================================
    # Image loading and display
    # =========================================================================
    
    def _get_image(self, frame_idx):
        """Get image from cache or load from ZIP (single mode)"""
        if frame_idx < 0 or frame_idx >= self.total_frames:
            return None
        
        if self.zip_file is None:
            return None
        
        with self._prefetch_lock:
            if frame_idx in self.cache:
                return self.cache[frame_idx]
        
        try:
            filename = self.image_files[frame_idx]
            with self.zip_file.open(filename) as f:
                img_data = f.read()
                img = Image.open(io.BytesIO(img_data))
                img_array = np.array(img)
            
            with self._prefetch_lock:
                self.cache[frame_idx] = img_array
                self._cleanup_cache(frame_idx)
            
            return img_array
            
        except Exception as e:
            print(f"TV: Error loading frame {frame_idx}: {str(e)}")
            return None
    
    def _get_image_from_tv(self, tv_num, frame_idx):
        """Get image from specific TV cache or load from ZIP"""
        if tv_num == 1:
            zip_file = self.tv1_zip
            images = self.tv1_images
            cache = self.tv1_cache
        else:
            zip_file = self.tv2_zip
            images = self.tv2_images
            cache = self.tv2_cache
        
        if frame_idx < 0 or frame_idx >= len(images):
            return None
        
        if zip_file is None:
            return None
        
        with self._prefetch_lock:
            if frame_idx in cache:
                return cache[frame_idx]
        
        try:
            filename = images[frame_idx]
            with zip_file.open(filename) as f:
                img_data = f.read()
                img = Image.open(io.BytesIO(img_data))
                img_array = np.array(img)
            
            with self._prefetch_lock:
                cache[frame_idx] = img_array
                # Limit cache size
                if len(cache) > self.cache_size // 2:
                    keys = list(cache.keys())
                    for k in keys[:len(cache) - self.cache_size // 2]:
                        if k != frame_idx:
                            del cache[k]
            
            return img_array
            
        except Exception as e:
            print(f"TV: Error loading TV{tv_num} frame {frame_idx}: {str(e)}")
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
        prefetch_count = min(30, self.total_frames)
        
        for i in range(1, prefetch_count + 1):
            if self._prefetch_stop:
                break
            
            frame_idx = center_frame + (i * direction)
            if 0 <= frame_idx < self.total_frames:
                if self.compare_mode:
                    # Prefetch both TVs
                    self._get_image_from_tv(1, frame_idx)
                    tv2_frame = self._find_matching_frame(
                        frame_idx, len(self.tv1_images), len(self.tv2_images))
                    self._get_image_from_tv(2, tv2_frame)
                else:
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
        """Display specified frame (single TV mode)"""
        if self.compare_mode:
            return self._display_compare_frame(frame_idx, update_ui)
        
        img_array = self._get_image(frame_idx)
        
        if img_array is None:
            return False
        
        filename = self.image_files[frame_idx]
        
        if self.im is None:
            self.ax.clear()
            self.im = self.ax.imshow(img_array, cmap='gray' if len(img_array.shape) == 2 else None)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            # Clear lines associated with this axis
            self.drawn_lines = [(l, ax) for l, ax in self.drawn_lines if ax != self.ax]
            if self.current_draw_ax == self.ax:
                self.current_line_obj = None
        else:
            self.im.set_data(img_array)
        
        time_sec = self._frame_to_time(frame_idx)
        self.ax.set_title(f"Frame {frame_idx + 1}/{self.total_frames} (t = {time_sec:.3f} s)")
        
        self.canvas.draw_idle()
        self.canvas.flush_events()

        # Update blitting background if draw mode is active
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

        self.current_frame = frame_idx

        if update_ui:
            self.frame_entry.delete(0, tk.END)
            self.frame_entry.insert(0, str(frame_idx + 1))
            self.time_entry.delete(0, tk.END)
            self.time_entry.insert(0, f'{time_sec:.3f}')
            self.filename_label.config(text=filename)

        return True

    def _display_compare_frame(self, frame_idx, update_ui=True):
        """Display frames in compare mode (TV01 and TV02 side by side)"""
        # TV01 frame (master)
        tv1_img = self._get_image_from_tv(1, frame_idx)
        tv1_time = self._frame_to_time(frame_idx)
        
        # Find matching TV02 frame by time (returns None if out of range)
        tv2_frame = self._find_matching_frame(
            frame_idx, len(self.tv1_images), len(self.tv2_images))
        
        # Get TV02 image only if frame is valid
        tv2_img = None
        tv2_time = tv1_time  # Use same time for display
        if tv2_frame is not None:
            tv2_img = self._get_image_from_tv(2, tv2_frame)
            tv2_time = self._frame_to_time(tv2_frame)
        
        # Check if at least one has data
        if tv1_img is None and tv2_img is None:
            return False
        
        # Display TV01 (right)
        if tv1_img is not None:
            if self.im1 is None:
                self.ax1.clear()
                self.im1 = self.ax1.imshow(tv1_img, cmap='gray' if len(tv1_img.shape) == 2 else None)
                self.ax1.set_xticks([])
                self.ax1.set_yticks([])
            else:
                self.im1.set_data(tv1_img)
            self.ax1.set_title(f"TV01 Frame {frame_idx + 1}/{len(self.tv1_images)} (t = {tv1_time:.3f} s)")
        else:
            # No data for TV01 at this time
            self.ax1.clear()
            self.ax1.set_facecolor('black')
            self.ax1.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                         fontsize=16, color='white', transform=self.ax1.transAxes)
            self.ax1.set_xticks([])
            self.ax1.set_yticks([])
            self.ax1.set_title(f"TV01 (t = {tv1_time:.3f} s) - No Data")
            self.im1 = None
        
        # Display TV02 (left)
        if tv2_img is not None and tv2_frame is not None:
            if self.im2 is None:
                self.ax2.clear()
                self.im2 = self.ax2.imshow(tv2_img, cmap='gray' if len(tv2_img.shape) == 2 else None)
                self.ax2.set_xticks([])
                self.ax2.set_yticks([])
            else:
                self.im2.set_data(tv2_img)
            self.ax2.set_title(f"TV02 Frame {tv2_frame + 1}/{len(self.tv2_images)} (t = {tv2_time:.3f} s)")
        else:
            # No data for TV02 at this time
            self.ax2.clear()
            self.ax2.set_facecolor('black')
            self.ax2.text(0.5, 0.5, 'No Data', ha='center', va='center', 
                         fontsize=16, color='white', transform=self.ax2.transAxes)
            self.ax2.set_xticks([])
            self.ax2.set_yticks([])
            self.ax2.set_title(f"TV02 (t = {tv1_time:.3f} s) - No Data")
            self.im2 = None
        
        self.canvas.draw_idle()
        self.canvas.flush_events()

        # Update blitting background if draw mode is active
        if self.draw_mode:
            self.canvas.draw()
            self._draw_background = self.canvas.copy_from_bbox(self.figure.bbox)

        self.current_frame = frame_idx

        # Always update filename display
        tv1_name = self.tv1_images[frame_idx] if frame_idx < len(self.tv1_images) else "N/A"
        tv2_name = self.tv2_images[tv2_frame] if tv2_frame is not None and tv2_frame < len(self.tv2_images) else "N/A"
        self.filename_label.config(text=f"TV01: {tv1_name} | TV02: {tv2_name}")
        
        if update_ui:
            self.frame_entry.delete(0, tk.END)
            self.frame_entry.insert(0, str(frame_idx + 1))
            self.time_entry.delete(0, tk.END)
            self.time_entry.insert(0, f'{tv1_time:.3f}')
        
        return True
    
    # =========================================================================
    # Frame navigation
    # =========================================================================
    
    def _on_slider_change(self, value):
        """Handle slider value change with debounce"""
        # Cancel any pending slider update
        if hasattr(self, '_slider_job') and self._slider_job is not None:
            self.frame.after_cancel(self._slider_job)
        
        # Schedule update after short delay (debounce)
        self._slider_job = self.frame.after(10, self._do_slider_update)
    
    def _do_slider_update(self):
        """Actually perform the slider update"""
        self._slider_job = None
        frame_idx = self.frame_var.get()
        self._display_frame(frame_idx, update_ui=True)
    
    def _on_frame_entry(self, event):
        """Handle frame entry input"""
        self._goto_frame()
    
    def _on_mouse_wheel(self, event):
        """Handle mouse wheel event for frame navigation"""
        if self.total_frames == 0:
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
    
    def _on_time_entry(self, event):
        """Handle time entry input"""
        self._goto_time()
    
    def _goto_time(self):
        """Go to frame at specified time"""
        if self.total_frames == 0:
            return
        
        try:
            time_sec = float(self.time_entry.get())
            frame_idx = self._time_to_frame(time_sec, self.total_frames)
            
            self.current_frame = frame_idx
            self.frame_var.set(frame_idx)
            self._display_frame(frame_idx)
            self._start_prefetch(frame_idx)
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid time in seconds")
    
    def _step_frame(self, delta):
        """Step forward/backward by delta frames"""
        if self.total_frames == 0:
            return
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
    
    # =========================================================================
    # Playback controls
    # =========================================================================
    
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
        self._fps_history = []
        
        self._start_prefetch(self.current_frame, direction=1)
        
        self._play_next_frame()
    
    def _pause_play(self):
        """Pause playback"""
        self.is_playing = False
        self.play_btn.config(text='Play')
        if self.play_job:
            self.frame.after_cancel(self.play_job)
            self.play_job = None
        
        self.actual_fps_label.config(text='Actual: -- FPS')
        
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
        
        current = self.frame_var.get()
        next_frame = current + 1
        
        if next_frame >= self.total_frames:
            if self.loop_var.get():
                next_frame = 0
            else:
                self._pause_play()
                return
        
        now = time.time()
        actual_frame_time = now - self._last_frame_time
        self._last_frame_time = now
        
        if actual_frame_time > 0:
            instant_fps = 1.0 / actual_frame_time
            if not hasattr(self, '_fps_history'):
                self._fps_history = []
            self._fps_history.append(instant_fps)
            if len(self._fps_history) > 10:
                self._fps_history.pop(0)
            avg_fps = sum(self._fps_history) / len(self._fps_history)
            self.actual_fps_label.config(text=f'Actual: {avg_fps:.1f} FPS')
        
        self.current_frame = next_frame
        self.frame_var.set(next_frame)
        self._display_frame(next_frame, update_ui=False)
        
        speed_str = self.speed_var.get()
        if speed_str == 'Max':
            target_delay = 1
        else:
            speed_mult = float(speed_str.replace('x', ''))
            base_delay = 100
            target_delay = max(1, int(base_delay / speed_mult))
        
        self.play_job = self.frame.after(target_delay, self._play_next_frame)
    
    # =========================================================================
    # Cleanup and settings
    # =========================================================================
    
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
        
        self._cleanup_compare_mode()
        self.cache.clear()
    
    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.get()
        }
        set_tab_settings("tv", settings)
    
    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("tv")
        
        if settings.get("shot"):
            self.shot_entry.delete(0, tk.END)
            self.shot_entry.insert(0, settings["shot"])