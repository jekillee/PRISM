#!/usr/bin/python3.8

"""
TV Startup Comparison tab for comparing startup sequences across multiple shots
Supports both TV01 and TV02 camera selection
Based on standalone TV01 startup montage viewer
"""

import tkinter as tk
from tkinter import ttk, messagebox
import zipfile
import os
import io
import numpy as np
from PIL import Image, ImageDraw, ImageFont

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from ui.ui_constants import (
    CONTROL_PANEL_WIDTH, PAD_X, PAD_Y,
    ENTRY_WIDTH_SHOT, BUTTON_WIDTH_MEDIUM
)
from ui.tv_utils import (
    TV_FPS, get_tv_startup_zip_path, find_available_tvs, frame_to_time_ms
)
from config.user_settings import get_tab_settings, set_tab_settings


class TVStartupTab:
    """TV startup comparison viewer tab - compare startup sequences across shots"""

    # Label column width for consistent alignment
    LABEL_COLUMN_WIDTH = 90

    # Fixed image dimensions
    FRAME_HEIGHT = 320
    FRAME_WIDTH = 135
    LABEL_WIDTH = 110

    def __init__(self, parent, app_config, diagnostic_config):
        self.parent = parent
        self.app_config = app_config
        self.diag_config = diagnostic_config

        self.frame = ttk.Frame(parent)
        self.toolbar = None

        # Shot list storage (shot numbers only, images generated on Plot)
        self.shot_list = []  # List of shot numbers
        self.current_tv = None  # TV type locked when first shot added
        self.last_canvas = None  # Combined stacked canvas
        self.available_tvs = []  # Available TVs for current shot

        # Plot references
        self.ax = None
        self.im = None

        # Track last fetched shot for validation
        self.last_fetched_shot = None

    def create_widgets(self):
        """Create TV Startup tab widgets"""
        # Right: Control panel (pack first so it doesn't resize)
        control_frame = ttk.Frame(self.frame, width=CONTROL_PANEL_WIDTH)
        control_frame.pack(side=tk.RIGHT, fill='y', expand=False)
        control_frame.pack_propagate(False)

        self._create_shot_controls(control_frame)
        self._create_display_controls(control_frame)
        self._create_plot_controls(control_frame)

        # Left: Image display area (no scroll)
        self._create_display()

        # Load saved settings
        self.load_settings()

    def _create_shot_controls(self, parent):
        """Create shot input section"""
        frame = ttk.LabelFrame(parent, text="1. Shot Selection", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)

        # Shot entry with up/down buttons and Fetch button
        ttk.Label(frame, text='Shot', anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')

        shot_frame = ttk.Frame(frame)
        shot_frame.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='ew')

        self.shot_entry = ttk.Entry(shot_frame, width=10)
        self.shot_entry.pack(side=tk.LEFT, fill='x', expand=True)
        self.shot_entry.bind('<Return>', lambda e: self._fetch_available_tvs())
        self.shot_entry.bind('<KeyRelease>', self._on_shot_changed)

        ttk.Button(shot_frame, text='\u25B2', width=2,
                   command=lambda: self._adjust_shot(1)).pack(side=tk.LEFT, padx=(2, 0))
        ttk.Button(shot_frame, text='\u25BC', width=2,
                   command=lambda: self._adjust_shot(-1)).pack(side=tk.LEFT)
        ttk.Button(shot_frame, text='Fetch', width=6,
                   command=self._fetch_available_tvs).pack(side=tk.LEFT, padx=(5, 0))

        # TV Selection (TV01 / TV02) with Add button
        ttk.Label(frame, text='TV', anchor='w').grid(
            row=1, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')

        tv_frame = ttk.Frame(frame)
        tv_frame.grid(row=1, column=1, padx=PAD_X, pady=PAD_Y, sticky='w')

        self.tv_var = tk.StringVar(value='TV01')
        self.tv01_radio = ttk.Radiobutton(tv_frame, text='TV01', variable=self.tv_var,
                                           value='TV01', state='disabled')
        self.tv01_radio.pack(side=tk.LEFT)
        self.tv02_radio = ttk.Radiobutton(tv_frame, text='TV02', variable=self.tv_var,
                                           value='TV02', state='disabled')
        self.tv02_radio.pack(side=tk.LEFT, padx=(10, 0))

        self.add_btn = ttk.Button(tv_frame, text='Add', command=self._add_shot,
                                   width=5, state='disabled')
        self.add_btn.pack(side=tk.LEFT, padx=(15, 0))

        # Status label
        self.status_label = tk.Label(frame, text="Ready", fg='gray',
                                      font=('TkDefaultFont', 9, 'bold'))
        self.status_label.grid(row=3, column=0, columnspan=2, padx=PAD_X, pady=2, sticky='w')

    def _create_display_controls(self, parent):
        """Create display options section"""
        frame = ttk.LabelFrame(parent, text="2. Added Shots", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        # Shots list
        list_frame = ttk.Frame(frame)
        list_frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        # Listbox with scrollbar
        listbox_frame = ttk.Frame(list_frame)
        listbox_frame.pack(fill='x', pady=(2, 0))

        self.shots_listbox = tk.Listbox(listbox_frame, height=6, width=25)
        self.shots_listbox.pack(side=tk.LEFT, fill='x', expand=True)

        # Bind keyboard delete keys
        self.shots_listbox.bind('<Delete>', lambda e: self._remove_selected())
        self.shots_listbox.bind('<BackSpace>', lambda e: self._remove_selected())

        scrollbar = ttk.Scrollbar(listbox_frame, orient='vertical',
                                   command=self.shots_listbox.yview)
        scrollbar.pack(side=tk.RIGHT, fill='y')
        self.shots_listbox.config(yscrollcommand=scrollbar.set)

        # Button frame for Remove Selected and Clear All
        btn_frame = ttk.Frame(list_frame)
        btn_frame.pack(fill='x', pady=(5, 0))

        ttk.Button(btn_frame, text='Remove Selected',
                   command=self._remove_selected, width=14).pack(side=tk.LEFT, padx=(0, 5))
        ttk.Button(btn_frame, text='Clear All',
                   command=self._clear_all, width=10).pack(side=tk.LEFT)

        # Hint for keyboard delete
        hint_label = ttk.Label(list_frame, text="(Delete/Backspace to remove)",
                               font=('TkDefaultFont', 8), foreground='gray')
        hint_label.pack(anchor='w', pady=(2, 0))

    def _create_plot_controls(self, parent):
        """Create plot section"""
        frame = ttk.LabelFrame(parent, text="3. Plot", labelanchor="n")
        frame.pack(fill='x', padx=PAD_X, pady=PAD_Y)

        frame.grid_columnconfigure(0, minsize=self.LABEL_COLUMN_WIDTH)
        frame.grid_columnconfigure(1, weight=1)

        # Frame range (start/end)
        ttk.Label(frame, text='Frame Range', anchor='w').grid(
            row=0, column=0, padx=PAD_X, pady=PAD_Y, sticky='w')

        range_frame = ttk.Frame(frame)
        range_frame.grid(row=0, column=1, padx=PAD_X, pady=PAD_Y, sticky='w')

        self.frame_start_var = tk.StringVar(value='1')
        ttk.Entry(range_frame, textvariable=self.frame_start_var, width=6).pack(side=tk.LEFT)

        ttk.Label(range_frame, text='~').pack(side=tk.LEFT, padx=5)

        self.frame_end_var = tk.StringVar(value='24')
        ttk.Entry(range_frame, textvariable=self.frame_end_var, width=6).pack(side=tk.LEFT)

        # Plot button
        ttk.Button(frame, text='Plot', command=self._plot,
                   width=15).grid(row=1, column=0, columnspan=2, padx=PAD_X, pady=PAD_Y, sticky='ew')

    def _create_display(self):
        """Create matplotlib display area (no scroll)"""
        # Container frame
        container = ttk.Frame(self.frame)
        container.pack(side=tk.LEFT, fill='both', expand=True)

        # Matplotlib figure
        self.figure = Figure(figsize=(10, 6), tight_layout=True)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_axis_off()
        self.ax.text(0.5, 0.5, 'Add shots and click Plot to compare startup sequences',
                     ha='center', va='center', fontsize=12,
                     transform=self.ax.transAxes, color='gray')

        # Matplotlib widget
        self.canvas = FigureCanvasTkAgg(self.figure, master=container)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        self.canvas.draw()

    def _adjust_shot(self, delta):
        """Adjust shot number by delta"""
        try:
            current = int(self.shot_entry.get())
            new_shot = max(1, current + delta)
            self.shot_entry.delete(0, tk.END)
            self.shot_entry.insert(0, str(new_shot))
            # Disable Add and TV buttons when shot changes
            self._disable_add_controls()
        except ValueError:
            pass

    def _on_shot_changed(self, event=None):
        """Handle shot entry value change - disable Add if shot differs from last fetched"""
        try:
            current_shot = int(self.shot_entry.get().strip())
        except ValueError:
            current_shot = None

        # If shot changed from last fetched, disable Add button and TV radios
        if current_shot != self.last_fetched_shot:
            self._disable_add_controls()

    def _disable_add_controls(self):
        """Disable Add button and TV radio buttons"""
        self.add_btn.config(state='disabled')
        self.tv01_radio.config(state='disabled')
        self.tv02_radio.config(state='disabled')

    def _fetch_available_tvs(self):
        """Fetch available TV startup files for current shot"""
        try:
            shot = int(self.shot_entry.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return

        # Store the fetched shot number
        self.last_fetched_shot = shot

        # Use tv_utils to find available TVs
        self.available_tvs = find_available_tvs(shot, startup=True)

        # Update radio buttons based on availability and lock status
        if self.current_tv is not None:
            # TV type is locked - only enable if matches locked type and available
            if 'TV01' in self.available_tvs and self.current_tv == 'TV01':
                self.tv01_radio.config(state='normal')
                self.tv_var.set('TV01')
            else:
                self.tv01_radio.config(state='disabled')

            if 'TV02' in self.available_tvs and self.current_tv == 'TV02':
                self.tv02_radio.config(state='normal')
                self.tv_var.set('TV02')
            else:
                self.tv02_radio.config(state='disabled')
        else:
            # No lock - enable based on availability
            if 'TV01' in self.available_tvs:
                self.tv01_radio.config(state='normal')
            else:
                self.tv01_radio.config(state='disabled')

            if 'TV02' in self.available_tvs:
                self.tv02_radio.config(state='normal')
            else:
                self.tv02_radio.config(state='disabled')

            # Auto-select first available
            if self.available_tvs:
                self.tv_var.set(self.available_tvs[0])

        # Enable/disable Add button
        if self.available_tvs:
            # Check if selected TV is available (considering lock)
            selected_tv = self.tv_var.get()
            if self.current_tv is not None:
                if selected_tv == self.current_tv and selected_tv in self.available_tvs:
                    self.add_btn.config(state='normal')
                else:
                    self.add_btn.config(state='disabled')
            else:
                if selected_tv in self.available_tvs:
                    self.add_btn.config(state='normal')
                else:
                    self.add_btn.config(state='disabled')

            self._set_status(f"Found: {', '.join(self.available_tvs)} for #{shot}", 'green')
            print(f"[TV Startup] Found {self.available_tvs} for shot #{shot}")
        else:
            self.add_btn.config(state='disabled')
            self._set_status(f"No TV startup data for #{shot}", 'red')
            print(f"[TV Startup] No TV startup data for shot #{shot}")

    def _add_time_label(self, canvas, frame_start, frame_end):
        """Add time labels to the top of canvas"""
        h, w, _ = canvas.shape
        ncol = frame_end - frame_start + 1

        try:
            font = ImageFont.truetype("DejaVuSans.ttf", size=25)
        except Exception:
            font = ImageFont.load_default()

        img = Image.fromarray((canvas * 255).astype(np.uint8), mode='RGB')
        draw = ImageDraw.Draw(img)

        # Add label every 4 frames
        for i in range(0, ncol, 4):
            frame_num = frame_start + i
            time_ms = frame_to_time_ms(frame_num)
            text = f'{time_ms:3.0f} ms'
            draw.text((i * self.FRAME_WIDTH + 10, 30), text, fill=(255, 255, 255), font=font)

        return np.asarray(img, dtype=np.float32) / 255.0

    def _add_shot_label(self, canvas, shot):
        """Add shot number label to left side of canvas"""
        h, w, _ = canvas.shape

        img = Image.new('RGB', (self.LABEL_WIDTH, h), (0, 0, 0))
        draw = ImageDraw.Draw(img)

        text_shot = str(shot)

        try:
            font_shot = ImageFont.truetype("DejaVuSans.ttf", size=max(16, h // 10))
        except Exception:
            font_shot = ImageFont.load_default()

        # Shot number position
        bbox_shot = draw.textbbox((0, 0), text_shot, font=font_shot)
        tw_shot = bbox_shot[2] - bbox_shot[0]
        th_shot = bbox_shot[3] - bbox_shot[1]

        x_shot = (self.LABEL_WIDTH - tw_shot) // 2
        y_shot = (h - th_shot) // 2

        draw.text((x_shot, y_shot), text_shot, fill=(255, 255, 255), font=font_shot)

        label_arr = np.asarray(img, dtype=np.float32) / 255.0
        return np.concatenate([label_arr, canvas], axis=1)

    def _load_startup_images(self, shot, tv_name, frame_start, frame_end):
        """Load startup images from ZIP file for specified frame range"""
        ncol = frame_end - frame_start + 1
        canvas = np.zeros((self.FRAME_HEIGHT, self.FRAME_WIDTH * ncol, 3), dtype=float)

        zip_name = get_tv_startup_zip_path(shot, tv_name)
        if not os.path.exists(zip_name):
            raise FileNotFoundError(f"ZIP not found: {zip_name}")

        print(f"[TV Startup] Loading {zip_name} (frames {frame_start}-{frame_end})")

        # Calculate BMP index from frame number
        # Original: index0 = int((time_i * 1000 + 100) * TV_FPS / 1000)
        # For frame 1 at time 0ms: index = (0 + 100) * 210 / 1000 = 21
        base_index = 21  # Frame 1 corresponds to index 21

        with zipfile.ZipFile(zip_name, 'r') as zf:
            for i, frame_num in enumerate(range(frame_start, frame_end + 1)):
                bmp_index = base_index + (frame_num - 1)
                bmp_name = f"{shot:06d}-{bmp_index:05d}.bmp"

                try:
                    data = zf.read(bmp_name)
                except KeyError:
                    continue

                # Read image from memory
                im = Image.open(io.BytesIO(data))
                out = im.resize((320, 320))
                rsize_arr = np.asarray(out, dtype=np.float32) / 255.0

                # Crop center portion
                crop = rsize_arr[:, 75:210, :]

                c = i * self.FRAME_WIDTH
                canvas[:, c:c + self.FRAME_WIDTH, :] = crop

        return canvas

    def _add_shot(self):
        """Add shot to list"""
        try:
            shot = int(self.shot_entry.get().strip())
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid shot number")
            return

        tv_name = self.tv_var.get()

        # Check if TV type is locked
        if self.current_tv is not None and self.current_tv != tv_name:
            messagebox.showwarning("Warning",
                f"Only {self.current_tv} shots can be added.\nClear all to change TV type.")
            return

        # Check if shot already added
        if shot in self.shot_list:
            messagebox.showwarning("Warning", f"Shot #{shot} is already added")
            return

        # Lock TV type on first shot
        if self.current_tv is None:
            self.current_tv = tv_name
            # Disable the other radio button
            if tv_name == 'TV01':
                self.tv02_radio.config(state='disabled')
            else:
                self.tv01_radio.config(state='disabled')

        # Add to list
        self.shot_list.append(shot)
        self.shots_listbox.insert(tk.END, f"#{shot}")

        # Disable Add and TV buttons after adding
        self._disable_add_controls()

        self._set_status(f"Added #{shot}. Total: {len(self.shot_list)} shots", 'green')
        print(f"[TV Startup] Added shot #{shot} ({tv_name})")

    def _remove_selected(self):
        """Remove selected shot from list"""
        selection = self.shots_listbox.curselection()
        if not selection:
            return

        index = selection[0]
        shot = self.shot_list[index]

        # Remove from storage
        del self.shot_list[index]
        self.shots_listbox.delete(index)

        # Unlock TV type if list is empty
        if len(self.shot_list) == 0:
            self.current_tv = None
            # Re-enable based on last fetched availability
            if 'TV01' in self.available_tvs:
                self.tv01_radio.config(state='normal')
            if 'TV02' in self.available_tvs:
                self.tv02_radio.config(state='normal')
            self._set_status("Ready", 'gray')
        else:
            self._set_status(f"Removed #{shot}. Total: {len(self.shot_list)} shots", 'green')

        print(f"[TV Startup] Removed shot #{shot}")

    def _clear_all(self):
        """Clear all shots"""
        self.shot_list = []
        self.current_tv = None
        self.last_canvas = None
        self.shots_listbox.delete(0, tk.END)

        # Re-enable radio buttons based on last fetched availability
        if 'TV01' in self.available_tvs:
            self.tv01_radio.config(state='normal')
        if 'TV02' in self.available_tvs:
            self.tv02_radio.config(state='normal')

        # Reset display
        self.ax.clear()
        self.ax.set_axis_off()
        self.ax.text(0.5, 0.5, 'Add shots and click Plot to compare startup sequences',
                     ha='center', va='center', fontsize=12,
                     transform=self.ax.transAxes, color='gray')

        self.canvas.draw()
        self._set_status("Cleared", 'gray')

        print("[TV Startup] Cleared all shots")

    def _plot(self):
        """Generate and display the comparison plot"""
        if not self.shot_list:
            messagebox.showwarning("Warning", "Please add at least one shot")
            return

        # Validate frame range
        try:
            frame_start = int(self.frame_start_var.get().strip())
            frame_end = int(self.frame_end_var.get().strip())
            if frame_start < 1:
                frame_start = 1
            if frame_end < frame_start:
                raise ValueError("End frame must be >= start frame")
            # Limit to max 50 frames
            if frame_end - frame_start + 1 > 50:
                frame_end = frame_start + 49
                self.frame_end_var.set(str(frame_end))
                messagebox.showinfo("Info", "Frame range limited to 50 frames maximum")
        except ValueError as e:
            messagebox.showerror("Error", f"Invalid frame range: {e}")
            return

        ncol = frame_end - frame_start + 1
        tv_name = self.current_tv

        self._set_status("Loading images...", 'blue')
        self.frame.update_idletasks()

        images = []
        try:
            for i, shot in enumerate(self.shot_list):
                # Load images for this shot
                canvas = self._load_startup_images(shot, tv_name, frame_start, frame_end)

                # Add time labels only for first shot
                if i == 0:
                    canvas = self._add_time_label(canvas, frame_start, frame_end)

                # Add shot label
                canvas = self._add_shot_label(canvas, shot)

                images.append(canvas)

            # Stack all canvases vertically
            stacked = np.vstack(images)
            self.last_canvas = stacked

            # Update plot
            self.ax.clear()
            self.ax.imshow(stacked)
            self.ax.set_axis_off()

            self.figure.tight_layout()
            self.canvas.draw()

            self._set_status(f"Plotted {len(self.shot_list)} shots", 'green')
            print(f"[TV Startup] Plotted {len(self.shot_list)} shots")

        except FileNotFoundError as e:
            self._set_status("File not found", 'red')
            messagebox.showerror("Error", str(e))
        except Exception as e:
            self._set_status("Error", 'red')
            messagebox.showerror("Error", f"Failed to plot:\n{str(e)}")

    def _set_status(self, text, color='gray'):
        """Update status label with color"""
        self.status_label.config(text=text, fg=color)
        self.status_label.update()

    def cleanup(self):
        """Cleanup resources when tab is closed"""
        self.shot_list = []
        self.last_canvas = None

    def save_settings(self):
        """Save current tab settings"""
        settings = {
            "shot": self.shot_entry.get(),
            "frame_start": self.frame_start_var.get(),
            "frame_end": self.frame_end_var.get(),
            "tv": self.tv_var.get()
        }
        set_tab_settings("tv_startup", settings)

    def load_settings(self):
        """Load and apply saved settings"""
        settings = get_tab_settings("tv_startup")

        if settings.get("shot"):
            self.shot_entry.delete(0, tk.END)
            self.shot_entry.insert(0, settings["shot"])

        if settings.get("frame_start"):
            self.frame_start_var.set(settings["frame_start"])

        if settings.get("frame_end"):
            self.frame_end_var.set(settings["frame_end"])

        if settings.get("tv"):
            self.tv_var.set(settings["tv"])
