#!/usr/bin/python3.8

"""
Custom navigation toolbar for matplotlib with axis control buttons
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk


class QuietNavigationToolbar(NavigationToolbar2Tk):
    """
    Custom navigation toolbar with 300 DPI PNG/EPS save functionality
    """

    def set_message(self, s):
        """Safely update the message string"""
        if hasattr(self, 'message'):
            self.message.set(s)

    def save_figure(self, *args):
        """Save figure"""
        # Get initial directory
        initial_dir = os.path.expanduser("~")

        # Hide hidden files in file dialog (Linux Tk)
        try:
            self.master.tk.call('tk_getOpenFile', '-foption')
        except:
            pass
        self.master.tk.call('set', '::tk::dialog::file::showHiddenVar', '0')

        # Open file dialog to select base filename
        filepath = filedialog.asksaveasfilename(
            parent=self.master,
            initialdir=initial_dir,
            title="Save figure",
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("All files", "*.*")]
        )

        if not filepath:
            return

        # Remove extension if user added one
        base_path = os.path.splitext(filepath)[0]

        # Define output paths
        png_path = f"{base_path}.png"

        try:
            # Get the figure from canvas
            fig = self.canvas.figure

            # Save as PNG
            fig.savefig(png_path, dpi=300, bbox_inches='tight',
                       facecolor='white', edgecolor='none')

            # Show success message
            messagebox.showinfo(
                "Save Complete",
                f"Saved:\n\n{png_path}"
            )

        except Exception as e:
            messagebox.showerror("Save Error", f"Failed to save figure:\n{str(e)}")


class AxisControlToolbar(QuietNavigationToolbar):
    """
    Navigation toolbar with axis control button (Qt-style)
    Single button opens a dialog with axis selector dropdown
    """

    def __init__(self, canvas, window, tab_instance=None):
        """
        Initialize toolbar with axis control button

        Args:
            canvas: The matplotlib canvas
            window: The parent window/frame
            tab_instance: The tab instance that has ax1, ax2, and canvas attributes
        """
        super().__init__(canvas, window)
        self.tab_instance = tab_instance
        self.axes_config = []  # List of (name, getter_func) tuples
        self.share_x = True

        # Add separator and axis control button
        self._add_axis_button()

    def _add_axis_button(self):
        """Add single axis control button to the toolbar"""
        # Add a separator
        separator = ttk.Separator(self, orient='vertical')
        separator.pack(side=tk.LEFT, fill='y', padx=5, pady=2)

        # Single Axes button
        self.axes_btn = tk.Button(self, text='Axes', width=5, relief='raised',
                                  command=self._show_axes_dialog)
        self.axes_btn.pack(side=tk.LEFT, padx=1)

    def configure_axes(self, has_y2=True, share_x=True, x_label='X', y1_label='Y1', y2_label='Y2'):
        """
        Configure available axes for the dialog

        Args:
            has_y2: Whether Y2 axis exists
            share_x: Whether X-axis is shared between ax1 and ax2
            x_label: Display label for X axis
            y1_label: Display label for Y1 axis
            y2_label: Display label for Y2 axis
        """
        self.share_x = share_x
        self.axes_config = [
            (f'X ({x_label})', 'x', 'ax1'),
            (f'Y1 ({y1_label})', 'y', 'ax1'),
        ]
        if has_y2:
            self.axes_config.append((f'Y2 ({y2_label})', 'y', 'ax2'))

    def _get_axis_and_limits(self, axis_type, ax_name):
        """Get axis object and current limits"""
        if self.tab_instance is None:
            return None, (0, 1)

        # Get the axis object
        ax = getattr(self.tab_instance, ax_name, None)
        if ax is None and ax_name == 'ax1':
            ax = getattr(self.tab_instance, 'ax', None)

        if ax is None:
            return None, (0, 1)

        # Get limits based on axis type
        if axis_type == 'x':
            return ax, ax.get_xlim()
        else:
            return ax, ax.get_ylim()

    def _show_axes_dialog(self):
        """Show Qt-style axes configuration dialog"""
        if self.tab_instance is None:
            return

        # Build axes config if not set
        if not self.axes_config:
            self.configure_axes()

        dialog = tk.Toplevel(self.master)
        dialog.title('Axis Limits')
        dialog.geometry('280x160')
        dialog.resizable(False, False)
        dialog.transient(self.master)
        dialog.grab_set()

        # Center the dialog
        dialog.update_idletasks()
        x = self.master.winfo_rootx() + 100
        y = self.master.winfo_rooty() + 100
        dialog.geometry(f'+{x}+{y}')

        # Frame for inputs
        frame = ttk.Frame(dialog, padding=10)
        frame.pack(fill='both', expand=True)

        # Axis selector dropdown
        ttk.Label(frame, text='Axis:').grid(row=0, column=0, padx=5, pady=5, sticky='e')
        axis_names = [cfg[0] for cfg in self.axes_config]
        axis_var = tk.StringVar(value=axis_names[0] if axis_names else '')
        axis_combo = ttk.Combobox(frame, textvariable=axis_var, values=axis_names,
                                  state='readonly', width=20)
        axis_combo.grid(row=0, column=1, padx=5, pady=5, sticky='w')

        # Min entry
        ttk.Label(frame, text='Min:').grid(row=1, column=0, padx=5, pady=5, sticky='e')
        min_entry = ttk.Entry(frame, width=18)
        min_entry.grid(row=1, column=1, padx=5, pady=5, sticky='w')

        # Max entry
        ttk.Label(frame, text='Max:').grid(row=2, column=0, padx=5, pady=5, sticky='e')
        max_entry = ttk.Entry(frame, width=18)
        max_entry.grid(row=2, column=1, padx=5, pady=5, sticky='w')

        def update_entries(*args):
            """Update min/max entries when axis selection changes"""
            selected = axis_var.get()
            for name, axis_type, ax_name in self.axes_config:
                if name == selected:
                    ax, limits = self._get_axis_and_limits(axis_type, ax_name)
                    min_entry.delete(0, tk.END)
                    max_entry.delete(0, tk.END)
                    min_entry.insert(0, f'{limits[0]:.6g}')
                    max_entry.insert(0, f'{limits[1]:.6g}')
                    break

        # Bind axis selection change
        axis_combo.bind('<<ComboboxSelected>>', update_entries)

        # Initialize with first axis
        update_entries()

        # Button frame
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=3, column=0, columnspan=2, pady=10)

        def on_apply():
            """Apply the current axis limits"""
            try:
                min_val = float(min_entry.get()) if min_entry.get().strip() else None
                max_val = float(max_entry.get()) if max_entry.get().strip() else None

                if min_val is not None and max_val is not None and min_val >= max_val:
                    messagebox.showerror("Error", "Min must be less than Max", parent=dialog)
                    return

                selected = axis_var.get()
                for name, axis_type, ax_name in self.axes_config:
                    if name == selected:
                        self._apply_limits(axis_type, ax_name, min_val, max_val)
                        break

                # Update entries to show applied values
                update_entries()

            except ValueError:
                messagebox.showerror("Error", "Please enter valid numbers", parent=dialog)

        def on_auto():
            """Auto-scale the selected axis"""
            selected = axis_var.get()
            for name, axis_type, ax_name in self.axes_config:
                if name == selected:
                    self._apply_limits(axis_type, ax_name, None, None, auto=True)
                    break
            # Update entries to show new values
            update_entries()

        ttk.Button(btn_frame, text='Apply', command=on_apply, width=8).pack(side=tk.LEFT, padx=3)
        ttk.Button(btn_frame, text='Auto', command=on_auto, width=8).pack(side=tk.LEFT, padx=3)
        ttk.Button(btn_frame, text='Close', command=dialog.destroy, width=8).pack(side=tk.LEFT, padx=3)

        # Bind Enter key
        min_entry.bind('<Return>', lambda e: on_apply())
        max_entry.bind('<Return>', lambda e: on_apply())

        # Focus on min entry
        min_entry.focus_set()
        min_entry.select_range(0, tk.END)

    def _apply_limits(self, axis_type, ax_name, min_val, max_val, auto=False):
        """Apply limits to specified axis"""
        if self.tab_instance is None:
            return

        canvas = getattr(self.tab_instance, 'canvas', None)
        if canvas is None:
            return

        # Get the axis object
        ax = getattr(self.tab_instance, ax_name, None)
        if ax is None and ax_name == 'ax1':
            ax = getattr(self.tab_instance, 'ax', None)

        if ax is None:
            return

        if axis_type == 'x':
            if auto:
                ax.autoscale(axis='x')
                # Also apply to ax2 if share_x is True
                if self.share_x:
                    ax2 = getattr(self.tab_instance, 'ax2', None)
                    if ax2 is not None:
                        ax2.autoscale(axis='x')
            else:
                ax.set_xlim(left=min_val, right=max_val)
                # Also apply to ax2 if share_x is True
                if self.share_x:
                    ax2 = getattr(self.tab_instance, 'ax2', None)
                    if ax2 is not None:
                        ax2.set_xlim(left=min_val, right=max_val)
        else:  # y axis
            if auto:
                ax.autoscale(axis='y')
            else:
                ax.set_ylim(bottom=min_val, top=max_val)

        canvas.draw_idle()
        self.push_current()
