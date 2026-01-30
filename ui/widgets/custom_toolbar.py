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
    Navigation toolbar with axis control buttons (X, Y1, Y2)
    Replaces the separate Axis Control Panel in the UI
    """

    def __init__(self, canvas, window, tab_instance=None):
        """
        Initialize toolbar with axis control buttons

        Args:
            canvas: The matplotlib canvas
            window: The parent window/frame
            tab_instance: The tab instance that has ax1, ax2, and canvas attributes
        """
        super().__init__(canvas, window)
        self.tab_instance = tab_instance

        # Add separator and axis control buttons
        self._add_axis_buttons()

    def _add_axis_buttons(self):
        """Add axis control buttons to the toolbar"""
        # Add a separator
        separator = ttk.Separator(self, orient='vertical')
        separator.pack(side=tk.LEFT, fill='y', padx=5, pady=2)

        # X-axis button
        self.x_btn = tk.Button(self, text='X', width=3, relief='raised',
                               command=self._show_x_dialog)
        self.x_btn.pack(side=tk.LEFT, padx=1)

        # Y1-axis button
        self.y1_btn = tk.Button(self, text='Y1', width=3, relief='raised',
                                command=self._show_y1_dialog)
        self.y1_btn.pack(side=tk.LEFT, padx=1)

        # Y2-axis button (may be hidden if no second axis)
        self.y2_btn = tk.Button(self, text='Y2', width=3, relief='raised',
                                command=self._show_y2_dialog)
        self.y2_btn.pack(side=tk.LEFT, padx=1)

    def configure_axes(self, has_y2=True, x_label='X', y1_label='Y1', y2_label='Y2'):
        """
        Configure which axis buttons to show and their labels

        Args:
            has_y2: Whether to show Y2 button
            x_label: Label for X button tooltip
            y1_label: Label for Y1 button tooltip
            y2_label: Label for Y2 button tooltip
        """
        self.x_label = x_label
        self.y1_label = y1_label
        self.y2_label = y2_label

        if has_y2:
            self.y2_btn.pack(side=tk.LEFT, padx=1)
        else:
            self.y2_btn.pack_forget()

    def _show_x_dialog(self):
        """Show dialog for X-axis limits"""
        if self.tab_instance is None:
            return

        ax = getattr(self.tab_instance, 'ax1', None) or getattr(self.tab_instance, 'ax', None)
        if ax is None:
            return

        current_xlim = ax.get_xlim()
        self._show_axis_dialog('X-axis', current_xlim, self._apply_x_limits)

    def _show_y1_dialog(self):
        """Show dialog for Y1-axis limits"""
        if self.tab_instance is None:
            return

        ax = getattr(self.tab_instance, 'ax1', None) or getattr(self.tab_instance, 'ax', None)
        if ax is None:
            return

        current_ylim = ax.get_ylim()
        self._show_axis_dialog('Y1-axis', current_ylim, self._apply_y1_limits)

    def _show_y2_dialog(self):
        """Show dialog for Y2-axis limits"""
        if self.tab_instance is None:
            return

        ax2 = getattr(self.tab_instance, 'ax2', None)
        if ax2 is None:
            return

        current_ylim = ax2.get_ylim()
        self._show_axis_dialog('Y2-axis', current_ylim, self._apply_y2_limits)

    def _show_axis_dialog(self, title, current_limits, apply_callback):
        """Show a dialog for setting axis limits"""
        dialog = tk.Toplevel(self.master)
        dialog.title(f'Set {title} Limits')
        dialog.geometry('250x120')
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

        # Min entry
        ttk.Label(frame, text='Min:').grid(row=0, column=0, padx=5, pady=5, sticky='e')
        min_entry = ttk.Entry(frame, width=15)
        min_entry.grid(row=0, column=1, padx=5, pady=5)
        min_entry.insert(0, f'{current_limits[0]:.6g}')

        # Max entry
        ttk.Label(frame, text='Max:').grid(row=1, column=0, padx=5, pady=5, sticky='e')
        max_entry = ttk.Entry(frame, width=15)
        max_entry.grid(row=1, column=1, padx=5, pady=5)
        max_entry.insert(0, f'{current_limits[1]:.6g}')

        # Button frame
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=2, column=0, columnspan=2, pady=10)

        def on_apply():
            try:
                min_val = float(min_entry.get()) if min_entry.get().strip() else None
                max_val = float(max_entry.get()) if max_entry.get().strip() else None

                if min_val is not None and max_val is not None and min_val >= max_val:
                    messagebox.showerror("Error", "Min must be less than Max", parent=dialog)
                    return

                apply_callback(min_val, max_val)
                dialog.destroy()
            except ValueError:
                messagebox.showerror("Error", "Please enter valid numbers", parent=dialog)

        def on_auto():
            """Auto-scale the axis"""
            apply_callback(None, None, auto=True)
            dialog.destroy()

        ttk.Button(btn_frame, text='Apply', command=on_apply, width=8).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Auto', command=on_auto, width=8).pack(side=tk.LEFT, padx=5)
        ttk.Button(btn_frame, text='Cancel', command=dialog.destroy, width=8).pack(side=tk.LEFT, padx=5)

        # Bind Enter key
        min_entry.bind('<Return>', lambda e: on_apply())
        max_entry.bind('<Return>', lambda e: on_apply())

        # Focus on min entry
        min_entry.focus_set()
        min_entry.select_range(0, tk.END)

    def _apply_x_limits(self, min_val, max_val, auto=False):
        """Apply X-axis limits"""
        if self.tab_instance is None:
            return

        ax1 = getattr(self.tab_instance, 'ax1', None) or getattr(self.tab_instance, 'ax', None)
        ax2 = getattr(self.tab_instance, 'ax2', None)
        canvas = getattr(self.tab_instance, 'canvas', None)

        if ax1 is None or canvas is None:
            return

        if auto:
            ax1.autoscale(axis='x')
            if ax2 is not None:
                ax2.autoscale(axis='x')
        else:
            ax1.set_xlim(left=min_val, right=max_val)
            if ax2 is not None:
                ax2.set_xlim(left=min_val, right=max_val)

        canvas.draw_idle()
        self.push_current()

    def _apply_y1_limits(self, min_val, max_val, auto=False):
        """Apply Y1-axis limits"""
        if self.tab_instance is None:
            return

        ax = getattr(self.tab_instance, 'ax1', None) or getattr(self.tab_instance, 'ax', None)
        canvas = getattr(self.tab_instance, 'canvas', None)

        if ax is None or canvas is None:
            return

        if auto:
            ax.autoscale(axis='y')
        else:
            ax.set_ylim(bottom=min_val, top=max_val)

        canvas.draw_idle()
        self.push_current()

    def _apply_y2_limits(self, min_val, max_val, auto=False):
        """Apply Y2-axis limits"""
        if self.tab_instance is None:
            return

        ax2 = getattr(self.tab_instance, 'ax2', None)
        canvas = getattr(self.tab_instance, 'canvas', None)

        if ax2 is None or canvas is None:
            return

        if auto:
            ax2.autoscale(axis='y')
        else:
            ax2.set_ylim(bottom=min_val, top=max_val)

        canvas.draw_idle()
        self.push_current()