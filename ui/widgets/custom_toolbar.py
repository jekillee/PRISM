#!/usr/bin/python3.8

"""
Custom navigation toolbar for matplotlib with axis control buttons
"""

import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox, TclError
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
        except TclError:
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
    Navigation toolbar with axis control button (PyQt5 Figure Options style)
    Single button opens a dialog to select figure and adjust X/Y axes
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
        self.figures_config = []  # List of (display_name, ax_attr_name) tuples
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

    def configure_axes(self, has_y2=True, share_x=True, ax1_label='Axes 1', ax2_label='Axes 2'):
        """
        Configure available figures/axes for the dialog

        Args:
            has_y2: Whether second axes (ax2) exists
            share_x: Whether X-axis is shared between ax1 and ax2
            ax1_label: Display label for first axes
            ax2_label: Display label for second axes
        """
        self.share_x = share_x
        self.figures_config = [
            (ax1_label, 'ax1'),
        ]
        if has_y2:
            self.figures_config.append((ax2_label, 'ax2'))

    def _get_axes(self, ax_name):
        """Get axis object"""
        if self.tab_instance is None:
            return None

        ax = getattr(self.tab_instance, ax_name, None)
        if ax is None and ax_name == 'ax1':
            ax = getattr(self.tab_instance, 'ax', None)
        return ax

    def _show_axes_dialog(self):
        """Show PyQt5-style Figure Options dialog"""
        if self.tab_instance is None:
            return

        # Build figures config if not set
        if not self.figures_config:
            self.configure_axes()

        dialog = tk.Toplevel(self.master)
        dialog.title('Axis Limits')
        dialog.geometry('340x180')
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

        # Figure/Axes selector dropdown
        ttk.Label(frame, text='Axes:').grid(row=0, column=0, padx=5, pady=5, sticky='e')
        figure_names = [cfg[0] for cfg in self.figures_config]
        figure_var = tk.StringVar(value=figure_names[0] if figure_names else '')
        figure_combo = ttk.Combobox(frame, textvariable=figure_var, values=figure_names,
                                    state='readonly', width=25)
        figure_combo.grid(row=0, column=1, columnspan=3, padx=5, pady=5, sticky='w')

        # X-axis row
        ttk.Label(frame, text='X-axis').grid(row=1, column=0, padx=5, pady=5, sticky='e')
        ttk.Label(frame, text='Min:').grid(row=1, column=1, padx=(5, 2), pady=5, sticky='e')
        x_min_entry = ttk.Entry(frame, width=10)
        x_min_entry.grid(row=1, column=2, padx=2, pady=5)
        ttk.Label(frame, text='Max:').grid(row=1, column=3, padx=(5, 2), pady=5, sticky='e')
        x_max_entry = ttk.Entry(frame, width=10)
        x_max_entry.grid(row=1, column=4, padx=2, pady=5)

        # Y-axis row
        ttk.Label(frame, text='Y-axis').grid(row=2, column=0, padx=5, pady=5, sticky='e')
        ttk.Label(frame, text='Min:').grid(row=2, column=1, padx=(5, 2), pady=5, sticky='e')
        y_min_entry = ttk.Entry(frame, width=10)
        y_min_entry.grid(row=2, column=2, padx=2, pady=5)
        ttk.Label(frame, text='Max:').grid(row=2, column=3, padx=(5, 2), pady=5, sticky='e')
        y_max_entry = ttk.Entry(frame, width=10)
        y_max_entry.grid(row=2, column=4, padx=2, pady=5)

        def get_current_ax_name():
            """Get the ax attribute name for current selection"""
            selected = figure_var.get()
            for name, ax_name in self.figures_config:
                if name == selected:
                    return ax_name
            return 'ax1'

        def update_entries(*args):
            """Update X/Y entries when figure selection changes"""
            ax_name = get_current_ax_name()
            ax = self._get_axes(ax_name)
            if ax is None:
                return

            xlim = ax.get_xlim()
            ylim = ax.get_ylim()

            x_min_entry.delete(0, tk.END)
            x_max_entry.delete(0, tk.END)
            y_min_entry.delete(0, tk.END)
            y_max_entry.delete(0, tk.END)

            x_min_entry.insert(0, f'{xlim[0]:.6g}')
            x_max_entry.insert(0, f'{xlim[1]:.6g}')
            y_min_entry.insert(0, f'{ylim[0]:.6g}')
            y_max_entry.insert(0, f'{ylim[1]:.6g}')

        # Bind figure selection change
        figure_combo.bind('<<ComboboxSelected>>', update_entries)

        # Initialize with first figure
        update_entries()

        # Button frame
        btn_frame = ttk.Frame(frame)
        btn_frame.grid(row=3, column=0, columnspan=5, pady=10)

        def on_apply():
            """Apply X and Y axis limits to current figure"""
            try:
                ax_name = get_current_ax_name()
                ax = self._get_axes(ax_name)
                if ax is None:
                    return

                canvas = getattr(self.tab_instance, 'canvas', None)
                if canvas is None:
                    return

                # Parse X limits
                x_min = float(x_min_entry.get()) if x_min_entry.get().strip() else None
                x_max = float(x_max_entry.get()) if x_max_entry.get().strip() else None

                # Parse Y limits
                y_min = float(y_min_entry.get()) if y_min_entry.get().strip() else None
                y_max = float(y_max_entry.get()) if y_max_entry.get().strip() else None

                # Validate
                if x_min is not None and x_max is not None and x_min >= x_max:
                    messagebox.showerror("Error", "X Min must be less than X Max", parent=dialog)
                    return
                if y_min is not None and y_max is not None and y_min >= y_max:
                    messagebox.showerror("Error", "Y Min must be less than Y Max", parent=dialog)
                    return

                # Apply X limits
                if x_min is not None and x_max is not None:
                    ax.set_xlim(left=x_min, right=x_max)
                    # Apply to other axes if share_x
                    if self.share_x:
                        for _, other_ax_name in self.figures_config:
                            if other_ax_name != ax_name:
                                other_ax = self._get_axes(other_ax_name)
                                if other_ax is not None:
                                    other_ax.set_xlim(left=x_min, right=x_max)

                # Apply Y limits
                if y_min is not None and y_max is not None:
                    ax.set_ylim(bottom=y_min, top=y_max)

                canvas.draw_idle()
                self.push_current()

                # Update entries to show applied values
                update_entries()

            except ValueError:
                messagebox.showerror("Error", "Please enter valid numbers", parent=dialog)

        def on_auto():
            """Auto-scale both axes for current figure"""
            ax_name = get_current_ax_name()
            ax = self._get_axes(ax_name)
            if ax is None:
                return

            canvas = getattr(self.tab_instance, 'canvas', None)
            if canvas is None:
                return

            ax.autoscale()

            # Apply X autoscale to other axes if share_x
            if self.share_x:
                xlim = ax.get_xlim()
                for _, other_ax_name in self.figures_config:
                    if other_ax_name != ax_name:
                        other_ax = self._get_axes(other_ax_name)
                        if other_ax is not None:
                            other_ax.set_xlim(xlim)

            canvas.draw_idle()
            self.push_current()

            # Update entries to show new values
            update_entries()

        ttk.Button(btn_frame, text='Apply', command=on_apply, width=8).pack(side=tk.LEFT, padx=3)
        ttk.Button(btn_frame, text='Auto', command=on_auto, width=8).pack(side=tk.LEFT, padx=3)
        ttk.Button(btn_frame, text='Close', command=dialog.destroy, width=8).pack(side=tk.LEFT, padx=3)

        # Bind Enter key to apply
        for entry in [x_min_entry, x_max_entry, y_min_entry, y_max_entry]:
            entry.bind('<Return>', lambda e: on_apply())

        # Focus on x_min entry
        x_min_entry.focus_set()
        x_min_entry.select_range(0, tk.END)
