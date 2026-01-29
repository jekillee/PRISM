#!/usr/bin/python3.8

"""
Custom navigation toolbar for matplotlib
"""

import os
from tkinter import filedialog, messagebox
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