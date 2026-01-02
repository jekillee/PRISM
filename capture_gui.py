#!/usr/bin/python3.8

"""
PRISM GUI Capture Tool for Paper Submission
Captures the initial PRISM window at 300 DPI for publication

Usage:
    python capture_gui.py [output_prefix]
    
Output:
    {output_prefix}_300dpi.png
    {output_prefix}_300dpi.eps
"""

import sys
import os
import tkinter as tk
from tkinter import ttk
from PIL import Image, ImageGrab

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from config.app_config import AppConfig, VERSION
from config.diagnostic_config import get_enabled_diagnostics
from data_loaders.efit_loader import EFITLoader
from plotting.plot_manager import PlotManager
from ui.tab_factory import TabFactory
from ui.widgets.custom_toolbar import QuietNavigationToolbar


class PRISMCaptureApp:
    """PRISM application modified for GUI capture"""
    
    def __init__(self, output_prefix='prism_gui'):
        self.output_prefix = output_prefix
        self.config = AppConfig()
        self.efit_loader = EFITLoader(self.config)
        self.plot_manager = PlotManager(self.config)
        
        self.root = tk.Tk()
        self.root.title(f'PRISM v{VERSION} - Plasma Research Integrated System for Multi-diagnostics')
        self.root.geometry("1500x700")
        
        self.tab_configs = []
        self.tab_cache = {}
        self.tab_frames = []
        
        self._create_widgets()
        self._create_bottom_bar()
        
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_changed)
        
        # Create first tab for initial display
        if self.tab_configs:
            self._create_tab_content(0)
            self._on_tab_changed(None)
        
        # Schedule capture after GUI is fully rendered
        self.root.after(1000, self._capture_and_save)
    
    def _create_widgets(self):
        """Create main application widgets with placeholder frames"""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(expand=True, fill='both')
        
        enabled_diagnostics = get_enabled_diagnostics()
        tab_types = ['profile', 'timetrace']
        
        for diag in enabled_diagnostics:
            for tab_type in tab_types:
                if not TabFactory.should_create_tab(diag, tab_type):
                    continue
                
                tab_name = TabFactory.get_tab_name(diag, tab_type)
                
                self.tab_configs.append({
                    'diagnostic': diag,
                    'tab_type': tab_type,
                    'tab_name': tab_name
                })
                
                placeholder = ttk.Frame(self.notebook)
                self.tab_frames.append(placeholder)
                self.notebook.add(placeholder, text=tab_name)
        
        # Spectrogram tab
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'spectrogram',
            'tab_name': TabFactory.get_tab_name(None, 'spectrogram')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='Spectrogram')

        # n-Mode Spectrum tab
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'nmode',
            'tab_name': TabFactory.get_tab_name(None, 'nmode')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='n-Mode Spectrum')
        
        # TV tab
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'tv',
            'tab_name': TabFactory.get_tab_name(None, 'tv')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='TV')
        
        # IRVB tab
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'irvb',
            'tab_name': TabFactory.get_tab_name(None, 'irvb')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='IRVB')    
        
    def _create_tab_content(self, tab_index):
        """Create actual tab content when first accessed"""
        if tab_index in self.tab_cache:
            return self.tab_cache[tab_index]
        
        config = self.tab_configs[tab_index]
        placeholder = self.tab_frames[tab_index]
        
        print(f"Creating tab: {config['tab_name']}...")
        
        tab = TabFactory.create_tab(
            placeholder,
            self.config,
            config['diagnostic'],
            config['tab_type'],
            self.efit_loader,
            self.plot_manager
        )
        
        if tab is not None:
            tab.frame.pack(expand=True, fill='both')
            self.tab_cache[tab_index] = tab
            print(f"  Tab created: {config['tab_name']}")
        
        return tab
    
    def _on_tab_changed(self, event):
        """Handle tab change - create tab if not exists, update toolbar"""
        current_tab_index = self.notebook.index(self.notebook.select())
        
        if current_tab_index not in self.tab_cache:
            self._create_tab_content(current_tab_index)
        
        current_tab = self.tab_cache.get(current_tab_index)
        
        # Destroy existing toolbar
        if hasattr(self, 'toolbar') and self.toolbar:
            try:
                if self.toolbar.mode != '':
                    self.toolbar.mode = ''
                if hasattr(self.toolbar.canvas, 'widgetlock'):
                    self.toolbar.canvas.widgetlock.release(self.toolbar)
            except:
                pass
            
            try:
                self.toolbar.destroy()
            except:
                pass
            
            self.toolbar = None
        
        # Create new toolbar for current tab
        if current_tab and hasattr(current_tab, 'canvas') and current_tab.canvas:
            try:
                self.toolbar = QuietNavigationToolbar(current_tab.canvas, self.toolbar_frame_container)
                self.toolbar.update()
                current_tab.toolbar = self.toolbar
            except Exception as e:
                print(f"Warning: Failed to create toolbar: {e}")
                self.toolbar = None
    
    def _create_bottom_bar(self):
        """Create bottom bar with toolbar and manual button"""
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        
        self.toolbar_frame_container = tk.Frame(self.bottom_frame)
        self.toolbar_frame_container.pack(side=tk.LEFT, padx=5)
        
        self.toolbar = None
        
        manual_button = ttk.Button(self.bottom_frame, text="View KSTAR Diagnostics Manual", 
                                   command=lambda: None)
        manual_button.pack(side=tk.RIGHT, padx=5)
        
        developer_label = tk.Label(
            self.bottom_frame, 
            text="Developed by Jekil Lee (jklee@kfe.re.kr)"
        )
        developer_label.pack(side=tk.RIGHT, padx=5)
    
    def _capture_and_save(self):
        """Capture window and save as 300 DPI PNG and EPS"""
        print("\n" + "=" * 50)
        print("Capturing GUI for paper submission...")
        print("=" * 50)
        
        # Force update to ensure all widgets are rendered
        self.root.update_idletasks()
        self.root.update()
        
        # Get window geometry
        x = self.root.winfo_rootx()
        y = self.root.winfo_rooty()
        w = self.root.winfo_width()
        h = self.root.winfo_height()
        
        print(f"Window position: ({x}, {y})")
        print(f"Window size: {w} x {h}")
        
        # Capture the window region
        bbox = (x, y, x + w, y + h)
        image = ImageGrab.grab(bbox)
        
        # Calculate physical size in inches for 300 DPI
        width_inch = w / 300.0
        height_inch = h / 300.0
        print(f"Physical size at 300 DPI: {width_inch:.2f} x {height_inch:.2f} inches")
        
        # Save as PNG with 300 DPI metadata
        png_path = f"{self.output_prefix}_300dpi.png"
        image.save(png_path, 'PNG', dpi=(300, 300))
        print(f"Saved: {png_path}")
        
        # Save as EPS
        eps_path = f"{self.output_prefix}_300dpi.eps"
        
        # Convert to RGB if necessary (EPS doesn't support RGBA)
        if image.mode == 'RGBA':
            rgb_image = Image.new('RGB', image.size, (255, 255, 255))
            rgb_image.paste(image, mask=image.split()[3])
            rgb_image.save(eps_path, 'EPS', dpi=(300, 300))
        else:
            image.save(eps_path, 'EPS', dpi=(300, 300))
        print(f"Saved: {eps_path}")
        
        print("\n" + "=" * 50)
        print("Capture complete!")
        print("=" * 50 + "\n")
        
        # Close the application
        self.root.after(500, self.root.destroy)
    
    def run(self):
        """Start the application"""
        self.root.mainloop()


def main():
    """Main entry point"""
    # Get output prefix from command line
    output_prefix = sys.argv[1] if len(sys.argv) > 1 else 'prism_gui'
    
    print("\n" + "=" * 50)
    print("PRISM GUI Capture Tool")
    print("For paper submission (300 DPI)")
    print("=" * 50)
    print(f"Output prefix: {output_prefix}")
    print("=" * 50 + "\n")
    
    app = PRISMCaptureApp(output_prefix)
    app.run()


if __name__ == "__main__":
    main()