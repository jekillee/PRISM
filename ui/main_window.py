#!/usr/bin/python3.8

"""
Main application window with on-demand tab initialization
"""

import tkinter as tk
from tkinter import ttk, messagebox
import subprocess
import os

from config.app_config import AppConfig, VERSION
from config.diagnostic_config import get_enabled_diagnostics
from data_loaders.efit_loader import EFITLoader
from plotting.plot_manager import PlotManager
from ui.tab_factory import TabFactory
from ui.widgets.custom_toolbar import QuietNavigationToolbar


class PRISMApp:
    """Main application class for PRISM with on-demand tab initialization"""
    
    def __init__(self):
        self.config = AppConfig()
        self.efit_loader = EFITLoader(self.config)
        self.plot_manager = PlotManager(self.config)
        
        self.root = tk.Tk()
        self.root.title(f'PRISM v{VERSION} - Plasma Research Integrated System for Multi-diagnostics')
        self.root.geometry("1500x700")
        
        # Tab management
        self.tab_configs = []  # Store tab configuration info
        self.tab_cache = {}    # Cache for created tabs {index: tab_instance}
        self.tab_frames = []   # Placeholder frames for each tab
        
        self._create_widgets()
        self._create_bottom_bar()
        
        self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_changed)
        
        # Print startup message BEFORE creating first tab
        self._print_startup_message()
        
        # Create first tab immediately for initial display
        if self.tab_configs:
            self._create_tab_content(0)
            self._on_tab_changed(None)
    
    def _create_widgets(self):
        """Create main application widgets with placeholder frames"""
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(expand=True, fill='both')
        
        enabled_diagnostics = get_enabled_diagnostics()
        tab_types = ['profile', 'timetrace']
        
        # Create diagnostic-specific tab configs
        for diag in enabled_diagnostics:
            for tab_type in tab_types:
                if not TabFactory.should_create_tab(diag, tab_type):
                    continue
                
                tab_name = TabFactory.get_tab_name(diag, tab_type)
                
                # Store config instead of creating tab
                self.tab_configs.append({
                    'diagnostic': diag,
                    'tab_type': tab_type,
                    'tab_name': tab_name
                })
                
                # Create placeholder frame
                placeholder = ttk.Frame(self.notebook)
                self.tab_frames.append(placeholder)
                self.notebook.add(placeholder, text=tab_name)
        
        # Spectrogram tab config
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'spectrogram',
            'tab_name': TabFactory.get_tab_name(None, 'spectrogram')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='Spectrogram')

	# n-Mode Spectrum tab config
        self.tab_configs.append({
	    'diagnostic': None,
	    'tab_type': 'nmode',
	    'tab_name': TabFactory.get_tab_name(None, 'nmode')
	})
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='n-Mode Spectrum')
        
        # TV tab config
        self.tab_configs.append({
            'diagnostic': None,
            'tab_type': 'tv',
            'tab_name': TabFactory.get_tab_name(None, 'tv')
        })
        placeholder = ttk.Frame(self.notebook)
        self.tab_frames.append(placeholder)
        self.notebook.add(placeholder, text='TV')
        
        # IRVB tab config
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
        
        # Create the actual tab with placeholder as parent
        tab = TabFactory.create_tab(
            placeholder,
            self.config,
            config['diagnostic'],
            config['tab_type'],
            self.efit_loader,
            self.plot_manager
        )
        
        if tab is not None:
            # Pack tab's frame inside placeholder
            tab.frame.pack(expand=True, fill='both')
            self.tab_cache[tab_index] = tab
            print(f"  Tab created: {config['tab_name']}")
        
        return tab
    
    def _on_tab_changed(self, event):
        """Handle tab change - create tab if not exists, update toolbar"""
        current_tab_index = self.notebook.index(self.notebook.select())
        
        # Create tab content if not already created
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
                                   command=self._show_manual)
        manual_button.pack(side=tk.RIGHT, padx=5)
        
        developer_label = tk.Label(
            self.bottom_frame, 
            text="Developed by Jekil Lee (jklee@kfe.re.kr)"
        )
        developer_label.pack(side=tk.RIGHT, padx=5)
    
    def _show_manual(self):
        """Open user manual PDF"""
        try:
            if not os.path.exists(self.config.MANUAL_PATH):
                messagebox.showinfo("Manual Not Found", 
                                  "User manual not found.\n\nPlease contact Jekil Lee (jklee@kfe.re.kr) for the manual.")
                return
            
            subprocess.run(["xdg-open", self.config.MANUAL_PATH], check=True)
            
        except Exception as e:
            messagebox.showerror("Error", f"Error opening manual: {str(e)}")
    
    def _print_startup_message(self):
        """Print startup message to console"""
        print("\n")
        print("+" + "=" * 62 + "+")
        print("|" + " " * 62 + "|")
        print("|" + f"PRISM v{VERSION}".center(62) + "|")
        print("|" + "Plasma Research Integrated System for Multi-diagnostics".center(62) + "|")
        print("|" + " " * 62 + "|")
        print("+" + "=" * 62 + "+")
        print("|" + " " * 62 + "|")
        print("|" + "Developed by Jekil Lee".center(62) + "|")
        print("|" + "jklee@kfe.re.kr".center(62) + "|")
        print("|" + " " * 62 + "|")
        print("+" + "=" * 62 + "+")
        print()
        print("Features:")
        print("  - On-demand tab initialization for faster startup")
        print("  - Modular diagnostic system architecture")
        print("  - Profile and time trace views")
        print("  - EFIT equilibrium mapping (efitrt1, efitrt2, efit01, efit02, efit04)")
        print("  - CES analysis type selection (mod/nn)")
        print("  - Spectrogram analysis (ECE, Mirnov, BES, TCI)")
        print("  - n-Mode spectrum analysis (toroidal mode numbers)")
        print("  - TV image sequence viewer with line drawing")
        print("  - IRVB 2D radiation profile viewer")
        print()
        print("Enabled diagnostics:")
        for diag in get_enabled_diagnostics():
            print(f"  - {diag}")
        print()
        print("Standalone tools:")
        print("  - Spectrogram & n-Mode (spec)")
        print("  - TV Viewer (tv)")
        print("  - IRVB Viewer (irvb)")
        print()
        print("=" * 64)
        print()
    
    def run(self):
        """Start the application main loop"""
        self.root.mainloop()