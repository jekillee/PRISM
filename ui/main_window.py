#!/usr/bin/python3.8

"""
Main application window with on-demand tab initialization
"""

import tkinter as tk
from tkinter import ttk, messagebox, TclError
import subprocess
import os

from config.app_config import AppConfig, VERSION
from config.diagnostic_config import get_enabled_diagnostics
from config.user_settings import (
    load_settings, save_settings, show_update_popup
)
from data_loaders.efit_loader import EFITLoader
from plotting.plot_manager import PlotManager
from ui.tab_factory import TabFactory
from ui.widgets.custom_toolbar import QuietNavigationToolbar, AxisControlToolbar


class PRISMApp:
    """Main application class for PRISM with on-demand tab initialization"""
    
    def __init__(self):
        self.config = AppConfig()
        self.efit_loader = EFITLoader(self.config)
        self.plot_manager = PlotManager(self.config)
        
        # Load user settings
        load_settings()
        
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
        
        # Show update popup if new version
        show_update_popup(self.root)
        
        # Bind window close event
        self.root.protocol("WM_DELETE_WINDOW", self._on_close)
    
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

        # Update toolbar for current tab
        self._update_toolbar(current_tab_index)

    def _update_toolbar(self, tab_index):
        """Update toolbar for the current tab"""
        # Destroy existing toolbar
        if hasattr(self, 'toolbar') and self.toolbar:
            try:
                self.toolbar.destroy()
            except TclError:
                pass
            self.toolbar = None

        current_tab = self.tab_cache.get(tab_index)
        if not current_tab or not hasattr(current_tab, 'canvas') or not current_tab.canvas:
            return

        config = self.tab_configs[tab_index]
        tab_type = config['tab_type']

        # TV: Use QuietNavigationToolbar (no Axes button)
        if tab_type == 'tv':
            self.toolbar = QuietNavigationToolbar(current_tab.canvas, self.toolbar_frame)
            self.toolbar.update()
            current_tab.toolbar = self.toolbar
            return

        # All other tabs: Use AxisControlToolbar
        self.toolbar = AxisControlToolbar(current_tab.canvas, self.toolbar_frame, tab_instance=current_tab)
        self.toolbar.update()
        current_tab.toolbar = self.toolbar

        # Configure axes based on tab type
        if tab_type == 'irvb':
            # IRVB: Configure dynamically based on loaded data
            self._configure_irvb_axes(current_tab)
        elif tab_type == 'spectrogram':
            self.toolbar.configure_axes(has_y2=False, ax1_label='Spectrogram')
        elif tab_type == 'nmode':
            self.toolbar.configure_axes(has_y2=True, ax1_label='n-mode', ax2_label='Amplitude')
        elif tab_type == 'profile':
            # Profile tabs
            if config['diagnostic'] == 'MSE':
                self.toolbar.configure_axes(has_y2=True, ax1_label='γ [rad]', ax2_label='q / j')
            else:
                self.toolbar.configure_axes(
                    has_y2=current_tab.param2 is not None,
                    ax1_label=current_tab.param1['label'],
                    ax2_label=current_tab.param2['label'] if current_tab.param2 else 'Axes 2'
                )
        elif tab_type == 'timetrace':
            # Timetrace tabs
            self.toolbar.configure_axes(
                has_y2=current_tab.param2 is not None,
                ax1_label=current_tab.param1['label'],
                ax2_label=current_tab.param2['label'] if current_tab.param2 else 'Axes 2'
            )

    def _configure_irvb_axes(self, tab):
        """Configure axes for IRVB tab with all time traces"""
        import numpy as np

        if not hasattr(tab, 'ax_traces') or not tab.ax_traces:
            # No data loaded yet, use default
            self.toolbar.configure_axes(has_y2=False, ax1_label='IRVB')
            return

        # Build figures_config for all time traces using psi_boundaries
        figures_config = []
        n_regions = len(tab.ax_traces)
        boundaries = [0] + getattr(tab, 'psi_boundaries', []) + [np.inf]

        for i, ax in enumerate(tab.ax_traces):
            setattr(tab, f'ax_trace_{i}', ax)
            psi_min = boundaries[i] if i < len(boundaries) else 0
            psi_max = boundaries[i + 1] if i + 1 < len(boundaries) else np.inf
            if psi_max == np.inf:
                label = f'ψ>{psi_min:.2f}'
            else:
                label = f'ψ={psi_min:.2f}-{psi_max:.2f}'
            figures_config.append((label, f'ax_trace_{i}'))

        self.toolbar.figures_config = figures_config
        self.toolbar.share_x = True
    
    def _create_bottom_bar(self):
        """Create bottom bar with toolbar and info"""
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)

        # Toolbar container (left side)
        self.toolbar_frame = tk.Frame(self.bottom_frame)
        self.toolbar_frame.pack(side=tk.LEFT, fill=tk.Y)
        self.toolbar = None

        # Right side: Manual button and developer info
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
    
    def _on_close(self):
        """Handle window close event"""
        # Save settings for all created tabs
        for tab_index, tab in self.tab_cache.items():
            if hasattr(tab, 'save_settings'):
                tab.save_settings()
        
        save_settings()
        self.root.destroy()
    
    def run(self):
        """Start the application main loop"""
        self.root.mainloop()