#!/usr/bin/python3.8

"""
Standalone launcher for individual PRISM tabs
"""

import tkinter as tk
from tkinter import ttk

from config.app_config import AppConfig, VERSION
from config.diagnostic_config import DIAGNOSTICS
from data_loaders.efit_loader import EFITLoader
from plotting.plot_manager import PlotManager
from ui.widgets.custom_toolbar import QuietNavigationToolbar


class StandaloneLauncher:
    """Launcher for standalone tab execution"""
    
    MODES = {
        'tivt': {
            'title': 'Ti, vT Viewer',
            'tabs': ['tivt_profile', 'tivt_timetrace']
        },
        'nete': {
            'title': 'ne, Te Viewer',
            'tabs': ['nete_profile', 'nete_timetrace']
        },
        'mse': {
            'title': 'MSE Viewer',
            'tabs': ['mse_profile', 'mse_timetrace']
        },
        'spec': {
            'title': 'Spectrogram Analyzer',
            'tabs': ['spectrogram', 'nmode']
        },
        'tv': {
            'title': 'TV Image Viewer',
            'tabs': ['tv']
        },
        'irvb': {
            'title': 'IRVB Viewer',
            'tabs': ['irvb']
        }
    }
    
    def __init__(self):
        self.config = AppConfig()
        self.efit_loader = EFITLoader(self.config)
        self.plot_manager = PlotManager(self.config)
    
    def run(self, mode):
        """Run standalone viewer for specified mode"""
        if mode not in self.MODES:
            print(f"Error: Unknown mode '{mode}'")
            print(f"Available modes: {', '.join(self.MODES.keys())}")
            return
        
        mode_config = self.MODES[mode]
        
        # Create root window
        self.root = tk.Tk()
        self.root.title(f"PRISM v{VERSION} - {mode_config['title']}")
        self.root.geometry("1400x700")
        
        # Print startup message
        self._print_startup_message(mode, mode_config['title'])
        
        # Create notebook if multiple tabs
        if len(mode_config['tabs']) > 1:
            self.notebook = ttk.Notebook(self.root)
            self.notebook.pack(expand=True, fill='both')
            
            self.tabs = []
            self.tab_cache = {}
            
            for i, tab_type in enumerate(mode_config['tabs']):
                placeholder = ttk.Frame(self.notebook)
                tab_name = self._get_tab_display_name(tab_type)
                self.notebook.add(placeholder, text=tab_name)
                self.tabs.append({'type': tab_type, 'frame': placeholder})
            
            # Create bottom bar
            self._create_bottom_bar()
            
            # Bind tab change event
            self.notebook.bind("<<NotebookTabChanged>>", self._on_tab_changed)
            
            # Create first tab
            self._create_tab_content(0)
            self._on_tab_changed(None)
        else:
            # Single tab mode
            tab_type = mode_config['tabs'][0]
            tab = self._create_single_tab(tab_type)
            
            if tab:
                tab.frame.pack(expand=True, fill='both')
                
                # Create bottom bar
                self._create_bottom_bar()
                
                # Setup toolbar
                if hasattr(tab, 'canvas'):
                    self.toolbar = QuietNavigationToolbar(tab.canvas, self.toolbar_frame_container)
                    self.toolbar.update()
                    tab.toolbar = self.toolbar
        
        self.root.mainloop()
    
    def _get_tab_display_name(self, tab_type):
        """Get display name for tab type"""
        names = {
            'tivt_profile': 'Ti, vT Profile',
            'tivt_timetrace': 'Ti, vT Time Trace',
            'nete_profile': 'ne, Te Profile',
            'nete_timetrace': 'ne, Te Time Trace',
            'mse_profile': 'MSE Profile',
            'mse_timetrace': 'MSE Time Trace',
            'spectrogram': 'Spectrogram',
            'nmode': 'n-Mode Spectrum',
            'tv': 'TV',
            'irvb': 'IRVB'
        }
        return names.get(tab_type, tab_type)
    
    def _create_single_tab(self, tab_type):
        """Create a single tab instance"""
        print(f"Creating tab: {self._get_tab_display_name(tab_type)}...")
        
        if tab_type == 'tivt_profile':
            from ui.tivt_profile_tab import TiVTProfileTab
            from data_loaders.ces_loader import CESLoader
            from core.file_parser import CESFileParser
            
            diag_config = DIAGNOSTICS['CES']
            data_loader = CESLoader(self.config, diag_config)
            file_parser = CESFileParser()
            
            tab = TiVTProfileTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='CES',
                tab_type='profile',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager,
                file_parser=file_parser
            )
            tab.create_widgets()
            
        elif tab_type == 'tivt_timetrace':
            from ui.tivt_timetrace_tab import TiVTTimeTraceTab
            from data_loaders.ces_loader import CESLoader
            from core.file_parser import CESFileParser
            
            diag_config = DIAGNOSTICS['CES']
            data_loader = CESLoader(self.config, diag_config)
            file_parser = CESFileParser()
            
            tab = TiVTTimeTraceTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='CES',
                tab_type='timetrace',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager,
                file_parser=file_parser
            )
            tab.create_widgets()
            
        elif tab_type == 'nete_profile':
            from ui.nete_profile_tab import NeTeProfileTab
            from data_loaders.thomson_loader import ThomsonLoader
            
            diag_config = DIAGNOSTICS['Thomson']
            data_loader = ThomsonLoader(self.config, diag_config)
            
            tab = NeTeProfileTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='Thomson',
                tab_type='profile',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager
            )
            tab.create_widgets()
            
        elif tab_type == 'nete_timetrace':
            from ui.nete_timetrace_tab import NeTeTimeTraceTab
            from data_loaders.thomson_loader import ThomsonLoader
            
            diag_config = DIAGNOSTICS['Thomson']
            data_loader = ThomsonLoader(self.config, diag_config)
            
            tab = NeTeTimeTraceTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='Thomson',
                tab_type='timetrace',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager
            )
            tab.create_widgets()
            
        elif tab_type == 'mse_profile':
            from ui.mse_profile_tab import MSEProfileTab
            from data_loaders.mse_loader import MSELoader
            
            diag_config = DIAGNOSTICS['MSE']
            data_loader = MSELoader(self.config, diag_config)
            
            tab = MSEProfileTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='MSE',
                tab_type='profile',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager
            )
            tab.create_widgets()
            
        elif tab_type == 'mse_timetrace':
            from ui.mse_timetrace_tab import MSETimeTraceTab
            from data_loaders.mse_loader import MSELoader
            
            diag_config = DIAGNOSTICS['MSE']
            data_loader = MSELoader(self.config, diag_config)
            
            tab = MSETimeTraceTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_name='MSE',
                tab_type='timetrace',
                data_loader=data_loader,
                efit_loader=self.efit_loader,
                plot_manager=self.plot_manager
            )
            tab.create_widgets()
            
        elif tab_type == 'spectrogram':
            from ui.spectrogram_tab import SpectrogramTab
            
            tab = SpectrogramTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            
        elif tab_type == 'nmode':
            from ui.nmode_spectrum_tab import NModeSpectrumTab
            
            tab = NModeSpectrumTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            
        elif tab_type == 'tv':
            from ui.tv_tab import TVTab
            
            tab = TVTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            
        elif tab_type == 'irvb':
            from ui.irvb_tab import IRVBTab
            
            tab = IRVBTab(
                parent=self.root,
                app_config=self.config,
                diagnostic_config=DIAGNOSTICS,
                efit_loader=self.efit_loader
            )
            tab.create_widgets()
            
        else:
            print(f"  Unknown tab type: {tab_type}")
            return None
        
        print(f"  Tab created: {self._get_tab_display_name(tab_type)}")
        return tab
    
    def _create_tab_content(self, tab_index):
        """Create tab content for notebook mode"""
        if tab_index in self.tab_cache:
            return self.tab_cache[tab_index]
        
        tab_info = self.tabs[tab_index]
        tab_type = tab_info['type']
        placeholder = tab_info['frame']
        
        # Temporarily set parent to placeholder for tab creation
        original_root = self.root
        self.root = placeholder
        
        tab = self._create_single_tab(tab_type)
        
        # Restore root
        self.root = original_root
        
        if tab:
            tab.frame.pack(expand=True, fill='both')
            self.tab_cache[tab_index] = tab
        
        return tab
    
    def _on_tab_changed(self, event):
        """Handle tab change event"""
        current_tab_index = self.notebook.index(self.notebook.select())
        
        # Create tab if not exists
        if current_tab_index not in self.tab_cache:
            self._create_tab_content(current_tab_index)
        
        current_tab = self.tab_cache.get(current_tab_index)
        
        # Update toolbar
        if hasattr(self, 'toolbar') and self.toolbar:
            try:
                self.toolbar.destroy()
            except:
                pass
        
        if current_tab and hasattr(current_tab, 'canvas'):
            self.toolbar = QuietNavigationToolbar(current_tab.canvas, self.toolbar_frame_container)
            self.toolbar.update()
            current_tab.toolbar = self.toolbar
    
    def _create_bottom_bar(self):
        """Create bottom bar with toolbar"""
        self.bottom_frame = tk.Frame(self.root)
        self.bottom_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=10, pady=5)
        
        self.toolbar_frame_container = tk.Frame(self.bottom_frame)
        self.toolbar_frame_container.pack(side=tk.LEFT, padx=5)
        
        self.toolbar = None
        
        developer_label = tk.Label(
            self.bottom_frame,
            text="Developed by Jekil Lee (jklee@kfe.re.kr)"
        )
        developer_label.pack(side=tk.RIGHT, padx=5)
    
    def _print_startup_message(self, mode, title):
        """Print startup message"""
        print("\n")
        print("+" + "=" * 50 + "+")
        print("|" + " " * 50 + "|")
        print("|" + f"PRISM v{VERSION} - {title}".center(50) + "|")
        print("|" + " " * 50 + "|")
        print("+" + "=" * 50 + "+")
        print("|" + " " * 50 + "|")
        print("|" + "Developed by Jekil Lee".center(50) + "|")
        print("|" + "jklee@kfe.re.kr".center(50) + "|")
        print("|" + " " * 50 + "|")
        print("+" + "=" * 50 + "+")
        print()