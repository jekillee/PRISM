#!/usr/bin/python3.8

"""
Factory for creating diagnostic tabs
Each diagnostic has its own tab class
"""

from config.diagnostic_config import DIAGNOSTICS

# Tab classes
from ui.nete_profile_tab import NeTeProfileTab
from ui.nete_timetrace_tab import NeTeTimeTraceTab
from ui.tivt_profile_tab import TiVTProfileTab
from ui.tivt_timetrace_tab import TiVTTimeTraceTab
from ui.mse_profile_tab import MSEProfileTab
from ui.mse_timetrace_tab import MSETimeTraceTab
from ui.spectrogram_tab import SpectrogramTab
from ui.nmode_spectrum_tab import NModeSpectrumTab
from ui.tv_tab import TVTab
from ui.irvb_tab import IRVBTab

# Data loaders
from data_loaders.ces_loader import CESLoader
from data_loaders.thomson_loader import ThomsonLoader
from data_loaders.ece_loader import ECELoader
from data_loaders.mse_loader import MSELoader
from data_loaders.xics_loader import XICSLoader

# File parsers
from core.file_parser import CESFileParser


class TabFactory:
    """Factory for creating diagnostic tabs"""
    
    # Map diagnostic names to their loader classes
    LOADER_MAP = {
        'CES': CESLoader,
        'XICS': XICSLoader,
        'Thomson': ThomsonLoader,
        'ECE': ECELoader,
        'MSE': MSELoader,
    }
    
    # Map diagnostic names to their file parser classes
    PARSER_MAP = {
        'CES': CESFileParser,
    }
    
    # Map diagnostics to profile tab classes
    PROFILE_TAB_MAP = {
        'Thomson': NeTeProfileTab,
        'CES': TiVTProfileTab,
        'MSE': MSEProfileTab,
    }
    
    # Map diagnostics to timetrace tab classes
    TIMETRACE_TAB_MAP = {
        'Thomson': NeTeTimeTraceTab,
        'CES': TiVTTimeTraceTab,
        'MSE': MSETimeTraceTab,
    }
    
    # Custom tab display names
    TAB_NAMES = {
        'Thomson': {
            'profile': 'ne, Te Profile',
            'timetrace': 'ne, Te Time Trace'
        },
        'CES': {
            'profile': 'Ti, vT Profile',
            'timetrace': 'Ti, vT Time Trace'
        },
        'MSE': {
            'profile': 'MSE Profile',
            'timetrace': 'MSE Time Trace'
        },
        'IRVB': {
            'viewer': 'IRVB'}
    }
    
    @staticmethod
    def get_tab_name(diagnostic_name, tab_type):
        """Get display name for tab"""
        if tab_type == 'spectrogram':
            return 'Spectrogram'
        
        if tab_type == 'nmode':
            return 'n-Mode Spectrum'
        
        if tab_type == 'tv':
            return 'TV'
        
        if tab_type == 'irvb':
            return 'IRVB'
        
        if diagnostic_name in TabFactory.TAB_NAMES:
            name = TabFactory.TAB_NAMES[diagnostic_name].get(tab_type)
            if name is not None:
                return name
        
        # Default naming
        if tab_type == 'profile':
            return f"{diagnostic_name} Profile"
        else:
            return f"{diagnostic_name} Time Trace"
    
    @staticmethod
    def should_create_tab(diagnostic_name, tab_type):
        """Check if tab should be created"""
        if tab_type == 'spectrogram':
            return True
        
        if tab_type == 'nmode':
            return True
        
        if tab_type == 'tv':
            return True
        
        if tab_type == 'irvb':
            return True
        
        if tab_type == 'profile':
            return diagnostic_name in TabFactory.PROFILE_TAB_MAP
        else:
            return diagnostic_name in TabFactory.TIMETRACE_TAB_MAP
    
    @staticmethod
    def create_tab(notebook, app_config, diagnostic_name, tab_type, 
                   efit_loader, plot_manager):
        """Create a diagnostic tab"""
        # Special case: Spectrogram tab
        if tab_type == 'spectrogram':
            tab = SpectrogramTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab
        
        # Special case: n-Mode Spectrum tab
        if tab_type == 'nmode':
            tab = NModeSpectrumTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab
        
        # Special case: TV tab
        if tab_type == 'tv':
            tab = TVTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab
        
        # Special case: IRVB tab (add after TV tab case)
        if tab_type == 'irvb':
            tab = IRVBTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS,
                efit_loader=efit_loader
            )
            tab.create_widgets()
            return tab
        
        # For other tab types, diagnostic_name must be valid
        if diagnostic_name is None or diagnostic_name not in DIAGNOSTICS:
            raise ValueError(f"Unknown diagnostic: {diagnostic_name}")
        
        if not TabFactory.should_create_tab(diagnostic_name, tab_type):
            return None
        
        diag_config = DIAGNOSTICS[diagnostic_name]
        
        # Get loader
        if diagnostic_name not in TabFactory.LOADER_MAP:
            raise ValueError(f"No loader implemented for: {diagnostic_name}")
        
        loader_class = TabFactory.LOADER_MAP[diagnostic_name]
        data_loader = loader_class(app_config, diag_config)
        
        # Get file parser (if available)
        file_parser = None
        if diag_config.get('file_loadable', False) and diagnostic_name in TabFactory.PARSER_MAP:
            parser_class = TabFactory.PARSER_MAP[diagnostic_name]
            file_parser = parser_class()
        
        # Get tab class
        if tab_type == 'profile':
            tab_class = TabFactory.PROFILE_TAB_MAP[diagnostic_name]
        else:
            tab_class = TabFactory.TIMETRACE_TAB_MAP[diagnostic_name]
        
        # Create tab instance
        if file_parser:
            tab = tab_class(
                parent=notebook,
                app_config=app_config,
                diagnostic_name=diagnostic_name,
                tab_type=tab_type,
                data_loader=data_loader,
                efit_loader=efit_loader,
                plot_manager=plot_manager,
                file_parser=file_parser
            )
        else:
            tab = tab_class(
                parent=notebook,
                app_config=app_config,
                diagnostic_name=diagnostic_name,
                tab_type=tab_type,
                data_loader=data_loader,
                efit_loader=efit_loader,
                plot_manager=plot_manager
            )
        
        tab.create_widgets()
        
        return tab