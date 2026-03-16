"""
Factory for creating diagnostic tabs
Each diagnostic has its own tab class

All tab/loader/parser imports are deferred to first use (lazy imports)
to minimize application startup time.
"""

from config.diagnostic_config import DIAGNOSTICS


def _get_loader_map():
    """Lazily import and return the loader class map"""
    from data_loaders.ces_loader import CESLoader
    from data_loaders.thomson_loader import ThomsonLoader
    from data_loaders.ece_loader import ECELoader
    from data_loaders.mse_loader import MSELoader
    from data_loaders.xics_loader import XICSLoader
    return {
        'CES': CESLoader,
        'XICS': XICSLoader,
        'Thomson': ThomsonLoader,
        'ECE': ECELoader,
        'MSE': MSELoader,
    }


def _get_parser_map():
    """Lazily import and return the file parser class map"""
    from core.file_parser import CESFileParser
    return {
        'CES': CESFileParser,
    }


def _get_profile_tab_map():
    """Lazily import and return the profile tab class map"""
    from ui.nete_profile_tab import NeTeProfileTab
    from ui.tivt_profile_tab import TiVTProfileTab
    from ui.mse_profile_tab import MSEProfileTab
    return {
        'Thomson': NeTeProfileTab,
        'CES': TiVTProfileTab,
        'MSE': MSEProfileTab,
    }


def _get_timetrace_tab_map():
    """Lazily import and return the timetrace tab class map"""
    from ui.nete_timetrace_tab import NeTeTimeTraceTab
    from ui.tivt_timetrace_tab import TiVTTimeTraceTab
    from ui.mse_timetrace_tab import MSETimeTraceTab
    return {
        'Thomson': NeTeTimeTraceTab,
        'CES': TiVTTimeTraceTab,
        'MSE': MSETimeTraceTab,
    }


class TabFactory:
    """Factory for creating diagnostic tabs"""

    # Static maps populated on first access
    LOADER_MAP = None
    PARSER_MAP = None
    PROFILE_TAB_MAP = None
    TIMETRACE_TAB_MAP = None

    # Diagnostics that have profile/timetrace tabs (no import needed for checks)
    _PROFILE_DIAGNOSTICS = {'Thomson', 'CES', 'MSE'}
    _TIMETRACE_DIAGNOSTICS = {'Thomson', 'CES', 'MSE'}
    
    # Custom tab display names
    TAB_NAMES = {
        'Thomson': {
            'profile': 'ne, Te',
            'timetrace': 'ne, Te'
        },
        'CES': {
            'profile': 'Ti, vT',
            'timetrace': 'Ti, vT'
        },
        'MSE': {
            'profile': 'MSE',
            'timetrace': 'MSE'
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
            return 'n-mode'

        if tab_type == 'tv':
            return 'TV'

        if tab_type == 'irvb':
            return 'IRVB'

        if tab_type == 'tv_startup':
            return 'TV Startup'

        if tab_type == 'neutron':
            return 'Neutron'

        if diagnostic_name in TabFactory.TAB_NAMES:
            name = TabFactory.TAB_NAMES[diagnostic_name].get(tab_type)
            if name is not None:
                return name
        
        # Default naming
        return diagnostic_name
    
    @staticmethod
    def should_create_tab(diagnostic_name, tab_type):
        """Check if tab should be created"""
        if tab_type in ('spectrogram', 'nmode', 'tv', 'irvb', 'tv_startup', 'neutron'):
            return True

        if tab_type == 'profile':
            return diagnostic_name in TabFactory._PROFILE_DIAGNOSTICS
        else:
            return diagnostic_name in TabFactory._TIMETRACE_DIAGNOSTICS
    
    @classmethod
    def _ensure_maps(cls):
        """Populate class-level maps on first use (lazy import)"""
        if cls.LOADER_MAP is None:
            cls.LOADER_MAP = _get_loader_map()
            cls.PARSER_MAP = _get_parser_map()
            cls.PROFILE_TAB_MAP = _get_profile_tab_map()
            cls.TIMETRACE_TAB_MAP = _get_timetrace_tab_map()

    @staticmethod
    def create_tab(notebook, app_config, diagnostic_name, tab_type,
                   efit_loader, plot_manager):
        """Create a diagnostic tab"""
        # Special case: Spectrogram tab
        if tab_type == 'spectrogram':
            from ui.spectrogram_tab import SpectrogramTab
            tab = SpectrogramTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab

        # Special case: n-Mode Spectrum tab
        if tab_type == 'nmode':
            from ui.nmode_spectrum_tab import NModeSpectrumTab
            tab = NModeSpectrumTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab

        # Special case: TV tab
        if tab_type == 'tv':
            from ui.tv_tab import TVTab
            tab = TVTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab

        # Special case: IRVB tab
        if tab_type == 'irvb':
            from ui.irvb_tab import IRVBTab
            tab = IRVBTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS,
                efit_loader=efit_loader
            )
            tab.create_widgets()
            return tab

        # Special case: TV Startup tab
        if tab_type == 'tv_startup':
            from ui.tv_startup_tab import TVStartupTab
            tab = TVStartupTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab

        # Special case: Neutron tab
        if tab_type == 'neutron':
            from ui.neutron_timetrace_tab import NeutronTimeTraceTab
            tab = NeutronTimeTraceTab(
                parent=notebook,
                app_config=app_config,
                diagnostic_config=DIAGNOSTICS
            )
            tab.create_widgets()
            return tab

        # For other tab types, diagnostic_name must be valid
        if diagnostic_name is None or diagnostic_name not in DIAGNOSTICS:
            raise ValueError(f"Unknown diagnostic: {diagnostic_name}")

        if not TabFactory.should_create_tab(diagnostic_name, tab_type):
            return None

        # Lazy-load all class maps on first tab creation
        TabFactory._ensure_maps()

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
        tab._restore_shot_from_settings()

        return tab