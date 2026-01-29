#!/usr/bin/python3.8

"""
Diagnostic system configuration
Each diagnostic defines its metadata for automatic tab generation
"""

DIAGNOSTICS = {
    'CES': {
        'name': 'CES',
        'full_name': 'Charge Exchange Spectroscopy',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': True,
        'analysis_types': {
            'mod': {'node_prefix': '\\CES_', 'display': 'mod'},
            'nn': {'node_prefix': '\\CESNN_', 'display': 'nn'}
        },
        'default_analysis_type': 'mod',
        'parameters': [
            {
                'name': 'Ti',
                'label': 'T$_i$ [keV]',
                'short_label': 'Ti[keV]',
                'mds_node_format': '%sTI%02d',
                'mds_error_format': '%sTI%02d:err_bar',
                'scale': 1e-3,
            },
            {
                'name': 'vT',
                'label': 'v$_T$ [km/s]',
                'short_label': 'vT[km/s]',
                'mds_node_format': '%sVT%02d',
                'mds_error_format': '%sVT%02d:err_bar',
                'scale': 1.0,
            }
        ],
        'radius': {
            'mds_node_format': '%sRT%02d',
            'scale': 1e-3,
        },
        'time': {
            'mds_node_format': 'dim_of(%sTI01)',
        },
        'channels': 32,
    },
    
    'XICS': {
        'name': 'XICS',
        'full_name': 'X-ray Imaging Crystal Spectrometer',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'Ti',
                'label': 'T$_i$ [keV]',
                'short_label': 'Ti[keV]',
                'mds_node': '\\TXCS_TI053',
                'scale': 1.0,  # Already in keV
            },
            {
                'name': 'vT',
                'label': 'v$_T$ [km/s]',
                'short_label': 'vT[km/s]',
                'mds_node': '\\TXCS_VR053',
                'scale': 1.0,  # Already in km/s
            }
        ],
        'radius': 1.8,  # Fixed R position [m]
        'channels': 1,
    },
    
    'Thomson': {
        'name': 'Thomson',
        'full_name': 'Thomson Scattering',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'Te',
                'label': 'T$_e$ [keV]',
                'short_label': 'Te[keV]',
                'scale': 1e-3,
            },
            {
                'name': 'ne',
                'label': 'n$_e$ [10$^{19}$/m$^3$]',
                'short_label': 'ne[1e19/m3]',
                'scale': 1e-19,
            }
        ],
        'channels': 31,
    },
    
    'ECE': {
        'name': 'ECE',
        'full_name': 'Electron Cyclotron Emission',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'single_plot': True,
        'r_limits': (1.3, 2.3),
        'efit_limits': (-1, 1),
        'parameters': [
            {
                'name': 'Te',
                'label': 'T$_e$ [keV]',
                'short_label': 'Te[keV]',
                'scale': 1e-3,
            }
        ],
        'channels': 76,
    },
    
    'MSE': {
        'name': 'MSE',
        'full_name': 'Motional Stark Effect',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'gamma',
                'label': r'$\gamma$ [rad]',
                'short_label': 'gamma[rad]',
                'scale': 1.0,
            },
            {
                'name': 'q',
                'label': 'q',
                'short_label': 'q',
                'scale': 1.0,
            }
        ],
        'raw_channels': 25,
        'profile_points': 20,
    },
    
    'Mirnov': {
        'name': 'Mirnov',
        'full_name': 'Mirnov Coil',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'dB',
                'label': 'dB/dt [T/s]',
                'short_label': 'dB/dt',
                'scale': 1.0,
            }
        ],
        # Toroidal array (MC1T)
        'toroidal_channels': [
            'MC1T02', 'MC1T03', 'MC1T04', 'MC1T05', 'MC1T06', 'MC1T07', 'MC1T08',
            'MC1T10', 'MC1T11', 'MC1T12', 'MC1T13', 'MC1T14', 'MC1T15', 'MC1T16',
            'MC1T19'
        ],
        # Poloidal array (MC1P + MC2P)
        'poloidal_channels': [
            'MC1P01', 'MC1P02', 'MC1P03', 'MC1P04', 'MC1P05', 'MC1P06',
            'MC1P08', 'MC1P10', 'MC1P11', 'MC1P12', 'MC1P13',
            'MC1P15', 'MC1P16', 'MC1P17', 'MC1P20', 'MC1P21', 'MC1P22',
            'MC2P10', 'MC2P11'
        ],
    },
    
    'BES': {
        'name': 'BES',
        'full_name': 'Beam Emission Spectroscopy',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'intensity',
                'label': 'Intensity [a.u.]',
                'short_label': 'Int[a.u.]',
                'scale': 1.0,
            }
        ],
        # 4 vertical (01-04) x 16 radial (01-16) = 64 channels
        'vertical_channels': 4,
        'radial_channels': 16,
        'node_format': '\\BES_%02d%02d:FOO',  # vertical, radial
        'rpos_format': '\\BES_%02d%02d:RPOS',
        'vpos_format': '\\BES_%02d%02d:VPOS',
    },
    
    'TCI': {
        'name': 'TCI',
        'full_name': 'Two-Color Interferometer',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'ne_line',
                'label': 'Line-averaged n$_e$ [a.u.]',
                'short_label': 'ne_line[a.u.]',
                'scale': 1.0,
            }
        ],
        # 5 channels (line-averaged density, no R/Z position)
        'channels': ['NE_TCI01', 'NE_TCI02', 'NE_TCI03', 'NE_TCI04', 'NE_TCI05'],
    },
    
    'TV': {
        'name': 'TV',
        'full_name': 'IVIS TV (Visible Camera)',
        'mds_tree': None,  # File-based, not MDS+
        'enabled': True,
        'file_loadable': True,
        'analysis_types': None,
        'standalone_viewer': True,  # Uses TV Viewer tab
        'supported_formats': ['.bmp', '.png', '.jpg', '.jpeg', '.gif', '.tiff'],
    },
    
    'ECEI': {
        'name': 'ECEI',
        'full_name': 'Electron Cyclotron Emission Imaging',
        'mds_tree': 'kstar',
        'enabled': True,
        'file_loadable': False,
        'analysis_types': None,
        'parameters': [
            {
                'name': 'Te_fluct',
                'label': 'T$_e$ fluctuation [a.u.]',
                'short_label': 'Te_fluct[a.u.]',
                'scale': 1.0,
            }
        ],
        # 24 vertical x 8 radial (frequency) channels per device
        'vertical_channels': 24,
        'radial_channels': 8,
        # Three devices
        'devices': ['GT', 'GR', 'HT'],
        # GR excluded channels (used by other diagnostics)
        'gr_excluded_vertical': [23, 24],
        # MDS node patterns
        'node_format': '\\ECEI_{device}{v:02d}{r:02d}:FOO',
        # Device parameter nodes
        'device_params': {
            'GT': {
                'mode': '\\GT_MODE',
                'lo_freq': '\\GT_LOFREQ',
                'lens_zoom': '\\GT_LENSZOOM',
                'lens_focus': '\\GT_LENSFOCUS',
            },
            'GR': {
                'mode': '\\GR_MODE',
                'lo_freq': '\\GR_LOFREQ',
                'lens_zoom': '\\GR_LENSZOOM',
                'lens_focus': '\\GR_LENSFOCUS',
            },
            'HT': {
                'mode': '\\HT_MODE',
                'lo_freq': '\\HT_LOFREQ',
                'lens_zoom': '\\HT_LENSZOOM',
                'lens_focus': '\\HT_LENSFOCUS',
            },
        },
        # Common parameter
        'tf_current_node': '\\ECEI_I_TF',
    },
}


def get_enabled_diagnostics():
    """Return list of enabled diagnostic names"""
    return [name for name, config in DIAGNOSTICS.items() if config['enabled']]