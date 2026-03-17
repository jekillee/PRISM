<p align="center">
  <img src="ui/icons/prism-logo-dark.svg" alt="PRISM Logo" width="128">
</p>

# PRISM

**Plasma Research Integrated System for Multi-diagnostics**

A modular diagnostic data visualization platform for KSTAR tokamak at Korea Institute of Fusion Energy (KFE).

## Features

- **Modular Architecture**: Easily extensible for new diagnostics
- **Sidebar Navigation**: Categorized tree view with lazy-loaded tabs
- **Profile & Time Trace Views**: For each diagnostic system
- **EFIT Mapping**: Support for multiple EFIT trees (efitrt1, efitrt2, efit01, efit02, efit04)
- **Multiple Data Sources**: MDS+ and file-based data
- **Unified Thomson/ECE Viewer**: Overlay ECE data on Thomson profiles
- **IP Fault Time Masking**: Automatic filtering of post-discharge data
- **Spectrogram Analysis**: FFT-based spectrogram for ECE, Mirnov, BES, TCI, ECEI
- **n-Mode Spectrum Analysis**: Toroidal mode number analysis from Mirnov coils
- **TV Image Viewer**: Sequential image viewer for visible camera data with line drawing
- **TV Startup Comparison**: Compare plasma startup sequences across multiple shots
- **IRVB Viewer**: 2D radiation profile with EFIT overlay and regional Prad analysis
- **Dark/Light Theme**: Runtime theme switching with persistence
- **Data Preview & Save**: Spreadsheet-style preview before exporting data
- **Per-Channel Visibility**: Show/hide individual channels in profile plots

## Supported Diagnostics

| Diagnostic | Description | Parameters |
|------------|-------------|------------|
| CES | Charge Exchange Spectroscopy | Ti, vT |
| XICS | X-ray Imaging Crystal Spectrometer | Ti, vT |
| Thomson | Thomson Scattering | Te, ne |
| ECE | Electron Cyclotron Emission | Te |
| ECEI | ECE Imaging | 2D Te fluctuation |
| MSE | Motional Stark Effect | gamma, q, j |
| BES | Beam Emission Spectroscopy | Spectrogram |
| TCI | Two-Color Interferometer | ne, Spectrogram |
| Mirnov | Magnetic Probes | Spectrogram, n-mode |
| Neutron | Fusion Neutron (near J-port) | Fission, He3, Diamond |
| TV | Visible Camera (IVIS) | Image sequence |
| IRVB | Infra-Red Video Bolometer | 2D Prad |

## Usage

| Command | Description |
|---------|-------------|
| `prism` | Full PRISM with sidebar navigation |
| `prism -s` | Select and launch individual diagnostic viewers |
| `prism -h` | Show help |

## Installation

### For nkstar Users (KFE Internal)

PRISM is pre-installed and available system-wide. No installation required.

```bash
# Full PRISM
prism

# Select a viewer
prism -s
```

### For External Users

PRISM requires access to KSTAR MDS+ server, which is only available from the nkstar server at KFE. External users cannot run PRISM directly.

This repository is provided for:
- Code reference and architecture review
- Adaptation to other tokamak facilities with their own MDS+ infrastructure

If you are interested in adapting PRISM for other fusion devices, please contact the author.

## Directory Structure

```
PRISM/
├── main.py                      # Main entry point (full & standalone)
├── run_prism.sh                 # Shell launcher script
├── config/
│   ├── app_config.py            # Global configuration
│   ├── diagnostic_config.py     # Diagnostic metadata
│   ├── user_settings.py         # User settings persistence
│   ├── thomson_positions.json   # Thomson channel R positions by shot range
│   └── mirnov_config.json       # Mirnov coil configurations by year
├── core/
│   ├── data_structures.py       # Data classes
│   └── file_parser.py           # File parser
├── data_loaders/
│   ├── base_loader.py           # Base loader class
│   ├── ces_loader.py            # CES loader
│   ├── xics_loader.py           # XICS loader
│   ├── thomson_loader.py        # Thomson loader
│   ├── ece_loader.py            # ECE loader
│   ├── ecei_loader.py           # ECEI loader
│   ├── mse_loader.py            # MSE loader
│   ├── bes_loader.py            # BES loader
│   ├── tci_loader.py            # TCI loader
│   ├── irvb_loader.py           # IRVB loader
│   ├── neutron_loader.py        # Neutron loader
│   └── efit_loader.py           # EFIT loader
├── ui/
│   ├── theme.py                 # Theme manager (dark/light QSS, palette, mpl)
│   ├── ui_constants.py          # UI constants and helpers
│   ├── main_window.py           # Main window with sidebar navigation
│   ├── base_tab.py              # Base tab class
│   ├── profile_base_tab.py      # Profile tab base class
│   ├── timetrace_base_tab.py    # Time trace tab base class
│   ├── tab_factory.py           # Tab creation factory
│   ├── tivt_profile_tab.py      # Ti,vT profile tab
│   ├── tivt_timetrace_tab.py    # Ti,vT time trace tab
│   ├── nete_profile_tab.py      # ne,Te profile tab
│   ├── nete_timetrace_tab.py    # ne,Te time trace tab
│   ├── mse_profile_tab.py       # MSE profile tab
│   ├── mse_timetrace_tab.py     # MSE time trace tab
│   ├── spectrogram_tab.py       # Spectrogram tab
│   ├── nmode_spectrum_tab.py    # n-Mode Spectrum tab
│   ├── tv_tab.py                # TV image viewer tab
│   ├── tv_startup_tab.py        # TV Startup Comparison tab
│   ├── tv_utils.py              # TV utility functions
│   ├── irvb_tab.py              # IRVB viewer tab
│   ├── neutron_timetrace_tab.py # Neutron time trace tab
│   ├── icons/                   # SVG icons (logo, themed widgets)
│   └── widgets/
│       └── custom_toolbar.py    # Custom matplotlib toolbar
├── plotting/
│   └── plot_manager.py          # Plot management
├── requirements.txt             # Python dependencies
├── README.md                    # This file
└── CHANGELOG.md                 # Version history
```

## Requirements

- Python 3.8+
- PySide6 (Qt 6) - GUI framework
- matplotlib - Plotting
- numpy, scipy - Numerical computation
- Pillow - Image processing (TV tab)
- MDSplus - KSTAR data access (KFE internal)

## Note for External Users

This software is developed and configured for KSTAR tokamak at KFE. The MDS+ server addresses, data paths, and tree structures in `config/` are KSTAR-specific and require access to the nkstar server.

To adapt for other facilities, modify:
- `config/app_config.py`: MDS+ server address, file paths
- `config/diagnostic_config.py`: MDS+ node paths and tree names

## Adding New Diagnostics

1. **Add configuration** in `diagnostic_config.py`
2. **Implement data loader** (inherit from `BaseDiagnosticLoader`)
3. **Create tab classes** (inherit from `BaseTab`)
4. **Register in TabFactory** (`tab_factory.py`)

See existing implementations for reference.

## Author

Jekil Lee (jklee@kfe.re.kr)

Korea Institute of Fusion Energy (KFE)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- KSTAR Diagnostics Team at KFE
- Contributors to the diagnostic data analysis routines