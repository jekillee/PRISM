<p align="center">
  <img src="ui/icons/prism-logo-dark.svg" alt="PRISM Logo" width="128">
</p>

# PRISM

**Plasma Research Integrated System for Multi-diagnostics**

A unified visualization and data analysis platform for KSTAR tokamak.

## Features

- **Modular Architecture**: Easily extensible for new diagnostics
- **Sidebar Navigation**: Categorized tree view with on-demand tabs
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
- **Profile Fitting**: mtanh, ptanh, EPED, spline, RBF, GPR with pedestal analysis
- **X-axis Selection**: R/ψN/ρpol/ρtor radio buttons in profile plot section
- **Interactive Channel Toggle**: Double-click data points to exclude/include channels from fitting
- **Dark/Light Theme**: Runtime theme switching with persistence
- **Data Preview & Save**: Spreadsheet-style preview before exporting data, with p-file (PEQDSK) format export for fitted profiles
- **Profile Browse**: Slider-based data browsing with playback before selecting time points
- **Time Averaging**: dt-based profile averaging for fitting and visualization
- **Collapsible Sidebar**: Expand/collapse category groups with state persistence
- **TRANSP CDF Viewer**: Load TRANSP output CDF files for profile and time trace visualization with variable filter/search and cross-run comparison
- **EFIT Viewer**: Time traces (AEQDSK scalars), profiles (GEQDSK), 2D equilibrium contours, and p-file viewer with MDS+ and file support
- **Grouped Sidebar**: Collapsible groups (Diagnostics/EFIT/BiProfile/TRANSP) with expand/collapse state persistence

## Supported Diagnostics

| Diagnostic | Description | Parameters |
|------------|-------------|------------|
| CES | Charge Exchange Spectroscopy | Ti, vT |
| XICS | X-ray Imaging Crystal Spectrometer | Ti, vT |
| Thomson | Thomson Scattering | Te, ne, ne_avg |
| ECE | Electron Cyclotron Emission | Te |
| ECEI | ECE Imaging | 2D Te fluctuation |
| MSE | Motional Stark Effect | gamma, q, j |
| BES | Beam Emission Spectroscopy | Spectrogram |
| TCI | Two-Color Interferometer | ne, Spectrogram |
| Mirnov | Magnetic Probes | Spectrogram, n-mode |
| Neutron | Fusion Neutron (near J-port) | Fission, He3, Diamond |
| TV | Visible Camera (IVIS) | Image sequence |
| IRVB | Infra-Red Video Bolometer | 2D Prad |
| EFIT | Equilibrium Fitting | Scalars, profiles, 2D ψ, p-file |
| TRANSP | Transport Analysis (BiProfile + CDF) | Ti, vT, Te, ne, CDF profiles/traces |

## Usage

| Command | Description |
|---------|-------------|
| `prism` | Full PRISM (Diagnostics + EFIT + BiProfile + TRANSP) |
| `prism -d`, `prism --diag` | Diagnostic data viewer only |
| `prism -e`, `prism --efit` | EFIT viewer |
| `prism -b`, `prism --biprofile` | BiProfile viewer |
| `prism -t`, `prism --transp` | TRANSP CDF viewer |
| `prism -s`, `prism --select` | Select and launch individual viewers |
| `prism -h`, `prism --help` | Show help |

## Installation

### For nkstar / ukstar Users (KFE Internal)

PRISM is pre-installed on both servers. No installation required.

```bash
# Full PRISM (Diagnostics + EFIT + BiProfile + TRANSP)
prism

# Diagnostic data viewer only
prism -d

# EFIT viewer
prism -e

# BiProfile viewer
prism -b

# TRANSP CDF viewer
prism -t

# Select and launch individual viewers
prism -s
```

| Server | PRISM Path | Python | Note |
|--------|-----------|--------|------|
| nkstar | `/home/users/jklee/PRISM` | 3.8 (`/usr/bin/python38`) | All tabs |
| ukstar | `/UKSTAR_HOME/jklee/PRISM` | 3.8 (`/usr/bin/python3.8`) | TV/IRVB excluded |

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
│   ├── file_parser.py           # File parser
│   └── fitting.py               # Profile fitting functions
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
│   ├── efit_loader.py           # EFIT loader
│   ├── biprofile_loader.py     # BiProfile loader (BIPROFILE + DIAG_PARAMS + raw)
│   ├── transp_cdf_loader.py   # TRANSP CDF (netCDF) loader
│   └── efit_viewer_loader.py  # EFIT viewer loader (MDS+, g-file, a-file)
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
│   ├── biprofile_profile_tab.py # BiProfile profile tab
│   ├── biprofile_timetrace_tab.py # BiProfile time trace tab
│   ├── transp_profile_tab.py   # TRANSP CDF profile tab
│   ├── transp_timetrace_tab.py # TRANSP CDF time trace tab
│   ├── efit_scalar_tab.py      # EFIT time traces (AEQDSK scalars)
│   ├── efit_profile_tab.py     # EFIT profiles (GEQDSK)
│   ├── efit_2d_tab.py          # EFIT 2D equilibrium contours
│   ├── efit_pfile_tab.py       # EFIT p-file viewer
│   ├── icons/                   # SVG icons (logo, themed widgets)
│   └── widgets/
│       ├── custom_toolbar.py    # Custom matplotlib toolbar
│       ├── toggle_switch.py    # iOS-style animated toggle
│       ├── preview_dialog.py   # Profile Browse dialog
│       └── biprofile_browse_dialog.py  # BiProfile Browse dialog
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
- scikit-learn - GPR fitting (optional)
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