# PRISM

**Plasma Research Integrated System for Multi-diagnostics**

A modular diagnostic data visualization platform for KSTAR tokamak at Korea Institute of Fusion Energy (KFE).

## Features

- **Modular Architecture**: Easily extensible for new diagnostics
- **Profile & Time Trace Views**: For each diagnostic system
- **EFIT Mapping**: Support for multiple EFIT trees (efitrt1, efitrt2, efit01, efit02, efit04)
- **Multiple Data Sources**: MDS+ and file-based data
- **Unified Thomson/ECE Viewer**: Overlay ECE data on Thomson profiles
- **IP Fault Time Masking**: Automatic filtering of post-discharge data
- **Spectrogram Analysis**: FFT-based spectrogram for ECE, Mirnov, BES, TCI, ECEI
- **n-Mode Spectrum Analysis**: Toroidal mode number analysis from Mirnov coils
- **TV Image Viewer**: Sequential image viewer for visible camera data with line drawing
- **IRVB Viewer**: 2D radiation profile with EFIT overlay and regional Prad analysis

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
| TV | Visible Camera (IVIS) | Image sequence |
| IRVB | Infra-Red Video Bolometer | 2D Prad |

## Standalone Viewers

| Command | Description |
|---------|-------------|
| `prism` | Full PRISM (all tabs) |
| `prism tivt` | Ti, vT (CES/XICS) viewer |
| `prism nete` | ne, Te (Thomson/ECE) viewer |
| `prism mse` | MSE viewer |
| `prism spec` | Spectrogram & n-Mode viewer |
| `prism tv` | TV image viewer |
| `prism irvb` | IRVB viewer |

## Installation

### For nkstar Users (KFE Internal)

PRISM is pre-installed and available system-wide. No installation required.

```bash
# Just run
prism

# Or standalone viewers
prism tivt
prism nete
prism mse
prism spec
prism tv
prism irvb
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
├── main.py                      # Main entry point
├── standalone_launcher.py       # Standalone mode launcher
├── run_prism.sh                 # Shell launcher script
├── config/
│   ├── app_config.py            # Global configuration
│   └── diagnostic_config.py     # Diagnostic metadata
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
│   └── efit_loader.py           # EFIT loader
├── ui/
│   ├── ui_constants.py          # UI constants
│   ├── main_window.py           # Main window
│   ├── base_tab.py              # Base tab class
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
│   ├── irvb_tab.py              # IRVB viewer tab
│   └── widgets/
│       └── custom_toolbar.py    # Custom matplotlib toolbar
├── plotting/
│   └── plot_manager.py          # Plot management
├── requirements.txt             # Python dependencies
├── README.md                    # This file
├── CHANGELOG.md                 # Version history
└── INSTALL.md                   # Installation guide
```

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