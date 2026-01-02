# Changelog

All notable changes to PRISM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-01-02

### Added
- Initial release of PRISM (Plasma Research Integrated System for Multi-diagnostics)
- **CES (Charge Exchange Spectroscopy)**
  - Profile and time trace views
  - mod/nn analysis type selection
  - File loading support (TGF format)
- **XICS (X-ray Imaging Crystal Spectrometer)**
  - Ti, vT profile and time trace views
  - Overlay capability on CES profiles
- **Thomson Scattering**
  - Te, ne profile and time trace views
  - IP fault time masking
- **ECE (Electron Cyclotron Emission)**
  - Te profile and time trace views
  - Configurable sampling rates (20Hz to 1000Hz)
  - Baseline correction (-1 to 0 sec)
  - Bad channel filtering
  - Overlay capability on Thomson profiles
- **ECEI (Electron Cyclotron Emission Imaging)**
  - 2D Te fluctuation imaging
  - Support for GT, GR, HT devices
  - ABCD matrix-based position calculation
  - Spectrogram analysis support
- **MSE (Motional Stark Effect)**
  - TGAMMA raw data visualization (25 channels)
  - q (safety factor) and j (current density) profiles
  - NB source auto-detection (NB-1A, NB-1B, NB-1C)
  - Bad channel detection and filtering
  - EFIT R_edge integration
- **Spectrogram Analysis**
  - ECE spectrogram support
  - ECEI spectrogram support
  - BES spectrogram support
  - TCI (Two-Color Interferometer) spectrogram support
  - Mirnov coil spectrogram (toroidal/poloidal arrays)
  - Dynamic range control (1-11 decades)
  - NFFT selection (256-4096)
- **n-Mode Spectrum Analysis**
  - Toroidal mode number analysis from Mirnov coil array
  - Parallel data loading from MDS+ with progress display
  - FFT-based spectral analysis with configurable time interval
  - Mode number calculation (n=1 to n=8)
  - Amplitude evolution tracking per mode
  - Sign selection (positive/negative/absolute/all modes)
  - Integrate option (dB/dt to B conversion)
  - Detrend option (8-segment linear detrend)
  - Contour and imshow plot types
  - Settings persistence across sessions
- **TV (IVIS Visible Camera) Viewer**
  - ZIP file loading for sequential images (BMP, PNG, JPG, etc.)
  - Frame navigation: slider, buttons, direct input
  - Mouse wheel support for frame navigation
  - Playback control: Play/Pause/Stop, FPS, frame skip, loop
  - Image caching for smooth playback
  - Auto-search TV01/TV02 files by shot number
  - Line drawing for paper figures:
    - Draw mode toggle
    - Multi-point curve drawing with spline smoothing
    - Customizable line color and width
    - Show/hide and clear options
- **IRVB (Infra-Red Video Bolometer) Viewer**
  - 2D radiation profile visualization (Prad in MW/m^3)
  - EFIT equilibrium overlay (LCFS, psi contours, X-point, limiter)
  - Regional Prad analysis with configurable psi boundaries
  - Time traces for each region
  - Frame-by-frame navigation with mouse wheel support
  - Data from HTTP server with local caching
- **Standalone Viewers**
  - Ti, vT standalone (tivt)
  - ne, Te standalone (nete)
  - MSE standalone (mse)
  - Spectrogram & n-Mode standalone (spec)
  - TV standalone (tv)
  - IRVB standalone (irvb)
- **EFIT Mapping**
  - Support for multiple EFIT trees (efitrt1, efitrt2, efit01, efit02, efit04)
  - Three coordinate systems: psi_N, rho_pol, rho_tor
- Modular architecture for easy diagnostic extension
- Unified Thomson/ECE overlay capability
- Symbolic link installation (prism command)
- Lazy tab loading for improved startup performance

### Technical Details
- Python 3.8+ compatible
- MDS+ database connectivity
- Tkinter-based GUI with matplotlib integration
- Pillow for image processing (TV tab)

---

## Version Numbering

- **MAJOR**: Incompatible API changes or major feature additions
- **MINOR**: New features in a backward-compatible manner
- **PATCH**: Backward-compatible bug fixes

## Maintainer

Jekil Lee (jklee@kfe.re.kr)
