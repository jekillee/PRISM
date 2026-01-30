# Changelog

All notable changes to PRISM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.8] - 2026-01-30

### Added
- **N-Mode Spectrum - Use Full Shot Length Option**
  - New checkbox "Use full shot length" in Parameters panel
  - When enabled, time input fields are disabled and calculation uses actual shot duration
  - Automatically detects shot's time range from MDS+ data
  - Setting persisted across sessions

## [1.1.7] - 2026-01-29

### Changed
- **MSE Profile - Improved q/j Visualization**
  - q/j profile now plotted as solid line with error band (fill_between)
  - Markers added at TGAMMA measurement positions (interpolated values)
  - Shows full pmse profile shape while highlighting actual measurement locations
  - Applied to both R plot and EFIT-mapped plot

- **MSE Text Export - Unified Format**
  - Profile and time trace text exports now save gamma, q, j in same row (like TiVT)
  - q and j are interpolated at TGAMMA measurement positions
  - Single table format: R, gamma, gamma_err, drr, q, q_err, j, j_err
  - Removed Section 1/Section 2 separation

- **Code Quality - Specific Exception Handling**
  - Replaced bare `except:` clauses with specific exceptions (`MdsException`, `TclError`, etc.)
  - Affected files: thomson_loader.py, bes_loader.py, ece_loader.py, efit_loader.py, xics_loader.py, base_tab.py
  - Improves debugging by preserving exception context

- **Code Quality - Type Hints**
  - Added Python type hints to core data structures (`core/data_structures.py`)
  - Added type hints to base loader class (`data_loaders/base_loader.py`)
  - Improves IDE support and code maintainability

### Added
- **ProfileBaseTab Base Class** (`ui/profile_base_tab.py`)
  - New abstract base class for profile tabs (TiVT, NeTe, MSE)
  - Applied to nete_profile_tab, tivt_profile_tab, mse_profile_tab
  - Extracts common profile plotting logic
  - Reduces code duplication in profile tab implementations

- **TimeTraceBaseTab Base Class** (`ui/timetrace_base_tab.py`)
  - New abstract base class for time trace tabs (TiVT, NeTe, MSE)
  - Applied to nete_timetrace_tab, tivt_timetrace_tab, mse_timetrace_tab
  - Extracts common time trace plotting logic
  - Reduces code duplication (~220 lines removed across 6 tab files)

## [1.1.6] - 2026-01-28

### Fixed
- **Thomson ne/Te - Asymmetric Error Bars**
  - Thomson error bars were incorrectly plotted as symmetric (averaged upper+lower)
  - Now properly displays asymmetric error bars using separate upper/lower bounds
  - Updated txt export to include `Te_err_u`, `Te_err_l`, `ne_err_u`, `ne_err_l` columns

### Changed
- **File Save Dialog - Hide Hidden Folders**
  - Save as txt/npz dialogs no longer show hidden folders (`.` prefixed)
  - Applied to all save dialogs (profile, time trace, IRVB, spectrogram, N-mode spectrum)
- **EFIT Source in TXT Export**
  - Profile tab txt files now include `# EFIT Source: {tree}` in header when EFIT mapping is applied
  - Applied to ne/Te, Ti/vT, and MSE profile tabs

## [1.1.5] - 2026-01-27

### Fixed
- **IRVB - Ptot/Recon Time Mismatch**
  - Fixed `ptot` not being sliced together with `time` and `recon` in `_slice_by_ip_fault` and `_slice_by_efit_time`
  - Previously, slicing only applied to `time` and `recon`, causing `ptot` indices to be misaligned with the actual time points

### Changed
- **IRVB - Example Script**
  - Added `Ptot` trace to regional Prad time trace plot in NPZ example script

## [1.1.4] - 2026-01-22

### Added
- **TV Tab - Compare Mode**
  - Side-by-side viewing of TV01 (right) and TV02 (left)
  - Time-synchronized display (matching frames by timestamp)
  - "No Data" display when one TV is out of range
  - Dropdown option: TV01, TV02, or TV01 + TV02
  - Default to TV01 + TV02 when both are available
- **TV Tab - Time Input**
  - Direct time input (seconds) in Frame Control panel
  - Navigate to specific time with Go button or Enter key
  - Time display auto-updates when navigating frames

### Changed
- **TV Tab UI Refactored**
  - Separated Search and Load buttons (previously combined as Fetch)
  - Search: finds available TVs and updates dropdown
  - Load: loads selected TV option
  - Removed "Load from File" row, replaced with "..." button

### Fixed
- **TV Tab Slider Sync**
  - Fixed frame/filename mismatch when using slider
  - Implemented debounce to prevent race condition in event queue

### Removed
- Frame skip option from Playback Control (simplified to 1-frame steps)

## [1.1.3] - 2026-01-19

### Fixed
- **N-Mode Spectrum**
  - Fixed "inhomogeneous array" error when Mirnov channels have different data points
  - Channels with different sampling rates are now interpolated to common time grid
  - Simplified plot titles: frequency plot shows only shot number, amplitude plot has no title

## [1.1.2] - 2026-01-19

### Fixed
- **MSE Time Trace**
  - Fixed `gamma_min` variable not initialized error
  - Changed q/j error display from errorbar to fill_between style (consistent with gamma)

## [1.1.1] - 2026-01-08

### Changed
- **N-Mode Spectrum UI Improvements**
  - Panel style unified (numbering, center-aligned labels)
  - Shot up/down button height matched with entry field
  - Tolerance/Fraction layout separated to two rows
- **TV Tab**
  - Removed Frame skip feature
- **TivT Profile/Timetrace**
  - Adjusted dropdown and button widths

## [1.1.0] - 2026-01-07

### Added
- **User Settings Persistence**
  - Automatic save/restore of user settings across sessions
  - Settings stored in ~/.config/prism/settings.json
- **Update Notification Popup**
  - Shows changelog on first launch after update
  - "Do not show again" option
- **NPZ Data Export with Example Scripts**
  - Spectrogram: Save data as NPZ with Python plotting example
  - N-Mode Spectrum: Save mode analysis results with example script
  - IRVB: Save 2D Prad data with equilibrium overlay example
  - Syntax-highlighted code display (Spyder dark theme)
- **TV Tab Enhancements**
  - 2026 campaign support (new file path structure)
  - Time display on frames: t = frame/210 - 0.1s
  - Shot number up/down buttons
- **Shot Navigation Buttons**
  - Added up/down buttons for Spectrogram, N-Mode Spectrum, IRVB, TV tabs
- **TCI IP Fault Masking**
  - Automatic data masking after IP fault time + 0.5s
  - Consistent with Thomson scattering masking

### Changed
- IP fault masking extended to fault_time + 0.5s (Thomson, TCI)
- Improved CES error messages for missing data

### Fixed
- Show Nodes annotation clipping (no longer affects figure size)
- neTe profile figure size consistency with other tabs
- TV tab status_label error on load

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