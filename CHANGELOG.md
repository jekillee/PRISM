# Changelog

All notable changes to PRISM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.2] - 2026-02-17

### Fixed
- **N-Mode Spectrum NPZ Export**: Fixed bug where exporting negative-mode-only data included positive mode amplitudes and unfiltered mode spectrum
  - `mode_spectrum` now correctly filtered by sign setting (pos/neg/abs)
  - `amplitude` array now contains only the selected sign's modes instead of all modes
  - Negative-mode export now preserves negative sign in `mode_spectrum` (consistent with `n_modes_list`)
  - Affected export modes: positive-only, negative-only, and absolute
- **N-Mode Spectrum Example Script**: Added missing `get_mode_color()` function definition that caused `NameError` when running the exported example script

### Changed
- **N-Mode Spectrum Channel Display**: Channel list now printed with numbered rows and toroidal angles instead of flat list

## [2.0.1] - 2026-02-15

### Added
- **PRISM Logo**: Added prism refraction logo (dark/light variants) across the application
  - Sidebar header (inline with PRISM title)
  - Select Viewer (`prism -s`) header
  - What's New dialog header
  - Window icon for all windows (main, selector, standalone)
  - README.md header

### Changed
- **Faster Startup**: Deferred first tab creation and What's New popup to after window is shown using `QTimer.singleShot()`; What's New shown first (lightweight), then first tab loaded after dialog closes
- **Select Viewer Compactness**: Reduced window size from 300x550 to 240x380 with compact button height (26px), tighter spacing, and inline logo+title header
- **README Updated**: Replaced outdated standalone viewer commands (`prism tivt`, `prism nete`, etc.) with current `prism -s` usage; added Sidebar Navigation to features list

## [2.0.0] - 2026-02-13

### Changed - GUI Framework Migration (tkinter → PySide6)
- **Complete GUI rewrite**: Migrated all 22 GUI files from tkinter/ttk to PySide6 (Qt 6)
- **Sidebar Navigation**: Replaced flat 11-tab `ttk.Notebook` with categorized `QTreeWidget` sidebar
  - 4 groups: Profiles, Time Traces, Spectral Analysis, Imaging
  - `QStackedWidget` for lazy-loaded tab content (preserved fast startup)
- **Widget Migration**:
  - `ttk.Frame` → `QWidget` with `QVBoxLayout`/`QHBoxLayout`/`QGridLayout`
  - `ttk.LabelFrame` → `QGroupBox`
  - `ttk.Entry` / `tk.StringVar` → `QLineEdit` (direct `.text()` access)
  - `ttk.Combobox` → `QComboBox`
  - `ttk.Checkbutton` / `tk.BooleanVar` → `QCheckBox` (`.isChecked()`)
  - `ttk.Radiobutton` / `tk.StringVar` → `QRadioButton` + `QButtonGroup`
  - `tk.Listbox` + `ttk.Scrollbar` → `QListWidget` (built-in scrollbars)
  - `ttk.Scale` → `QSlider`
  - `tk.Toplevel` → `QDialog`
  - `messagebox` → `QMessageBox`
  - `filedialog` → `QFileDialog`
- **Matplotlib Backend**: `FigureCanvasTkAgg` → `FigureCanvasQTAgg`, `NavigationToolbar2Tk` → `NavigationToolbar2QT`
- **Layout System**: `pack()`/`grid()` → `QVBoxLayout`/`QHBoxLayout`/`QGridLayout`/`QSplitter`
- **Timer Pattern**: `frame.after()` → `QTimer.singleShot()`
- **Entry Point**: `tk.Tk()` + `root.mainloop()` → `QApplication` + `app.exec()`

### Added
- **Dark/Light Theme Toggle**: Runtime theme switching via sidebar button
  - `ThemeManager` class in new `ui/theme.py` module (~650 lines)
  - Applies to all widgets, matplotlib canvases, toolbar, sidebar, and scrollbars
  - Theme preference saved to `settings.json` and restored on startup
- **Shot Number Persistence**: All Profile and TimeTrace tabs save/restore shot number to `settings.json`
  - Implemented in `BaseTab._restore_shot_from_settings()` and `BaseTab.save_settings()`
  - Works in both main window and standalone mode
- **Data Preview & Save**: Spreadsheet-style preview dialog before saving data
  - Replaced separate "Save as .txt" / "Save as .npz" buttons with unified "Preview & Save" button
  - Tabular data preview with copy and save functionality
- **Per-Channel Visibility Toggle**: Profile tabs (TiVT, NeTe, MSE) now have "Select Channels" button
  - Select/deselect individual channels to dim (gray) on the profile plot
  - Channel list built from selected listbox entries (not from cache)
  - ECE uses channel-number-based keys for correct cross-shot behavior
- **Custom Combobox Arrow Icons**: SVG arrows for dark/light themes (`ui/icons/`)

### Changed - Standalone Mode
- **Unified Launcher**: Removed per-tab CLI commands (`prism tivt`, `prism spec`, etc.) and `standalone_launcher.py`
  - `prism` launches full PRISM with sidebar
  - `prism -s` (or `--select`) opens a tab selector with buttons for each viewer
  - Button click launches that single tab using the same `PRISMApp` class (consistent theme, toolbar, settings)

### Changed - Startup Message
- Removed "Enabled diagnostics" list from console output; now shows only the banner box
- Added update date next to version in banner (e.g., `PRISM v2.0.0 (2026-02-13)`)
- `UPDATE_DATE` constant added to `config/app_config.py`

### Changed - UI Polish
- **Shot Up/Down Buttons**: Unified to compact vertical stack layout across all tabs
- **EFIT Labels**: Unicode subscript notation for ψ_N, ρ_pol, ρ_tor
- **Sidebar**: Adjusted title font size (PRISM 25px, version 12px), improved spacing
- **Control Panels**: Unified panel names, box widths, and section numbering across all tabs
- **TV Tab**: Improved toolbar and dropdown styling
- **TV Startup Tab**: Improved Add button and panel layout
- **Playback Controls**: Redesigned layout for TV and IRVB tabs
- **Combobox Styling**: Fixed text clipping near dropdown arrow, unified height (`min-height: 20px`) across all dropdowns
- **RadioButton Styling**: Square indicator (`border-radius: 3px`) with blue fill on checked; added disabled state styling
- **MSE q/j RadioButtons**: Stretch-based layout for even spacing with Plot button
- **Widget Widths**: Adjusted shot entry, analysis combo, and diagnostic combo widths for better fit; removed unnecessary `setMinimumWidth` constraints from shot entry fields
- **Scrollable Control Panels**: Profile and time trace control areas now scroll when window is small
- **Matplotlib Integration**: Dark-themed axes, grid, and tick labels; save-to-white for exported figures
- **What's New Dialog**: Tab labels changed from "v2.x"/"v1.x" to "v2"/"v1"

### Fixed
- **PYTHONPATH in Launcher**: `run_prism.sh` now exports PySide6 `site-packages` path so all users can access PySide6 without local installation
- **X11 Forwarding Compatibility**: Environment variables set before QApplication init (`main.py` and `run_prism.sh`)
  - `QT_XCB_GL_INTEGRATION=none`, `QT_XCB_NO_XI2=1`, `QT_X11_NO_MITSHM=1`
  - `LIBGL_ALWAYS_SOFTWARE=1`, `QT_QUICK_BACKEND=software`
  - `unset WAYLAND_DISPLAY` (NoMachine warning), `NO_AT_BRIDGE=1`
- **ne/Te Profile ECE Listbox**: Fixed ECE-only mode not populating listbox (`'ECE only'` → `'ECE'` match)
- **ECE Channel Disable Bug**: Changed channel keys from index-based (`ECE_0`) to channel-number-based (`ECE_ch62`) to prevent cross-shot mismatch when disabling channels with different BT
- **Select Channels Scope**: Channel list now based on selected listbox entries, not entire data cache (applied to CES, Thomson, MSE, ECE)
- **Dependency Versions**: Pinned `PySide6>=6.2.4`, `matplotlib>=3.5.0` in requirements.txt

### Removed
- `#!/usr/bin/python3.8` shebang from all 34 Python files
- Unused tkinter-era widget size constants (`ENTRY_WIDTH_*`, `BUTTON_WIDTH_*`, `PAD_*`, etc.)
- `AxisControlToolbar` class from `custom_toolbar.py` (unused since v1.1.9 toolbar redesign)
- Duplicate `get_mode_color()` function in `nmode_spectrum_tab.py`
- `DARK_STYLE` backward-compatibility shim from `ui_constants.py`
- `examples/sidebar_demo.py` prototype file
- `standalone_launcher.py` — standalone mode now handled by `PRISMApp(mode='select')` in `main_window.py`

### Technical
- Version bumped from 1.1.11 to 2.0.0 (major version for framework change)
- New `ui/theme.py` module with `ThemeManager`, dark/light QSS, QPalette, matplotlib rcParams
- SVG icon files: `arrow-down-dark.svg`, `arrow-down-light.svg`, `check-white.svg`
- `SidebarNav` class in `ui/main_window.py` for categorized navigation
- `BaseTab._settings_key` mapping for per-tab settings persistence

## [1.1.11] - 2026-02-05

### Added
- **TV Startup Comparison Tab**: New tab for comparing plasma startup sequences across multiple shots
  - Support for TV01 or TV02 camera selection (locked after first shot added)
  - Fetch button to check TV availability before adding shot
  - Configurable frame range (start/end, max 50 frames)
  - Stack multiple shots vertically for visual comparison
  - Shot number with up/down buttons for easy navigation
  - Add button next to TV radio buttons (disabled until Fetch pressed)
  - Time labels displayed at top (every 4 frames)
  - Shot number labels on left side
  - Plot button to generate comparison image
  - Keyboard delete support (Delete/Backspace) for removing shots from list
  - Settings persistence (shot, frame range, TV type)
  - Standalone mode support: `prism startup`

### Changed
- **Standalone Commands**: Split `prism spec` into two separate commands
  - `prism spec` - Spectrogram only
  - `prism nmode` - n-Mode Spectrum only

### Refactored
- **TV Utilities Module** (`ui/tv_utils.py`): Extracted common TV functions
  - `TV_FPS`, `TV_OFFSET` constants
  - `get_year_from_shot()`, `get_campaign_from_shot()`, `get_campaign_from_year()`
  - `get_tv_zip_path()`, `get_tv_startup_zip_path()`
  - `find_available_tvs()` with startup option
  - `frame_to_time()`, `time_to_frame()`, `frame_to_time_ms()`
  - Applied to both TV tab and TV Startup tab (reduced code duplication)

## [1.1.10] - 2026-02-04

### Fixed
- **n-Mode Spectrum**: Fixed "too many values to unpack" error when loading Mirnov data
  - `_comment` key in mirnov_config.json was causing unpacking failure in `get_shot_year()`

### Changed
- **Console Output Consistency**: Unified print message format across all tabs and loaders
  - All messages now use `[DiagnosticName]` prefix (e.g., `[CES]`, `[Thomson]`, `[ECE]`, `[IRVB]`)
  - Sub-messages indented with `[DiagnosticName]   ` format
  - Diagnostic-specific prefixes: `[CES]`, `[XICS]`, `[Thomson]`, `[ECE]`, `[TCI]`, `[MSE]`, `[IRVB]`, `[TV]`, `[n-Mode]`, `[Spectrogram]`, `[EFIT]`
  - Removed redundant "Creating tab:" messages
  - Capitalized first letter of messages for consistency (e.g., "Data loaded", "Not available")

- **Status Label**: Added bold font and color-coded status labels to Spectrogram, TV, IRVB tabs
  - Gray: Ready, Blue: Running, Green: Done, Red: Error

## [1.1.9] - 2026-02-02

### Changed
- **Axis Control - Toolbar Integration**
  - Removed separate Axis Control Panel from all tabs
  - Added X, Y1, Y2 axis control buttons to navigation toolbar
  - Click buttons to open dialog for setting axis min/max limits
  - "Auto" button for automatic axis scaling
  - Toolbar moved to main window bottom bar (same row as "Developed by" label and Manual button)
  - TV tab uses QuietNavigationToolbar (no Axes button)
  - IRVB tab: All psi time traces selectable in Axes dialog (dynamic configuration based on psi boundaries), 2D plot excluded

- **UI Consistency Improvements**
  - Fetch button width unified (width=8) across all tabs
  - TV tab: "Search" button renamed to "Fetch"
  - Ti/vT, ne/Te, MSE tabs (profile & timetrace): Added shot up/down buttons (▲/▼) next to shot input
  - ne/Te profile dropdown: Changed options to 'TS+ECE', 'TS', 'ECE (100Hz)', 'ECE (1kHz)' (removed "only" suffix), reduced width to 10
  - ne/Te timetrace dropdown: Reduced width to 10
  - Save buttons: Changed "Save as NPZ" to "Save as .npz" for consistency with ".txt" format
  - IRVB tab: Changed limiter line color from black to white for better visibility

- **Save Data Section**
  - Ti/vT, ne/Te, MSE timetrace/profile tabs: Added numbered LabelFrame for Save Data section (consistent with standalone tabs)

### Technical
- New `AxisControlToolbar` class in `ui/widgets/custom_toolbar.py`
- New `QuietNavigationToolbar` class for tabs without Axes control (TV)
- Toolbar management centralized in `main_window.py` with `_update_toolbar()` method
- Dynamic axes configuration for IRVB with `_configure_irvb_axes()` method
- Shot adjustment implemented via `_adjust_shot(delta)` method in profile/timetrace tabs

### Refactoring
- **Configuration Externalization**
  - Thomson channel positions moved to `config/thomson_positions.json` (8 shot ranges)
  - Mirnov coil configurations moved to `config/mirnov_config.json` (yearly signs, angles, channel names)
  - Shot-to-year mapping now configurable via JSON

- **Code Quality Improvements**
  - Replaced all bare `except:` with specific exceptions (`TclError`, `Exception`)
  - Moved `get_ip_fault_time()` to `BaseDiagnosticLoader` base class (removed from thomson, ece, mse, tci loaders)
  - Added `IP_THRESHOLD_KA` constant for IP fault detection threshold
  - Fixed division by zero risk in `EFITData2D.get_normalized_psi()` with denominator check
  - Fixed IndexError risk in `CESFileParser` with array bounds validation
  - Added None checks for time arrays in Thomson and TCI loaders
  - Centralized hardcoded values: `CONTACT_EMAIL`, `AUTHOR_NAME`, `IRVB_SERVER` in `app_config.py`
  - Fixed mixed indentation (tabs → spaces) in `main_window.py`

- **IP Fault Masking Unification**
  - Added `IP_FAULT_OFFSET` constant (0.5s) to `base_loader.py`
  - Added `get_valid_time_mask()` method to `BaseDiagnosticLoader` for centralized masking
  - All loaders (Thomson, ECE, TCI) now use unified masking logic
  - Fixed ECE loader which was missing the +0.5s offset (now consistent with Thomson/TCI)

## [1.1.8] - 2026-01-29

### Added
- **N-Mode Spectrum - Use Full Shot Length Option**
  - New checkbox "Use full shot length" in Parameters panel
  - When enabled, time input fields are disabled and calculation uses actual shot duration
  - Automatically detects shot's time range from MDS+ data
  - Setting persisted across sessions

### Changed
- **UI Consistency Improvements**
  - Spectrogram, TV, IRVB tabs: Added consistent label column width (LABEL_COLUMN_WIDTH)
  - All tabs: Labels left-aligned with consistent spacing to input fields
  - Spectrogram: Removed colons from labels, changed "Shot" to "Loaded Shot" in Select Signal
  - TV, IRVB: Simplified Frame Control - removed navigation buttons, kept < > next to slider
  - TV, IRVB: Frame Control uses grid layout for consistent Frame/Time alignment
  - TV: Added shot up/down buttons (▲/▼)
  - IRVB: Added Time [s] input field for direct time navigation

- **N-Mode Spectrum UI Reorganization**
  - Separated Parameters panel (inputs only) and Plot panel (Calculate and Plot button, options)
  - Plot type changed from dropdown to radio buttons
  - Contour levels disabled when imshow selected
  - Time/Freq entry widths unified (width=10)

- **Save Data Panel**
  - Spectrogram, N-Mode Spectrum, IRVB: "Save as NPZ" and "Example Script" buttons now in single row
  - Changed "Show Example Script" to "Example Script"

- **TV Draw Line Improvements**
  - Draw Mode button now shows green background when ON for better visibility
  - Right-click finishes current line; next left-click starts a new separate line
  - Support for drawing multiple lines on both TV01 and TV02 in compare mode
  - Frame/Time inputs use fixed-width Entry for total frames display
  - Added Linestyle dropdown (dashed/solid/dotted), default is dashed
  - Default line color changed to white
  - Fixed TV02 image flickering in compare mode when drawing with blitting optimization

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
- Tkinter-based GUI with matplotlib integration (migrated to PySide6 in v2.0.0)
- Pillow for image processing (TV tab)

---

## Version Numbering

- **MAJOR**: Incompatible API changes or major feature additions
- **MINOR**: New features in a backward-compatible manner
- **PATCH**: Backward-compatible bug fixes

## Maintainer

Jekil Lee (jklee@kfe.re.kr)