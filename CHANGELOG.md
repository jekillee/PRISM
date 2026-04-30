# Changelog

All notable changes to PRISM will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.5.4] - 2026-04-30

### Added — TRANSP UFILE Tab (single, dispatches by class)
- **TRANSP UFILE Tab**: New tab under `TRANSP > Input > UFILE`. Open any TRANSP run directory; all UFILEs appear in one dropdown, sorted by physical category (equilibrium → profile → heating → history). Plot dispatches by class:
  - 1D time traces (BDI, CUR, VSF, ...) → single line
  - profile×time (NER/TER/TIO/OME/VTR/QRA/ZFC/GRB/PRS/QPR/...) → all-times overlay with Time colorbar
  - multi-channel heating (NBI/NB2/ECP/ECA/ECB) → per-source lines
  - geometry (LIM/MMX/...) → EFIT 2D-style (R, Z) layout
- **U-File parser** (`core/ufile_parser.py`): dependency-free TRANSP U-File reader; TIME may be on X0 or X1 (name-first detection, data-pattern fallback)
- **Directory loader** (`data_loaders/transp_ufile_loader.py`): regex `^([A-Za-z]+)(\d+)\.([A-Za-z0-9]{2,4})$` accepts any prefix (F / OMF / OMFXP / ...). Auto-classifies via NDIM + axis name/units
- **NBI 4-quant rendering**: 2-panel Pinj + Vacc. Source count resolved from NBEAM/NGYRO/NSOURCE scalar, or auto-detected from data pattern (last quarter of columns in [0, 1] = Cfull/Chalf fractions). Pinj legend shows mean (Cfull, Chalf, Cthird) per source
- **MMX flux surfaces**: Fourier reconstruction `R(θ, x), Z(θ, x)` from 4-block scaled moments (`[Rcm | Rsm | Zcm | Zsm]` with x^m radial scaling). ~20 nested surfaces from magnetic axis to boundary, viridis-colored by `x = √(φ/φ_lim)`, with right-side colorbar
- **LIM auto-overlay**: Limiter contour drawn in gray on every geometry plot as a static reference frame; not separately selectable
- **Profile 2D / 3D radio**: 2D = lines + Time colorbar + legend (ellipsis for >20 items); 3D = `plot_surface(x, time, value)`
- **Time slider**: Activated only for 3D UFILEs (MMX). Live-updates on drag
- **UFILE catalog** (`config/transp_ufiles.json`): label/category mapping (history/equilibrium/profile/heating); unknown extensions fall back to the file's own `F.name` header

### Added — Browse Dialog enhancements
- **2D / 3D View toggle**: Added to ProfilePreviewDialog, TimeTracePreviewDialog, TranspProfilePreviewDialog, BiProfilePreviewDialog. 3D uses `mpl_toolkits.mplot3d` surface plot with cube-shifted layout so the z-axis label doesn't clip
- **ELM detection**: New `Find ELM peaks` toggle. Difference-of-Gaussians (mild=0.15ms, heavy=15× scale) + scipy `find_peaks` with prominence/distance filtering catches grassy ELMs. Adaptive D-α zoom (5% of Ip duration, clamped [0.1, 1.0] s)
- **Fix Axes** option for stable axis ranges across slider drags
- **Browse Time Trace channel list**: Now also includes ECE and TCI channels (was: Thomson only)

### Added — EFIT enhancements
- **EFIT 2D rational surfaces**: New toggle showing q = 1, 3/2, 2, 5/2, 3, 4, 5 surfaces as dashed contours over the equilibrium
- **EFIT Profiles X-axis = R**: R [m] is now selectable alongside ψ_N / ρ_pol / ρ_tor
- **EFIT multi-tree accumulation**: Profiles and 2D tabs accumulate fetched trees (efitrt1, efit01, ...) for direct cross-tree comparison instead of replacing the cache
- **g-File auto-detect**: Opening g-files no longer mixes them with previously fetched MDS+ trees
- **EFIT 2D limiter** rendered uniformly in gray
- **EFIT tree saved to settings**, persisted across sessions

### Added — Profile fitting (PRISM-style pedestal output)
- **Pedestal output (PRISM style)**: After fitting, the console now prints a unified pedestal summary — Height, Top, Foot, Width, Location — with ± uncertainties and both normalized (ψ_N) and physical (cm/m via EFIT) units. mtanh, ptanh, and EPED share the same definitions through `_resolve_ped_position` and `_compute_ped_height` (e.g., EPED implicitly treats position as `1 − a3/2`)
- **Fit axis tracking** (`fit_x_axis`): Fitting in ψ_N then plotting in ρ_pol no longer overlays the wrong curve — the overlay auto-hides on coordinate mismatch
- **Fit-dt tracking** (`fit_dt_s`): When fitting with dt averaging, markers now reflect the same averaged data used in the fit (black-edged markers + `±50ms` legend tag); single-slice when only mapping
- **dt averaging RMS error propagation**: `mean(σ)/√N` → `√(mean(σ²)/N)`

### Added — Toolbar / theme / TV
- **EPS save button**: Custom matplotlib toolbar now exposes PNG / SVG / EPS save buttons
- **TV campaign mapping** (`config/tv_campaigns.json`): Per-shot-range campaign directory lookup. JSON-config-first with directory-scan fallback
- **QSpinBox / QDoubleSpinBox theme styling**: Explicit dark/light theme rules added so spin controls render correctly under both themes (was: time-slice spinbox showed system colors)

### Changed — Reorganization
- **TRANSP sidebar**: `TRANSP > Input > UFILE` and `TRANSP > Output > Profiles / Time Traces` (was: flat `TRANSP CDF Profiles / Time Traces`)
- **Directory restructure**: `ui/*.py` flat layout (~30 files) → `ui/tabs/` hierarchy with `base_tab.py`/`profile_base_tab.py`/`timetrace_base_tab.py` at the top, and grouped subfolders `diagnostics/{profiles, timetraces, spectral, imaging}`, `biprofile`, `transp`, `efit`. All tab files moved with `git mv`; `tab_factory.py` lazy-imports updated to new paths
- **Plot style**: All UFILE plots use `'.-'` (lw=1.5, markersize=3) matching CDF Time Trace
- **2D plot layout (UFILE)**: EFIT 2D-style — `Figure((10, 6))`, `set_aspect('equal', adjustable='box')`, fixed `xlim=(1.0, 2.5)`, `ylim=(-1.5, 1.5)` for KSTAR (R, Z) plots. Colorbar placed in an explicit cax that auto-snaps to the right edge of the compressed data box (resize-aware) instead of stealing axes width
- **Status messages**: TRANSP UFILE / CDF tabs now use the bottom status bar (`show_status`) instead of inline labels — consistent with diagnostic tabs
- **CDF Browse ELM group**: Removed from TRANSP CDF preview dialogs (CDF has no Dα signal)
- **vT ymin via robust z-range**: 1–99 percentile + edge skipping so noise no longer drives ylim too negative
- **Time Trace 3D ordering**: 3D Browse for time trace tabs sorts by R value (when available) and maps slider to sorted position, instead of raw channel index
- **EFIT 2D limiter**: Now rendered all-gray
- **ECE tqdm**: Uses `sys.stdout` so progress shows even when stderr is redirected (`run_prism.sh` workflow)

### Fixed
- **Profile×Time UFILE classification**: GRB/PRS/QPR/NER/TER/TIO/OME/VTR/QRA/ZFC where X0 is spatial and X1 is TIME are now correctly classified (name-based time detection wins over the data-pattern heuristic — X0 in [0, 1] no longer false-matches as time)
- **CDF Browse 3D variable change**: Surface re-renders for the new variable in 3D mode (was: only the highlight line updated, surface stayed stale)
- **CDF Browse Variable dropdown**: Always-visible scrollbar + max 20 visible items (matches main combo)
- **U-File parser leading-space**: Some TRANSP U-File writers prepend a space to every line; column-based fallback was reading units `'rm'` instead of `'m'` for LIM. Fixed by `lstrip` before column splitting
- **AEQDSK length units**: RCEN, RMAXIS, AMINOR, GAPIN, GAPOUT, ... cm → m; AREA / AREAO cm² → m²; VOLUME cm³ → m³. Labels uniformly SI. X-point (RXPT1/2, ZXPT1/2) added
- **EFIT MDS+ scalar units**: Auto-detect cgs storage (e.g., ROUT in EFIT01) and convert to SI based on value magnitude
- **EFIT Time Traces shot persistence**: Shot saved/restored across sessions (matches Profiles, 2D)
- **Toolbar after tab switch**: Active pan/zoom mode held the canvas widgetlock, blocking the new tab's toolbar. Now deactivated before destruction
- **EPED PED HEIGHT**: Was reading `a2` directly (incorrect for some functions). Now evaluates the function at the resolved position. EPED also infers position implicitly as `1 - a3/2`
- **Mock vs Prod migration tests**: Time-averaged plot was being drawn even when only EFIT mapping was computed (no fit yet). Single-slice for mapping plots; averaged data only when fit was performed with dt
- **3D zlabel clipping**: Multi-iteration `_layout_3d_axes(ax, shift_frac=0.30, shrink_frac=0.10)` shifts the cube left so the right-side z-label has room
- **ECEI partial-node fallback**: Device parameter nodes (`\GT_MODE`, ...) not always uploaded → fallback to channel list without R/Z positions; missing nodes logged but loader doesn't crash
- **NameError actual_time**: Fixed in nete profile fit raw-data preview path
- **BiProfile 3D View**: BiProfilePreviewDialog now has its own View toggle and `_render_3d` (didn't inherit from `_PreviewBase`)
- **Browse dialog ECE/TCI access**: `_get_timetrace_preview_channels` now adds ECE and TCI channels even when only those are loaded (was: Thomson only)
- **Theme-aware spinbox styling**: Time-slice QSpinBox / QDoubleSpinBox now follow dark/light themes
- **Settings restore order** (UFILE tab): `_restore_settings()` now runs before `_create_widgets()` so radio buttons / combos reflect persisted values (was: widgets initialized to defaults, then settings overwrote the model only — UI fell out of sync, e.g. 2D radio shown but plot rendered 3D)

## [2.5.3] - 2026-04-27

### Added
- **EFIT Multi-tree Selection**: EFIT Profiles and 2D tabs now accumulate fetched trees instead of replacing them. Multiple EFIT trees (efitrt1, efit01, etc.) can be loaded together for cross-comparison
- **"Or" Label**: EFIT Profiles, Time Traces, 2D tabs and Ti/vT Profiles, Time Traces tabs show "Or" prefix before Open File buttons for clarity

### Changed
- **EFIT Profiles/2D Layout**: EFIT tree dropdown moved between up/down buttons and Fetch button (matches Time Traces layout). Dropdown shows plain tree names with consistent width across all EFIT tabs
- **EFIT Profiles Browse Button**: "Load EFIT data to browse" → "Load data to browse" (consistent with other EFIT tabs)
- **File Rename**: `ui/efit_scalar_tab.py` → `ui/efit_timetrace_tab.py`

### Fixed
- **Toolbar After Tab Switch**: Pan/zoom now works after switching tabs and returning. Previously, an active pan/zoom mode left the canvas widgetlock held by the destroyed toolbar, preventing the new toolbar from operating
- **g-File Available List**: Opening g-files no longer mixes them with previously fetched MDS+ trees in the Available list; only the just-opened files are shown while cache is preserved for already-selected entries
- **EFIT a-File Units**: Length quantities (RCEN, RMAXIS, AMINOR, GAPIN, GAPOUT, etc.) converted from cm to m, area (AREA, AREAO) cm² to m², volume (VOLUME) cm³ to m³. Labels updated to be uniformly SI. X-point labels (RXPT1/2, ZXPT1/2) added
- **EFIT MDS+ Scalar Units**: Auto-detect cgs storage in MDS+ AEQDSK scalars (e.g., ROUT in EFIT01) and convert to SI based on value magnitude
- **EFIT Time Traces Shot Persistence**: Shot number now saved/restored across sessions (matches Profiles and 2D tabs)

## [2.5.2] - 2026-04-20

### Added
- **p-File Export**: Preview & Save dialog now includes "p-File Format" mode after profile fitting. Exports fitted profiles in PEQDSK format (401-point psi_N grid) with Copy to Clipboard and Save as .txt
- **Startup Banner**: Rainbow-colored ASCII art PRISM logo with version info and usage guide printed at launch
- **Rainbow PRISM Title**: Sidebar, Select Viewer, and What's New dialog display PRISM in rainbow colors

### Changed
- **Startup Banner**: Removed old box-style banner from main_window, replaced with ASCII art in main.py

### Fixed
- **ne/Te Raw Data Preview**: Fixed `NameError: actual_time` when previewing Thomson raw data after fitting

## [2.5.1] - 2026-04-17

### Added
- **TS_NE_AVG**: Line-averaged electron density along TCI02 chord (`\TS_NE_AVG`) in Thomson loader and ne/Te Time Trace tab

### Fixed
- **ECEI Spectrogram**: Graceful fallback when ECEI device parameter nodes are missing (TREE-E-NODATA). Channel list loads without R, Z positions
- **ECEI Position Diagnostics**: Log specific missing MDS+ nodes for ECEI position calculation

## [2.5.0] - 2026-04-15

### Added
- **EFIT Viewer** (`prism -e`): New standalone EFIT viewer with four tabs
  - Time Traces: AEQDSK scalar variables from MDS+ or a-files (multi-select). EFITVIEWER_KFE mapping standard
  - Profiles: GEQDSK 1D profiles from MDS+ or g-files. X-axis: ψ_N/ρ_pol/ρ_tor
  - 2D: Contour plots (ψ_N/ρ_pol/ρ_tor) with separatrix, boundary, magnetic axis, X-points, limiter. Configurable SOL contours
  - p-File: PEQDSK profile viewer with Browse
- **Diagnostic-only Mode** (`prism -d`): Launches with Diagnostics sidebar only (no EFIT/BiProfile/TRANSP)
- **BiProfile Standalone** (`prism -b`): Dedicated BiProfile viewer mode
- **Grouped Sidebar**: Full PRISM sidebar organized into collapsible groups (Diagnostics → EFIT → BiProfile → TRANSP) with text arrows inheriting font style/color. Divider lines between groups
- **a-file Parser**: AEQDSK parser using EFITVIEWER_KFE standard (Block 1: 24 scalars, Block 2: 41 scalars). BETAN computed from BETAT/(Ip[MA]/(a[m]×BT[T]))
- **g-file Parser**: GEQDSK parser with Fortran format support and limiter data (XLIM/YLIM)
- **p-file Parser**: PEQDSK parser (block headers with psinorm variables)
- **Time Trace Browse**: Slider-based channel browsing in all Time Trace tabs (Ti/vT, ne/Te, MSE)
- **TRANSP Browse Dialogs**: Profile Browse (variable selector + time slider) and Time Trace Browse (variable selector, all runs plotted)
- **TRANSP X-axis Selection**: R/ψ_N/ρ_pol/ρ_tor radio buttons in Plot section. R = LFS major radius [m]
- **TRANSP/BiProfile Preview & Save**: Spreadsheet-style data preview with Copy to Clipboard and Save as .csv
- **EFIT/BiProfile/TRANSP in prism -s**: Tab selector with grouped QGroupBox layout

### Changed
- **Sidebar Restructuring**: Four collapsible groups with Unicode arrows (▶/▼). Sub-categories indented, leaf tabs further indented. QFrame dividers between groups
- **Selector UI** (`prism -s`): QGroupBox per group (Diagnostics/EFIT/BiProfile/TRANSP) with blue header, gray sub-labels. Version/date in header, author in footer
- **Browse Rename**: All "Preview" buttons renamed to "Browse". "Preview" reserved for Save Data section
- **Preview Dialog Consolidation**: All preview dialogs merged into single `preview_dialog.py`
- **TRANSP Profile Workflow**: Open CDF keeps cache for cross-CDF comparison; Selected list preserved across CDF opens
- **Variable Settings**: TRANSP/EFIT variable selection saved/restored across sessions
- **Non-native File Dialogs**: All non-modal QFileDialog use DontUseNativeDialog to suppress GTK transient parent warnings
- **Startup Banner**: Updated with all launch modes (prism, -d, -e, -b, -t, -s, -h)
- **run_prism.sh**: All modes stderr-redirected (`2>/dev/null`) to suppress MDSplus buffer_free messages

### Fixed
- **BiProfile fetch in full PRISM**: Fixed `AttributeError: '_bi_shot_data'` (missing initialization in full mode)
- **EFIT Scalar Tab**: Fixed `QGridLayout` not imported, variable dropdown empty (indentation bug), Browse KeyError
- **EFIT a-file Mapping**: Corrected RMAXIS/ZMAXIS positions, QSTAR vs Q0 distinction, header skip
- **EFIT g-file/p-file Cache**: Added time to cache key to avoid collision for same shot different times
- **EFIT 2D Tab**: Full rewrite (removed references to non-existent self.app_config/self.efit_loader)
- **prism -s Crash**: _SingleTabWindow now handles BiProfile/TRANSP/EFIT tab types
- **Profile Browse Line Style**: p-file Browse now uses line+marker (was marker-only)

## [2.4.2] - 2026-04-14

### Added
- **Time Trace Browse**: Slider-based channel browsing in all Time Trace tabs (Ti/vT, ne/Te, MSE). Browse button in Select Data section
- **TRANSP Browse Dialogs**: Profile Browse (variable selector + time slider) and Time Trace Browse (variable selector, all runs plotted)
- **TRANSP X-axis Selection**: R/ψ_N/ρ_pol/ρ_tor radio buttons in Plot section (same layout as PRISM Profiles). R = LFS major radius [m]
- **TRANSP Preview & Save**: Spreadsheet-style data preview with Copy to Clipboard and Save as .csv (PRISM-style format with all coordinate columns)
- **BiProfile Preview & Save**: Same Preview & Save for both BiProfile Profile and Time Trace tabs
- **BiProfile/TRANSP in prism -s**: Tab selector now includes BiProfile and TRANSP tabs with grouped layout
- **BiProfile/TRANSP in full PRISM sidebar**: Hierarchical categories (BiProfile > Profiles/Time Traces, TRANSP flat)

### Changed
- **Browse Rename**: All "Preview" buttons renamed to "Browse" across all tabs (PRISM, BiProfile, TRANSP). "Preview" reserved for Save Data section
- **Preview Dialog Consolidation**: All preview dialogs merged into single `preview_dialog.py` (ProfilePreviewDialog, TimeTracePreviewDialog, TranspProfilePreviewDialog, TranspTimeTracePreviewDialog, BiProfilePreviewDialog)
- **TRANSP Profile Workflow**: Open CDF keeps cache for cross-CDF comparison; Selected list preserved across CDF opens; re-selecting same CDF switches without reload
- **TRANSP CDF Filter**: `.CDF` only by default in file dialog
- **TRANSP Status Label**: Moved to Load CDF section; colored (green) instead of gray
- **Plot Option Auto-apply**: Font size changes in Option dialog now immediately re-plot (all custom tabs)
- **Selector UI**: Category headers in blue (#0d6efd); TRANSP shown as single button row without sub-labels
- **Variable Settings**: TRANSP variable selection saved/restored across sessions

### Fixed
- **TRANSP Filter**: Filtering to non-matching text no longer clears Available listbox
- **TRANSP CDF Reopen**: Previously opened CDF can be re-selected (switches via cache instead of being silently ignored)
- **prism -b/--bi removed**: Only `-t`/`--transp`/`transp` supported

## [2.4.1] - 2026-04-14

### Added
- **TRANSP CDF Viewer**: New TRANSP category in `prism transp` sidebar with Profile and Time Trace tabs for TRANSP output CDF files
  - Profile tab: Open CDF → filter/select variable → select time points → plot radial profiles. Cross-CDF comparison via Selected list
  - Time Trace tab: Open multiple CDFs → filter/select variable → compare runs. Remove Selected Run button for run management
  - Auto-detect default CDF directory (nkstar: `~/`, ukstar: `/UKSTAR_HOME/data/transp/{user}/`)
  - Duplicate CDF detection (same file cannot be loaded twice)
  - netCDF4 backend with scipy fallback; auto-categorization of profile (2D) vs time trace (1D) variables
  - Non-modal file dialog for CDF open
- **ECE tqdm Progress Bar**: Terminal progress bar during parallel ECE channel loading (fallback to log if tqdm unavailable)
- **BiProfile Settings Persistence**: Shot number, color mode, font sizes, show nodes, TS scale state saved/restored across sessions

### Changed
- **`prism transp` Command**: Renamed from `prism bi` (old aliases removed). Launch via `-t`/`--transp`/`transp`. Window title: "TRANSP Viewer"
- **Startup Banner & What's New**: Now shown in all modes (full, select, transp)
- **BiProfile Default Colormap**: Changed from Fixed(tab10) to Gradient(viridis) for all BiProfile tabs
- **BiProfile TS ne Scale**: Scale applied by default; toggle renamed to "Unapply TS ne Scale"
- **TRANSP Sidebar**: BiProfile and TRANSP as independent collapsible root nodes with ▼/▶ arrows and state persistence

### Fixed
- **ne/Te Profile ECE+TS**: Fixed ECE data not showing in TS+ECE mode (undefined `dt_s` variable in `plot_data()` caused silent NameError)
- **BiProfile Preview**: Fixed import error (`BiProfilePreviewDialog` → `BiProfileBrowseDialog`)
- **BiProfile Color Bug**: Fixed all entries showing same color with Gradient colormaps (viridis etc.) by using mode prefix instead of `hasattr(cmap, 'colors')`

## [2.4.0] - 2026-04-13

### Added
- **Profile Browse Dialog**: Slider-based data browsing before Select Data in all Profile tabs (Ti/vT, ne/Te, MSE). Browse button opens non-modal dialog with playback, frame/time navigation, Fix Axes toggle, and Add to Selected
- **Time Averaging (dt)**: Toggle + dt input in Fitting section. Averages profiles over [t-dt, t+dt] for fitting and plotting. Marker distinction: filled (raw) vs filled+black edge (averaged). Validation for invalid/negative/too-small dt values
- **Preview & Save Averaged Data**: New "Averaged Data" tab in Preview & Save when dt is active, alongside Raw Data showing all individual time slices in the dt window
- **ECE Parallel Loading**: Channel-level parallel MDS+ loading (1 worker per channel). ~3-5x speedup across all BT values
- **BiProfile Select Data**: BiProfile tabs now use Available/Selected listbox workflow with Browse dialog, matching PRISM Profiles pattern. Colormap selection in Plot Options
- **Sidebar Collapse/Expand**: Click category headers (▼/▶) to collapse/expand. State saved to user settings
- **Dropdown Settings Persistence**: CES analysis type, Thomson diagnostic selection saved/restored across sessions

### Changed
- **ukstar Support**: Auto-detect nkstar/ukstar via hostname. Single `run_prism.sh` and `main.py` for both servers (Python 3.8). TV/TV Startup/IRVB tabs hidden on ukstar (nkstar-local resources)
- **Sidebar**: "Spectral Analysis" renamed to "Spectral". Neutron accent color removed
- **Non-modal Dialogs**: All popup dialogs (Plot Options, Preview & Save, Channel Selector, Fitting Options, Example Script, What's New) converted to non-modal for independent window interaction
- **R-shift in Export**: Preview & Save Raw Data now includes R-shift correction
- **Disabled Widget Styling**: QLineEdit:disabled and QComboBox:disabled now visually distinct in both dark and light themes

### Fixed
- **Toolbar RuntimeError**: Fixed crash when hovering canvas after tab switching (disconnects matplotlib callbacks before toolbar destruction)

## [2.3.4] - 2026-04-10

### Added
- **BiProfile Viewer**: `prism bi` mode for BIPROFILE Bayesian inference fitting results (Ti/vT, ne/Te profiles & time traces)

### Fixed
- **Toolbar RuntimeError**: Fixed crash when hovering canvas after tab switching (deleted C++ QLabel)

## [2.3.3] - 2026-04-07

### Changed
- **Thomson 2024 R Positions**: Updated 2024 campaign channel R positions confirmed by Dr. Jongha Lee (2026-04-07)

## [2.3.2] - 2026-03-30

### Fixed
- **Profile Plot Stale Data**: Fixed bug where adding/removing selected data after EFIT mapping or fitting would still plot the old selection. Now EFIT mapping, fit results, and flux x-axis are automatically reset when the selected data list changes.
- **TCI Validation HFS Fix**: Fixed incorrect synthetic line-averaged density for inner channels (especially TCI01). HFS psi_N sign convention (negative) caused ne=0 along half the beam path. Now uses abs(psi_N) for symmetric ne evaluation.
- **Fitting Option Dialog**: Apply button now closes the dialog after applying parameters
- **Fix Checkbox Visibility**: Fix checkboxes in Fitting Option dialog are now visible in dark theme (replaced QTableWidgetItem checkbox with QCheckBox widget)

### Added
- **TCI Validation**: ne/Te profile tab now has a TCI Validation toggle in Fitting section. After fitting ne, computes synthetic line-averaged density along each TCI chord (TCI01~05) from the fitted ne profile and compares with measured TCI values in the terminal. Setting is saved to user preferences.
- **Toggle Switch Widget**: iOS-style animated toggle switch replaces checkboxes for Show Nodes and TCI Validation for a cleaner, more consistent UI
- **ECE TF Coil Current Log**: ECE loader now prints TF coil current (kA) during data loading
- **Preview Copy to Clipboard**: Added "Copy to Clipboard" button in Preview & Save dialog for all modes (Raw Data, Fitted Profile, Fit Parameters)
- **Toolbar Coordinate Display**: Restored mouse coordinate display in navigation toolbar (was suppressed)

### Changed
- **Fitting Section Disabled Until EFIT**: "5. Fitting" section is disabled until EFIT mapping is computed, with label indicating requirement
- **Show Nodes Toggle**: Replaced QCheckBox with animated toggle switch in profile tabs (ne/Te, Ti/vT, MSE)
- **Preview Save by Mode**: Save button in Preview & Save dialog now saves the currently viewed mode (Raw Data / Fitted Profile / Fit Parameters) instead of always saving raw data. Fit Parameters saves as .txt instead of .csv. File extension is auto-appended if missing.

## [2.3.1] - 2026-03-26

### Fixed
- **Clipboard Copy**: Multi-cell selection in Preview & Save spreadsheet now copies all selected cells (Ctrl+C)
- **ne/Te Profile Preview**: Raw Data tab now shows all columns instead of only 2 (header parsing fix)
- **ne/Te Time Trace Preview**: Same header parsing fix applied to time trace preview
- **MSE gamma Y-axis**: Fixed y-axis clipping by using nanmin/nanmax instead of nanpercentile for enabled channels
- **Formula Display**: Fixed mtanh formula — separated pedestal (`y_ped`) and core terms to avoid `y` on both sides; removed incorrect `(x < 1-a3)` condition from EPED formula

### Added
- **Fitted Profile Source Column**: Added source column (CES/Thomson/ECE) to Fitted Profile preview in Ti/vT and ne/Te tabs

### Changed
- **Fitting X-axis Auto-expand**: Plot x-axis automatically expands when fitting X Range exceeds default limits (±0.1 margin)

## [2.3.0] - 2026-03-24

### Added
- **Profile Fitting**: Added profile fitting feature to Ti/vT and ne/Te profile tabs
  - Fitting functions: mtanh (8-param), ptanh (7-param), EPED (6-param), spline, RBF, GPR
  - RBF: Adjustable number of Gaussian bases (default 5); fewer = smoother, more = detailed
  - GPR: Squared-exponential kernel with ±2σ uncertainty bands
  - Option dialog with editable parameter table (Value, Min, Max, Fix) and LaTeX formula
  - Bold function title and single-line description with word wrap in all fit dialogs
  - Console output: parameter table with values, errors, and bounds after fitting
  - Pedestal summary: HEIGHT, FOOT, WIDTH (with cm conversion), LOC (with R conversion)
  - Fitted curves overlaid on profile plots (both R-space and flux-space)
  - Fitting restricted to flux coordinates (ψₙ, ρₚₒₗ, ρₜₒᵣ); R-space shows warning
  - Preview & Save extended with Fitted Profile and Fit Parameters views
  - Fit function selection and RBF n_bases saved to user settings
- **X-axis Selection**: Unified plot dispatch with R/ψₙ/ρₚₒₗ/ρₜₒᵣ radio buttons
  - "X-axis" label left-aligned, radio buttons right-aligned with equal spacing
  - Flux radios disabled until EFIT mapping is computed
  - Single Plot button dispatches to R-space or flux-space plot
  - Show Nodes checkbox toggles node labels immediately without re-plotting

### Changed
- **UI Restructuring (Profile Tabs)**: Reordered sections — EFIT Mapping (3) → Plot (4) → Fitting (5) → Save (6)
  - Removed separate EFIT Plot button; Plot button now handles both R and flux coordinates
  - All internal re-plot calls (style dialog, fitting, channel toggle) route through unified dispatcher
- **MSE q/j Selection**: Changed from radio buttons to dropdown (QComboBox) in both profile and timetrace tabs
- **Select Channels Dialog**: Single-line hint text; added note that unchecked channels are excluded from fitting
- **Show Nodes**: Checkbox now toggles node annotations immediately (unchecking removes without re-plot)
- **Fitting Descriptions**: Cleaned up mtanh, ptanh, EPED descriptions to remove external-tool references
- **Pedestal Coordinate**: Width/position conversion now uses actual fitting coordinate instead of hardcoded ψₙ
- **Formula Color**: LaTeX formula in fit dialogs uses system text color instead of hardcoded gray
- **Formula Display**: mtanh definition shown as `where,` line; ptanh formula merged to single line; a1 placed first in both
- **Channel Toggle**: Double-click message now shows channel label matching Select Channels dialog
- **Channel Toggle Axis**: Double-click preserves axis limits; Plot button auto-scales y-axis to enabled channels only
- **Flux Plot Shading**: Replaced vertical dashed lines at x=0,1 with gray shaded regions outside plasma (x≤0, x≥1)
- **Flux X-axis Range**: Default x-axis range changed from [0, 1.05] to [-0.1, 1.1]

## [2.2.3] - 2026-03-23

### Fixed
- **Multi-user Python Environment Isolation**: Prevented package version conflicts when other users run PRISM by filtering other users' `~/.local` paths from `sys.path` at startup
- **Runtime Warnings**: Suppressed numpy version mismatch warnings on startup

## [2.2.2] - 2026-03-20

### Changed
- **Sidebar Version Display**: Show update date next to version — `v2.2.2 (2026-03-20)`, reduced margin between PRISM title and version
- **Global Shot Placeholder**: Changed placeholder text from `Shot #` to `Enter Shot`

### Fixed
- **Ti/vT TimeTrace File Data Not Plotted**: Fixed CES file data (tgf) not rendering — `_parse_entry` returned `tgf01` as source but cache key used `tgf`, causing cache miss

## [2.2.1] - 2026-03-18

### Changed
- **Apply All Button Icon**: Changed from `SP_DialogApplyButton` to `SP_DialogOkButton`
- **Neutron Shared X-Axis**: Linked x-axes across 3 subplots with `sharex`
- **Neutron Backspace Shortcut**: Added Backspace key shortcut for removing shots (in addition to Delete)
- **Help Message**: Enhanced `prism -h` with Available Tabs, Keyboard Shortcuts, Data Sources, Settings path

### Fixed
- **Thomson Negative Errorbar Crash**: Fixed `'yerr' must not contain negative values` crash — negative error values now clipped to 0 (errorbar not drawn for those points)
- **Toolbar QAction RuntimeError**: Suppressed `Internal C++ object already deleted` error on zoom/pan — PySide6/matplotlib compatibility issue

## [2.2.0] - 2026-03-16

### Added
- **SVG Save Button**: Added "SVG" button to toolbar next to PNG save — saves figure as SVG with white background
- **N-Mode Contour Line Width**: Added contour line width (0.1–5.0, default 0.8) control to N-Mode Plot Options dialog
- **N-Mode Amplitude Line Width**: Added amplitude line width (0.5–5.0, default 1.5) control to N-Mode Plot Options dialog
- **N-Mode Selective Mode Checkboxes**: Replaced n-modes slider (1–8) with 5 checkboxes (n=1~5) in Plot section — always calculates nmodes=5, filters at plot time
- **N-Mode Title Font Size**: Added title font size (6–24, default 12) control to N-Mode Plot Options dialog
- **Spectrogram Title Font Size**: Added title font size (6–24, default 12) control to Spectrogram Plot Options dialog
- **Global Shot Input**: Added shot number entry with "Apply All" icon button to sidebar — sets shot number across all tabs at once
- **Neutron Time Trace Tab**: New tab for fusion neutron diagnostics (near J-port) — 3x1 subplot layout (Fission Chamber, He3 Counter, Diamond Detector) with shot overplot
- **CES nn Paper**: Added `docs/CES/Lee_2026_Fusion_Eng._Des._222_115518.pdf`

### Changed
- **Unified Initial Canvas Styling**: All tabs now apply consistent font sizes (label 12, tick 10) and grid style (`--`, lw=0.3) on initial empty canvas, matching post-data-load styling
  - Affected: Profile tabs, Time Trace tabs, Spectrogram, N-Mode Spectrum, MSE Profile
- **Profile Top Margin**: Changed from `top=0.93` to `top=0.95` (unified with Time Trace)
- **PNG Save Background**: Changed from figure facecolor (dark gray) to white for clean exports
- **PNG/SVG Save**: Removed `bbox_inches='tight'` to preserve manual margin settings
- **PNG Save Button**: Replaced default matplotlib save icon with "PNG" text button matching SVG button style
- **Toolbar Icon Colors**: All toolbar icons (except theme toggle) now colorized per theme — black in light theme, white in dark theme
- **TV Tab Control Panel**: Added `ScrollBarAlwaysOn` and consistent margins (`9, 9, 9, 9`) — panel width now matches all other tabs
- **TV Startup Control Panel**: Reduced Shot List listbox stretch ratio (3:1→1:1) and input field widths (100→60px) to match other tabs' control panel width
- **TV Compare Mode Filename**: Changed filename display from `|` separator to newline for better readability
- **Edit Axis Dropdown Labels**: All tabs now show meaningful axis names in matplotlib's "Customize" → "Select Axes" dropdown instead of empty/coordinate strings
  - Profile: parameter labels (e.g., "Ti [keV]", "vT [km/s]")
  - Time Trace: parameter labels
  - MSE: "γ [rad]", "q"
  - Spectrogram: "Spectrogram"
  - N-Mode: "Frequency [kHz]", "Amplitude [Gauss]"
  - TV: "TV Image", "TV01", "TV02"
  - TV Startup: "TV Startup"
  - IRVB: "P_rad (psi=...)", "2D Profile"
- **N-Mode Time [s] Layout**: Changed from 2-row span grid to HBoxLayout with stretch, matching Freq [kHz] layout
- **N-Mode Imshow Alpha**: Fixed alpha to 1.0 (removed adjustable alpha option)
- **Spectrogram Colorbar**: Removed colorbar from spectrogram plot for cleaner display
- **Spectrogram Top Margin**: Changed from `top=0.95` to `top=0.92`
- **N-Mode Top Margin**: Changed from `top=0.95` to `top=0.92`
- **IRVB Margins**: Changed from `top=0.95, right=0.95` to `top=0.92, right=0.93`
- **Toolbar Icon Size**: Fixed inconsistent icon sizes between dark/light themes by setting `setIconSize(QSize(24, 24))`
- **Toolbar Button Colors**: PNG/SVG button colors now match toolbar icon colors (dark: #cccccc, light: #555555)
- **MSE Gamma Conversion**: Changed TGAMMA display from raw to degrees — `arctan(tgamma) * 180/π + 90`, reference line moved from 0 to 90 deg
- **MSE Labels**: Updated axis labels and export headers from `γ [rad]` to `γ [deg]`
- **IRVB Legend Colors**: Legend text and background now theme-aware via `rcParams` (previously hardcoded white/gray)
- **IRVB Time Entry Width**: Increased from 60px to 75px for full timestamp display
- **Sidebar New Tab Highlight**: New tabs (Neutron) shown with accent color in sidebar

### Fixed
- **N-Mode Imshow Colors**: Fixed black/wrong colors when using "Default" palette with imshow mode — changed from `LinearSegmentedColormap` (which created gradient artifacts for dark colors like #000000) to `ListedColormap` for flat single-color rendering
- **N-Mode Initial Canvas**: Fixed missing units on amplitude axis label (`'Amplitude'` → `'Amplitude [Gauss]'`)
- **N-Mode Contour Levels Disable**: Added visual feedback (grayed background) when contour levels input is disabled in imshow mode
- **MSE Time Trace ylim**: Fixed gamma axis bottom limit too low — initialization changed from `(0, 0)` to `(inf, -inf)` for proper auto-ranging
- **TV/TV Startup tight_layout Warning**: Suppressed matplotlib tight_layout warnings by switching to explicit `subplots_adjust`

### Removed
- **IRVB IP Fault Masking**: Removed automatic IP fault time filtering from IRVB data loading — data now shows full time range without `t_ip_fault` cutoff

## [2.1.0] - 2026-02-24

### Added
- **Plot Options Dialog**: Unified "Plot Options" button across all visualization tabs
  - Profile and Time Trace tabs: Color mode (9 colormap options), label font size, legend font size, tick font size; separator repositioned above button
  - Spectrogram tab: Colorbar selection (8 colormaps), label font size, tick font size
  - N-Mode Spectrum tab: Color mode (fixed/gradient), label font size, legend font size, tick font size; ColorManager integration
  - IRVB tab: Two-section dialog — Time Trace (color mode, label/legend/tick font sizes) and 2D Plot (colorbar, LCFS/mag. axis/flux contour colors, label/tick font sizes)
  - All settings persist across sessions
- **TV Startup Time Range**: Added Time Range row (ms) below Frame Range with bidirectional sync — editing frame values auto-updates time values and vice versa
- **N-Mode Default Color Palette**: Added "Default" option to N-Mode color dropdown — restores the original 15-color palette used in previous versions
- **Docs Folder**: Created `docs/` folder for documentation; "View Manual" button replaced with "View Docs" that opens the docs folder

### Changed
- **TV Startup Mixed TV Support**: TV01 and TV02 shots can now be freely mixed in the same comparison (no TV type locking)
- **TV Startup Shot Labels**: Left-side labels now show shot number and TV name (e.g., "38000\nTV01") instead of just shot number
- **TV Startup Frame Numbering**: Frame 1 now corresponds to the first frame in the ZIP file (consistent with TV tab), instead of the 0ms frame
- **TV Startup Section Rename**: "2. Added Shots" renamed to "2. Shot List" for consistency with other tabs
- **Select Viewer UI**: Widened selector to 360px with horizontal button rows per category; dynamic window height based on content; theme follows full PRISM theme setting automatically
- **Select Data Listbox Height**: Increased listbox height from 120px to 180px (1.5x) for better visibility
- **Option Button Width**: Option button now proportionally sized (3:1 ratio with Plot button) to match EFIT Mapping button width
- **Scrollable Control Panels**: Fixed scroll behavior — control panels now properly scroll when window is resized smaller than content
- **Spectrogram/N-Mode Margins**: Replaced `tight_layout` with explicit `subplots_adjust` so toolbar margin adjustments apply correctly
- **N-Mode Figure Margins**: Default margins now match Time Traces layout (left=0.10, right=0.97, top=0.95, bottom=0.10, hspace=0.15)
- **IRVB Plot Options**: Added Limiter color option to 2D Plot section

### Fixed
- **View Docs Path**: Docs folder path now resolves relative to application root directory instead of hardcoded absolute path
- **N-Mode RGBA Error**: Fixed `Invalid RGBA argument` error caused by incorrect discrete colormap indexing in ColorManager
- **Spectrogram Colorbar Error**: Fixed `KeyError` on re-plot caused by `colorbar.remove()` — now clears entire figure and recreates subplot
- **What's New Dialog**: v2 tab now shows all v2.x entries (2.1.0 + 2.0.x) instead of only the current minor version
- **Control Panel Stretch**: Fixed Select Data box stretching when window height exceeds control panel content
- **Zero Reference Lines**: vT=0 (Ti/vT), tgamma=0, q=1, j=0 (MSE) dashed lines now theme-aware (dark: white, light: gray); previously hardcoded silver/gray with inconsistent alpha
- **Profile Left Margins**: Increased left margin from 0.08 to 0.10 for all profile tabs (Ti/vT, ne/Te, MSE) to prevent y-axis label clipping
- **Grid Styling Unified**: All tabs (except TV, TV Startup, Spectrogram) now use consistent grid style (`lw=0.3`, theme-aware color); fixed IRVB hardcoded `lightgray` and N-Mode Spectrum `lw=0.5` + hardcoded `gray`
- **Theme Switch Updates**: `apply_theme_to_figure()` now refreshes grid colors and zero reference lines (`gid='zero_ref'`) on theme change
- **Spectrogram Shot Label**: Fixed loaded shot number always showing white text; now follows theme text color

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