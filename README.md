<p align="center">
  <img src="ui/icons/prism-logo-dark.svg" alt="PRISM Logo" width="128">
</p>

# PRISM

**Plasma Research Integrated System for Multi-diagnostics**

A unified visualization and data analysis platform for KSTAR tokamak.

## Citation

If you use PRISM in your work, please cite:

> J.K. Lee, *An integrated multi-diagnostic visualization platform for KSTAR tokamak*,
> **Fusion Engineering and Design**, Volume 228 (2026), 115786.
> ISSN 0920-3796.
> https://doi.org/10.1016/j.fusengdes.2026.115786

BibTeX:
```bibtex
@article{Lee2026PRISM,
  author  = {J. K. Lee},
  title   = {An integrated multi-diagnostic visualization platform for KSTAR tokamak},
  journal = {Fusion Engineering and Design},
  volume  = {228},
  pages   = {115786},
  year    = {2026},
  issn    = {0920-3796},
  doi     = {10.1016/j.fusengdes.2026.115786}
}
```

## Features

- **Modular GUI** — sidebar viewer (Diagnostics / EFIT / BiProfile / TRANSP) with on-demand tabs, dark/light themes, and persistent session state.
- **Profiles & time traces** — every diagnostic on a selectable R / ψ_N / ρ_pol / ρ_tor x-axis with multi-tree EFIT mapping.
- **Diagnostics** — ne/Te (Thomson + ECE), Ti/vT (CES + XICS), line-averaged density (TCI / mm-wave / FIR), and Neutron.
- **Profile fitting** — mtanh / spline / GPR and more, with pedestal parameters and uncertainties.
- **Spectral** — FFT spectrogram (ECE / Mirnov / BES / TCI / ECEI) and toroidal n-mode spectrum from Mirnov coils.
- **Imaging & EFIT** — TV camera viewer, IRVB 2D radiation with EFIT overlay + regional P_rad, and EFIT viewer (scalars / profiles / 2D contours / p-file).
- **BiProfile & TRANSP** — bundled/derived-quantity profiles and TRANSP UFILE input + CDF output.
- **Headless batch + Python SDK** — compute n-mode and IRVB without the GUI (`prism nmode` / `prism irvb`, or `from batch import …`), cached as one `.npz` per shot.
- **Export** — spreadsheet-style data preview, text / p-file output, and PNG / SVG / EPS figure save.

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
| Interferometer | mm-wave (midplane) / FIR (R=1.8 m) | ne (line-averaged) |
| Mirnov | Magnetic Probes | Spectrogram, n-mode |
| Neutron | Fusion Neutron (near J-port) | Fission, He3, Diamond |
| TV | Visible Camera (IVIS) | Image sequence |
| IRVB | Infra-Red Video Bolometer | 2D Prad |
| EFIT | Equilibrium Fitting | Scalars, profiles, 2D ψ, p-file |
| BiProfile Derived | Derived plasma quantities from BiProfile + EFIT | p_e/p_i/p_tot, dT/dR, R/L_T, R/L_n, β, ν*, ω_pe, ω_ci, c_s, ρ_i, E_r (v_θ=0 / neoclassical), neoclassical v_θ for D and C⁶⁺, ω_imp−ω_main, σ_neo, η_neo, j_BS |
| TRANSP | Transport Analysis | Ti, vT, Te, ne, UFILE profiles/traces/2D, CDF profiles/traces |

## Usage

| Command | Description |
|---------|-------------|
| `prism` | Full PRISM (Diagnostics + EFIT + BiProfile + TRANSP) |
| `prism -d`, `prism --diag` | Diagnostic data viewer only |
| `prism -e`, `prism --efit` | EFIT viewer |
| `prism -b`, `prism --biprofile` | BiProfile viewer |
| `prism -t`, `prism --transp` | TRANSP CDF viewer |
| `prism -s`, `prism --select` | Select and launch individual viewers |
| `prism nmode [args]` | Headless n-mode compute + cache, no GUI (see `prism nmode -h`) |
| `prism irvb [args]` | Headless IRVB radiation + regional Prad compute + cache, no GUI (see `prism irvb -h`) |
| `prism python <script>` | Run a script in PRISM's environment (bundled Python + `import batch`), no path setup |
| `prism examples` | List the runnable batch-API example scripts as ready-to-run commands |
| `prism -h`, `prism --help` | Show help |

### Headless batch API (Python SDK)

Compute without the GUI and cache one `.npz` per shot (same layout as the GUI's saved
NPZ). Run on nkstar/ukstar, inside the KSTAR network.

```python
# n-mode (toroidal mode-number) spectrum
from batch import NModeJobSpec, run
spec = NModeJobSpec(shot=40848, tmin=2.0, tmax=8.0, fmin=0, fmax=100)
result, status = run(spec)                 # status: "computed" | "cached"
result.frequency, result.mode_spectrum, result.amplitude_max, result.amplitude_sum

# IRVB 2D radiation + regional Prad
from batch import IRVBJobSpec, run
spec = IRVBJobSpec(shot=40085, efit_tree="efit01", psi_boundaries=(0.7, 1.0))
result, status = run(spec)
result.prad_2d      # (n_frames, nZ, nR) radiated power density [MW/m^3]
result.psi_n        # (n_frames, nZ, nR) psi_N on the IRVB grid — mask any region
result.region_prad  # (n_regions, n_frames) per-region Prad [MW]

# Many shots (either spec type); one failure never stops the batch
from batch import run_many
specs = [IRVBJobSpec(shot=s) for s in (40085, 40090, 40095)]
outcomes = run_many(specs, progress=lambda o: print(o.status, o.spec.shot))
```

Run your script with the **`prism python`** helper — it uses PRISM's bundled Python
(scipy / MDSplus / h5py) and puts `batch` on the path, so you don't need to know where
PRISM is installed or set any environment yourself:

```bash
prism python my_analysis.py
prism python -c "import batch; help(batch)"
```

**Importing `batch` from your own Python session (on nkstar/ukstar) instead of using
`prism python`?** `batch` is not on `sys.path`, so add the PRISM install directory
yourself, and run under PRISM's bundled interpreter
(`<PRISM>/vendor/cpython-3.8/bin/python3.8`), which has the dependencies (scipy,
MDSplus, h5py) — a bare system Python usually lacks MDSplus:

```python
import sys
sys.path.insert(0, "/home/users/jklee/PRISM")      # nkstar
# sys.path.insert(0, "/UKSTAR_HOME/jklee/PRISM")    # ukstar
from batch import IRVBJobSpec, run
```

**`run()` returns the result as a variable** (`result, status = run(spec)`) — the
`.npz` under `~/prism_results/<subsystem>/` is just a cache (override with
`$PRISM_ARCHIVE_ROOT` or the `archive_root=` argument; `overwrite=True` recomputes).
To compute **in-memory only, with no file written**, call the compute function
directly:

```python
from batch.compute import compute_irvb, compute_nmode
result = compute_irvb(IRVBJobSpec(shot=40085))   # returns the result, writes nothing
```

**Not sure what a result contains?** Just `print(result)` — it lists every field with
its shape and units. Runnable, self-documenting examples live in
[`examples/`](examples/) (`batch_nmode.py`, `batch_irvb.py`) — run **`prism examples`**
to list them as ready-to-run `prism python <path>` commands (no need to know the
install path).

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
| nkstar | `/home/users/jklee/PRISM` | bundled `vendor/cpython-3.8/bin/python3.8` (3.8.20) | All tabs |
| ukstar | `/UKSTAR_HOME/jklee/PRISM` | bundled `vendor/cpython-3.8/bin/python3.8` (3.8.20) | All tabs except TV viewers (nkstar-only); IRVB available |

### For External Users

PRISM requires access to KSTAR MDS+ server, which is only available from KFE internal servers (nkstar / ukstar). External users cannot run PRISM directly.

This repository is provided for:
- Code reference and architecture review
- Adaptation to other tokamak facilities with their own MDS+ infrastructure

If you are interested in adapting PRISM for other fusion devices, please contact the author.

### Bundled Python (`vendor/`, v2.6.0+)

Since v2.6.0 PRISM bundles **both the Python 3.8 interpreter and every third-party package** under `vendor/`, so the server's system Python and `~/.local/lib/...` are never touched. The layout is the standard CPython tree:

```
vendor/
└── cpython-3.8/                       # portable Linux Python interpreter
    ├── bin/python3.8
    └── lib/python3.8/
        ├── (stdlib)
        └── site-packages/             # PySide6, numpy, scipy, matplotlib, Pillow,
                                       # netCDF4, scikit-learn, MDSplus, h5py, hdf5plugin
```

`vendor/` is git-ignored and distributed via the existing rsync deploy flow rather than git.

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
│   ├── mirnov_config.json       # Mirnov coil configurations by year
│   ├── irvb_limiter.json        # IRVB divertor limiter polygons (carbon / tungsten)
│   └── radius_tces.csv          # CES channel radii by shot range (RT-node fallback)
├── core/
│   ├── data_structures.py       # Data classes
│   ├── derived_quantities.py    # Derived plasma quantities (neoclassical, Sauter, ...)
│   ├── file_parser.py           # File parser
│   ├── ufile_parser.py          # TRANSP UFILE parser
│   ├── fitting.py               # Profile fitting functions
│   └── nmode.py                 # Qt-free n-mode compute core (shared by GUI/SDK/CLI)
├── batch/                       # Headless batch SDK + CLI (n-mode, IRVB)
│   ├── jobspec.py               # NModeJobSpec / IRVBJobSpec (cache-key inputs)
│   ├── compute.py               # compute_nmode / compute_irvb (wrap core + loaders)
│   ├── results.py               # NModeResult / IRVBResult (+ .npz save/load)
│   ├── archive.py               # per-shot .npz cache
│   ├── runner.py                # run / run_many (subsystem dispatch)
│   └── cli.py                   # prism nmode / prism irvb
├── examples/                    # runnable batch-API examples
│   ├── batch_nmode.py           # n-mode: run + inspect result
│   └── batch_irvb.py            # IRVB: run + inspect + custom psi region
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
│   ├── irvb_loader.py           # IRVB loaders (2D recon + \IRVB1_PRAD time trace)
│   ├── neutron_loader.py        # Neutron loader
│   ├── efit_loader.py           # EFIT loader (profile-mapping ψ_N)
│   ├── biprofile_loader.py      # BiProfile loader (BIPROFILE + DIAG_PARAMS + raw)
│   ├── transp_loader.py         # TRANSP loaders (CDF output + U-File input)
│   └── efit_viewer_loader.py    # EFIT viewer loader (MDS+, g-file, a-file)
├── ui/
│   ├── theme.py                       # Theme manager (dark/light QSS, palette, mpl)
│   ├── ui_constants.py                # UI constants and helpers
│   ├── main_window.py                 # Main window with sidebar navigation
│   ├── tab_factory.py                 # Tab creation factory
│   ├── icons/                         # SVG icons (logo, themed widgets)
│   ├── widgets/
│   │   ├── custom_toolbar.py          # Custom matplotlib toolbar
│   │   ├── toggle_switch.py           # iOS-style animated toggle
│   │   └── preview_dialog.py          # Browse dialogs
│   └── tabs/
│       ├── base_tab.py                # Base tab class
│       ├── profile_base_tab.py        # Profile tab base class
│       ├── timetrace_base_tab.py      # Time trace tab base class
│       ├── diagnostics/
│       │   ├── profiles/
│       │   │   ├── ion_profile_tab.py        # Ion (Ti, vT) profile  (CES + XICS)
│       │   │   ├── electron_profile_tab.py   # Electron (ne, Te) profile  (Thomson + ECE)
│       │   │   └── mse_profile_tab.py        # MSE q/j profile
│       │   ├── timetraces/
│       │   │   ├── ion_timetrace_tab.py        # Ion Ti/vT time trace  (CES + XICS)
│       │   │   ├── electron_timetrace_tab.py   # Electron ne/Te time trace  (Thomson/ECE/Interferometer)
│       │   │   ├── mse_timetrace_tab.py        # MSE q/j time trace
│       │   │   ├── neutron_timetrace_tab.py    # Neutron rate time trace
│       │   │   └── irvb_timetrace_tab.py       # IRVB \IRVB1_PRAD time trace
│       │   ├── spectral/
│       │   │   ├── spectrogram_tab.py          # FFT spectrogram (ECE/Mirnov/BES/TCI/ECEI)
│       │   │   └── nmode_spectrum_tab.py       # Toroidal n-mode spectrum (Mirnov)
│       │   └── imaging/
│       │       ├── tv_tab.py                   # TV visible-camera image viewer
│       │       ├── tv_startup_tab.py           # TV startup comparison
│       │       ├── tv_utils.py                 # TV image helpers
│       │       └── irvb_tab.py                 # IRVB 2D Prad viewer
│       ├── biprofile/
│       │   ├── biprofile_profile_tab.py        # BiProfile Ti/vT, ne/Te profiles
│       │   ├── biprofile_timetrace_tab.py      # BiProfile time traces
│       │   └── biprofile_derived_profile_tab.py # Derived quantities (neoclassical, transport)
│       ├── transp/
│       │   ├── transp_profile_tab.py         # CDF Output - Profiles
│       │   ├── transp_timetrace_tab.py       # CDF Output - Time Traces
│       │   └── transp_ufile_tab.py           # UFILE Input (single dispatch tab)
│       └── efit/
│           ├── efit_profile_tab.py           # GEQDSK profiles
│           ├── efit_timetrace_tab.py         # AEQDSK scalars
│           ├── efit_2d_tab.py                # 2D equilibrium contours
│           └── efit_pfile_tab.py             # p-file viewer
├── plotting/
│   └── plot_manager.py                # Plot management
├── setup_vendor.ps1                   # Office-PC bundler (Python + pip wheels into vendor/)
├── requirements.txt                   # Documented dependencies (vendor/ is authoritative)
├── README.md
└── CHANGELOG.md
```

## Requirements

- Python 3.8+
- PySide6 (Qt 6) - GUI framework
- matplotlib - Plotting
- numpy, scipy - Numerical computation
- scikit-learn - GPR fitting (optional)
- Pillow - Image processing (TV tab)
- netCDF4 - TRANSP CDF loading
- h5py, hdf5plugin - raw-Mirnov HDF5 archive (blosc-zstd)
- MDSplus - KSTAR data access (KFE internal)

## Note for External Users

This software is developed and configured for KSTAR tokamak at KFE. The MDS+ server addresses, data paths, and tree structures in `config/` are KSTAR-specific and require access to KFE internal servers (nkstar / ukstar).

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

## Copyright

Copyright © 2026 Korea Institute of Fusion Energy. All rights reserved.

Korea Copyright Commission Registration No. C-2026-023829 (18 May 2026).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

This research was supported by the R&D Programs "High Performance Tokamak Plasma Research & Development (EN2601-17)" and "Korea-US Collaboration Research for High Performance Plasma on Tungsten Divertor (EN2603-02)" through the Korea Institute of Fusion Energy (KFE), and by the National Research Foundation (NRF) (No. RS-2026-25545529), funded by the Korean government (MSIT), Republic of Korea.