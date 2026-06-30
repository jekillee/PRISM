#!/usr/bin/python3

"""
PRISM v2.6.2
Main entry point

Plasma Research Integrated System for Multi-diagnostics
A modular viewer for KSTAR diagnostic data including:
- CES (Charge Exchange Spectroscopy)
- Thomson Scattering
- ECE (Electron Cyclotron Emission)
- MSE (Motional Stark Effect)
- Neutron (Fission Chamber, He3 Counter, Diamond Detector)
- Spectrogram Analysis
- TV Image Viewer
- IRVB (Infra-Red Video Bolometer)
- TRANSP (Transport Analysis Profile Fitting)

Usage:
    python main.py              # Full PRISM (all tabs with sidebar)
    python main.py --select     # Tab selector (choose one viewer)
    python main.py -s           # Tab selector (short form)
    python main.py --transp     # TRANSP viewer
    python main.py transp       # TRANSP viewer (short form)

Author: Jekil Lee (jklee@kfe.re.kr)
"""

import os
import sys
import warnings

# Sanity-check that the bundled Python is the one running us — i.e. that
# run_prism.sh launched us via vendor/cpython-3.8/bin/python3.8. If someone
# invokes main.py with a system Python, the bundled third-party packages
# (PySide6, numpy, netCDF4, ...) won't be on sys.path and PRISM will silently
# fall back to whatever the system Python provides. Fail fast instead.
_PRISM_DIR = os.path.dirname(os.path.abspath(__file__))
_BUNDLED_PY = os.path.join(_PRISM_DIR, 'vendor', 'cpython-3.8', 'bin', 'python3.8')
if not os.path.exists(_BUNDLED_PY):
    sys.stderr.write(
        f"[PRISM] FATAL: bundled Python missing.\n"
        f"[PRISM]   expected: {_BUNDLED_PY}\n"
        f"[PRISM]   fix:      run setup_vendor.ps1 on the office PC, then\n"
        f"[PRISM]             rsync the PRISM directory (with vendor/) to the server.\n"
    )
    sys.exit(1)

# Suppress numpy compiletime version mismatch warnings (defensive)
warnings.filterwarnings("ignore", message=".*compiletime version.*", category=RuntimeWarning)
warnings.filterwarnings("ignore", message=".*binary incompatibility.*", category=RuntimeWarning)

# Suppress matplotlib tight_layout warnings (IRVB colorbar axes are incompatible)
warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)

# X11 forwarding compatibility (e.g. Bitvise SSH + Xming, NoMachine)
# - Disable GLX (OpenGL) to prevent GLX initialization crash on remote X servers
# - Disable XInput 2 to prevent XCB BadRequest errors (major code 130)
# - Disable MIT-SHM shared memory (XCB error 170) — Xming does not support it
# - Force software OpenGL rendering to avoid GPU-related X11 crashes
# - Unset WAYLAND_DISPLAY to suppress "Ignoring WAYLAND_DISPLAY on Gnome" warning
os.environ.setdefault('QT_XCB_GL_INTEGRATION', 'none')
os.environ.setdefault('QT_XCB_NO_XI2', '1')
os.environ.setdefault('QT_X11_NO_MITSHM', '1')
os.environ.setdefault('LIBGL_ALWAYS_SOFTWARE', '1')
os.environ.setdefault('QT_QUICK_BACKEND', 'software')
if 'WAYLAND_DISPLAY' in os.environ:
    del os.environ['WAYLAND_DISPLAY']

# Suppress MDSplus debug messages (buffer_free, etc.)
os.environ.setdefault('MDSPLUS_DEBUG', '0')

# Force Qt to use a fixed 96 DPI for font metrics so PRISM looks identical
# across access paths whose X servers report different DPIs (e.g. NoMachine's
# server-side Xorg = 100 dpi vs MobaXterm SSH X11 forward = 96 dpi).
os.environ.setdefault('QT_FONT_DPI', '96')

from PySide6.QtWidgets import QApplication
from ui.theme import ThemeManager


def create_app():
    """Create and configure QApplication with saved theme"""
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    # Pin glyph rendering: load PRISM-bundled DejaVu Sans (shipped with
    # matplotlib) and force it as the default app font with hinting/AA
    # hints chosen to eliminate path-dependent variation between NoMachine
    # direct and SSH X11 forwarded sessions.
    from PySide6.QtGui import QFontDatabase, QFont
    _font_path = os.path.join(
        _PRISM_DIR, 'vendor', 'cpython-3.8', 'lib', 'python3.8',
        'site-packages', 'matplotlib', 'mpl-data', 'fonts', 'ttf',
        'DejaVuSans.ttf')
    if os.path.isfile(_font_path):
        QFontDatabase.addApplicationFont(_font_path)
        f = QFont('DejaVu Sans', 10)
        f.setHintingPreference(QFont.PreferNoHinting)
        f.setStyleStrategy(QFont.PreferAntialias)
        app.setFont(f)

    # Load saved theme from settings, default to 'dark'
    from config.user_settings import load_settings, get_theme
    load_settings()
    saved_theme = get_theme()
    ThemeManager.apply_theme(saved_theme)

    return app


def _print_startup():
    from config.app_config import VERSION, UPDATE_DATE
    v = f"v{VERSION} ({UPDATE_DATE})"
    print(f"""
\033[91m  ██████╗ \033[33m██████╗ \033[93m██╗\033[92m███████╗\033[94m███╗   ███╗\033[0m
\033[91m  ██╔══██╗\033[33m██╔══██╗\033[93m██║\033[92m██╔════╝\033[94m████╗ ████║\033[0m
\033[91m  ██████╔╝\033[33m██████╔╝\033[93m██║\033[92m███████╗\033[94m██╔████╔██║\033[0m
\033[91m  ██╔═══╝ \033[33m██╔══██╗\033[93m██║\033[92m╚════██║\033[94m██║╚██╔╝██║\033[0m
\033[91m  ██║     \033[33m██║  ██║\033[93m██║\033[92m███████║\033[94m██║ ╚═╝ ██║\033[0m
\033[91m  ╚═╝     \033[33m╚═╝  ╚═╝\033[93m╚═╝\033[92m╚══════╝\033[94m╚═╝     ╚═╝\033[0m
  {v}

  \033[91mP\033[0mlasma \033[33mR\033[0mesearch \033[93mI\033[0mntegrated \033[92mS\033[0mystem for \033[94mM\033[0multi-diagnostics:
  A unified visualization and data analysis platform for KSTAR tokamak.

\033[37m  Developed by Jekil Lee (jklee@kfe.re.kr)
  Copyright © 2026 Korea Institute of Fusion Energy. All rights reserved.
  Korea Copyright Commission Registration No. C-2026-023829 (18 May 2026).\033[0m

\033[37m  If you use PRISM in your work, please cite:
    J.K. Lee, "An integrated multi-diagnostic visualization platform for
    KSTAR tokamak," Fusion Engineering and Design, Vol. 228, 2026, 115786.
    https://doi.org/10.1016/j.fusengdes.2026.115786\033[0m

\033[37m  For usage and command-line options, run `prism -h`.\033[0m
""")


def main():
    """Main entry point with mode selection"""
    # Parse command line argument
    arg = sys.argv[1] if len(sys.argv) > 1 else ''

    _print_startup()

    app = create_app()

    from ui.main_window import PRISMApp

    if arg in ('-d', '--diag', 'diag'):
        window = PRISMApp(mode='diag')
    elif arg in ('-b', '--biprofile', 'biprofile'):
        window = PRISMApp(mode='biprofile')
    elif arg in ('-t', '--transp', 'transp'):
        window = PRISMApp(mode='transp')
    elif arg in ('-e', '--efit', 'efit'):
        window = PRISMApp(mode='efit')
    elif arg in ('-s', '--select'):
        window = PRISMApp(mode='select')
    else:
        window = PRISMApp(mode='')

    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
