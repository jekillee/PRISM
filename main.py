#!/usr/bin/python3

"""
PRISM v2.0
Main entry point

Plasma Research Integrated System for Multi-diagnostics
A modular viewer for KSTAR diagnostic data including:
- CES (Charge Exchange Spectroscopy)
- Thomson Scattering
- ECE (Electron Cyclotron Emission)
- MSE (Motional Stark Effect)
- Spectrogram Analysis
- TV Image Viewer
- IRVB (Infra-Red Video Bolometer)

Usage:
    python main.py              # Full PRISM (all tabs with sidebar)
    python main.py --select     # Tab selector (choose one viewer)
    python main.py -s           # Tab selector (short form)

Author: Jekil Lee (jklee@kfe.re.kr)
"""

import os
import sys
import warnings

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

from PySide6.QtWidgets import QApplication
from ui.theme import ThemeManager


def create_app():
    """Create and configure QApplication with saved theme"""
    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    # Load saved theme from settings, default to 'dark'
    from config.user_settings import load_settings, get_theme
    load_settings()
    saved_theme = get_theme()
    ThemeManager.apply_theme(saved_theme)

    return app


def main():
    """Main entry point with mode selection"""
    # Parse command line argument
    arg = sys.argv[1] if len(sys.argv) > 1 else ''

    app = create_app()

    from ui.main_window import PRISMApp

    if arg in ('-s', '--select'):
        window = PRISMApp(mode='select')
    else:
        window = PRISMApp(mode='')

    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
