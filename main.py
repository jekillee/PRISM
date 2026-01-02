#!/usr/bin/python3.8

"""
PRISM v1.0
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
    python main.py          # Full PRISM
    python main.py ces      # CES standalone
    python main.py nete     # ne,Te standalone
    python main.py mse      # MSE standalone
    python main.py spec     # Spectrogram standalone
    python main.py tv       # TV standalone
    python main.py irvb     # IRVB standalone

Author: Jekil Lee (jklee@kfe.re.kr)
"""

import sys


def main():
    """Main entry point with mode selection"""
    # Get mode from command line argument
    mode = sys.argv[1] if len(sys.argv) > 1 else ''
    
    if mode == '' or mode == 'all':
        # Full PRISM
        from ui.main_window import PRISMApp
        app = PRISMApp()
        app.run()
    else:
        # Standalone mode
        from standalone_launcher import StandaloneLauncher
        launcher = StandaloneLauncher()
        launcher.run(mode)


if __name__ == "__main__":
    main()