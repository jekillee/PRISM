#!/bin/bash
#
# PRISM Launcher Script
# Plasma Research Integrated System for Multi-diagnostics
#
# After symbolic link installation, use 'prism' command:
#   prism              # Full PRISM
#   prism tivt         # Ti, vT Profile + TimeTrace (CES + XICS)
#   prism nete         # ne, Te Profile + TimeTrace
#   prism mse          # MSE Profile + TimeTrace
#   prism spec         # Spectrogram
#   prism tv           # TV Viewer
#   prism irvb         # IRVB Viewer
#   prism --help       # Show help
#
# Author: Jekil Lee (jklee@kfe.re.kr)
#

# PRISM installation directory
PRISM_HOME="/home/users/jklee/PRISM"

# Python path
PYTHON_PATH="/usr/bin/python38"

# Check if directory exists
if [ ! -d "$PRISM_HOME" ]; then
    echo "Error: PRISM directory not found at $PRISM_HOME"
    exit 1
fi

# Check if main.py exists
if [ ! -f "$PRISM_HOME/main.py" ]; then
    echo "Error: main.py not found in $PRISM_HOME"
    exit 1
fi

# Set PYTHONPATH to include PRISM directory
export PYTHONPATH="$PRISM_HOME:$PYTHONPATH"

# Show help
show_help() {
    echo ""
    echo "PRISM - Plasma Research Integrated System for Multi-diagnostics"
    echo ""
    echo "Usage: prism [mode]"
    echo ""
    echo "Modes:"
    echo "  (none)    Launch full PRISM with all tabs"
    echo "  tivt      Ti, vT Profile + TimeTrace viewer (CES + XICS)"
    echo "  nete      ne, Te Profile + TimeTrace viewer (Thomson/ECE)"
    echo "  mse       MSE Profile + TimeTrace viewer"
    echo "  spec      Spectrogram analyzer"
    echo "  tv        TV image sequence viewer"
    echo "  irvb      IRVB 2D radiation profile viewer"
    echo "  --help    Show this help message"
    echo ""
    echo "Examples:"
    echo "  prism              # Launch full PRISM"
    echo "  prism tivt         # Launch Ti, vT standalone"
    echo "  prism nete         # Launch ne, Te standalone"
    echo "  prism mse          # Launch MSE standalone"
    echo "  prism spec         # Launch Spectrogram standalone"
    echo "  prism tv           # Launch TV Viewer standalone"
    echo "  prism irvb         # Launch IRVB Viewer standalone"
    echo ""
    echo "Author: Jekil Lee (jklee@kfe.re.kr)"
    echo ""
}

# Parse arguments
case "$1" in
    --help|-h)
        show_help
        exit 0
        ;;
    tivt|nete|mse|spec|tv|irvb|"")
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py "$1"
        ;;
    *)
        echo "Error: Unknown mode '$1'"
        echo "Use 'prism --help' for available options."
        exit 1
        ;;
esac
