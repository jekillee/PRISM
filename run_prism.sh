#!/bin/bash
#
# PRISM Launcher Script
# Plasma Research Integrated System for Multi-diagnostics
#
# Usage:
#   prism              # Full PRISM (all tabs with sidebar)
#   prism -s           # Tab selector (choose one viewer to launch)
#   prism --help       # Show help
#
# Author: Jekil Lee (jklee@kfe.re.kr)
#

# Detect server and set paths accordingly
SERVER_HOST=$(hostname)

if [[ "$SERVER_HOST" == nkstar* ]]; then
    PRISM_HOME="/home/users/jklee/PRISM"
    PYTHON_PATH="/usr/bin/python38"
    PYSIDE6_SITE="/home/users/jklee/.local/lib/python3.8/site-packages"
elif [[ "$SERVER_HOST" == ukstar* ]]; then
    PRISM_HOME="/UKSTAR_HOME/jklee/PRISM"
    PYTHON_PATH="/usr/bin/python3.8"
    PYSIDE6_SITE="/UKSTAR_HOME/jklee/.local/lib/python3.8/site-packages"
else
    echo "========================================================"
    echo "  PRISM: Unknown server '$SERVER_HOST'"
    echo "  PRISM is configured for nkstar and ukstar only."
    echo "  Please check your server or contact Jekil Lee (jklee@kfe.re.kr)"
    echo "========================================================"
    exit 1
fi

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

# Set PYTHONPATH to PRISM directory and jklee's site-packages only
export PYTHONPATH="$PRISM_HOME:$PYSIDE6_SITE"

# Suppress WAYLAND_DISPLAY warning on Gnome (NoMachine)
unset WAYLAND_DISPLAY

# Fix X11 errors when using remote X servers (Xming, VcXsrv, etc.)
# - Disable MIT-SHM shared memory (XCB error 170)
# - Force software OpenGL rendering
# - Disable accessibility bus warnings
export QT_X11_NO_MITSHM=1
export QT_XCB_GL_INTEGRATION=none
export LIBGL_ALWAYS_SOFTWARE=1
export QT_QUICK_BACKEND=software
export NO_AT_BRIDGE=1

# Show help
show_help() {
    # Extract version from app_config.py
    local version
    version=$(grep -oP 'VERSION\s*=\s*"\K[^"]+' "$PRISM_HOME/config/app_config.py" 2>/dev/null || echo "unknown")

    echo ""
    echo "PRISM v${version} - Plasma Research Integrated System for Multi-diagnostics"
    echo ""
    echo "Usage: prism [option]"
    echo ""
    echo "Options:"
    echo "  (none)                Full PRISM (Diagnostics + EFIT + BiProfile + TRANSP)"
    echo "  -d, --diag            Diagnostic data viewer"
    echo "  -e, --efit            EFIT viewer"
    echo "  -b, --biprofile       BiProfile viewer"
    echo "  -t, --transp          TRANSP CDF viewer"
    echo "  -s, --select          Select and launch individual viewers"
    echo "  -h, --help            Show this help message"
    echo ""
    echo "Available Tabs:"
    echo "  Diagnostics"
    echo "    Profiles      Ti/vT (CES), ne/Te (Thomson+ECE), MSE (gamma/q/j)"
    echo "    Time Traces   Ti/vT, ne/Te, MSE, Neutron (Fission/He3/Diamond)"
    echo "    Spectral      Spectrogram (ECE/ECEI/BES/Mirnov), n-Mode Spectrum"
    echo "    Imaging       TV Viewer, TV Startup, IRVB"
    echo "  EFIT            Time Traces, Profiles, 2D, p-File"
    echo "  BiProfile       Ti/vT, ne/Te Profiles & Time Traces"
    echo "  TRANSP          Ti/vT, ne/Te Profiles & Time Traces (CDF)"
    echo ""
    echo "Keyboard Shortcuts:"
    echo "  Delete/Backspace  Remove selected items from list"
    echo "  Enter             Fetch data (in shot input field)"
    echo ""
    echo "Data Sources:"
    echo "  MDS+    mdsr.kstar.kfe.re.kr:8005 (CES, Thomson, ECE, MSE, Neutron, ...)"
    echo "  HTTP    172.17.112.125 (IRVB)"
    echo "  File    TV images, TRANSP CDF, EFIT g/a/p-files"
    echo ""
    echo "Settings: ~/.config/prism/settings.json"
    echo ""
    echo "Examples:"
    echo "  prism              # Full PRISM (all groups)"
    echo "  prism -d           # Diagnostic data viewer"
    echo "  prism -e           # EFIT viewer"
    echo "  prism -b           # BiProfile viewer"
    echo "  prism -t           # TRANSP CDF viewer"
    echo "  prism -s           # Select and launch individual viewers"
    echo ""
    echo "Jekil Lee (jklee@kfe.re.kr)"
    echo ""
}

# Parse arguments
case "$1" in
    --help|-h)
        show_help
        exit 0
        ;;
    -s|--select)
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py --select 2>/dev/null
        ;;
    -d|--diag|diag|"")
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py 2>/dev/null
        ;;
    -b|--biprofile|biprofile)
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py --biprofile 2>/dev/null
        ;;
    -t|--transp|transp)
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py --transp 2>/dev/null
        ;;
    -e|--efit|efit)
        cd "$PRISM_HOME"
        $PYTHON_PATH main.py --efit 2>/dev/null
        ;;
    *)
        echo "Error: Unknown option '$1'"
        echo "Use 'prism --help' for available options."
        exit 1
        ;;
esac
