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
elif [[ "$SERVER_HOST" == ukstar* ]]; then
    PRISM_HOME="/UKSTAR_HOME/jklee/PRISM"
else
    echo "========================================================"
    echo "  PRISM: Unknown server '$SERVER_HOST'"
    echo "  PRISM is configured for nkstar and ukstar only."
    echo "  Please check your server or contact Jekil Lee (jklee@kfe.re.kr)"
    echo "========================================================"
    exit 1
fi

# Bundled Python interpreter + third-party packages (no system dependency).
PYTHON_PATH="$PRISM_HOME/vendor/cpython-3.8/bin/python3.8"

# Check if directory exists
if [ ! -d "$PRISM_HOME" ]; then
    echo "Error: PRISM directory not found at $PRISM_HOME"
    exit 1
fi

# Check if bundled Python exists
if [ ! -x "$PYTHON_PATH" ]; then
    echo "Error: bundled Python not found at $PYTHON_PATH"
    echo "       run setup_vendor.ps1 on the office PC and redeploy."
    exit 1
fi

# Check if main.py exists
if [ ! -f "$PRISM_HOME/main.py" ]; then
    echo "Error: main.py not found in $PRISM_HOME"
    exit 1
fi

# PYTHONPATH only needs PRISM_HOME (third-party packages live inside the
# bundled Python's own site-packages, which Python finds automatically).
# PYTHONNOUSERSITE: hard-disable ~/.local/lib/python3.8/site-packages so the
# bundled stack is the only one loaded, regardless of who runs PRISM.
export PYTHONPATH="$PRISM_HOME"
export PYTHONNOUSERSITE=1

# Headless batch archive root (n-mode result cache). Lives under the shared
# /tmp/prism tree (alongside /tmp/prism/irvb). Override by exporting
# PRISM_ARCHIVE_ROOT before launching. (/tmp is ephemeral — a cache, not storage.)
export PRISM_ARCHIVE_ROOT="${PRISM_ARCHIVE_ROOT:-/tmp/prism}"

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
    echo "A multi-diagnostics viewer for the KSTAR tokamak. Launches a GUI viewer,"
    echo "or runs a headless n-mode batch job. Usage: prism [option | nmode args]"
    echo ""
    echo "Quick start:"
    echo "  prism                                      # Full GUI (all viewers)"
    echo "  prism -s                                   # Pick one viewer to launch"
    echo "  prism nmode --shot 40848 --tmin 2 --tmax 8 # Headless n-mode -> cache"
    echo "  prism -h                                   # This help"
    echo ""
    echo "GUI viewers (opens the PySide6 app):"
    echo "  (none)            Full PRISM (Diagnostics + EFIT + BiProfile + TRANSP)"
    echo "  -d, --diag        Diagnostic data viewer"
    echo "  -e, --efit        EFIT viewer"
    echo "  -b, --biprofile   BiProfile viewer"
    echo "  -t, --transp      TRANSP CDF viewer"
    echo "  -s, --select      Select and launch individual viewers"
    echo "  -h, --help        Show this help"
    echo ""
    echo "  Available tabs:"
    echo "    Diagnostics"
    echo "      Profiles      Ti/vT (CES), ne/Te (Thomson+ECE), MSE (gamma/q/j)"
    echo "      Time Traces   Ti/vT, ne/Te, MSE, Neutron (Fission/He3/Diamond)"
    echo "      Spectral      Spectrogram (ECE/ECEI/BES/Mirnov), n-Mode Spectrum"
    echo "      Imaging       TV Viewer, TV Startup, IRVB"
    echo "    EFIT            Time Traces, Profiles, 2D, p-File"
    echo "    BiProfile       Ti/vT, ne/Te Profiles & Time Traces"
    echo "    TRANSP"
    echo "      Input         UFILE viewer (single tab)"
    echo "      Output        CDF Profiles & Time Traces"
    echo ""
    echo "  Keyboard:  Delete/Backspace = remove selected list items"
    echo "             Enter            = fetch data (in shot input field)"
    echo ""
    echo "Headless batch (no GUI, no X11):  prism nmode [args]"
    echo "  Compute a toroidal n-mode spectrum and save it as .npz (one file per shot"
    echo "  under the archive root). If a shot's file exists you're asked to"
    echo "  overwrite or skip ([a]ll / [s]kip-all / [q]uit for batches)."
    echo "    --shot N | --shots N N    single | multiple shots"
    echo "    --shot-range A B          inclusive range, e.g. 10000 10100"
    echo "    --tmin / --tmax           time window [s]       (default: full shot)"
    echo "    --fmin / --fmax           frequency range [kHz] (default 0 / 100)"
    echo "    --t-interval S            FFT window length [s]  (default 0.01)"
    echo "    --msign 0|1|-1|2          abs | pos | neg | all  (default 1=pos)"
    echo "    --nmodes N                number of modes        (default 5)"
    echo "    --tol T                   mode-fit residual tolerance (default 0.8)"
    echo "    --frac F                  amplitude threshold fraction (default 0.01)"
    echo "    --integrate               integrate dB/dt -> B"
    echo "    --no-detrend              disable per-window detrend"
    echo "  Output: \$PRISM_ARCHIVE_ROOT/nmode/nmode_<shot>.npz  (default /tmp/prism)"
    echo "  Full per-flag help: prism nmode -h"
    echo "  Examples:"
    echo "    prism nmode --shot 40848                          # full shot, defaults"
    echo "    prism nmode --shots 40848 40850 --tmin 2 --tmax 8 # multi-shot window"
    echo "    prism nmode --shot-range 10000 10100              # 101 shots (full shot)"
    echo "    prism nmode --shot 40848 --msign -1 --integrate   # neg-n, dB/dt -> B"
    echo ""
    echo "Data sources:"
    echo "  MDS+   mdsr.kstar.kfe.re.kr:8005 (CES, Thomson, ECE, MSE, Neutron, ...)"
    echo "  HTTP   172.17.112.125 (IRVB)"
    echo "  File   TV images, TRANSP UFILE dirs, TRANSP CDF, EFIT g/a/p-files"
    echo ""
    echo "Settings:  ~/.config/prism/settings.json"
    echo "Env:       PRISM_ARCHIVE_ROOT = n-mode cache root (default /tmp/prism)"
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
    nmode)
        # Headless n-mode compute + cache (no GUI). Pass remaining args to the
        # batch CLI; do NOT suppress stderr so progress/errors are visible.
        cd "$PRISM_HOME"
        $PYTHON_PATH -m batch.cli "$@"
        ;;
    *)
        echo "Error: Unknown option '$1'"
        echo "Use 'prism --help' for available options."
        exit 1
        ;;
esac
