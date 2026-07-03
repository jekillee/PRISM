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

# n-mode batch results default to ~/prism_results/nmode/ (config.app_config
# PRISM_RESULTS_ROOT). We intentionally do NOT force PRISM_ARCHIVE_ROOT here so
# that per-user default applies; export PRISM_ARCHIVE_ROOT before launching to
# override. (Raw Mirnov archive: /PRISM/mirnov_archive; IRVB cache: /PRISM/irvb.)

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
    echo "or runs a headless batch job. Usage: prism [option | nmode args | irvb args]"
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
    echo "Headless batch (no GUI / X11) - compute + cache one .npz per shot:"
    echo "  prism nmode [args]   toroidal n-mode spectrum (Mirnov coils)"
    echo "  prism irvb  [args]   IRVB 2D radiation + regional Prad"
    echo "  Shots:    --shot N  |  --shots N N  |  --shot-range A B"
    echo "  Flags:    prism nmode -h  /  prism irvb -h    (full per-command list)"
    echo "  Output:   \$PRISM_ARCHIVE_ROOT/{nmode,irvb}/<sub>_<shot>.npz  (default ~/prism_results)"
    echo "  Existing: prompt to overwrite / skip ([a]ll / [s]kip-all / [q]uit)"
    echo "  Examples: prism nmode --shot 40848 --tmin 2 --tmax 8"
    echo "            prism irvb  --shot 40085 --psi-boundaries 0.3 0.7 1.0"
    echo ""
    echo "Python SDK (same compute; run() also returns the result as a variable):"
    echo "  from batch import NModeJobSpec, IRVBJobSpec, run, run_many"
    echo "  result, status = run(IRVBJobSpec(shot=40085))    # print(result) lists all fields"
    echo "  Run a script in PRISM's env (no path setup):  prism python my_analysis.py"
    echo "  Runnable examples (lists ready-to-run commands):  prism examples"
    echo "  API reference:     prism python -c 'import batch; help(batch)'"
    echo ""
    echo "Data sources:"
    echo "  MDS+   mdsr.kstar.kfe.re.kr:8005 (CES, Thomson, ECE, MSE, Neutron, ...)"
    echo "  HTTP   172.17.112.125 (IRVB)"
    echo "  File   TV images, TRANSP UFILE dirs, TRANSP CDF, EFIT g/a/p-files"
    echo ""
    echo "Settings:  ~/.config/prism/settings.json"
    echo "Archives:"
    echo "    n-mode results    ~/prism_results/nmode/      (per-user; \$PRISM_ARCHIVE_ROOT overrides)"
    echo "    IRVB results      ~/prism_results/irvb/       (per-user; \$PRISM_ARCHIVE_ROOT overrides)"
    echo "    raw Mirnov        /PRISM/mirnov_archive/     (shared, auto-saved on first load)"
    echo "    IRVB cache        /PRISM/irvb/               (shared, GUI reconstruction cache)"
    echo "Env:       PRISM_ARCHIVE_ROOT = batch result root (default ~/prism_results)"
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
    nmode|irvb)
        # Headless n-mode / IRVB compute + cache (no GUI). Pass all args (incl. the
        # subcommand) to the batch CLI; do NOT suppress stderr so progress/errors
        # are visible.
        cd "$PRISM_HOME"
        $PYTHON_PATH -m batch.cli "$@"
        ;;
    examples)
        # List the bundled batch-API example scripts as ready-to-run commands
        # (absolute paths, so they work from any directory).
        echo "PRISM batch-API examples ($PRISM_HOME/examples):"
        found=0
        for f in "$PRISM_HOME"/examples/*.py; do
            [ -e "$f" ] || continue
            echo "  prism python $f"
            found=1
        done
        [ "$found" -eq 0 ] && echo "  (none found)"
        echo ""
        echo "Run one, or copy it and edit: prism python /path/to/your/copy.py"
        ;;
    python|exec)
        # Run a user script (or -c "...") in PRISM's environment: bundled Python
        # + PYTHONPATH so `import batch` and all bundled deps resolve, without the
        # user needing to know where PRISM is installed. Runs from the current
        # directory, so relative paths in the script still work.
        #   prism python my_analysis.py
        #   prism python -c "import batch; help(batch)"
        shift
        exec "$PYTHON_PATH" "$@"
        ;;
    *)
        echo "Error: Unknown option '$1'"
        echo "Use 'prism --help' for available options."
        exit 1
        ;;
esac
