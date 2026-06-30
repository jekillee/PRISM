# Populate PRISM's vendor/ directory with all third-party Python dependencies
# pinned to Linux x86_64 + Python 3.8 (nkstar / ukstar target environment).
#
# Run this once on the office PC after first checkout, or whenever
# requirements.txt changes. The downloaded wheels are cross-platform-targeted
# (manylinux2014_x86_64 + Python 3.8) so they install on the Linux servers
# even though we ran pip on Windows.
#
# Usage (Windows PowerShell 5.1 — the default Windows PowerShell):
#   Run as Administrator (right-click PowerShell -> Run as administrator)
#   so the Linux symlinks inside python-build-standalone's tarball can be
#   extracted onto NTFS. Then in the elevated session:
#
#   cd D:\Research\Code\Code Scripts\PRISM\v2.6.0
#   .\setup_vendor.ps1
#
# If you hit an execution-policy error, run once per session:
#   Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
#
# After it finishes, transfer the entire PRISM directory (including vendor/)
# to nkstar via scp/rsync.

$ErrorActionPreference = 'Stop'

$PRISM_DIR  = Split-Path -Parent $MyInvocation.MyCommand.Path
$VENDOR_ROOT = Join-Path $PRISM_DIR 'vendor'
$PY_DIR     = Join-Path $VENDOR_ROOT 'cpython-3.8'
$VENDOR_DIR = Join-Path $PY_DIR 'lib\python3.8\site-packages'

Write-Output "[setup_vendor] PRISM dir:  $PRISM_DIR"
Write-Output "[setup_vendor] Python dir: $PY_DIR"
Write-Output "[setup_vendor] Site-pkgs:  $VENDOR_DIR"
New-Item -ItemType Directory -Force -Path $VENDOR_ROOT | Out-Null

# Bundle Python 3.8 interpreter for Linux x86_64 (python-build-standalone).
# REQUIREMENT: PowerShell must be running as Administrator (or Windows Developer
# Mode enabled). The tarball contains Linux symlinks (python3 -> python3.8,
# libpython3.8.so -> ...so.1.0); NTFS rejects symlink creation without
# elevated privileges. Without admin, tar will print "Invalid argument" errors.
if (-not (Test-Path $PY_DIR)) {
    # Last python-build-standalone release including Python 3.8.20 for
    # Linux x86_64 (Python 3.8 reached EOL Oct 2024 and was dropped after this).
    $PY_URL = 'https://github.com/astral-sh/python-build-standalone/releases/download/20241002/cpython-3.8.20+20241002-x86_64-unknown-linux-gnu-install_only.tar.gz'
    $Tarball = Join-Path $env:TEMP 'prism_cpython38.tar.gz'

    Write-Output "[setup_vendor] Downloading portable Python 3.8.20 (Linux x86_64)..."
    Invoke-WebRequest -Uri $PY_URL -OutFile $Tarball

    # Clean up any partial extract from a previous failed run.
    $PartialDir = Join-Path $VENDOR_ROOT 'python'
    if (Test-Path $PartialDir) {
        Write-Output "[setup_vendor] Removing partial extract: $PartialDir"
        Remove-Item -Recurse -Force $PartialDir
    }

    # Exclude `python/share/terminfo` — ncurses terminfo contains filenames
    # that NTFS rejects (Windows reserved names like 'aux', 'con', ...).
    # PRISM is a GUI app and never needs terminfo, so skipping it is safe.
    Write-Output "[setup_vendor] Extracting Python (admin required for symlinks; terminfo excluded) ..."
    tar -xzf $Tarball -C $VENDOR_ROOT --exclude='python/share/terminfo' --exclude='python/share/terminfo/*'
    if ($LASTEXITCODE -ne 0) {
        throw "tar extract failed -- run PowerShell as Administrator (or enable Windows Developer Mode) so symlinks can be created."
    }

    # tarball expands to vendor\python\ — rename for clarity
    Rename-Item -Path $PartialDir -NewName 'cpython-3.8'
    Remove-Item $Tarball

    Write-Output "[setup_vendor] Python bundled at: $PY_DIR"
} else {
    Write-Output "[setup_vendor] Python already bundled, skipping download."
}

# Strict platform pinning so the wheels run on both nkstar and ukstar
# (CentOS / RHEL 7+, glibc >= 2.17). manylinux2014 is the safe baseline.
$PY_VER   = '3.8'
$PLATFORM = 'manylinux2014_x86_64'

$PipFlags = @(
    "--target=$VENDOR_DIR",
    "--python-version", $PY_VER,
    "--platform",       $PLATFORM,
    "--only-binary=:all:",
    "--upgrade"
)

$CorePackages = @(
    'PySide6',
    'numpy',
    'scipy',
    'matplotlib',
    'Pillow',
    'netCDF4',
    'scikit-learn',
    'h5py',                  # raw-Mirnov NAS archive (manylinux wheel bundles libHDF5; gzip built-in)
    'hdf5plugin',            # blosc-zstd filter for the archive (default compression)
    # Conditional deps that pip cross-install (`--python-version 3.8` on
    # a Windows 3.10 host) skips because environment markers are evaluated
    # against the host Python, not the target. Explicitly listing them
    # forces inclusion in vendor/.
    'importlib_resources',   # matplotlib requires for py<3.10
    'importlib_metadata',    # various packages on py<3.8 path
    'zipp',                  # dep of importlib_metadata
    'typing_extensions'      # various packages, py<3.10 conditional
)

Write-Output "[setup_vendor] Installing core packages..."
& pip install @PipFlags @CorePackages
if ($LASTEXITCODE -ne 0) {
    throw "pip install failed for core packages (exit $LASTEXITCODE)"
}

# Note: MDSplus sits in vendor/cpython-3.8/lib/python3.8/site-packages/MDSplus
# as-is. This script does not manage it.

Write-Output "[setup_vendor] Done."
$size = (Get-ChildItem -Recurse $VENDOR_DIR | Measure-Object -Property Length -Sum).Sum / 1MB
Write-Output ("[setup_vendor] Vendor size: {0:N1} MB" -f $size)
