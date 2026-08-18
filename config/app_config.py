"""
PRISM - Plasma Research Integrated System for Multi-diagnostics
Global application configuration
"""

import os

VERSION = "2.6.5"
UPDATE_DATE = "2026-08-18"
APP_NAME = "PRISM"
APP_FULL_NAME = "Plasma Research Integrated System for Multi-diagnostics"

# Contact information
CONTACT_EMAIL = "jklee@kfe.re.kr"
AUTHOR_NAME = "Jekil Lee"

# Application root directory (where config/ lives → parent)
_APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Shared scratch root on the servers (nkstar/ukstar). Everything PRISM caches or
# computes outside the GUI lives here under a per-subsystem subdir:
#   <PRISM_TMP_ROOT>/nmode/<shot>/...   headless n-mode archive (batch SDK)
#   <PRISM_TMP_ROOT>/irvb/...           IRVB .mat download cache
# Override with the PRISM_TMP_ROOT env var. World-writable so multiple users share
# one tree. (Note: /tmp is ephemeral — cleared on reboot; this is a cache, not
# durable storage.)
PRISM_TMP_ROOT = os.environ.get('PRISM_TMP_ROOT', '/tmp/prism')

# Durable results root in the user's HOME directory (~/prism_results). Kept
# separate from the PRISM install dir (~/PRISM) so results never mix with code.
# Batch/API n-mode outputs are written under a per-subsystem subdir here:
#   <PRISM_RESULTS_ROOT>/nmode/nmode_<shot>.npz   ->  ~/prism_results/nmode/nmode_<shot>.npz
# Override the whole root with the PRISM_RESULTS_ROOT env var.
PRISM_RESULTS_ROOT = os.environ.get('PRISM_RESULTS_ROOT', os.path.expanduser('~/prism_results'))

# Shared PRISM NAS root (mounted at /PRISM on nkstar & ukstar). Shared caches that
# are not per-user live here, e.g. the IRVB .mat cache (<PRISM_NAS_ROOT>/irvb) and
# the raw-Mirnov archive. Override the mount point with the PRISM_NAS_ROOT env var.
PRISM_NAS_ROOT = os.environ.get('PRISM_NAS_ROOT', '/PRISM')


def ensure_shared_dir(path, mode=0o777):
    """Create `path` (and parents) and best-effort chmod to `mode` so multiple
    users on a shared server can read/write one /tmp/prism. chmod failures (e.g.
    not the directory owner) are ignored. Returns the path."""
    os.makedirs(path, exist_ok=True)
    try:
        os.chmod(path, mode)
    except OSError:
        pass
    return path


class AppConfig:
    """Global application settings"""
    def __init__(self):
        # MDS+ connection
        self.MDS_IP = 'mdsr.kstar.kfe.re.kr:8005'

        # IRVB data server (Flask app since 2025; static .mat files live under
        # /data/ for recent shots and /data/afterCampaign/ for post-campaign shots).
        self.IRVB_SERVER = 'http://172.17.112.125/data'

        # EFIT tree options
        self.EFIT_TREES = {
            "efitrt1 (RT for PCS)": "efitrt1",
            "efitrt2 (RT with MSE)": "efitrt2",
            "efit01 (MAG)": "efit01",
            "efit02 (MSE)": "efit02",
            "efit04 (MAG+drift)": "efit04"
        }
        self.DEFAULT_EFIT_TREE = "efitrt1"

        # Paths
        self.DOCS_PATH = os.path.join(_APP_ROOT, 'docs')
        self.CES_RESULT_PATH = '/home/users/jklee/CESresults/'

        # Plot configuration
        self.FIGURE_SIZE = (10, 6)
        self.R_LIMITS = (1.8, 2.3)
        self.R_EDGE = 2.2

        # Unit conversion factors
        self.R_SCALE = 1e-3   # mm to m
        self.TI_SCALE = 1e-3  # eV to keV
        self.TE_SCALE = 1e-3  # eV to keV
        self.NE_SCALE = 1e-19 # m-3 to 1e19m-3
