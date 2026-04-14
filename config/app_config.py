"""
PRISM - Plasma Research Integrated System for Multi-diagnostics
Global application configuration
"""

import os

VERSION = "2.4.1"
UPDATE_DATE = "2026-04-14"
APP_NAME = "PRISM"
APP_FULL_NAME = "Plasma Research Integrated System for Multi-diagnostics"

# Contact information
CONTACT_EMAIL = "jklee@kfe.re.kr"
AUTHOR_NAME = "Jekil Lee"

# Application root directory (where config/ lives → parent)
_APP_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class AppConfig:
    """Global application settings"""
    def __init__(self):
        # MDS+ connection
        self.MDS_IP = 'mdsr.kstar.kfe.re.kr:8005'

        # IRVB data server
        self.IRVB_SERVER = 'http://172.17.112.125/data_ana'

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
