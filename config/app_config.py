#!/usr/bin/python3.8

"""
PRISM - Plasma Research Integrated System for Multi-diagnostics
Global application configuration
"""

VERSION = "1.0.0"
APP_NAME = "PRISM"
APP_FULL_NAME = "Plasma Research Integrated System for Multi-diagnostics"


class AppConfig:
    """Global application settings"""
    def __init__(self):
        # MDS+ connection
        self.MDS_IP = 'mdsr.kstar.kfe.re.kr:8005'
        
        # EFIT tree options
        self.EFIT_TREES = {
            "efitrt1 (RT for PCS)": "efitrt1",
            "efitrt2 (RT with MSE)": "efitrt2",
            "efit01 (MAG)": "efit01",
            "efit02 (MSE)": "efit02",
            "efit04 (MAG+)": "efit04"
        }
        self.DEFAULT_EFIT_TREE = "efitrt1"
        
        # Paths
        self.MANUAL_PATH = "/home/users/jklee/PRISM/DiagnosticsManual(250814).pdf"
        self.CES_RESULT_PATH = '/home/users/jklee/CESresults/'
        
        # Plot configuration
        self.FIGURE_SIZE = (10, 5)
        self.R_LIMITS = (1.8, 2.3)
        self.R_EDGE = 2.2
        
        # Unit conversion factors
        self.R_SCALE = 1e-3   # mm to m
        self.TI_SCALE = 1e-3  # eV to keV
        self.TE_SCALE = 1e-3  # eV to keV
        self.NE_SCALE = 1e-19 # particles/m³ to 10¹⁹ particles/m³