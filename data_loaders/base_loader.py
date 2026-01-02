#!/usr/bin/python3.8

"""
Abstract base class for diagnostic data loaders
"""

from abc import ABC, abstractmethod
from MDSplus import Connection


class BaseDiagnosticLoader(ABC):
    """Abstract base class for all diagnostic loaders"""
    
    def __init__(self, config, diagnostic_config):
        self.app_config = config
        self.diag_config = diagnostic_config
        self.mds_ip = config.MDS_IP
    
    @abstractmethod
    def load_data(self, shot_number, analysis_type=None):
        """Load diagnostic data from MDS+
        
        analysis_type: Optional analysis method (e.g., 'mod', 'nn' for CES)
        """
        pass
    
    def _connect_mds(self, shot_number):
        """Helper to establish MDS+ connection"""
        mds = Connection(self.mds_ip)
        mds.openTree(self.diag_config['mds_tree'], shot_number)
        return mds
    
    def _close_mds(self, mds, shot_number):
        """Helper to close MDS+ connection"""
        mds.closeTree(self.diag_config['mds_tree'], shot_number)