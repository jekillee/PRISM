#!/usr/bin/python3.8

"""
Abstract base class for diagnostic data loaders
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
from MDSplus import Connection

from core.data_structures import DiagnosticData


class BaseDiagnosticLoader(ABC):
    """Abstract base class for all diagnostic loaders

    Provides common MDS+ connection handling and defines the interface
    for all diagnostic-specific loaders.
    """

    def __init__(
        self,
        config: Any,
        diagnostic_config: Dict[str, Any]
    ) -> None:
        """Initialize the loader with configuration

        Args:
            config: Application configuration object (AppConfig)
            diagnostic_config: Diagnostic-specific configuration dictionary
        """
        self.app_config = config
        self.diag_config = diagnostic_config
        self.mds_ip: str = config.MDS_IP

    @abstractmethod
    def load_data(
        self,
        shot_number: int,
        analysis_type: Optional[str] = None
    ) -> DiagnosticData:
        """Load diagnostic data from MDS+

        Args:
            shot_number: KSTAR shot number
            analysis_type: Optional analysis method (e.g., 'mod', 'nn' for CES)

        Returns:
            DiagnosticData object containing the loaded measurements

        Raises:
            RuntimeError: If data loading fails
        """
        pass

    def _connect_mds(self, shot_number: int) -> Connection:
        """Establish MDS+ connection and open the diagnostic tree

        Args:
            shot_number: KSTAR shot number

        Returns:
            MDSplus Connection object with tree opened
        """
        mds = Connection(self.mds_ip)
        mds.openTree(self.diag_config['mds_tree'], shot_number)
        return mds

    def _close_mds(self, mds: Connection, shot_number: int) -> None:
        """Close MDS+ connection

        Args:
            mds: MDSplus Connection object
            shot_number: KSTAR shot number
        """
        mds.closeTree(self.diag_config['mds_tree'], shot_number)
