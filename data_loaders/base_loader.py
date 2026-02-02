#!/usr/bin/python3.8

"""
Abstract base class for diagnostic data loaders
"""

from abc import ABC, abstractmethod
from typing import Dict, Any, Optional
import numpy as np
from MDSplus import Connection

from core.data_structures import DiagnosticData

# IP fault detection threshold in kA
IP_THRESHOLD_KA = 100

# Time offset after IP fault to include in valid data (seconds)
IP_FAULT_OFFSET = 0.5


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

    def get_ip_fault_time(self, shot_number: int) -> Optional[float]:
        """Get IP fault time (last time when Ip > threshold)

        Used to mask invalid data after plasma current drops.

        Args:
            shot_number: KSTAR shot number

        Returns:
            IP fault time in seconds, or None if not available
        """
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            ip_time = mds.get('dim_of(\\pcrc03)').data()
            ip_data = mds.get('\\pcrc03/-1e3').data()  # Convert to kA
            mds.closeTree('kstar', shot_number)

            # Find indices where Ip > threshold
            valid_indices = np.where(ip_data > IP_THRESHOLD_KA)[0]
            if len(valid_indices) > 0:
                return ip_time[valid_indices[-1]]  # Last valid time
            return None
        except Exception as e:
            print(f"Warning: Could not get IP fault time: {str(e)}")
            return None

    def get_valid_time_mask(
        self,
        time_array: np.ndarray,
        ip_fault_time: Optional[float]
    ) -> np.ndarray:
        """Get boolean mask for valid time points

        Returns mask where time > 0 and time <= ip_fault_time + IP_FAULT_OFFSET.
        This provides consistent IP fault masking across all diagnostic loaders.

        Args:
            time_array: Array of time values in seconds
            ip_fault_time: IP fault time from get_ip_fault_time(), or None

        Returns:
            Boolean numpy array (True = valid time point)
        """
        if ip_fault_time is None:
            return time_array > 0
        return (time_array > 0) & (time_array <= ip_fault_time + IP_FAULT_OFFSET)
