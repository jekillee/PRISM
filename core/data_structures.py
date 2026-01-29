#!/usr/bin/python3.8

"""
Data structures for diagnostic measurements
"""

from typing import Dict, Any, Optional, Tuple
import numpy as np
from numpy.typing import NDArray


class DiagnosticData:
    """Generic diagnostic data structure"""

    def __init__(
        self,
        time: NDArray[np.float64],
        radius: NDArray[np.float64],
        measurements: Dict[str, Dict[str, Any]],
        source: str = 'mdsplus',
        analysis_type: Optional[str] = None
    ) -> None:
        """
        Args:
            time: 1D array of time points
            radius: 1D array of radial positions
            measurements: dict of {param_name: {'data': 2D array, 'error': 2D array, ...}}
            source: data source identifier ('mdsplus' or file name)
            analysis_type: analysis method (e.g., 'mod', 'nn' for CES)
        """
        self.time = time
        self.radius = radius
        self.measurements = measurements
        self.source = source
        self.analysis_type = analysis_type

    def get_parameter(
        self, param_name: str
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Get data and error for a specific parameter

        Args:
            param_name: Name of the parameter (e.g., 'Ti', 'Te', 'ne')

        Returns:
            Tuple of (data, error) arrays

        Raises:
            KeyError: If parameter not found in measurements
        """
        if param_name not in self.measurements:
            raise KeyError(f"Parameter {param_name} not found in measurements")

        param_dict = self.measurements[param_name]
        data = param_dict['data']

        # Check if asymmetric errors exist (Thomson)
        if 'error_upper' in param_dict and 'error_lower' in param_dict:
            error_upper = param_dict['error_upper']
            error_lower = param_dict['error_lower']
            error = (error_upper + error_lower) / 2
        else:
            error = param_dict['error']

        return data, error

    def get_parameter_asymmetric(
        self, param_name: str
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
        """Get data with asymmetric errors (upper and lower separately)

        Args:
            param_name: Name of the parameter

        Returns:
            Tuple of (data, error_upper, error_lower) arrays

        Raises:
            KeyError: If parameter not found in measurements
        """
        if param_name not in self.measurements:
            raise KeyError(f"Parameter {param_name} not found in measurements")

        param_dict = self.measurements[param_name]
        data = param_dict['data']

        if 'error_upper' in param_dict and 'error_lower' in param_dict:
            return data, param_dict['error_upper'], param_dict['error_lower']

        error = param_dict['error']
        return data, error, error


class EFITData:
    """EFIT equilibrium data structure for 1D profile mapping"""

    def __init__(
        self,
        time: NDArray[np.float64],
        radius: NDArray[np.float64],
        psi_n: NDArray[np.float64],
        rho_pol: NDArray[np.float64],
        rho_tor: NDArray[np.float64]
    ) -> None:
        self.time = time
        self.radius = radius
        self.psi_n = psi_n
        self.rho_pol = rho_pol
        self.rho_tor = rho_tor


class EFITData2D:
    """EFIT 2D equilibrium data structure for poloidal cross-section visualization"""

    def __init__(
        self,
        time: NDArray[np.float64],
        r_grid: NDArray[np.float64],
        z_grid: NDArray[np.float64],
        psirz: NDArray[np.float64],
        simag: NDArray[np.float64],
        sibry: NDArray[np.float64],
        bdry_r: NDArray[np.float64],
        bdry_z: NDArray[np.float64],
        nbdry: NDArray[np.int_],
        limiter_r: Optional[NDArray[np.float64]] = None,
        limiter_z: Optional[NDArray[np.float64]] = None,
        rmaxis: Optional[NDArray[np.float64]] = None,
        zmaxis: Optional[NDArray[np.float64]] = None
    ) -> None:
        """
        Args:
            time: 1D array [nt]
            r_grid: 1D array [nr]
            z_grid: 1D array [nz]
            psirz: 3D array [nt, nz, nr]
            simag: 1D array [nt] - psi at magnetic axis
            sibry: 1D array [nt] - psi at boundary
            bdry_r: 2D array [nt, max_bdry]
            bdry_z: 2D array [nt, max_bdry]
            nbdry: 1D array [nt] - number of boundary points
            limiter_r: 1D array - limiter R coordinates
            limiter_z: 1D array - limiter Z coordinates
            rmaxis: 1D array [nt] - magnetic axis R
            zmaxis: 1D array [nt] - magnetic axis Z
        """
        self.time = time
        self.r_grid = r_grid
        self.z_grid = z_grid
        self.psirz = psirz
        self.simag = simag
        self.sibry = sibry
        self.bdry_r = bdry_r
        self.bdry_z = bdry_z
        self.nbdry = nbdry
        self.limiter_r = limiter_r
        self.limiter_z = limiter_z
        self.rmaxis = rmaxis
        self.zmaxis = zmaxis

    def get_psi_normalized(self, time_idx: int) -> NDArray[np.float64]:
        """Get normalized psi at given time index

        Args:
            time_idx: Time index

        Returns:
            2D array of normalized psi values
        """
        psirz = self.psirz[time_idx]
        simag = self.simag[time_idx]
        sibry = self.sibry[time_idx]
        return (psirz - simag) / (sibry - simag)

    def get_boundary(
        self, time_idx: int
    ) -> Tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Get boundary R, Z coordinates at given time index

        Args:
            time_idx: Time index

        Returns:
            Tuple of (R, Z) boundary coordinate arrays
        """
        nbdry = int(self.nbdry[time_idx])
        return self.bdry_r[time_idx, :nbdry], self.bdry_z[time_idx, :nbdry]

    def get_magnetic_axis(self, time_idx: int) -> Tuple[float, float]:
        """Get magnetic axis R, Z coordinates at given time index

        Args:
            time_idx: Time index

        Returns:
            Tuple of (R, Z) magnetic axis coordinates
        """
        return self.rmaxis[time_idx], self.zmaxis[time_idx]

    def find_time_index(self, time: float) -> int:
        """Find closest time index for given time value

        Args:
            time: Time value to search for

        Returns:
            Index of closest time point
        """
        return int(np.argmin(np.abs(self.time - time)))
