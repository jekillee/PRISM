#!/usr/bin/python3.8

"""
Data structures for diagnostic measurements
"""


class DiagnosticData:
    """Generic diagnostic data structure"""
    def __init__(self, time, radius, measurements, source='mdsplus', analysis_type=None):
        """
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
    
    def get_parameter(self, param_name):
        """Get data and error for a specific parameter"""
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


class EFITData:
    """EFIT equilibrium data structure for 1D profile mapping"""
    def __init__(self, time, radius, psi_n, rho_pol, rho_tor):
        self.time = time
        self.radius = radius
        self.psi_n = psi_n
        self.rho_pol = rho_pol
        self.rho_tor = rho_tor


class EFITData2D:
    """EFIT 2D equilibrium data structure for poloidal cross-section visualization"""
    def __init__(self, time, r_grid, z_grid, psirz, simag, sibry,
                 bdry_r, bdry_z, nbdry, limiter_r=None, limiter_z=None,
                 rmaxis=None, zmaxis=None):
        self.time = time          # 1D array [nt]
        self.r_grid = r_grid      # 1D array [nr]
        self.z_grid = z_grid      # 1D array [nz]
        self.psirz = psirz        # 3D array [nt, nz, nr]
        self.simag = simag        # 1D array [nt] - psi at magnetic axis
        self.sibry = sibry        # 1D array [nt] - psi at boundary
        self.bdry_r = bdry_r      # 2D array [nt, max_bdry]
        self.bdry_z = bdry_z      # 2D array [nt, max_bdry]
        self.nbdry = nbdry        # 1D array [nt] - number of boundary points
        self.limiter_r = limiter_r  # 1D array - limiter R coordinates
        self.limiter_z = limiter_z  # 1D array - limiter Z coordinates
        self.rmaxis = rmaxis      # 1D array [nt] - magnetic axis R
        self.zmaxis = zmaxis      # 1D array [nt] - magnetic axis Z
    
    def get_psi_normalized(self, time_idx):
        """Get normalized psi at given time index"""
        psirz = self.psirz[time_idx]
        simag = self.simag[time_idx]
        sibry = self.sibry[time_idx]
        return (psirz - simag) / (sibry - simag)
    
    def get_boundary(self, time_idx):
        """Get boundary R, Z coordinates at given time index"""
        nbdry = int(self.nbdry[time_idx])
        return self.bdry_r[time_idx, :nbdry], self.bdry_z[time_idx, :nbdry]
    
    def get_magnetic_axis(self, time_idx):
        """Get magnetic axis R, Z coordinates at given time index"""
        return self.rmaxis[time_idx], self.zmaxis[time_idx]
    
    def find_time_index(self, time):
        """Find closest time index for given time value"""
        import numpy as np
        return np.argmin(np.abs(self.time - time))