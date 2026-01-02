#!/usr/bin/python3.8

"""
XICS (X-ray Imaging Crystal Spectrometer) data loader
Provides Ti and vT at R=1.8m with offset correction
"""

import numpy as np
from MDSplus import Connection
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class XICSLoader(BaseDiagnosticLoader):
    """Loader for XICS diagnostic data"""
    
    # XICS measurement position
    R_POSITION = 1.8  # [m]
    
    # MDS+ node names
    VT_NODE = "\\TXCS_VR053"
    TI_NODE = "\\TXCS_TI053"
    
    def __init__(self, config, diagnostic_config):
        super().__init__(config, diagnostic_config)
        self.ces_loader = None
    
    def _get_ces_loader(self):
        """Lazy initialization of CES loader for offset calculation"""
        if self.ces_loader is None:
            from data_loaders.ces_loader import CESLoader
            from config.diagnostic_config import DIAGNOSTICS
            self.ces_loader = CESLoader(self.app_config, DIAGNOSTICS['CES'])
        return self.ces_loader
    
    def _load_ces_vt_ch01(self, shot_number):
        """Load CES vT CH01 data with fallback: mod -> nn -> raw"""
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        # Try mod first, then nn
        for prefix, analysis_type in [('\\CES_', 'mod'), ('\\CESNN_', 'nn')]:
            try:
                time_node = f"dim_of({prefix}TI01)"
                vt_node = f"{prefix}VT01"
                
                time_arr = mds.get(time_node).data()
                vt_data = mds.get(vt_node).data()
                
                # Check if data is valid
                if len(time_arr) > 0 and not np.all(np.isnan(vt_data)):
                    mds.closeTree('kstar', shot_number)
                    return time_arr, vt_data, analysis_type
            except:
                continue
        
        mds.closeTree('kstar', shot_number)
        return None, None, None
    
    def _calculate_vt_offset(self, xics_time, xics_vt, ces_time, ces_vt):
        """Calculate vT offset using first overlapping time point"""
        if ces_time is None or ces_vt is None:
            return 0.0
        
        # Find overlapping time range
        overlap_start = max(xics_time[0], ces_time[0])
        overlap_end = min(xics_time[-1], ces_time[-1])
        
        if overlap_start >= overlap_end:
            print(f"XICS: No time overlap with CES (XICS: {xics_time[0]:.2f}-{xics_time[-1]:.2f}s, "
                  f"CES: {ces_time[0]:.2f}-{ces_time[-1]:.2f}s)")
            return 0.0
        
        # Find first XICS time point within overlap
        overlap_mask = (xics_time >= overlap_start) & (xics_time <= overlap_end)
        if not np.any(overlap_mask):
            return 0.0
        
        first_overlap_idx = np.where(overlap_mask)[0][0]
        first_overlap_time = xics_time[first_overlap_idx]
        first_xics_vt = xics_vt[first_overlap_idx]
        
        # Interpolate CES vT at this time
        ces_vt_interp = np.interp(first_overlap_time, ces_time, ces_vt)
        
        # Offset = XICS - CES (to be subtracted from XICS)
        offset = first_xics_vt - ces_vt_interp
        
        print(f"XICS: Offset calc at t={first_overlap_time:.3f}s: "
              f"XICS_vT={first_xics_vt:.1f}, CES_vT={ces_vt_interp:.1f}")
        
        return offset
    
    def load_data(self, shot_number, analysis_type=None):
        """Load XICS data from MDS+"""
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            
            # Load XICS time, Ti, vT
            time_arr = mds.get(f"dim_of({self.VT_NODE})").data()
            vt_data = mds.get(self.VT_NODE).data()
            ti_data = mds.get(self.TI_NODE).data()
            
            mds.closeTree('kstar', shot_number)
            
            # Convert to numpy arrays
            time_arr = np.array(time_arr)
            vt_data = np.array(vt_data)
            ti_data = np.array(ti_data)
            
            # Load CES vT for offset calculation
            ces_time, ces_vt, ces_type = self._load_ces_vt_ch01(shot_number)
            
            # Calculate and apply offset
            offset = self._calculate_vt_offset(time_arr, vt_data, ces_time, ces_vt)
            vt_data_corrected = vt_data - offset
            
            if ces_type:
                print(f"XICS vT offset: {offset:.3f} km/s (based on CES {ces_type})")
            else:
                print(f"XICS vT offset: not applied (no CES data)")
            
            # XICS Ti is already in keV, vT is already in km/s
            # No unit conversion needed
            
            # Single radius point at R=1.8m
            radius = np.array([self.R_POSITION])
            
            # Reshape data to (1, n_time) for consistency with other loaders
            ti_data = ti_data.reshape(1, -1)
            vt_data_corrected = vt_data_corrected.reshape(1, -1)
            
            # No error data available for XICS
            ti_err = np.zeros_like(ti_data)
            vt_err = np.zeros_like(vt_data_corrected)
            
            # Package measurements
            measurements = {
                'Ti': {'data': ti_data, 'error': ti_err},
                'vT': {'data': vt_data_corrected, 'error': vt_err},
                'vT_offset': offset
            }
            
            return DiagnosticData(time_arr, radius, measurements, source='XICS')
            
        except Exception as e:
            # Return None instead of raising error (silent skip)
            print(f"XICS data not available for shot {shot_number}: {str(e)}")
            return None
    
    def get_data_at_time(self, shot_number, time_point, dt_avg=0.01):
        """Get XICS Ti, vT at specific time point"""
        data = self.load_data(shot_number)
        
        if data is None:
            return None
        
        # Find time index
        time_idx = np.argmin(np.abs(data.time - time_point))
        
        # Check if time is within reasonable range
        if np.abs(data.time[time_idx] - time_point) > dt_avg:
            return None
        
        Ti_data, Ti_err = data.get_parameter('Ti')
        vT_data, vT_err = data.get_parameter('vT')
        
        return {
            'R': self.R_POSITION,
            'Ti': Ti_data[0, time_idx],
            'Ti_err': Ti_err[0, time_idx],
            'vT': vT_data[0, time_idx],
            'vT_err': vT_err[0, time_idx],
            'time_actual': data.time[time_idx],
            'vT_offset': data.measurements.get('vT_offset', 0.0)
        }
