"""
Thomson Scattering data loader implementation
"""

import json
from pathlib import Path
import numpy as np
from MDSplus import MdsException
from data_loaders.base_loader import BaseDiagnosticLoader, IP_FAULT_OFFSET
from core.data_structures import DiagnosticData

# Load Thomson positions configuration
_CONFIG_PATH = Path(__file__).parent.parent / 'config' / 'thomson_positions.json'
with open(_CONFIG_PATH, 'r') as f:
    _THOMSON_CONFIG = json.load(f)


class ThomsonLoader(BaseDiagnosticLoader):
    """Loader for Thomson Scattering diagnostic data"""

    def load_data(self, shot_number, analysis_type=None):
        """Load Thomson scattering data from MDS+
        
        analysis_type: Not used for Thomson (kept for interface consistency)
        """
        try:
            mds = self._connect_mds(shot_number)
            
            nbchanc = 14  # core channels
            nbchane = 17  # edge channels
            
            core_positions, edge_positions = self._get_hardcoded_positions(shot_number)
            
            te_data, ne_data = [], []
            te_errh_data, te_errl_data = [], []
            ne_errh_data, ne_errl_data = [], []
            ts_position = []
            times = None
            
            # Load core channels
            for i in range(nbchanc):
                try:
                    pos = core_positions[i]
                    if pos is None:
                        continue
                    te_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_TE').data()
                    ne_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_NE').data()
                    te_errh_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_TERRH').data()
                    te_errl_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_TERRL').data()
                    ne_errh_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_NERRH').data()
                    ne_errl_raw = mds.get(f'\\TS_CORE{i+1}:CORE{i+1}_NERRL').data()
                    
                    if times is None:
                        times = mds.get(f'dim_of(\\TS_CORE{i+1}:CORE{i+1}_TE)').data()
                        times = self._squeeze_if_needed(times)
                    
                    te_raw = self._squeeze_if_needed(te_raw)
                    ne_raw = self._squeeze_if_needed(ne_raw)
                    te_errh_raw = self._squeeze_if_needed(te_errh_raw)
                    te_errl_raw = self._squeeze_if_needed(te_errl_raw)
                    ne_errh_raw = self._squeeze_if_needed(ne_errh_raw)
                    ne_errl_raw = self._squeeze_if_needed(ne_errl_raw)
                    
                    te_raw[np.isnan(te_raw)] = 0.
                    ne_raw[np.isnan(ne_raw)] = 0.
                    te_errh_raw[np.isnan(te_errh_raw)] = 0.
                    te_errl_raw[np.isnan(te_errl_raw)] = 0.
                    ne_errh_raw[np.isnan(ne_errh_raw)] = 0.
                    ne_errl_raw[np.isnan(ne_errl_raw)] = 0.
                    
                    if np.max(te_raw) > 0. and np.max(ne_raw) > 0.:
                        te_data.append(te_raw)
                        ne_data.append(ne_raw)
                        te_errh_data.append(te_errh_raw)
                        te_errl_data.append(te_errl_raw)
                        ne_errh_data.append(ne_errh_raw)
                        ne_errl_data.append(ne_errl_raw)
                        ts_position.append(pos / 1000.)
                        
                except MdsException:
                    print(f'TS_CORE{i+1} not available')
            
            # Load edge channels
            for i in range(nbchane):
                try:
                    pos = edge_positions[i]
                    if pos is None:
                        continue
                    te_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_TE').data()
                    ne_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_NE').data()
                    te_errh_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_TERRH').data()
                    te_errl_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_TERRL').data()
                    ne_errh_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_NERRH').data()
                    ne_errl_raw = mds.get(f'\\TS_EDGE{i+1}:EDGE{i+1}_NERRL').data()
                    
                    if times is None:
                        times = mds.get(f'dim_of(\\TS_EDGE{i+1}:EDGE{i+1}_TE)').data()
                        times = self._squeeze_if_needed(times)
                    
                    te_raw = self._squeeze_if_needed(te_raw)
                    ne_raw = self._squeeze_if_needed(ne_raw)
                    te_errh_raw = self._squeeze_if_needed(te_errh_raw)
                    te_errl_raw = self._squeeze_if_needed(te_errl_raw)
                    ne_errh_raw = self._squeeze_if_needed(ne_errh_raw)
                    ne_errl_raw = self._squeeze_if_needed(ne_errl_raw)
                    
                    te_raw[np.isnan(te_raw)] = 0.
                    ne_raw[np.isnan(ne_raw)] = 0.
                    te_errh_raw[np.isnan(te_errh_raw)] = 0.
                    te_errl_raw[np.isnan(te_errl_raw)] = 0.
                    ne_errh_raw[np.isnan(ne_errh_raw)] = 0.
                    ne_errl_raw[np.isnan(ne_errl_raw)] = 0.
                    
                    if np.max(te_raw) > 0. and np.max(ne_raw) > 0.:
                        te_data.append(te_raw)
                        ne_data.append(ne_raw)
                        te_errh_data.append(te_errh_raw)
                        te_errl_data.append(te_errl_raw)
                        ne_errh_data.append(ne_errh_raw)
                        ne_errl_data.append(ne_errl_raw)
                        ts_position.append(pos / 1000.)
                        
                except MdsException:
                    print(f'TS_EDGE{i+1} not available')
            
            self._close_mds(mds, shot_number)
            
            if not te_data or times is None:
                raise RuntimeError("No valid Thomson data found")
            
            # Sort by radial position
            sorted_indices = np.argsort(ts_position)
            ts_position = np.array(ts_position)[sorted_indices]
            te_data = np.array(te_data)[sorted_indices]
            ne_data = np.array(ne_data)[sorted_indices]
            te_errh_data = np.array(te_errh_data)[sorted_indices]
            te_errl_data = np.array(te_errl_data)[sorted_indices]
            ne_errh_data = np.array(ne_errh_data)[sorted_indices]
            ne_errl_data = np.array(ne_errl_data)[sorted_indices]
            
            # Unit conversion
            np.seterr(invalid='ignore')
            temperature = te_data * self.app_config.TE_SCALE
            density = ne_data * self.app_config.NE_SCALE
            temp_error_upper = te_errh_data * self.app_config.TE_SCALE
            temp_error_lower = te_errl_data * self.app_config.TE_SCALE
            density_error_upper = ne_errh_data * self.app_config.NE_SCALE
            density_error_lower = ne_errl_data * self.app_config.NE_SCALE
            
            # Set minimum error values
            temp_error_upper[temp_error_upper == 0] = 0.01
            temp_error_lower[temp_error_lower == 0] = 0.01
            density_error_upper[density_error_upper == 0] = 0.01
            density_error_lower[density_error_lower == 0] = 0.01
            
            # Get ip_fault_time and mask data
            ip_fault_time = self.get_ip_fault_time(shot_number)
            if ip_fault_time is None:
                raise RuntimeError(f"Shot {shot_number}: Ip did not exceed 100 kA (failed shot)")

            # Mask time range using centralized method
            valid_time_mask = self.get_valid_time_mask(times, ip_fault_time)
            times = times[valid_time_mask]
            temperature = temperature[:, valid_time_mask]
            density = density[:, valid_time_mask]
            temp_error_upper = temp_error_upper[:, valid_time_mask]
            temp_error_lower = temp_error_lower[:, valid_time_mask]
            density_error_upper = density_error_upper[:, valid_time_mask]
            density_error_lower = density_error_lower[:, valid_time_mask]

            print(f"[Thomson]   Data masked to IP fault time + {IP_FAULT_OFFSET}s: {ip_fault_time + IP_FAULT_OFFSET:.3f} s")
            
            measurements = {
                'Te': {
                    'data': temperature, 
                    'error': (temp_error_upper + temp_error_lower) / 2,
                    'error_upper': temp_error_upper,
                    'error_lower': temp_error_lower
                },
                'ne': {
                    'data': density, 
                    'error': (density_error_upper + density_error_lower) / 2,
                    'error_upper': density_error_upper,
                    'error_lower': density_error_lower
                }
            }
            
            return DiagnosticData(times, ts_position, measurements, 
                                 source='mdsplus', analysis_type=None)
            
        except Exception as e:
            raise RuntimeError(f"Failed to load Thomson data for shot {shot_number}: {str(e)}")
    
    def _squeeze_if_needed(self, data):
        """Remove unnecessary outer dimensions"""
        data = np.array(data)
        while data.ndim > 1 and data.shape[0] == 1:
            data = np.squeeze(data, axis=0)
        return data
    
    def _get_hardcoded_positions(self, shot_num):
        """Get channel positions from config file based on shot number"""
        for entry in _THOMSON_CONFIG['positions']:
            shot_min, shot_max = entry['shot_range']
            if shot_max is None:
                shot_max = float('inf')
            if shot_min <= shot_num <= shot_max:
                return entry['core'], entry['edge']

        # Fallback to last entry if no match (shouldn't happen)
        last_entry = _THOMSON_CONFIG['positions'][-1]
        return last_entry['core'], last_entry['edge']
