#!/usr/bin/python3.8

"""
CES data loader implementation with analysis type support
"""

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class CESLoader(BaseDiagnosticLoader):
    """Loader for CES diagnostic data"""
    
    def load_data(self, shot_number, analysis_type='mod'):
        """Load CES data from MDS+
        
        analysis_type: 'mod' for beam modulation or 'nn' for neural network
        """
        try:
            mds = self._connect_mds(shot_number)
            
            # Get node prefix based on analysis type
            analysis_config = self.diag_config['analysis_types'].get(analysis_type, 
                             self.diag_config['analysis_types']['mod'])
            node_prefix = analysis_config['node_prefix']
            
            # Load time
            time_node = self.diag_config['time']['mds_node_format'] % node_prefix
            time_arr = mds.get(time_node).data()
            
            # Load radius, temperature, and velocity for all channels
            radius_list, ti_list, ti_err_list, vt_list, vt_err_list = [], [], [], [], []
            
            n_channels = self.diag_config['channels']
            
            for i in range(1, n_channels + 1):
                # Radius
                radius_node = self.diag_config['radius']['mds_node_format'] % (node_prefix, i)
                radius_list.append(mds.get(radius_node).data()[0])
                
                # Temperature
                ti_node = self.diag_config['parameters'][0]['mds_node_format'] % (node_prefix, i)
                ti_err_node = self.diag_config['parameters'][0]['mds_error_format'] % (node_prefix, i)
                ti_list.append(mds.get(ti_node).data())
                ti_err_list.append(mds.get(ti_err_node).data())
                
                # Velocity
                vt_node = self.diag_config['parameters'][1]['mds_node_format'] % (node_prefix, i)
                vt_err_node = self.diag_config['parameters'][1]['mds_error_format'] % (node_prefix, i)
                vt_list.append(mds.get(vt_node).data())
                vt_err_list.append(mds.get(vt_err_node).data())
            
            self._close_mds(mds, shot_number)
            
            # Convert to numpy arrays and apply scaling
            np.seterr(invalid='ignore')
            
            radius = np.array(radius_list) * self.diag_config['radius']['scale']
            ti_data = np.array(ti_list) * self.diag_config['parameters'][0]['scale']
            ti_err = np.array(ti_err_list) * self.diag_config['parameters'][0]['scale']
            vt_data = np.array(vt_list) * self.diag_config['parameters'][1]['scale']
            vt_err = np.array(vt_err_list) * self.diag_config['parameters'][1]['scale']
            
            # Clean negative error values
            ti_err[ti_err < 0] = np.nan
            vt_err[vt_err < 0] = np.nan
            
            # Package measurements
            measurements = {
                'Ti': {'data': ti_data, 'error': ti_err},
                'vT': {'data': vt_data, 'error': vt_err}
            }
            
            return DiagnosticData(time_arr, radius, measurements, 
                                 source='mdsplus', analysis_type=analysis_type)
            
        except Exception as e:
            raise RuntimeError(f"Failed to load CES data for shot {shot_number}: {str(e)}")