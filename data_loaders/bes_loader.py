#!/usr/bin/python3.8

"""
BES (Beam Emission Spectroscopy) data loader
"""

import numpy as np
from MDSplus import Connection
from data_loaders.base_loader import BaseDiagnosticLoader


class BESLoader(BaseDiagnosticLoader):
    """Loader for BES diagnostic data"""
    
    # Channel structure: 4 vertical x 16 radial
    VERTICAL_CHANNELS = 4
    RADIAL_CHANNELS = 16
    
    def __init__(self, config, diagnostic_config):
        super().__init__(config, diagnostic_config)
        self._position_cache = {}  # {shot: {'channels': [], 'R': [], 'Z': []}}
    
    def get_channel_positions(self, shot_number):
        """Get R, Z positions for all BES channels"""
        if shot_number in self._position_cache:
            return self._position_cache[shot_number]
        
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        try:
            channels = []
            R_values = []
            Z_values = []
            
            for v in range(1, self.VERTICAL_CHANNELS + 1):
                for r in range(1, self.RADIAL_CHANNELS + 1):
                    ch_name = f'{v:02d}{r:02d}'
                    try:
                        R_pos = mds.get(f'\\BES_{ch_name}:RPOS').data()
                        V_pos = mds.get(f'\\BES_{ch_name}:VPOS').data()
                        
                        # Handle array or scalar
                        if hasattr(R_pos, '__len__'):
                            R_pos = float(R_pos[0])
                        else:
                            R_pos = float(R_pos)
                        if hasattr(V_pos, '__len__'):
                            V_pos = float(V_pos[0])
                        else:
                            V_pos = float(V_pos)
                        
                        # Convert mm to m
                        R_pos = R_pos * 1e-3
                        V_pos = V_pos * 1e-3
                        
                        channels.append(ch_name)
                        R_values.append(R_pos)
                        Z_values.append(V_pos)
                        
                    except:
                        pass
            
            mds.closeTree('kstar', shot_number)
            
            result = {
                'channels': channels,
                'R': np.array(R_values),
                'Z': np.array(Z_values)
            }
            
            self._position_cache[shot_number] = result
            return result
            
        except Exception as e:
            mds.closeTree('kstar', shot_number)
            raise e
    
    def load_channel_data(self, shot_number, channel_name):
        """Load single BES channel data
        
        channel_name: '0108' format (vertical + radial)
        """
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        try:
            node_name = f'\\BES_{channel_name}:FOO'
            
            data = mds.get(node_name).data()
            time = mds.get(f'dim_of({node_name})').data()
            
            mds.closeTree('kstar', shot_number)
            
            # Get position info
            positions = self.get_channel_positions(shot_number)
            R_pos = None
            Z_pos = None
            
            if channel_name in positions['channels']:
                idx = positions['channels'].index(channel_name)
                R_pos = positions['R'][idx]
                Z_pos = positions['Z'][idx]
            
            return {
                'time': np.array(time),
                'data': np.array(data),
                'channel': channel_name,
                'R': R_pos,
                'Z': Z_pos
            }
            
        except Exception as e:
            mds.closeTree('kstar', shot_number)
            raise RuntimeError(f"Failed to load BES {channel_name}: {str(e)}")
    
    def load_data(self, shot_number, analysis_type=None):
        """Load BES channel info (interface compatibility)"""
        return self.get_channel_positions(shot_number)