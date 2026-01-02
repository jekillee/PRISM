#!/usr/bin/python3.8

"""
TCI (Two-Color Interferometer) data loader for line-averaged density
"""

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class TCILoader(BaseDiagnosticLoader):
    """Loader for TCI line-averaged density data"""
    
    N_CHANNELS = 5
    
    # Sampling rate options (same as ECE for consistency)
    SAMPLING_RATES = {
        '100Hz': 0.01,
        '1kHz': 0.001
    }
    
    def load_data(self, shot_number, analysis_type=None, sampling_rate=None):
        """Load TCI data from MDS+
        
        sampling_rate: sampling period in seconds (None=original, 0.01=100Hz, 0.001=1kHz)
        """
        mds = self._connect_mds(shot_number)
        
        try:
            channels = []
            ne_list = []
            time_data = None
            
            for ch in range(1, self.N_CHANNELS + 1):
                node_name = f'\\ne_tci{ch:02d}'
                try:
                    if sampling_rate is not None:
                        # Use resample for downsampled data
                        data = mds.get(f'resample({node_name}, 0, *, {sampling_rate})').data()
                        time = np.linspace(0, len(data) * sampling_rate, len(data))
                    else:
                        # Original data without resampling
                        data = mds.get(node_name).data()
                        time = mds.get(f'dim_of({node_name})').data()
                    
                    channels.append(ch)
                    ne_list.append(data)
                    
                    if time_data is None:
                        time_data = time
                        
                except Exception as e:
                    print(f"TCI channel {ch} not available: {str(e)}")
            
            if not channels:
                raise ValueError("No TCI channels available")
            
            # Create data array (n_channels x n_time)
            n_time = len(time_data)
            ne_data = np.full((len(channels), n_time), np.nan)
            
            for i, ne in enumerate(ne_list):
                if len(ne) == n_time:
                    ne_data[i, :] = ne
                else:
                    # Interpolate if time arrays differ
                    ne_data[i, :] = np.interp(time_data, 
                                              np.linspace(time_data[0], time_data[-1], len(ne)), 
                                              ne)
            
            # Use channel numbers as radius placeholder
            radius = np.array(channels, dtype=float)
            
            measurements = {
                'ne': {
                    'data': ne_data,
                    'error': np.zeros_like(ne_data),  # No error data
                    'channels': channels
                }
            }
            
            return DiagnosticData(
                time=time_data,
                radius=radius,
                measurements=measurements,
                source='mdsplus'
            )
            
        finally:
            self._close_mds(mds, shot_number)
    
    def load_raw_data(self, shot_number, channel):
        """Load TCI raw signal data from MDS+
        
        channel: channel number (1-5)
        Returns: dict with 'time' and 'data' arrays
        """
        mds = self._connect_mds(shot_number)
        
        try:
            node_name = f'\\tci{channel:02d}_ls:foo'
            
            data = mds.get(node_name).data()
            time = mds.get(f'dim_of({node_name})').data()
            
            # Raw signal has one extra time frame at the beginning
            time = time[1:]
            
            return {
                'time': np.array(time),
                'data': np.array(data)
            }
            
        finally:
            self._close_mds(mds, shot_number)