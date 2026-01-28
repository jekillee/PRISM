#!/usr/bin/python3.8

"""
ECE (Electron Cyclotron Emission) data loader
"""

import numpy as np
from MDSplus import Connection
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class ECELoader(BaseDiagnosticLoader):
    """Loader for ECE diagnostic data"""
    
    # ECE R range (overlapping with Thomson: 1.8-2.3m)
    R_MIN = 1.8
    R_MAX = 2.3
    
    # Spectrogram R range (wider for spectrogram use)
    R_MIN_SPEC = 1.3
    R_MAX_SPEC = 2.3
    
    # Known bad ECE channels (always excluded)
    BAD_CHANNELS = [1, 2, 3, 27, 28, 48, 49, 52, 59, 72]
    
    # Sampling rate options
    SAMPLING_RATES = {
        '20Hz': 0.05,
        '50Hz': 0.02,
        '100Hz': 0.01,
        '500Hz': 0.002,
        '1000Hz': 0.001
    }
    DEFAULT_SAMPLING_RATE = 0.01  # 100Hz
    
    def __init__(self, app_config, diag_config=None):
        super().__init__(app_config, diag_config)
        self._position_cache = {}  # {shot: {'channels': [], 'R': [], ...}}
    
    def get_channel_positions(self, shot_number, r_min=None, r_max=None):
        """Get R, Z positions for ECE channels
        
        Returns dict with:
            channels: array of channel numbers
            R: array in [m]
            Z: array in [m] (always 0 for ECE)
            I_TF: TF coil current [kA]
        """
        cache_key = shot_number
        if cache_key in self._position_cache:
            return self._position_cache[cache_key]
        
        # Use spectrogram range if not specified
        if r_min is None:
            r_min = self.R_MIN_SPEC
        if r_max is None:
            r_max = self.R_MAX_SPEC
        
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        try:
            # Get TF coil current
            I_TF = np.mean(mds.get('\\PCITFMSRD').data()) / 1e3  # [kA]
            
            # Get ECE channel frequencies
            freq = []
            for i in range(76):
                freq.append(int(mds.get(f'\\ECE{i+1:02d}:FREQ').data()))
            freq = np.array(freq)
            
            mds.closeTree('kstar', shot_number)
            
            # Calculate R positions (2nd harmonic)
            R_2nd = 10.0415 * I_TF / freq
            
            # Filter channels within R range
            ch_all = np.arange(1, 77)
            mask_in_range = (R_2nd >= r_min) & (R_2nd <= r_max)
            
            # Exclude known bad channels
            mask_good = ~np.isin(ch_all, self.BAD_CHANNELS)
            mask_valid = mask_in_range & mask_good
            
            channels = ch_all[mask_valid]
            R_values = R_2nd[mask_valid]
            Z_values = np.zeros_like(R_values)  # ECE Z position is always 0
            
            # Sort by R
            sort_idx = np.argsort(R_values)
            channels = channels[sort_idx]
            R_values = R_values[sort_idx]
            Z_values = Z_values[sort_idx]
            
            result = {
                'channels': channels,
                'R': R_values,
                'Z': Z_values,
                'I_TF': I_TF
            }
            
            self._position_cache[cache_key] = result
            return result
            
        except Exception as e:
            mds.closeTree('kstar', shot_number)
            raise e
    
    def get_ip_fault_time(self, shot_number):
        """Get IP fault time (last time when Ip > 100 kA)"""
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            ip_time = mds.get('dim_of(\\pcrc03)').data()
            ip_data = mds.get('\\pcrc03/-1e3').data()  # Convert to kA
            mds.closeTree('kstar', shot_number)
            
            # Find indices where Ip > 100 kA
            valid_indices = np.where(ip_data > 100)[0]
            if len(valid_indices) > 0:
                return ip_time[valid_indices[-1]]  # Last valid time
            return None
        except Exception as e:
            print(f"Warning: Could not get IP fault time: {str(e)}")
            return None
    
    def load_data(self, shot_number, analysis_type=None, sampling_rate=None):
        """Load ECE data from MDS+
        
        sampling_rate: sampling period in seconds (default 0.01 = 100Hz)
        """
        if sampling_rate is None:
            sampling_rate = self.DEFAULT_SAMPLING_RATE
            
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            
            # Get TF coil current
            I_TF = np.mean(mds.get('\\PCITFMSRD').data()) / 1e3  # [kA]
            
            # Get ECE channel frequencies
            ch_all = np.arange(1, 77)
            freq = []
            for i in range(76):
                freq.append(int(mds.get(f'\\ECE{i+1:02d}:FREQ').data()))
            freq = np.array(freq)
            
            # Calculate resonance positions (2nd and 3rd harmonic)
            R_2nd_all = 10.0415 * I_TF / freq
            R_3rd_all = 15.0623 * I_TF / freq
            
            # Create masks for channels within R range
            mask_2nd_inside = (R_2nd_all >= self.R_MIN) & (R_2nd_all <= self.R_MAX)
            mask_3rd_inside = (R_3rd_all >= self.R_MIN) & (R_3rd_all <= self.R_MAX)
            mask_overlap = mask_2nd_inside & mask_3rd_inside
            mask_valid = mask_2nd_inside & (~mask_overlap)
            
            # Load all channels where 2nd harmonic is inside boundary
            mask_load = mask_2nd_inside
            ch = ch_all[mask_load]
            R = R_2nd_all[mask_load]
            
            # Store masks for loaded channels
            mask_overlap_load = mask_overlap[mask_load]
            mask_valid_load = mask_valid[mask_load]
            
            if len(ch) == 0:
                mds.closeTree('kstar', shot_number)
                raise RuntimeError("No ECE channels in valid R range")
            
            # Get time base with specified sampling rate (start from t=0)
            first_ch = ch[0]
            mds.get(f'SetTimeContext(0,*,{sampling_rate})').data()
            time = mds.get(f'dim_of(\\ECE{first_ch:02d})').data()
            
            # Load ECE data for all channels in mask_load
            Te_list = []
            loaded_indices = []
            
            print(f"  Loading ECE data for {len(ch)} channels...")
            import time as time_module
            start_time = time_module.time()
            
            for idx, channel in enumerate(ch):
                try:
                    data = mds.get(f'\\ECE{channel:02d}').data()
                    Te_list.append(data)
                    loaded_indices.append(idx)
                except:
                    print(f'    ECE{channel:02d} not available')
                
                # Progress bar
                progress = (idx + 1) / len(ch)
                bar_length = 40
                filled = int(bar_length * progress)
                bar = '█' * filled + '░' * (bar_length - filled)
                print(f'\r    [{bar}] {idx+1}/{len(ch)} ({progress*100:.1f}%)', end='', flush=True)
            
            elapsed = time_module.time() - start_time
            print(f'\n  Completed in {elapsed:.2f} sec')
            
            mds.get('SetTimeContext(,,)').data()
            
            if len(Te_list) == 0:
                mds.closeTree('kstar', shot_number)
                raise RuntimeError("No valid ECE data found")
            
            # Convert to arrays (only successfully loaded channels)
            loaded_indices = np.array(loaded_indices)
            Te_data = np.array(Te_list)
            R_arr = R[loaded_indices]
            ch_arr = ch[loaded_indices]
            overlap_arr = mask_overlap_load[loaded_indices]
            valid_arr = mask_valid_load[loaded_indices]
            
            # Load baseline for correction (-1 to 0 sec)
            mds.get(f'SetTimeContext(-1,0,{sampling_rate})').data()
            baseline_list = []
            
            print(f"  Loading baseline data (-1 to 0 s)...")
            for channel in ch_arr:
                try:
                    baseline_data = mds.get(f'\\ECE{channel:02d}').data()
                    baseline_avg = np.mean(baseline_data)
                    baseline_list.append(baseline_avg)
                except:
                    baseline_list.append(0.0)
            
            # Subtract baseline from main data
            baseline_arr = np.array(baseline_list).reshape(-1, 1)
            Te_data = Te_data - baseline_arr
            print(f"  Baseline correction applied")
            
            mds.get('SetTimeContext(,,)').data()
            mds.closeTree('kstar', shot_number)
            
            # Filter out bad channels (Te_mean == 0)
            Te_mean = np.mean(Te_data, axis=1)
            mask_good = (Te_mean != 0)
            
            if np.sum(~mask_good) > 0:
                bad_channels = ch_arr[~mask_good]
                print(f"  Removing bad channels: {bad_channels}")
            
            Te_data = Te_data[mask_good]
            R_arr = R_arr[mask_good]
            ch_arr = ch_arr[mask_good]
            overlap_arr = overlap_arr[mask_good]
            valid_arr = valid_arr[mask_good]
            
            # Sort by R
            sort_idx = np.argsort(R_arr)
            R_arr = R_arr[sort_idx]
            Te_data = Te_data[sort_idx]
            ch_arr = ch_arr[sort_idx]
            overlap_arr = overlap_arr[sort_idx]
            valid_arr = valid_arr[sort_idx]
            
            # Get ip_fault_time and mask data
            ip_fault_time = self.get_ip_fault_time(shot_number)
            if ip_fault_time is None:
                raise RuntimeError(f"Shot {shot_number}: Ip did not exceed 100 kA (failed shot)")
            
            # Mask time range: 0 < time <= ip_fault_time
            valid_time_mask = (time > 0) & (time <= ip_fault_time)
            time = time[valid_time_mask]
            Te_data = Te_data[:, valid_time_mask]
            
            print(f"  Data masked to IP fault time: {ip_fault_time:.3f} s")
            
            # Apply unit conversion (eV to keV)
            np.seterr(invalid='ignore')
            Te_data = Te_data * self.app_config.TE_SCALE
            
            # No error bars for ECE
            Te_err = np.zeros_like(Te_data)
            
            # Package measurements
            measurements = {
                'Te': {
                    'data': Te_data, 
                    'error': Te_err,
                    'overlap_mask': overlap_arr,
                    'valid_mask': valid_arr,
                    'channels': ch_arr
                }
            }
            
            # Store additional info
            result = DiagnosticData(time, R_arr, measurements, source='mdsplus')
            result.I_TF = I_TF
            
            return result
            
        except Exception as e:
            raise RuntimeError(f"Failed to load ECE data for shot {shot_number}: {str(e)}")
    
    def get_profile_at_time(self, shot_number, time_point, dt_avg=0.005, sampling_rate=None):
        """Get ECE Te profile at specific time
        
        time_point: target time [s]
        dt_avg: averaging window half-width [s] (default Â±5ms)
        sampling_rate: sampling period in seconds (default 0.01 = 100Hz)
        """
        try:
            # Load full data
            data = self.load_data(shot_number, sampling_rate=sampling_rate)
            
            # Find time indices for averaging
            idx_avg = np.where((data.time >= (time_point - dt_avg)) & 
                              (data.time <= (time_point + dt_avg)))[0]
            
            if len(idx_avg) == 0:
                return None
            
            time_actual = np.mean(data.time[idx_avg])
            
            # Get Te profile averaged over time window
            Te_data = data.measurements['Te']['data']
            Te_profile = np.mean(Te_data[:, idx_avg], axis=1)
            
            return {
                'R': data.radius,
                'Te': Te_profile,
                'channels': data.measurements['Te']['channels'],
                'overlap_mask': data.measurements['Te']['overlap_mask'],
                'valid_mask': data.measurements['Te']['valid_mask'],
                'time_actual': time_actual,
                'I_TF': data.I_TF
            }
            
        except Exception as e:
            print(f"Error getting ECE profile: {str(e)}")
            return None