#!/usr/bin/python3.8

"""
ECEI (Electron Cyclotron Emission Imaging) data loader
Supports GT, GR, HT devices with ABCD matrix position calculation
"""

import numpy as np
from MDSplus import Connection


class ECEILoader:
    """ECEI data loader with channel position calculation"""
    
    # Device configurations
    DEVICES = ['GT', 'GR', 'HT']
    
    # Channel structure: 24 vertical x 8 radial (frequency)
    VERTICAL_CHANNELS = 24
    RADIAL_CHANNELS = 8
    
    # GR device excluded channels (used by other diagnostics)
    GR_EXCLUDED_VERTICAL = [23, 24]
    
    def __init__(self, mds_ip):
        """Initialize with MDS+ server IP"""
        self.mds_ip = mds_ip
        self._position_cache = {}  # {(shot, device): {'R': 2D, 'Z': 2D, ...}}
    
    def load_data(self, shot_number, device='GT', v_ch=1, r_ch=1):
        """Load ECEI raw data for a specific channel"""
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        try:
            node_name = f'\\ECEI_{device}{v_ch:02d}{r_ch:02d}:FOO'
            data = mds.get(node_name).data()
            time = mds.get(f'dim_of({node_name})').data()
            
            return {
                'time': np.array(time),
                'data': np.array(data),
                'shot': shot_number,
                'device': device,
                'v_ch': v_ch,
                'r_ch': r_ch,
                'channel': f'{device}{v_ch:02d}{r_ch:02d}'
            }
        finally:
            mds.closeTree('kstar', shot_number)
    
    def get_channel_positions(self, shot_number, device):
        """Get R, Z, Angle positions for all channels of a device
        
        Returns dict with:
            R: (24, 8) array in [m]
            Z: (24, 8) array in [m]
            A: (24, 8) array in [rad]
            channels: list of channel names
            Bt, mode, LO, sz, sf: device parameters
        """
        cache_key = (shot_number, device)
        if cache_key in self._position_cache:
            return self._position_cache[cache_key]
        
        mds = Connection(self.mds_ip)
        mds.openTree('kstar', shot_number)
        
        try:
            # Load device parameters from MDS+
            Bt_raw = mds.get('\\ECEI_I_TF').data() * 0.0995556
            Bt = round(Bt_raw, 1)  # Round to 1 decimal place
            
            mode = mds.get(f'\\{device}_MODE').data()
            LO = mds.get(f'\\{device}_LOFREQ').data()
            sz = mds.get(f'\\{device}_LENSZOOM').data()
            sf = mds.get(f'\\{device}_LENSFOCUS').data()
            
            mds.closeTree('kstar', shot_number)
            
            # Calculate R positions for 8 radial channels [cm]
            # MATLAB: rr = 180*27.99*m*Bt./((7:-1:0)*0.9 + 2.6 + LO)
            rr_cm = 180 * 27.99 * mode * Bt / ((np.arange(7, -1, -1)) * 0.9 + 2.6 + LO)
            
            # R position array (24, 8) in [m]
            R_position = np.tile(rr_cm, (self.VERTICAL_CHANNELS, 1)) / 100  # cm to m
            
            # Calculate Z and Angle positions using ABCD matrix
            Z_position = np.zeros((self.VERTICAL_CHANNELS, self.RADIAL_CHANNELS))
            A_position = np.zeros((self.VERTICAL_CHANNELS, self.RADIAL_CHANNELS))
            
            for r_idx in range(self.RADIAL_CHANNELS):
                # Convert R from cm to mm for ABCD calculation
                R_mm = rr_cm[r_idx] * 10
                Z_col, A_col = self._calculate_z_angle(device, R_mm, sz, sf)
                Z_position[:, r_idx] = Z_col
                A_position[:, r_idx] = A_col
            
            # Build channel name list
            channels = []
            for v in range(1, self.VERTICAL_CHANNELS + 1):
                for r in range(1, self.RADIAL_CHANNELS + 1):
                    channels.append(f'{device}{v:02d}{r:02d}')
            
            result = {
                'R': R_position,      # [m], (24, 8)
                'Z': Z_position,      # [m], (24, 8)
                'A': A_position,      # [rad], (24, 8)
                'channels': channels,
                'Bt': Bt,
                'mode': mode,
                'LO': LO,
                'sz': sz,
                'sf': sf
            }
            
            self._position_cache[cache_key] = result
            return result
            
        except Exception as e:
            mds.closeTree('kstar', shot_number)
            raise e
    
    def _calculate_z_angle(self, device, R_mm, sz, sf):
        """Calculate Z and Angle positions using ABCD matrix
        
        R_mm: R position in mm
        Returns: (Z array [m], Angle array [rad]) for 24 vertical channels
        """
        sp = 2300 - R_mm
        
        # Vertical position from reference axis [mm]
        # MATLAB: zECEI = ((24:-1:1) - 12.5)*14
        z_ecei = (np.arange(24, 0, -1) - 12.5) * 14
        # Angle against reference axis at ECEI array box
        a_ecei = np.zeros(24)
        
        # Build ABCD matrix based on device
        if device in ['GT', 'GR']:
            ABCD = self._build_abcd_gt_gr(sp, sz, sf)
        else:  # HT
            ABCD = self._build_abcd_ht(sp, sz, sf)
        
        # Calculate z and angle at R
        # za = ABCD * [zECEI; aECEI]
        za = ABCD @ np.vstack([z_ecei, a_ecei])
        z_R = za[0, :] / 1000  # mm to m
        a_R = za[1, :]         # rad
        
        return z_R, a_R
    
    def _build_abcd_gt_gr(self, sp, sz, sf):
        """Build ABCD matrix for GT and GR devices (MATLAB formula)"""
        n = 1.52
        
        # Matrix multiplication in MATLAB order: M1 * M2 * M3 * ...
        matrices = [
            np.array([[1, sp + (1954 - sz)], [0, 1]]),
            np.array([[1, 0], [(n - 1) / (-1000), n]]),
            np.array([[1, 160], [0, 1]]),
            np.array([[1, 0], [(1 - n) / (1000 * n), 1 / n]]),
            np.array([[1, 2280 - (1954 + 160 - sz)], [0, 1]]),
            np.array([[1, 0], [(n - 1) / 1000, n]]),
            np.array([[1, 20], [0, 1]]),
            np.array([[1, 0], [0, 1 / n]]),
            np.array([[1, 4288 - (2280 + 20) - sf], [0, 1]]),
            np.array([[1, 0], [(n - 1) / (-1200), n]]),
            np.array([[1, 140], [0, 1]]),
            np.array([[1, 0], [(1 - n) / (1200 * n), 1 / n]]),
            np.array([[1, 4520 - (4288 + 140 - sf)], [0, 1]]),
            np.array([[1, 0], [0, n]]),
            np.array([[1, 30], [0, 1]]),
            np.array([[1, 0], [0, 1 / n]]),
            np.array([[1, 4940 - (4520 + 30)], [0, 1]]),
        ]
        
        # MATLAB order: M1 @ M2 @ M3 @ ...
        M = matrices[0]
        for m in matrices[1:]:
            M = M @ m
        
        return M
    
    def _build_abcd_ht(self, sp, sz, sf):
        """Build ABCD matrix for HT device (MATLAB formula)"""
        n = 1.526
        
        # Matrix multiplication in MATLAB order: M1 * M2 * M3 * ...
        matrices = [
            np.array([[1, sp + 2553.01], [0, 1]]),
            np.array([[1, 0], [(n - 1) / (-695), n]]),
            np.array([[1, 150], [0, 1]]),
            np.array([[1, 0], [0, 1 / n]]),
            np.array([[1, 4500.41 - (2553.01 + 150) - sz], [0, 1]]),
            np.array([[1, 0], [0, n]]),
            np.array([[1, 40], [0, 1]]),
            np.array([[1, 0], [(1 - n) / (-515 * n), 1 / n]]),
            np.array([[1, 6122.41 - (4500.41 + 40) - sf + sz], [0, 1]]),
            np.array([[1, 0], [0, n]]),
            np.array([[1, 150], [0, 1]]),
            np.array([[1, 0], [(1 - n) / (630 * n), 1 / n]]),
            np.array([[1, 6478.41 - (6122.41 + 150 - sf)], [0, 1]]),
            np.array([[1, 0], [0, n]]),
            np.array([[1, 40], [0, 1]]),
            np.array([[1, 0], [0, 1 / n]]),
            np.array([[1, 7161.01 - (6478.41 + 40)], [0, 1]]),
        ]
        
        # MATLAB order: M1 @ M2 @ M3 @ ...
        M = matrices[0]
        for m in matrices[1:]:
            M = M @ m
        
        return M
    
    def get_valid_channels(self, device):
        """Get list of valid channel tuples (v, r) for a device"""
        channels = []
        
        for v in range(1, self.VERTICAL_CHANNELS + 1):
            # Skip excluded channels for GR
            if device == 'GR' and v in self.GR_EXCLUDED_VERTICAL:
                continue
            
            for r in range(1, self.RADIAL_CHANNELS + 1):
                channels.append((v, r))
        
        return channels
    
    def get_channel_list(self, shot_number, device):
        """Get formatted channel list with positions for dropdown display
        
        Returns list of tuples: (display_string, channel_name, R, Z)
        """
        positions = self.get_channel_positions(shot_number, device)
        
        result = []
        for v in range(1, self.VERTICAL_CHANNELS + 1):
            # Skip excluded channels for GR
            if device == 'GR' and v in self.GR_EXCLUDED_VERTICAL:
                continue
            
            for r in range(1, self.RADIAL_CHANNELS + 1):
                ch_name = f'{device}{v:02d}{r:02d}'
                R = positions['R'][v - 1, r - 1]
                Z = positions['Z'][v - 1, r - 1]
                display = f'{ch_name} (R={R:.3f}m, Z={Z:.3f}m)'
                result.append((display, ch_name, R, Z))
        
        return result