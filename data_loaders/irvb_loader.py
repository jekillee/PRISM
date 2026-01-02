#!/usr/bin/python3.8

"""
IRVB (Infra-Red Video Bolometer) data loader
Loads reconstructed radiation data from HTTP server
"""

import os
import numpy as np
import scipy.io
import urllib.request


class IRVBData:
    """IRVB data structure"""
    def __init__(self, time, recon, ptot, x_grid, y_grid):
        self.time = time          # 1D array [s]
        self.recon = recon        # 3D array [time, y, x] [MW/m^3]
        self.ptot = ptot          # 1D array [MW]
        self.x_grid = x_grid      # 1D array [m] (R direction)
        self.y_grid = y_grid      # 1D array [m] (Z direction)


class IRVBLoader:
    """Loader for IRVB radiation data"""
    
    IRVB_SERVER = 'http://172.17.112.125/data_ana'
    LOCAL_CACHE_DIR = '/tmp/prism_irvb_cache'
    
    # IRVB grid parameters (fixed)
    NX = 100
    NY = 100
    X_START = 1.2555   # R start [m]
    X_STEP = 0.011     # R step [m]
    Y_START = -1.4355  # Z start [m]
    Y_STEP = 0.029     # Z step [m]
    
    def __init__(self, config):
        self.config = config
        self._ensure_cache_dir()
    
    def _ensure_cache_dir(self):
        """Create cache directory with full permissions for all users"""
        if not os.path.isdir(self.LOCAL_CACHE_DIR):
            os.makedirs(self.LOCAL_CACHE_DIR, mode=0o777, exist_ok=True)
            try:
                os.chmod(self.LOCAL_CACHE_DIR, 0o777)
            except PermissionError:
                pass
    
    def load_data(self, shot_number):
        """Load IRVB data from server or cache"""
        mat_file = f'shot-{shot_number:06d}_data.mat'
        local_path = os.path.join(self.LOCAL_CACHE_DIR, mat_file)
        url = f'{self.IRVB_SERVER}/{mat_file}'
        
        # Download if not cached
        if not os.path.isfile(local_path):
            print(f"IRVB: Downloading {mat_file} from server...")
            try:
                urllib.request.urlretrieve(url, local_path)
                # Set file permissions so other users can overwrite
                os.chmod(local_path, 0o666)
                print(f"IRVB: Download complete")
            except Exception as e:
                raise RuntimeError(f"Failed to download IRVB data: {str(e)}")
        else:
            print(f"IRVB: Using cached file {mat_file}")
        
        # Load .mat file
        try:
            data = scipy.io.loadmat(local_path)
            
            recon = np.array(data['recon_MW'])         # [time, y, x]
            time = np.array(data['time_Sec'][0])       # [s]
            ptot = np.array(data['Ptot_MW'])[0]        # [MW]
            
            # Generate grid arrays
            x_grid = np.linspace(
                self.X_START, 
                self.X_START + self.X_STEP * (self.NX - 1), 
                self.NX
            )
            y_grid = np.linspace(
                self.Y_START, 
                self.Y_START + self.Y_STEP * (self.NY - 1), 
                self.NY
            )
            
            print(f"IRVB: Loaded {len(time)} timepoints, grid {self.NX}x{self.NY}")
            
            return IRVBData(time, recon, ptot, x_grid, y_grid)
            
        except Exception as e:
            raise RuntimeError(f"Failed to load IRVB data: {str(e)}")