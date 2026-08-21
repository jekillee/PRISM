"""
IRVB (Infra-Red Video Bolometer) data loaders

- IRVBLoader:     reconstructed 2D radiation (.mat) from the IRVB HTTP server
- IRVBPradLoader: total radiated power node \\IRVB1_PRAD from MDS+ (time trace)
"""

import os
import numpy as np
import scipy.io
import urllib.request
import urllib.error

from config.app_config import PRISM_NAS_ROOT, ensure_shared_dir
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


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

    LOCAL_CACHE_DIR = os.path.join(PRISM_NAS_ROOT, 'irvb')

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
        """Create the shared IRVB cache dir on the PRISM NAS (/PRISM/irvb),
        best-effort world-writable so all users can refresh it."""
        ensure_shared_dir(self.LOCAL_CACHE_DIR)      # IRVB cache: /PRISM/irvb
    
    def load_data(self, shot_number):
        """Load IRVB data from server or cache"""
        mat_file = f'shot-{shot_number:06d}_data.mat'
        local_path = os.path.join(self.LOCAL_CACHE_DIR, mat_file)

        # Always re-download a fresh copy (IRVB reconstructions can be updated), even
        # if a file already exists. The server serves recent shots at
        # <IRVB_SERVER>/<file> and post-campaign (older) shots at
        # <IRVB_SERVER>/afterCampaign/<file>, so try both. Download to a temp file
        # and atomically replace, so an existing good file survives a failed fetch.
        base = self.config.IRVB_SERVER.rstrip('/')
        candidates = [f'{base}/{mat_file}', f'{base}/afterCampaign/{mat_file}']
        tmp_path = f'{local_path}.{os.getpid()}.tmp'
        downloaded = False
        for url in candidates:
            print(f"[IRVB] Downloading {mat_file} from {url} ...")
            try:
                urllib.request.urlretrieve(url, tmp_path)
                os.replace(tmp_path, local_path)   # atomic overwrite of any prior copy
                try:
                    os.chmod(local_path, 0o666)    # let other users overwrite
                except OSError:
                    pass
                downloaded = True
                break
            except urllib.error.HTTPError as e:
                if os.path.isfile(tmp_path):
                    os.remove(tmp_path)
                if e.code == 404:
                    continue   # not at this location — try the next
                raise RuntimeError(
                    f"IRVB server returned HTTP {e.code} for shot #{shot_number} ({url})"
                )
            except Exception as e:
                if os.path.isfile(tmp_path):
                    os.remove(tmp_path)
                raise RuntimeError(f"Failed to download IRVB data: {str(e)}")
        if not downloaded:
            raise RuntimeError(
                f"No IRVB data on the server for shot #{shot_number} (404 under "
                f"{base}/ and {base}/afterCampaign/). It may have no IRVB reconstruction."
            )

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
            
            print(f"[IRVB] Loaded {len(time)} timepoints, grid {self.NX}x{self.NY}")
            
            return IRVBData(time, recon, ptot, x_grid, y_grid)

        except Exception as e:
            raise RuntimeError(f"Failed to load IRVB data: {str(e)}")


class IRVBPradLoader(BaseDiagnosticLoader):
    """Loader for the IRVB total radiated power node ``\\IRVB1_PRAD`` (MDS+).

    Used by the IRVB Time Trace tab to plot / compare ``\\IRVB1_PRAD`` (the value
    the IRVB analysis publishes to MDS) across shots.
    """

    NODE = '\\IRVB1_PRAD'

    def load_data(self, shot_number, analysis_type=None, apply_mask=True):
        """Load ``\\IRVB1_PRAD`` (MW) vs time from MDS+.

        apply_mask=True trims to the IP-fault window (like the other time traces);
        pass False to get the full-range node as-is (used by the IRVB imaging tab).
        """
        mds = self._connect_mds(shot_number)
        try:
            data = np.asarray(mds.get(self.NODE).data()).flatten()
            time = np.asarray(mds.get(f'dim_of({self.NODE})').data()).flatten()

            n = min(len(data), len(time))
            data = data[:n]
            time = time[:n]

            # IP fault masking (consistent with other time traces)
            if apply_mask:
                ip_fault_time = self.get_ip_fault_time(shot_number)
                if ip_fault_time is not None:
                    mask = self.get_valid_time_mask(time, ip_fault_time)
                    time = time[mask]
                    data = data[mask]

            if len(time) == 0:
                raise ValueError(f"No {self.NODE} data for shot #{shot_number}")

            print(f"[IRVB1_PRAD] {len(data)} points, "
                  f"t = {time[0]:.3f} ~ {time[-1]:.3f} s")

            measurements = {
                'prad': {'time': time, 'data': data, 'name': 'IRVB1 Prad'},
            }
            return DiagnosticData(
                time=time,
                radius=np.array([1.0]),
                measurements=measurements,
                source='mdsplus',
            )
        finally:
            self._close_mds(mds, shot_number)