#!/usr/bin/python3.8

"""
EFIT equilibrium data loader
Supports both 1D profile mapping and 2D poloidal cross-section data
"""

import numpy as np
from MDSplus import Connection
from scipy.interpolate import interp1d, interp2d
from core.data_structures import EFITData, EFITData2D


class EFITLoader:
    """Loader for EFIT equilibrium data"""
    
    def __init__(self, config):
        self.config = config
        self.mds_ip = config.MDS_IP
    
    def load_efit_data(self, shot_number, efit_tree=None):
        """Load EFIT equilibrium data for 1D profile mapping (R -> psi_n, rho)"""
        if efit_tree is None:
            efit_tree = self.config.DEFAULT_EFIT_TREE
            
        try:
            mds = Connection(self.mds_ip)
            mds.openTree(efit_tree, shot_number)
            
            time_arr = mds.get('\\gtime').data()
            # efitrt1, efitrt2 use seconds; others use milliseconds
            if efit_tree not in ['efitrt1', 'efitrt2']:
                time_arr = time_arr / 1e3
            
            nw = int(mds.get('\\mw').data()[0])
            r_grid = mds.get('\\r').data()
            z_grid = mds.get('\\z').data()
            
            simag = mds.get('\\ssimag').data()
            sibry = mds.get('\\ssibry').data()
            psi_rz = mds.get('\\psirz').data()
            qpsi = mds.get('\\qpsi').data()
            
            mds.closeTree(efit_tree, shot_number)
            
            psi_for_qpsi = np.linspace(0, 1.0, nw)
            psi_n_frames, rho_pol_frames, rho_tor_frames = [], [], []
            
            for t in range(len(time_arr)):
                psi_n, rho_pol, rho_tor = self._process_efit_frame(
                    r_grid, z_grid, psi_rz[t], simag[t], sibry[t], 
                    qpsi[t], psi_for_qpsi, nw
                )
                psi_n_frames.append(psi_n)
                rho_pol_frames.append(rho_pol)
                rho_tor_frames.append(rho_tor)
            
            return EFITData(time_arr, r_grid, psi_n_frames, rho_pol_frames, rho_tor_frames)
            
        except Exception as e:
            raise RuntimeError(f"Failed to load EFIT data for shot {shot_number} from {efit_tree}: {str(e)}")
    
    def load_efit_2d(self, shot_number, efit_tree=None):
        """Load 2D EFIT data for poloidal cross-section visualization"""
        if efit_tree is None:
            efit_tree = self.config.DEFAULT_EFIT_TREE
        
        try:
            mds = Connection(self.mds_ip)
            mds.openTree(efit_tree, shot_number)
            
            time_arr = mds.get('\\gtime').data()
            if efit_tree not in ['efitrt1', 'efitrt2']:
                time_arr = time_arr / 1e3
            
            r_grid = mds.get('\\r').data()
            z_grid = mds.get('\\z').data()
            simag = mds.get('\\ssimag').data()
            sibry = mds.get('\\ssibry').data()
            psirz = mds.get('\\psirz').data()
            
            # Boundary data - bdry_data shape: [nt, npts, 2]
            bdry_data = mds.get('\\bdry').data()
            bdry_r = bdry_data[:, :, 0]
            bdry_z = bdry_data[:, :, 1]
            nbdry = mds.get('\\nbdry').data()
            
            # Magnetic axis
            rmaxis = mds.get('\\rmaxis').data()
            zmaxis = mds.get('\\zmaxis').data()
            
            # Limiter data (optional)
            try:
                lim_data = mds.get('\\lim').data()
                limiter_r = lim_data[:, 0]
                limiter_z = lim_data[:, 1]
            except:
                limiter_r = None
                limiter_z = None
            
            mds.closeTree(efit_tree, shot_number)
            
            return EFITData2D(
                time_arr, r_grid, z_grid, psirz, simag, sibry,
                bdry_r, bdry_z, nbdry, limiter_r, limiter_z,
                rmaxis, zmaxis
            )
            
        except Exception as e:
            raise RuntimeError(f"Failed to load 2D EFIT data for shot {shot_number} from {efit_tree}: {str(e)}")
    
    def _process_efit_frame(self, r_grid, z_grid, psi_rz, simag, sibry, qpsi, psi_for_qpsi, nw):
        """Process a single EFIT time frame for 1D profile mapping"""
        interp_psi_rz = interp2d(r_grid, z_grid, psi_rz)
        psi_rz_z0 = interp_psi_rz(r_grid, 0)
        psi_n = (psi_rz_z0 - simag) / (sibry - simag)
        
        rho_pol = np.sqrt(psi_n)
        
        interp_qpsi = interp1d(psi_for_qpsi[:-1], qpsi[:-1], fill_value='extrapolate')
        
        phi = np.zeros(nw)
        for i in range(nw-1):
            phi[i+1] = np.trapz(interp_qpsi(psi_for_qpsi)[0:i+2], x=psi_for_qpsi[0:i+2])
        
        rho_tor_from_qpsi = np.sqrt(phi / phi[-1])
        interp_rho_tor = interp1d(psi_for_qpsi, rho_tor_from_qpsi, fill_value='extrapolate')
        rho_tor = interp_rho_tor(psi_n)
        
        # Handle HFS (High Field Side)
        hfs_indices = [i for i in range(len(psi_n)) if i > 0 and psi_n[i] < psi_n[i-1]]
        psi_n[hfs_indices] = -psi_n[hfs_indices]
        rho_pol[hfs_indices] = -rho_pol[hfs_indices]
        rho_tor[hfs_indices] = -rho_tor[hfs_indices]
        
        return psi_n, rho_pol, rho_tor