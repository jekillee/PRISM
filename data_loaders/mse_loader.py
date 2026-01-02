#!/usr/bin/python3.8

"""
MSE (Motional Stark Effect) data loader
"""

import numpy as np
from MDSplus import Connection
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class MSELoader(BaseDiagnosticLoader):
    """Loader for MSE diagnostic data"""
    
    # Number of raw channels and processed points
    N_RAW_CHANNELS = 25
    N_PROF_POINTS = 20
    
    # NB source drr table (mm -> m conversion in load_data)
    DRR_TABLE = {
        'NB-1A': [32.1913, 27.7642, 23.6813, 20.2123, 17.0620, 14.5040, 12.2766, 10.4646, 
                  8.85010, 7.63049, 6.62952, 5.89612, 5.40161, 5.11340, 5.02441, 5.11572, 
                  5.39307, 5.83228, 6.41235, 7.13867, 8.00610, 8.94043, 10.0667, 11.2344, 12.5549],
        'NB-1B': [124.427, 114.534, 105.092, 96.7568, 88.8531, 82.1125, 75.9073, 70.5181, 
                  65.3152, 60.9685, 56.9159, 53.3745, 50.2976, 47.5640, 45.1602, 43.1921, 
                  41.3933, 39.9307, 38.7666, 37.8420, 37.1472, 36.6963, 36.4189, 36.3428, 36.4434],
        'NB-1C': [32.1259, 32.8662, 33.3674, 33.6091, 33.6267, 33.4414, 33.0685, 32.5482, 
                  31.8264, 31.0122, 30.0271, 28.9324, 27.7430, 26.4329, 25.0051, 23.5601, 
                  21.9141, 20.2148, 18.4749, 16.6453, 14.7280, 12.8577, 10.7783, 8.75977, 6.59985]
    }
    
    def _determine_nb_source(self, aa1_value):
        """Determine NB source from AA1GAM value"""
        if aa1_value > 0.84:
            return 'NB-1B'
        elif aa1_value > 0.80:
            return 'NB-1A'
        else:
            return 'NB-1C'
    
    def _detect_bad_channels(self, tgamma):
        """Detect bad channels (mean TGAMMA ~ -1)"""
        tgamma_mean = np.mean(tgamma, axis=1)
        bad_mask = (tgamma_mean > -1.05) & (tgamma_mean < -0.95)
        return ~bad_mask  # Return good mask
    
    def get_ip_fault_time(self, shot_number):
        """Get IP fault time (last time when Ip > 100 kA)"""
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            ip_time = mds.get('dim_of(\\pcrc03)').data()
            ip_data = mds.get('\\pcrc03/-1e3').data()
            mds.closeTree('kstar', shot_number)
            
            valid_indices = np.where(ip_data > 100)[0]
            if len(valid_indices) > 0:
                return ip_time[valid_indices[-1]]
            return None
        except Exception as e:
            print(f"Warning: Could not get IP fault time: {str(e)}")
            return None
    
    def load_data(self, shot_number, analysis_type=None):
        """Load MSE data from MDS+"""
        try:
            mds = Connection(self.mds_ip)
            mds.openTree('kstar', shot_number)
            
            # =================================================================
            # Load raw MSE data: TGAMMA, SGAMMA (25 channels)
            # =================================================================
            time_raw = mds.get('dim_of(\\TGAMMA01)').data()
            nt_raw = len(time_raw)
            
            tgamma = np.zeros((self.N_RAW_CHANNELS, nt_raw))
            sgamma = np.zeros((self.N_RAW_CHANNELS, nt_raw))
            
            for i in range(self.N_RAW_CHANNELS):
                ch = f'{i+1:02d}'
                tgamma[i, :] = mds.get(f'\\TGAMMA{ch}').data()
                sgamma[i, :] = mds.get(f'\\SGAMMA{ch}').data()
            
            # =================================================================
            # Load processed MSE data: j, q, r/a profiles (20 points)
            # =================================================================
            time_prof = mds.get('dim_of(\\pmse_jv01)').data()
            nt_prof = len(time_prof)
            
            j_data = np.zeros((self.N_PROF_POINTS, nt_prof))
            j_err = np.zeros((self.N_PROF_POINTS, nt_prof))
            q_data = np.zeros((self.N_PROF_POINTS, nt_prof))
            q_err = np.zeros((self.N_PROF_POINTS, nt_prof))
            roa_data = np.zeros((self.N_PROF_POINTS, nt_prof))
            
            for i in range(self.N_PROF_POINTS):
                ch = f'{i+1:02d}'
                j_data[i, :] = mds.get(f'\\pmse_jv{ch}').data()
                j_err[i, :] = mds.get(f'\\pmse_je{ch}').data()
                q_data[i, :] = mds.get(f'\\pmse_qv{ch}').data()
                q_err[i, :] = mds.get(f'\\pmse_qe{ch}').data()
                roa_data[i, :] = mds.get(f'\\pmse_av{ch}').data()
            
            # q0, magnetic axis
            q0 = mds.get('\\pmse_q0v').data()
            q0_err = mds.get('\\pmse_q0e').data()
            magx = mds.get('\\pmse_magxv').data()
            magx_err = mds.get('\\pmse_magxe').data()
            
            # =================================================================
            # Load static data: R positions, NB source
            # =================================================================
            rr = mds.get('\\RRRGAM').data()
            aa1 = mds.get('\\AA1GAM').data()
            
            mds.closeTree('kstar', shot_number)
            
            # =================================================================
            # Load R_edge from efitrt1: rsurf + aminor
            # =================================================================
            mds.openTree('efitrt1', shot_number)
            time_efit = mds.get('\\gtime').data()
            rsurf = mds.get('\\rsurf').data()
            aminor = mds.get('\\aminor').data()
            mds.closeTree('efitrt1', shot_number)
            
            # R_edge = rsurf + aminor (time-dependent)
            r_edge_efit = rsurf + aminor
            
            # =================================================================
            # Determine NB source and drr
            # =================================================================
            nb_source = self._determine_nb_source(aa1[0])
            drr = np.array(self.DRR_TABLE[nb_source]) / 1000.0
            
            print(f"MSE data loaded: NB source = {nb_source}")
            print(f"  Raw: {tgamma.shape}, time = {time_raw[0]:.3f} ~ {time_raw[-1]:.3f} s")
            print(f"  Profile: {q_data.shape}, time = {time_prof[0]:.3f} ~ {time_prof[-1]:.3f} s")
            
            # =================================================================
            # Detect bad channels
            # =================================================================
            good_mask = self._detect_bad_channels(tgamma)
            bad_channels = np.where(~good_mask)[0] + 1
            if len(bad_channels) > 0:
                print(f"  Bad channels: {list(bad_channels)}")
            
            # =================================================================
            # Get IP fault time and mask data
            # =================================================================
            ip_fault_time = self.get_ip_fault_time(shot_number)
            if ip_fault_time is None:
                raise RuntimeError(f"Shot {shot_number}: Ip did not exceed 100 kA")
            
            # Mask raw data
            valid_raw_mask = (time_raw > 0) & (time_raw <= ip_fault_time)
            time_raw = time_raw[valid_raw_mask]
            tgamma = tgamma[:, valid_raw_mask]
            sgamma = sgamma[:, valid_raw_mask]
            
            # Mask profile data
            valid_prof_mask = (time_prof > 0) & (time_prof <= ip_fault_time)
            time_prof = time_prof[valid_prof_mask]
            j_data = j_data[:, valid_prof_mask]
            j_err = j_err[:, valid_prof_mask]
            q_data = q_data[:, valid_prof_mask]
            q_err = q_err[:, valid_prof_mask]
            roa_data = roa_data[:, valid_prof_mask]
            q0 = q0[valid_prof_mask]
            q0_err = q0_err[valid_prof_mask]
            magx = magx[valid_prof_mask]
            magx_err = magx_err[valid_prof_mask]
            
            print(f"  Data masked to IP fault time: {ip_fault_time:.3f} s")
            
            # Clean negative error values
            np.seterr(invalid='ignore')
            sgamma = np.clip(sgamma, 0, None)
            j_err = np.clip(j_err, 0, None)
            q_err = np.clip(q_err, 0, None)
            
            # =================================================================
            # Package measurements
            # =================================================================
            measurements = {
                'tgamma': {
                    'data': tgamma,
                    'error': sgamma,
                    'good_mask': good_mask,
                    'R': rr,
                    'drr': drr
                },
                'j': {
                    'data': j_data,
                    'error': j_err,
                    'roa': roa_data
                },
                'q': {
                    'data': q_data,
                    'error': q_err,
                    'roa': roa_data
                },
                'q0': {
                    'data': q0,
                    'error': q0_err
                },
                'magx': {
                    'data': magx,
                    'error': magx_err
                },
                'r_edge': {
                    'time': time_efit,
                    'data': r_edge_efit
                }
            }
            
            # Create result with raw time as primary
            result = DiagnosticData(time_raw, rr, measurements, source='mdsplus')
            result.time_prof = time_prof
            result.nb_source = nb_source
            
            return result
            
        except Exception as e:
            raise RuntimeError(f"Failed to load MSE data for shot {shot_number}: {str(e)}")
    
    def get_R_from_roa(self, roa, magx, r_edge):
        """Convert r/a to R coordinate
        
        R = R_axis + r/a * (R_edge - R_axis)
        """
        return magx + roa * (r_edge - magx)
