#!/usr/bin/python3.8

"""
Thomson Scattering data loader implementation
"""

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData


class ThomsonLoader(BaseDiagnosticLoader):
    """Loader for Thomson Scattering diagnostic data"""
    
    def get_ip_fault_time(self, shot_number):
        """Get IP fault time (last time when Ip > 100 kA)"""
        try:
            mds = self._connect_mds(shot_number)
            ip_time = mds.get('dim_of(\\pcrc03)').data()
            ip_data = mds.get('\\pcrc03/-1e3').data()  # Convert to kA
            self._close_mds(mds, shot_number)
            
            # Find indices where Ip > 100 kA
            valid_indices = np.where(ip_data > 100)[0]
            if len(valid_indices) > 0:
                return ip_time[valid_indices[-1]]  # Last valid time
            return None
        except Exception as e:
            print(f"Warning: Could not get IP fault time: {str(e)}")
            return None
    
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
                        
                except:
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
                        
                except:
                    print(f'TS_EDGE{i+1} not available')
            
            self._close_mds(mds, shot_number)
            
            if not te_data:
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
            
            # Mask time range: 0 < time <= ip_fault_time
            valid_time_mask = (times > 0) & (times <= ip_fault_time)
            times = times[valid_time_mask]
            temperature = temperature[:, valid_time_mask]
            density = density[:, valid_time_mask]
            temp_error_upper = temp_error_upper[:, valid_time_mask]
            temp_error_lower = temp_error_lower[:, valid_time_mask]
            density_error_upper = density_error_upper[:, valid_time_mask]
            density_error_lower = density_error_lower[:, valid_time_mask]
            
            print(f"  Data masked to IP fault time: {ip_fault_time:.3f} s")
            
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
        """Get hardcoded channel positions based on shot number"""
        if shot_num < 21779:
            core = [1806, 1826, 1848, 1871, 1894, 1917, 1942, 1966, 1991, 2016, 2041, 2068, 2093, 2119]
            edge = [2124, 2137, 2143, 2149, 2156, 2162, 2177, 2191, 2202, 2216, 2229, 2242, 2257, 2271, 2285, 2297, 2311]
        elif 21779 <= shot_num <= 24081:
            core = [1806, 1826, 1848, 1871, 1894, 1917, 1942, 1966, 1991, 2016, 2041, 2068, 2093, 2119]
            edge = [2124, 2137, 2143, 2149, 2156, 2162, 2177, 2191, 2202, 2216, 2229, 2242, 2257, 2271, 2285, 2297, 2311]
        elif 24082 <= shot_num <= 27400:
            core = [1802, 1823, 1845, 1866, 1889, 1912, 1937, 1961, 1986, 2011, 2036, 2061, 2086, 2114]
            edge = [2114, 2126, 2138, 2150, 2164, 2177, 2189, 2196, 2202, 2208, 2214, 2221, 2235, 2248, 2261, 2275, 2289]
        elif 27401 <= shot_num <= 30445:
            core = [1797, 1818, 1841, 1862, 1884, 1908, 1931, 1954, 1979, 2004, 2030, 2056, 2082, 2108]
            edge = [2108, 2120, 2133, 2146, 2153, 2171, 2183, 2190, 2197, 2203, 2209, 2216, 2229, 2243, 2254, 2269, 2282]
        elif 30446 <= shot_num <= 32768:
            core = [1800, 1820, 1842, 1864, 1887, 1909, 1933, 1956, 1983, 2008, 2033, 2058, 2083, 2109]
            edge = [2111, 2124, 2136, 2149, 2161, 2174, 2187, 2194, 2200, 2207, 2212, 2219, 2232, 2245, 2259, 2272, 2285]
        elif 32769 <= shot_num <= 34836:
            core = [1800, 1820, 1842, 1864, 1887, 1911, 1933, 1958, 1981, 2006, 2030, 2058, 2083, 2110]
            edge = [2111, 2124, 2137, 2147, 2162, 2175, 2187, 2195, 2200, 2205, 2210, 2217, 2233, 2245, 2259, 2271, 2290]
        elif 34837 <= shot_num <= 37741:
            core = [1795, 1817, 1839, 1861, 1884, 1908, 1930, 1954, 1978, 2003, 2029, 2054, 2081, 2108]
            edge = [2108, 2116, 2125, 2136, 2145, 2155, 2163, 2171, 2180, 2191, 2200, 2213, 2222, 2232, 2253, 2274, 2296]
        else:
            core = [1785, 1807, 1828, 1851, 1873, 1896, 1920, 1944, 1968, 1993, 2018, 2043, 2070, 2095]
            edge = [2096, 2105, 2113, 2123, 2131, 2141, 2149, 2158, 2168, 2178, 2197, 2217, 2237, 2258, 2280, None, None]
        
        return core, edge