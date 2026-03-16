"""
Neutron diagnostics data loader
Micro-Fission Chamber, He3 Counter, Diamond Detector near J-port
"""

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader, IP_FAULT_OFFSET
from core.data_structures import DiagnosticData


class NeutronLoader(BaseDiagnosticLoader):
    """Loader for fusion neutron diagnostics data from MDS+"""

    DETECTORS = {
        'fission': {
            'name': 'Fission Chamber',
            'node': '\\NTRN_INT05:FOO',
        },
        'he3': {
            'name': 'He3 Counter',
            'node': '\\NTRN_INT08:FOO',
        },
        'diamond': {
            'name': 'Diamond Detector',
            'node': '\\NTRN_INT03:FOO',
        },
    }

    def load_data(self, shot_number, analysis_type=None):
        """Load neutron detector data from MDS+"""
        mds = self._connect_mds(shot_number)

        try:
            measurements = {}

            for key, det in self.DETECTORS.items():
                node = det['node']
                try:
                    raw_data = mds.get(node).data()
                    raw_time = mds.get(f'dim_of({node})').data()

                    # Flatten to 1D (MDS+ may return multi-dimensional)
                    data = np.array(raw_data).flatten()
                    time = np.array(raw_time).flatten()

                    # Ensure matching lengths
                    n = min(len(data), len(time))
                    data = data[:n]
                    time = time[:n]

                    measurements[key] = {
                        'data': data,
                        'time': time,
                        'name': det['name'],
                    }
                    print(f"[Neutron] {det['name']}: {len(data)} points, "
                          f"t = {time[0]:.3f} ~ {time[-1]:.3f} s")

                except Exception as e:
                    print(f"[Neutron] {det['name']} not available: {str(e)}")

            if not measurements:
                raise ValueError("No neutron detector data available")

            # Apply IP fault masking
            ip_fault_time = self.get_ip_fault_time(shot_number)
            if ip_fault_time is not None:
                for key in measurements:
                    t = measurements[key]['time']
                    mask = self.get_valid_time_mask(t, ip_fault_time)
                    measurements[key]['time'] = t[mask]
                    measurements[key]['data'] = measurements[key]['data'][mask]
                print(f"[Neutron] Data masked to IP fault time + "
                      f"{IP_FAULT_OFFSET}s: {ip_fault_time + IP_FAULT_OFFSET:.3f} s")

            # Use first available detector's time as reference (after masking)
            ref_key = next(iter(measurements))
            ref_time = measurements[ref_key]['time']

            return DiagnosticData(
                time=ref_time,
                radius=np.array([1, 2, 3], dtype=float),
                measurements=measurements,
                source='mdsplus'
            )

        finally:
            self._close_mds(mds, shot_number)
