"""
CES data loader implementation with analysis type support
"""

import csv
from pathlib import Path

import numpy as np
from data_loaders.base_loader import BaseDiagnosticLoader
from core.data_structures import DiagnosticData

# Fallback radius calibration table (channel positions in mm) keyed by shot range.
# Used when the \CES_RTxx position nodes are missing from MDS+ (typically old shots).
_RADIUS_CSV_PATH = Path(__file__).parent.parent / 'config' / 'radius_tces.csv'


class CESLoader(BaseDiagnosticLoader):
    """Loader for CES diagnostic data"""

    def _load_radius_from_csv(self, shot_number, n_channels):
        """Load channel radii from config/radius_tces.csv as a fallback.

        The CSV holds calibration rows keyed by shotNo; for a given shot the row
        with the largest shotNo <= shot_number is used. Returns raw values (mm),
        so the caller still applies the configured radius scale.
        """
        rows = []
        with open(_RADIUS_CSV_PATH, newline='', encoding='utf-8-sig') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    rows.append((int(row['shotNo']), row))
                except (KeyError, ValueError):
                    continue

        if not rows:
            raise RuntimeError(f"No calibration rows in {_RADIUS_CSV_PATH.name}")

        rows.sort(key=lambda r: r[0])
        # Largest shotNo <= shot_number, else the earliest row (baseline)
        selected = rows[0]
        for shot_no, row in rows:
            if shot_no <= shot_number:
                selected = (shot_no, row)
            else:
                break

        cal_shot, cal_row = selected
        print(f"[CES] Using radius_tces.csv fallback for shot {shot_number} "
              f"(calibration row shotNo={cal_shot})")

        radius_list = []
        for i in range(1, n_channels + 1):
            radius_list.append(float(cal_row[f'ch{i:02d}']))
        return radius_list

    def load_data(self, shot_number, analysis_type='mod'):
        """Load CES data from MDS+

        analysis_type: 'mod' for beam modulation or 'nn' for neural network
        """
        try:
            mds = self._connect_mds(shot_number)

            # Get node prefix based on analysis type
            analysis_config = self.diag_config['analysis_types'].get(analysis_type,
                             self.diag_config['analysis_types']['mod'])
            node_prefix = analysis_config['node_prefix']

            # Load time
            time_node = self.diag_config['time']['mds_node_format'] % node_prefix
            time_arr = mds.get(time_node).data()

            n_channels = self.diag_config['channels']

            # Load radius for all channels; fall back to CSV if RT nodes are missing
            radius_list = []
            try:
                for i in range(1, n_channels + 1):
                    radius_node = self.diag_config['radius']['mds_node_format'] % (node_prefix, i)
                    radius_list.append(mds.get(radius_node).data()[0])
            except Exception as e:
                print(f"[CES] RT position nodes not available for shot {shot_number} "
                      f"({str(e).strip()}); falling back to config/radius_tces.csv")
                radius_list = self._load_radius_from_csv(shot_number, n_channels)

            # Load temperature and velocity for all channels
            ti_list, ti_err_list, vt_list, vt_err_list = [], [], [], []
            for i in range(1, n_channels + 1):
                # Temperature
                ti_node = self.diag_config['parameters'][0]['mds_node_format'] % (node_prefix, i)
                ti_err_node = self.diag_config['parameters'][0]['mds_error_format'] % (node_prefix, i)
                ti_list.append(mds.get(ti_node).data())
                ti_err_list.append(mds.get(ti_err_node).data())

                # Velocity
                vt_node = self.diag_config['parameters'][1]['mds_node_format'] % (node_prefix, i)
                vt_err_node = self.diag_config['parameters'][1]['mds_error_format'] % (node_prefix, i)
                vt_list.append(mds.get(vt_node).data())
                vt_err_list.append(mds.get(vt_err_node).data())

            self._close_mds(mds, shot_number)

            # Convert to numpy arrays and apply scaling
            np.seterr(invalid='ignore')

            radius = np.array(radius_list) * self.diag_config['radius']['scale']
            ti_data = np.array(ti_list) * self.diag_config['parameters'][0]['scale']
            ti_err = np.array(ti_err_list) * self.diag_config['parameters'][0]['scale']
            vt_data = np.array(vt_list) * self.diag_config['parameters'][1]['scale']
            vt_err = np.array(vt_err_list) * self.diag_config['parameters'][1]['scale']

            # Clean negative error values
            ti_err[ti_err < 0] = np.nan
            vt_err[vt_err < 0] = np.nan

            # Package measurements
            measurements = {
                'Ti': {'data': ti_data, 'error': ti_err},
                'vT': {'data': vt_data, 'error': vt_err}
            }

            return DiagnosticData(time_arr, radius, measurements,
                                 source='mdsplus', analysis_type=analysis_type)

        except Exception as e:
            raise RuntimeError(f"Failed to load CES data for shot {shot_number}: {str(e)}")
