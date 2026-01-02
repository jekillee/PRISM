#!/usr/bin/python3.8

"""
File parser for CES diagnostic result files
"""

import numpy as np
from core.data_structures import DiagnosticData


class CESFileParser:
    """Parser for CES result files (TGF format)"""
    
    @staticmethod
    def parse_file(filepath):
        """Parse CES result file and return DiagnosticData"""
        try:
            with open(filepath, 'r') as f:
                lines = f.readlines()
        except Exception as e:
            raise ValueError("Cannot read file. Invalid file format.\nPlease contact Jekil Lee (jklee@kfe.re.kr)")
        
        # Parse header
        metadata = {}
        data_start_idx = None
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            if line.startswith('# Name ='):
                full_name = line.split('=')[1].strip().strip("'")
                if '_ces_' in full_name:
                    metadata['name'] = full_name.split('_ces_')[1]
                elif 'ces_' in full_name:
                    metadata['name'] = full_name.split('ces_')[1]
                else:
                    metadata['name'] = full_name
            
            elif line.startswith('# ShotNo ='):
                metadata['shot'] = int(line.split('=')[1].strip())
            
            elif line.startswith('# DimSize ='):
                dims = line.split('=')[1].strip().split(',')
                metadata['n_time'] = int(dims[0])
                metadata['n_radius'] = int(dims[1])
            
            elif line == '# [data]':
                data_start_idx = i + 1
                break
        
        if data_start_idx is None or 'name' not in metadata or 'shot' not in metadata:
            raise ValueError("Invalid file format.\nPlease contact Jekil Lee (jklee@kfe.re.kr)")
        
        # Parse data section
        data_lines = []
        for line in lines[data_start_idx:]:
            line = line.strip()
            if line and not line.startswith('#'):
                data_lines.append(line)
        
        # Convert to numpy array
        data_array = []
        for line in data_lines:
            values = [float(x.strip()) for x in line.split(',')]
            data_array.append(values)
        
        data_array = np.array(data_array)
        
        # Extract columns: Time, R, Ti, Tierr, vT, vTerr
        time_col = data_array[:, 0]
        radius_col = data_array[:, 1]
        ti_col = data_array[:, 2]
        tierr_col = data_array[:, 3]
        vt_col = data_array[:, 4]
        vterr_col = data_array[:, 5]
        
        # Get unique values
        unique_times = np.unique(time_col)
        unique_radii = np.unique(radius_col)
        
        n_time = len(unique_times)
        n_radius = len(unique_radii)
        
        # Reshape into 2D arrays (radius x time)
        temperature = np.full((n_radius, n_time), np.nan)
        temp_error = np.full((n_radius, n_time), np.nan)
        velocity = np.full((n_radius, n_time), np.nan)
        vel_error = np.full((n_radius, n_time), np.nan)
        
        for i, row in enumerate(data_array):
            t = row[0]
            r = row[1]
            
            t_idx = np.where(unique_times == t)[0][0]
            r_idx = np.where(unique_radii == r)[0][0]
            
            temperature[r_idx, t_idx] = row[2]
            temp_error[r_idx, t_idx] = row[3]
            velocity[r_idx, t_idx] = row[4]
            vel_error[r_idx, t_idx] = row[5]
        
        # Convert units: eV to keV
        temperature = temperature * 1e-3
        temp_error = temp_error * 1e-3
        
        # Clean invalid data
        temp_error[temp_error < 0] = np.nan
        vel_error[vel_error < 0] = np.nan
        temperature[temperature == 0] = np.nan
        
        # Package into DiagnosticData structure
        measurements = {
            'Ti': {'data': temperature, 'error': temp_error},
            'vT': {'data': velocity, 'error': vel_error}
        }
        
        return DiagnosticData(
            time=unique_times,
            radius=unique_radii,
            measurements=measurements,
            source=metadata.get('name', 'file')
        ), metadata.get('shot', 0)