import csv
import math
import requests
from .N_B_variables import N_index_normalized, B_index
import urllib3
import argparse
import sys
import os
from datetime import datetime
from . import MiddleMan as OTSOLib

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

def convert_to_gsw_coordinates(datetime_str, vx_gse, vy_gse, vz_gse, bx_gse, by_gse, bz_gse, g, h):
    try:
        dt = datetime.strptime(datetime_str, '%Y-%m-%d %H:%M:%S')
        year = dt.year
        day_of_year = dt.timetuple().tm_yday
        hour = dt.hour
        minute = dt.minute
        second = dt.second
        
        date = [year, day_of_year, hour, minute, second,0]
        position_gse = [bx_gse, by_gse, bz_gse]
        Wind = [vx_gse, vy_gse, vz_gse]

        
        position_gsw = OTSOLib.gse2gswtsy15(date, position_gse, Wind, g, h)
        
        bt_gsw = math.sqrt(position_gsw[1]**2 + position_gsw[2]**2)

        return float(position_gsw[0]), float(position_gsw[1]), float(position_gsw[2]), float(bt_gsw)

    except Exception as e:
        print(f"Warning: GSW conversion failed for {datetime_str}: {e}")
        return None, None, None, None

def is_valid_value(value_str):
    if not value_str or value_str.strip() == '':
        return False
    
    try:
        val = float(value_str)
        # Check for common missing data indicators
        if val in [9999.99, 99999.9, 99.99, 999.99, -9999.99, -99999.9, -99.99, -999.99]:
            return False
        if math.isnan(val) or math.isinf(val):
            return False
        return True
    except (ValueError, TypeError):
        return False

def compute_corrected_symh(symh_values, pdyn_values):
    try:
        import numpy as np
        
        symh_array = np.array(symh_values, dtype=float)
        pdyn_array = np.array(pdyn_values, dtype=float)
    
        symh_clean = np.where(
            (symh_array == -9223372036854775808) | 
            (symh_array >= 99999) | 
            (symh_array <= -99999) | 
            (~np.isfinite(symh_array)), 
            np.nan, 
            symh_array
        )
        
        pdyn_clean = np.where(
            (pdyn_array < 0) | 
            (pdyn_array >= 99999) | 
            (~np.isfinite(pdyn_array)), 
            np.nan, 
            pdyn_array
        )
        
        symh_30min = compute_rolling_average_centered(symh_clean, window_size=6)
        pdyn_30min = compute_rolling_average_centered(pdyn_clean, window_size=6)
        
        valid_mask = ~np.isnan(symh_30min) & ~np.isnan(pdyn_30min)
        symh_corrected = np.full_like(symh_array, np.nan)
        symh_corrected[valid_mask] = 0.8 * symh_30min[valid_mask] - 13 * np.sqrt(pdyn_30min[valid_mask])
        
        return symh_corrected.tolist()
        
    except ImportError:
        print("Warning: numpy not available, using basic corrected SYM-H computation")
        return compute_corrected_symh_basic(symh_values, pdyn_values)
    except Exception as e:
        print(f"Error computing corrected SYM-H: {e}")
        return [None] * len(symh_values)

def compute_corrected_symh_basic(symh_values, pdyn_values):
    result = []
    
    for i, (symh, pdyn) in enumerate(zip(symh_values, pdyn_values)):
        try:
            # Check for valid values
            if (symh is None or pdyn is None or 
                not is_valid_symh_value(symh) or not is_valid_pdyn_value(pdyn)):
                result.append(None)
                continue
            
            # Simple 30-minute centered average (6 points)
            start_idx = max(0, i - 3)
            end_idx = min(len(symh_values), i + 4)
            
            # Get valid values in window
            symh_window = []
            pdyn_window = []
            
            for j in range(start_idx, end_idx):
                if (j < len(symh_values) and j < len(pdyn_values) and 
                    is_valid_symh_value(symh_values[j]) and is_valid_pdyn_value(pdyn_values[j])):
                    symh_window.append(float(symh_values[j]))
                    pdyn_window.append(float(pdyn_values[j]))
            
            if len(symh_window) >= 3 and len(pdyn_window) >= 3:
                avg_symh = sum(symh_window) / len(symh_window)
                avg_pdyn = sum(pdyn_window) / len(pdyn_window)
                corrected = 0.8 * avg_symh - 13 * math.sqrt(avg_pdyn)
                result.append(corrected)
            else:
                result.append(None)
                
        except Exception:
            result.append(None)
    
    return result

def compute_rolling_average_centered(values, window_size=6):
    try:
        import numpy as np
        n = len(values)
        rolling_avg = np.full(n, np.nan)
        half_window = window_size // 2
        
        input_values = np.array(values)
        
        for i in range(n):
            start_idx = max(0, i - half_window)
            end_idx = min(n, i + half_window + 1)
            window_values = input_values[start_idx:end_idx]
            
            valid_values = window_values[~np.isnan(window_values)]
            if len(valid_values) > 0:
                rolling_avg[i] = np.mean(valid_values)
        return rolling_avg
    except ImportError:
        return compute_rolling_average_centered_basic(values, window_size)

def compute_rolling_average_centered_basic(values, window_size=6):
    n = len(values)
    result = []
    half_window = window_size // 2
    
    for i in range(n):
        start_idx = max(0, i - half_window)
        end_idx = min(n, i + half_window + 1)
        
        valid_values = []
        for j in range(start_idx, end_idx):
            if j < len(values) and values[j] is not None:
                try:
                    val = float(values[j])
                    if not math.isnan(val) and math.isfinite(val):
                        valid_values.append(val)
                except (ValueError, TypeError):
                    pass
        
        if len(valid_values) > 0:
            result.append(sum(valid_values) / len(valid_values))
        else:
            result.append(None)
    
    return result

def is_valid_symh_value(value):
    try:
        val = float(value)
        return not (val >= 99999 or val <= -99999 or math.isnan(val) or not math.isfinite(val))
    except (ValueError, TypeError):
        return False

def is_valid_pdyn_value(value):
    try:
        val = float(value)
        return not (val < 0 or val >= 99999 or math.isnan(val) or not math.isfinite(val))
    except (ValueError, TypeError):
        return False

def compute_instantaneous_indices(row, g, h):
    try:
        if not (is_valid_value(row['BX_nT_GSE_GSM']) and
                is_valid_value(row['BY_nT_GSE']) and 
                is_valid_value(row['BZ_nT_GSE']) and
                is_valid_value(row['Vx_Velocity_km_s']) and
                is_valid_value(row['Vy_Velocity_km_s']) and
                is_valid_value(row['Vz_Velocity_km_s']) and
                is_valid_value(row['Proton_Density_n_cc'])):
            return (None, None, None, None, None, None, None, None, None, None, None, None)
        
        # Extract GSE magnetic field components
        BX_GSE = float(row['BX_nT_GSE_GSM'])
        BY_GSE = float(row['BY_nT_GSE'])
        BZ_GSE = float(row['BZ_nT_GSE'])
        
        # Extract velocity components
        Vx = float(row['Vx_Velocity_km_s'])
        Vy = float(row['Vy_Velocity_km_s'])
        Vz = float(row['Vz_Velocity_km_s'])
        
        # Extract proton density
        Np = float(row['Proton_Density_n_cc'])
        
        # Convert this instantaneous point from GSE to GSW
        bx_gsw, by_gsw, bz_gsw, bt_gsw = convert_to_gsw_coordinates(
            row['DateTime'], Vx, Vy, Vz, BX_GSE, BY_GSE, BZ_GSE, g, h
        )
        
        if any(b is None for b in [bx_gsw, by_gsw, bz_gsw]):
            return (None, None, None, None, None, None, None, None, None, None, None, None)
        
        # Calculate solar wind speed
        V = math.sqrt(Vx**2 + Vy**2 + Vz**2)
        
        # Calculate clock angle in GSW coordinates
        thetac_gsw = math.atan2(by_gsw, bz_gsw)
        if thetac_gsw < 0:
            thetac_gsw += 2 * math.pi
        
        # Compute N and B indices for this instantaneous point
        N_norm = N_index_normalized(V, bt_gsw, thetac_gsw)
        B_idx = B_index(Np, V, bt_gsw, thetac_gsw)
        
        return (N_norm, B_idx, bx_gsw, by_gsw, bz_gsw, bt_gsw, V, thetac_gsw, Vx, Vy, Vz, Np)
        
    except Exception as e:
        print(f"Warning: Failed to compute indices for {row.get('DateTime', 'unknown')}: {e}")
        return (None, None, None, None, None, None, None, None, None, None, None, None)

def calculate_rolling_average_indices(all_rows, window_size=7):
    results = []
    
    for i, current_row in enumerate(all_rows):
        result = {
            'DateTime': current_row['DateTime'],
            'BX_nT_GSW': current_row.get('BX_nT_GSW', ''),
            'BY_nT_GSW': current_row.get('BY_nT_GSW', ''),
            'BZ_nT_GSW': current_row.get('BZ_nT_GSW', ''),
            'Bt_nT_GSW': current_row.get('Bt_nT_GSW', ''),
            'V_km_s': current_row.get('V_km_s', ''),
            'Vx_km_s': current_row.get('Vx_km_s', ''),
            'Vy_km_s': current_row.get('Vy_km_s', ''),
            'Vz_km_s': current_row.get('Vz_km_s', ''),
            'Np_cm3': current_row.get('Np_cm3', ''),
            'Clock_angle_rad': current_row.get('Clock_angle_rad', ''),
            'Instantaneous_N_index': current_row.get('N_index_inst', ''),
            'Instantaneous_B_index': current_row.get('B_index_inst', ''),
            'N_index': '',
            'B_index': '',
            'SYM_H': current_row.get('SYM_H', ''),
            'Valid_points_in_window': ''
        }
        
        # Need enough prior data for a full window
        if i < (window_size - 1):
            results.append(result)
            continue
        
        # Extract the window: current point + (window_size-1) previous points
        start_idx = i - (window_size - 1)
        end_idx = i + 1
        window_data = all_rows[start_idx:end_idx]
        
        # Collect valid N and B index values from the window
        valid_N_indices = []
        valid_B_indices = []
        
        for row in window_data:
            if row.get('N_index_inst') and row['N_index_inst'] != '':
                try:
                    valid_N_indices.append(float(row['N_index_inst']))
                except (ValueError, TypeError):
                    pass
            
            if row.get('B_index_inst') and row['B_index_inst'] != '':
                try:
                    valid_B_indices.append(float(row['B_index_inst']))
                except (ValueError, TypeError):
                    pass
        
        # Compute rolling averages if we have enough valid points
        min_points = 3
        
        if len(valid_N_indices) >= min_points:
            result['N_index'] = f"{sum(valid_N_indices) / len(valid_N_indices):.6f}"
        
        if len(valid_B_indices) >= min_points:
            result['B_index'] = f"{sum(valid_B_indices) / len(valid_B_indices):.6f}"
        
        result['Valid_points_in_window'] = str(min(len(valid_N_indices), len(valid_B_indices)))
        
        results.append(result)
    
    return results

def convert_omni_to_tempfile(omni_filename, output_filename='tempfile.csv'):
    # print(f"\nConverting {omni_filename} to {output_filename}...")
    
    # Column specifications for OMNI 5-minute fixed-width format
    colspecs = [
        (0, 4),    # Year
        (4, 8),    # Day
        (8, 11),   # Hour
        (11, 14),  # Minute
        (67, 75),  # Bx GSE
        (75, 83),  # By GSE
        (83, 91),  # Bz GSE
        (91, 99),  # By GSM
        (99, 107), # Bz GSM
        (131, 139), # Vx GSE
        (139, 147), # Vy GSE
        (147, 155), # Vz GSE
        (155, 162), # Proton density
        (171, 177), # Flow pressure
        (269, 275), # SYM/H index
    ]
    
    column_names = [
        'Year', 'Day', 'Hour', 'Minute',
        'Bx_GSE_nT', 'By_GSE_nT', 'Bz_GSE_nT',
        'By_GSM_nT', 'Bz_GSM_nT',
        'Vx_GSE_km_s', 'Vy_GSE_km_s', 'Vz_GSE_km_s',
        'Proton_density_n_cc', 'Flow_pressure_nPa',
        'SYM_H_nT'
    ]
    
    try:
        import pandas as pd
        import numpy as np
        
        df = pd.read_fwf(omni_filename, colspecs=colspecs, names=column_names, header=None)
        
        # Replace fill values with NaN
        fill_values = [9999.9, 99.99, 999999, 9999999, 99999.0, 999.9, 99999.99,
                      999.99, 9999.99, 99999.9, -9999.9, -99.99]
        for fill_val in fill_values:
            df = df.replace(fill_val, np.nan)
        
        # Interpolate missing solar wind parameters between valid values
        solar_wind_columns = ['Bx_GSE_nT', 'By_GSE_nT', 'Bz_GSE_nT', 
                             'By_GSM_nT', 'Bz_GSM_nT',
                             'Vx_GSE_km_s', 'Vy_GSE_km_s', 'Vz_GSE_km_s',
                             'Proton_density_n_cc', 'Flow_pressure_nPa']
        
        # Count NaNs before interpolation
        nans_before = df[solar_wind_columns].isna().sum().sum()
        
        # Interpolate linearly between valid points (does not extrapolate)
        df[solar_wind_columns] = df[solar_wind_columns].interpolate(method='linear', limit_direction='both', limit_area='inside')
        
        # Count NaNs after interpolation
        nans_after = df[solar_wind_columns].isna().sum().sum()
        # print(f"  Interpolated {nans_before - nans_after} missing values")
        # print(f"  Remaining NaN values: {nans_after} (gaps at boundaries)")
        
        df['DateTime'] = pd.to_datetime(
            df['Year'].astype(int).astype(str) + '-' + 
            df['Day'].astype(int).astype(str) + '-' + 
            df['Hour'].astype(int).astype(str) + '-' + 
            df['Minute'].astype(int).astype(str),
            format='%Y-%j-%H-%M'
        )
        
        df['DateTime'] = df['DateTime'].dt.strftime('%Y-%m-%d %H:%M:%S')
        

        output_df = pd.DataFrame({
            'DateTime': df['DateTime'],
            'BX_nT_GSE_GSM': df['Bx_GSE_nT'],
            'BY_nT_GSE': df['By_GSE_nT'],
            'BZ_nT_GSE': df['Bz_GSE_nT'],
            'BY_nT_GSM': df['By_GSM_nT'],
            'BZ_nT_GSM': df['Bz_GSM_nT'],
            'Vx_Velocity_km_s': df['Vx_GSE_km_s'],
            'Vy_Velocity_km_s': df['Vy_GSE_km_s'],
            'Vz_Velocity_km_s': df['Vz_GSE_km_s'],
            'Proton_Density_n_cc': df['Proton_density_n_cc'],
            'Flow_pressure_nPa': df['Flow_pressure_nPa'],
            'SYM_H': df['SYM_H_nT']
        })
        
        # Compute corrected SYM-H using TSY15 methodology
        # print("  Computing corrected SYM-H using TSY15 methodology...")
        symh_corrected = compute_corrected_symh(
            output_df['SYM_H'].tolist(), 
            output_df['Flow_pressure_nPa'].tolist()
        )
        output_df['SYM_H'] = symh_corrected
        
        # print(f"  Added corrected SYM-H column")
        
        script_dir = os.path.join(os.path.dirname(__file__), "")
        output_path = os.path.join(script_dir, output_filename)
        output_df.to_csv(output_path, index=False)
        # print(f"Successfully converted {len(output_df)} rows to {output_path}")
        return True
        
    except ImportError as e:
        print(f"Error: Required package not available: {e}")
        print("Please install pandas and numpy: pip install pandas numpy")
        return False
    except Exception as e:
        print(f"Error converting OMNI data: {e}")
        import traceback
        traceback.print_exc()
        return False

def convert_tsy15_to_tempfile(tsy15_filename, output_filename='tempfile.csv'):
    # print(f"\nConverting {tsy15_filename} to {output_filename}...")
    
    # Column specifications for TSY15 fixed-width format
    colspecs = [
        (0, 4),    # IYEAR
        (4, 8),    # IDAY
        (8, 11),   # IHOUR
        (11, 14),  # MIN
        (14, 22),  # Bx IMF (GSW)
        (22, 30),  # By IMF (GSW)
        (30, 38),  # Bz IMF (GSW)
        (38, 46),  # Vx
        (46, 54),  # Vy
        (54, 62),  # Vz
        (62, 69),  # Np
        (78, 85),  # Sym-H
        (103, 110), # Pdyn
        (110, 118), # N-index
        (118, 126), # B-index
    ]
    
    column_names = [
        'IYEAR', 'IDAY', 'IHOUR', 'MIN',
        'Bx_IMF_nT', 'By_IMF_nT', 'Bz_IMF_nT',
        'Vx_km_s', 'Vy_km_s', 'Vz_km_s',
        'Np_cm3', 'SymH_nT', 'Pdyn_nPa',
        'N_index', 'B_index'
    ]
    
    try:
        import pandas as pd
        import numpy as np
        
        # Read the fixed-width file
        df = pd.read_fwf(tsy15_filename, colspecs=colspecs, names=column_names, header=None)
        
        # Replace fill values with NaN
        fill_values = [999999, 9999.99, 99.99]
        for fill_val in fill_values:
            df = df.replace(fill_val, np.nan)
        
        # Create DateTime column
        df['DateTime'] = pd.to_datetime(
            df['IYEAR'].astype(int).astype(str) + '-' + 
            df['IDAY'].astype(int).astype(str) + '-' + 
            df['IHOUR'].astype(int).astype(str) + '-' + 
            df['MIN'].astype(int).astype(str),
            format='%Y-%j-%H-%M'
        )
        
        # Format DateTime as string
        df['DateTime'] = df['DateTime'].dt.strftime('%Y-%m-%d %H:%M:%S')
        
        # Create output DataFrame with required columns including pre-computed indices
        output_df = pd.DataFrame({
            'DateTime': df['DateTime'],
            'BX_nT_GSE_GSM': df['Bx_IMF_nT'],
            'BY_nT_GSE': df['By_IMF_nT'],
            'BZ_nT_GSE': df['Bz_IMF_nT'],
            'BY_nT_GSM': df['By_IMF_nT'],
            'BZ_nT_GSM': df['Bz_IMF_nT'],
            'Vx_Velocity_km_s': df['Vx_km_s'],
            'Vy_Velocity_km_s': df['Vy_km_s'],
            'Vz_Velocity_km_s': df['Vz_km_s'],
            'Proton_Density_n_cc': df['Np_cm3'],
            'Flow_pressure_nPa': df['Pdyn_nPa'],
            'SYM_H': df['SymH_nT'],
            'N_index_normalised': df['N_index'],
            'B_index': df['B_index']
        })
        
        # Compute corrected SYM-H for TSY15 data (TSY15 provides raw SYM-H, needs correction)
        # print("  Computing corrected SYM-H for TSY15 data using TSY15 methodology...")
        symh_corrected = compute_corrected_symh(
            output_df['SYM_H'].tolist(), 
            output_df['Flow_pressure_nPa'].tolist()
        )
        output_df['SYM_H'] = symh_corrected
        
        # print(f"  Added corrected SYM-H column")
        
        # Save to CSV in OTSO package directory
        script_dir = os.path.join(os.path.dirname(__file__), "")
        output_path = os.path.join(script_dir, output_filename)
        output_df.to_csv(output_path, index=False)
        # print(f"Successfully converted {len(output_df)} rows to {output_path}")
        return True
        
    except ImportError as e:
        print(f"Error: Required package not available: {e}")
        print("Please install pandas and numpy: pip install pandas numpy")
        return False
    except Exception as e:
        print(f"Error converting TSY15 data: {e}")
        import traceback
        traceback.print_exc()
        return False

def check_and_download_tsy15_data(year):
    tsy15_url = f"https://geo.phys.spbu.ru/~tsyganenko/models/ta15/{year}_OMNI_5m_with_TA15_drivers.dat"
    
    # Save to OTSO package directory
    script_dir = os.path.join(os.path.dirname(__file__), "")
    local_filename = os.path.join(script_dir, f"{year}_OMNI_5m_with_TA15_drivers.dat")
    
    # print(f"Checking for TSY15 pre-computed data for year {year}...")
    # print(f"URL: {tsy15_url}")
    
    try:
        response = requests.get(tsy15_url, stream=True, verify=False, timeout=30)
        response.raise_for_status()
        
        with open(local_filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        # print(f"TSY15 data downloaded successfully: {local_filename}")
        return True, local_filename
            
    except requests.exceptions.RequestException as e:
        print(f"TSY15 data not available: {e}")
        return False, None
    except Exception as e:
        print(f"Error processing TSY15 data: {e}")
        return False, None

def download_omni_file(year):
    url = f"https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/omni_5min{year}.asc"
    
    # Save to OTSO package directory
    script_dir = os.path.join(os.path.dirname(__file__), "")
    local_filename = os.path.join(script_dir, f"omni_5min{year}.asc")
    
    try:
        total_size = 0
        try:
            response = requests.head(url, timeout=10)
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))
        except:
            print("Could not determine file size, showing download amount only...")
        
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        downloaded_size = 0
        last_progress_update = 0
        
        with open(local_filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    downloaded_size += len(chunk)
                    
                    # Update progress every 1MB to avoid too frequent updates
                    if downloaded_size - last_progress_update >= 1024 * 1024:  # 1MB
                        last_progress_update = downloaded_size
                        
                        # Display progress
                        if total_size > 0:
                            progress = (downloaded_size / total_size) * 100
                            bar_length = 50
                            filled_length = int(bar_length * downloaded_size // total_size)
                            bar = '█' * filled_length + '-' * (bar_length - filled_length)
                            
                            # Convert bytes to MB for display
                            downloaded_mb = downloaded_size / (1024 * 1024)
                            total_mb = total_size / (1024 * 1024)
                            
                            #print(f'\rProgress: |{bar}| {progress:.1f}% ({downloaded_mb:.1f}/{total_mb:.1f} MB)', end='', flush=True)
                        else:
                            # If we don't know total size, just show downloaded amount
                            downloaded_mb = downloaded_size / (1024 * 1024)
                            #print(f'\rDownloaded: {downloaded_mb:.1f} MB', end='', flush=True)
        
        #print()  # New line after progress bar
        
        # Final size check
        final_mb = downloaded_size / (1024 * 1024)
        #print(f"Downloaded successfully: {local_filename} ({final_mb:.1f} MB)")
        return True, local_filename
        
    except requests.exceptions.Timeout:
        print(f"Error: Download timed out")
        return False, None
    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return False, None
    except Exception as e:
        print(f"Unexpected error during download: {e}")
        return False, None

def check_and_download_data(year, check_tsy15=False):
    if year is None:
        #print("Error: No year specified")
        #print("Please specify a year to download data")
        return False, None, None
    
    #print("="*60)
    #print(f"DATA DOWNLOAD FOR YEAR {year}")
    #print("="*60)
    
    # Check TSY15 only if explicitly requested
    if check_tsy15:
        tsy15_success, tsy15_file = check_and_download_tsy15_data(year)
    
    # Use OMNI (default behavior)
    omni_success, omni_file = download_omni_file(year)
    
    if omni_success:
        # print(f"\nUsing OMNI 5-minute data for year {year}")
        return True, omni_file, 'omni'
    else:
        # print(f"Failed to download data for year {year}")
        return False, None, None

def process_omni_with_instantaneous_indices(g, h, year):
    # Get the directory of the current script (OTSO package location)
    script_dir = os.path.join(os.path.dirname(__file__), "")
    input_file = os.path.join(script_dir, 'tempfile.csv')
    output_file = os.path.join(script_dir, f'TSY15_{year}.csv')
    window_size = 7  # 30 minutes = 7 points for 5-minute resolution
    
    #print(f"\nProcessing tempfile.csv with CORRECTED methodology...")
    #print(f"Note: Solar wind parameters were interpolated during conversion")
    #print(f"Step 1: Computing instantaneous indices for each point")
    #print(f"Step 2: Applying {window_size}-point rolling average to indices")
    
    # Load all data
    all_rows = []
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            all_rows.append(row)
    
    #print(f"Loaded {len(all_rows)} rows of data")
    
    # Step 1: Compute instantaneous indices for each point
    #print("\nComputing instantaneous indices...")
    processed_rows = []
    
    for i, row in enumerate(all_rows):
        (N_norm, B_idx, bx_gsw, by_gsw, bz_gsw, bt_gsw, V, thetac_gsw, Vx, Vy, Vz, Np) = compute_instantaneous_indices(row, g, h)
        
        processed_row = row.copy()
        processed_row['BX_nT_GSW'] = f"{bx_gsw:.3f}" if bx_gsw is not None else ""
        processed_row['BY_nT_GSW'] = f"{by_gsw:.3f}" if by_gsw is not None else ""
        processed_row['BZ_nT_GSW'] = f"{bz_gsw:.3f}" if bz_gsw is not None else ""
        processed_row['Bt_nT_GSW'] = f"{bt_gsw:.3f}" if bt_gsw is not None else ""
        processed_row['V_km_s'] = f"{V:.2f}" if V is not None else ""
        processed_row['Vx_km_s'] = f"{Vx:.2f}" if Vx is not None else ""
        processed_row['Vy_km_s'] = f"{Vy:.2f}" if Vy is not None else ""
        processed_row['Vz_km_s'] = f"{Vz:.2f}" if Vz is not None else ""
        processed_row['Np_cm3'] = f"{Np:.3f}" if Np is not None else ""
        processed_row['Clock_angle_rad'] = f"{thetac_gsw:.4f}" if thetac_gsw is not None else ""
        processed_row['N_index_inst'] = f"{N_norm:.6f}" if N_norm is not None else ""
        processed_row['B_index_inst'] = f"{B_idx:.6f}" if B_idx is not None else ""
        
        processed_rows.append(processed_row)
        
        #if (i + 1) % 1000 == 0:
            #print(f"  Processed {i + 1} rows...")
    
    #print(f"Instantaneous indices computed for {len(processed_rows)} rows")
    
    # Step 2: Apply rolling average to the indices
    #print(f"\nApplying {window_size}-point rolling average to indices...")
    results = calculate_rolling_average_indices(processed_rows, window_size)
    
    # Write results
    with open(output_file, 'w', newline='') as outfile:
        fieldnames = ['DateTime', 'BX_nT_GSW', 'BY_nT_GSW', 'BZ_nT_GSW', 'Bt_nT_GSW',
                     'V_km_s', 'Vx_km_s', 'Vy_km_s', 'Vz_km_s', 'Np_cm3', 'Clock_angle_rad',
                     'Instantaneous_N_index', 'Instantaneous_B_index',
                     'N_index', 'B_index', 'SYM_H', 'Valid_points_in_window']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)
    
    #print(f"\nProcessing complete!")
    #print(f"Output saved as: {output_file}")
    
    # Show preview
    #print(f"\nPreview of results (first 10 rows with rolling indices):")
    #print(f"{'DateTime':<20} {'N_inst':>12} {'N_roll':>12} {'B_inst':>12} {'B_roll':>12}")
    #print(f"{'-'*70}")
    
    count = 0
    for row in results:
        if row['N_index'] and row['B_index']:
            #print(f"{row['DateTime']:<20} {row['Instantaneous_N_index']:>12} {row['N_index']:>12} "
            #      f"{row['Instantaneous_B_index']:>12} {row['B_index']:>12}")
            count += 1
            if count >= 10:
                break

def process_year(year, g, h, check_tsy15=False):
    # Validate year
    current_year = datetime.now().year
    #if year < 1995 or year > current_year:
    #    print(f"Error: Year must be between 1995 and {current_year}")
    #    return False, None, None
    
    # Download data
    success, filename, data_source = check_and_download_data(year, check_tsy15)
    
    if not success:
        print("Failed to download data file.")
        return False, None, None
    
    # Process data based on source
    if data_source == 'tsy15':
        #print("\n" + "="*60)
        #print("DOWNLOADED TSY15 PRE-COMPUTED DATA")
        #print("="*60)
        #print(f"File: {filename}")
        #print("Note: TSY15 data already contains computed N and B indices")
        #print("Converting to CSV format...")
        
        conversion_success = convert_tsy15_to_tempfile(filename)
        
        
        if not conversion_success:
            #print("\nCONVERSION FAILED")
            return False, None, data_source
            
        script_dir = os.path.join(os.path.dirname(__file__), "")
        final_output_file = os.path.join(script_dir, f'TSY15_{year}.csv')
        
        #print("\n" + "="*60)
        #print("PROCESSING COMPLETED SUCCESSFULLY!")
        #print("="*60)
        #print(f"Output file: {final_output_file}")
        return True, final_output_file, data_source
        
    elif data_source == 'omni':
        #print("\n" + "="*60)
        #print("DOWNLOADED OMNI 5-MINUTE DATA")
        #print("="*60)
        #print(f"File: {filename}")
        
        conversion_success = convert_omni_to_tempfile(filename)
        
        if not conversion_success:
            #print("\nCONVERSION FAILED")
            return False, None, data_source
        
        # Check if output file was created in OTSO directory
        script_dir = os.path.join(os.path.dirname(__file__), "")
        expected_tempfile = os.path.join(script_dir, 'tempfile.csv')
        if not os.path.exists(expected_tempfile):
            print(f"\ntempfile.csv was not created at {expected_tempfile}")
            return False, None, data_source
        
        #print("\n" + "="*60)
        #print("COMPUTING INDICES WITH CORRECTED METHODOLOGY")
        #print("="*60)
        
        try:
            final_output_file = process_omni_with_instantaneous_indices(g, h, year)
            #print("\n" + "="*60)
            #print("PROCESSING COMPLETED SUCCESSFULLY!")
            #print("="*60)
            #print(f"Output file: {final_output_file}")
            return True, final_output_file, data_source
        except FileNotFoundError:
            script_dir = os.path.join(os.path.dirname(__file__), "")
            #print(f"Error: tempfile.csv not found during processing at {script_dir}")
            return False, None, data_source
        except Exception as e:
            #print(f"An error occurred during processing: {e}")
            import traceback
            traceback.print_exc()
            return False, None, data_source
    
    return False, None, None

def TSY15_OTSO_Download_v2(OMNI_year,g,h,check_tsy15=False):
    success, output_file, data_source = process_year(OMNI_year, g, h, check_tsy15)
    
    if not success:
        sys.exit(1)
    
    sys.exit(0)
