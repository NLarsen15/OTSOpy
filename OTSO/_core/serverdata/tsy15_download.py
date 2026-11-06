import csv
import math
import requests
import urllib3
import os
from datetime import datetime
import numpy as np
import pandas as pd

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils.tsy_params_utils import N_index_normalized, B_index
from .scratch_dir import SCRATCH_DIR
from .prior_year import load_prior_year_tail

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

        position_gsw = np.zeros(3, dtype=np.float64)

        OTSOLib.gse2gswtsy15(date, position_gse, Wind, g, h, len(g), position_gsw)
        
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
    """Mean of the valid values in [i - w//2, i + w//2] (truncated at the ends; NaN where a
    window has none). Vectorised: NaN contributes 0 to the sum and not to the count, which
    adds the same values in the same order as averaging the valid values row by row."""
    import numpy as np
    from numpy.lib.stride_tricks import sliding_window_view
    values = np.asarray(values, dtype=float)
    half = window_size // 2
    padded = np.concatenate([np.full(half, np.nan), values, np.full(half, np.nan)])
    windows = sliding_window_view(padded, 2 * half + 1)
    valid = ~np.isnan(windows)
    count = valid.sum(axis=1)
    total = np.where(valid, windows, 0.0).sum(axis=1)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(count > 0, total / np.maximum(count, 1), np.nan)

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
        
        script_dir = SCRATCH_DIR
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
        script_dir = SCRATCH_DIR
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
    script_dir = SCRATCH_DIR
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
    
    script_dir = SCRATCH_DIR
    local_filename = os.path.join(script_dir, f"omni_5min{year}.asc")

    try:
        total_size = 0
        try:
            response = requests.head(url, timeout=10)
            response.raise_for_status()
            total_size = int(response.headers.get('content-length', 0))
        except requests.exceptions.RequestException:
            #print("Could not determine file size, showing download amount only...")
            pass

        response = requests.get(url, stream=True, timeout=30)
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

    # Prefer Tsyganenko's own pre-computed driver file when requested and available -
    # it already contains computed N/B indices, so it skips the OMNI recompute below.
    if check_tsy15:
        tsy15_success, tsy15_file = check_and_download_tsy15_data(year)
        if tsy15_success:
            return True, tsy15_file, 'tsy15'

    # Use OMNI (default behavior, and the fallback if TSY15 wasn't requested/available)
    omni_success, omni_file = download_omni_file(year)

    if omni_success:
        # print(f"\nUsing OMNI 5-minute data for year {year}")
        return True, omni_file, 'omni'
    else:
        # print(f"Failed to download data for year {year}")
        return False, None, None

def build_warmup_rows(year, window_size):
    """
    Builds synthetic 'instantaneous'-index rows from the tail of the previous year's
    cached data, so the rolling average has a real trailing window for the first points
    of the year instead of starting cold.

    This is an approximation: raw OMNI data for prior years isn't kept once a year
    finishes downloading (see prior_year.load_prior_year_tail), only that year's own
    final, already-rolling-averaged N_index/B_index. Using those directly as the
    instantaneous seed values isn't exact, but it's far closer to reality than the
    zero-history cold start it replaces, and it's the only history actually available
    without re-downloading and reprocessing the previous year's raw OMNI feed.
    """
    prior_tail = load_prior_year_tail(year, hours=1)
    if prior_tail is None or prior_tail.empty:
        return []
    tail = prior_tail.tail(window_size - 1)
    warmup_rows = []
    for _, r in tail.iterrows():
        warmup_rows.append({
            'DateTime': r['Date'].strftime('%Y-%m-%d %H:%M:%S'),
            'N_index_inst': f"{r['N_index']:.6f}" if pd.notna(r.get('N_index')) else '',
            'B_index_inst': f"{r['B_index']:.6f}" if pd.notna(r.get('B_index')) else '',
            'SYM_H': f"{r['SYM_H']:.3f}" if pd.notna(r.get('SYM_H')) else '',
        })
    return warmup_rows

_GSE_COLUMNS = ['BX_nT_GSE_GSM', 'BY_nT_GSE', 'BZ_nT_GSE',
                'Vx_Velocity_km_s', 'Vy_Velocity_km_s', 'Vz_Velocity_km_s', 'Proton_Density_n_cc']


def _to_float(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return np.nan


def _fmt(values, spec):
    return [format(v, spec) if v == v else "" for v in values]


def process_omni_with_instantaneous_indices(g, h, year):
    """Instantaneous N-index / B-index for every valid row (IMF rotated to GSW with the
    Fortran gse2gswTSY15), then their 7-point trailing mean (at least 3 valid points),
    warm-started from the previous year's tail. Writes TSY15_{year}.csv.

    Same rules and output as the former row-by-row version (see compute_instantaneous_indices
    and calculate_rolling_average_indices), with the bookkeeping vectorised."""
    script_dir = SCRATCH_DIR
    input_file = os.path.join(script_dir, 'tempfile.csv')
    output_file = os.path.join(script_dir, f'TSY15_{year}.csv')
    window_size = 7  # 30 minutes = 7 points for 5-minute resolution

    text = pd.read_csv(input_file, dtype=str, keep_default_na=False)
    n = len(text)

    invalid_markers = {9999.99, 99999.9, 99.99, 999.99, -9999.99, -99999.9, -99.99, -999.99}
    numeric = {}
    valid = np.ones(n, dtype=bool)
    for col in _GSE_COLUMNS:
        values = np.array([_to_float(v) for v in text[col]], dtype=float)  # float(), as before
        numeric[col] = values
        valid &= np.isfinite(values) & ~np.isin(values, list(invalid_markers))

    bx = np.full(n, np.nan); by = np.full(n, np.nan); bz = np.full(n, np.nan); bt = np.full(n, np.nan)
    speed = np.full(n, np.nan); theta = np.full(n, np.nan)
    n_inst = np.full(n, np.nan); b_inst = np.full(n, np.nan)
    vx, vy, vz = numeric['Vx_Velocity_km_s'], numeric['Vy_Velocity_km_s'], numeric['Vz_Velocity_km_s']
    np_cm3 = numeric['Proton_Density_n_cc']
    dates = text['DateTime'].to_numpy()

    for i in np.flatnonzero(valid):
        x, y, z, t = convert_to_gsw_coordinates(dates[i], vx[i], vy[i], vz[i],
                                                numeric['BX_nT_GSE_GSM'][i], numeric['BY_nT_GSE'][i],
                                                numeric['BZ_nT_GSE'][i], g, h)
        if x is None:
            valid[i] = False
            continue
        V = math.sqrt(vx[i]**2 + vy[i]**2 + vz[i]**2)
        thetac = math.atan2(y, z)
        if thetac < 0:
            thetac += 2 * math.pi
        bx[i], by[i], bz[i], bt[i], speed[i], theta[i] = x, y, z, t, V, thetac
        n_inst[i] = N_index_normalized(V, t, thetac)
        b_inst[i] = B_index(np_cm3[i], V, t, thetac)

    def keep_valid(values):
        return np.where(valid, values, np.nan)

    out = pd.DataFrame({
        'DateTime': dates,
        'BX_nT_GSW': _fmt(bx, '.3f'), 'BY_nT_GSW': _fmt(by, '.3f'), 'BZ_nT_GSW': _fmt(bz, '.3f'),
        'Bt_nT_GSW': _fmt(bt, '.3f'), 'V_km_s': _fmt(speed, '.2f'),
        'Vx_km_s': _fmt(keep_valid(vx), '.2f'), 'Vy_km_s': _fmt(keep_valid(vy), '.2f'),
        'Vz_km_s': _fmt(keep_valid(vz), '.2f'), 'Np_cm3': _fmt(keep_valid(np_cm3), '.3f'),
        'Clock_angle_rad': _fmt(theta, '.4f'),
        'Instantaneous_N_index': _fmt(n_inst, '.6f'), 'Instantaneous_B_index': _fmt(b_inst, '.6f'),
    })

    # The rolling mean works on the values as written (6 decimals), as before
    n_written = np.array([float(v) if v else np.nan for v in out['Instantaneous_N_index']])
    b_written = np.array([float(v) if v else np.nan for v in out['Instantaneous_B_index']])
    warmup_rows = build_warmup_rows(year, window_size)
    n_warm = len(warmup_rows)
    if n_warm:
        n_written = np.concatenate([[float(r['N_index_inst']) if r['N_index_inst'] else np.nan
                                     for r in warmup_rows], n_written])
        b_written = np.concatenate([[float(r['B_index_inst']) if r['B_index_inst'] else np.nan
                                     for r in warmup_rows], b_written])

    def trailing(values):
        from numpy.lib.stride_tricks import sliding_window_view
        m = len(values)
        mean = np.full(m, np.nan)
        count = np.zeros(m, dtype=int)
        if m >= window_size:
            windows = sliding_window_view(values, window_size)
            ok = ~np.isnan(windows)
            count[window_size - 1:] = ok.sum(axis=1)
            total = np.where(ok, windows, 0.0).sum(axis=1)
            with np.errstate(invalid='ignore', divide='ignore'):
                mean[window_size - 1:] = np.where(count[window_size - 1:] >= 3,
                                                  total / np.maximum(count[window_size - 1:], 1), np.nan)
        return mean[n_warm:], count[n_warm:]

    n_mean, n_count = trailing(n_written)
    b_mean, b_count = trailing(b_written)
    full_window = (np.arange(n) + n_warm) >= window_size - 1
    out['N_index'] = _fmt(n_mean, '.6f')
    out['B_index'] = _fmt(b_mean, '.6f')
    out['SYM_H'] = text['SYM_H'].to_numpy() if 'SYM_H' in text.columns else ""
    out['Valid_points_in_window'] = np.where(full_window, np.minimum(n_count, b_count).astype(str), "")

    out.to_csv(output_file, index=False)

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
            
        script_dir = SCRATCH_DIR
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
        script_dir = SCRATCH_DIR
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
            script_dir = SCRATCH_DIR
            #print(f"Error: tempfile.csv not found during processing at {script_dir}")
            return False, None, data_source
        except Exception as e:
            #print(f"An error occurred during processing: {e}")
            import traceback
            traceback.print_exc()
            return False, None, data_source
    
    return False, None, None
