import pandas as pd
import numpy as np
import os
from typing import Tuple

def decimal_year_from_date_array(date_array) -> float:
    year, day_of_year, hour, minute, secs, seconds = date_array
    
    # Total seconds in the day
    total_seconds = hour * 3600 + minute * 60 + secs + seconds
    
    # Fraction of the day
    day_fraction = total_seconds / 86400.0  # 86400 seconds in a day
    
    # Adjust day of year by the fraction
    adjusted_day = day_of_year + day_fraction
    
    # Convert to decimal year (assuming 365.25 days per year)
    decimal_year = year + (adjusted_day - 1) / 365.25
    
    return decimal_year

def load_igrf_data(csv_file: str) -> pd.DataFrame:
    script_dir = os.path.dirname(os.path.realpath(__file__))
    subfolder_name = ''
    csv_file_name = csv_file
    
    csv_file_path = os.path.join(script_dir, subfolder_name, csv_file_name)
    return pd.read_csv(csv_file_path)

def get_coefficient_value(df: pd.DataFrame, coeff_type: str, degree: int, order: int, year_col: str) -> float:
    mask = (df['coeff'] == coeff_type) & (df['SH_degree'] == degree) & (df['SH_order'] == order)
    rows = df[mask]
    
    if len(rows) == 0:
        return 0.0
    
    if year_col not in rows.columns:
        return 0.0
    
    value = rows[year_col].iloc[0]
    return float(value) if not pd.isna(value) else 0.0

def geopack_index(degree: int, order: int) -> int:
    if degree == 0:
        return 1  # g(0,0) or h(0,0) goes to index 1
    
    # Count coefficients before this (degree, order)
    count = 1  # Start at 1 to skip index 1
    
    for n in range(1, degree + 1):
        for m in range(0, n + 1):
            count += 1
            if n == degree and m == order:
                return count
    
    return count

def create_coefficient_arrays(df: pd.DataFrame, year_col: str) -> Tuple[np.ndarray, np.ndarray]:
    G = np.zeros(137)  # Index 0-136 for 15th order
    H = np.zeros(137)
    
    for _, row in df.iterrows():
        coeff_type = row['coeff']
        degree = int(row['SH_degree'])
        order = int(row['SH_order'])
        if year_col not in df.columns:
            continue
        value = row[year_col]
        if pd.isna(value):
            value = 0.0
        else:
            value = float(value)
        index = geopack_index(degree, order)
        if 1 <= index <= 136:
            if coeff_type == 'g':
                G[index] = value
            elif coeff_type == 'h':
                H[index] = value
    
    return G, H

def interpolate_coefficients(df: pd.DataFrame, target_year: float) -> Tuple[np.ndarray, np.ndarray, float]:
    year_cols = [col for col in df.columns if col not in ['coeff', 'SH_degree', 'SH_order']]
    
    regular_years = []
    secular_variation_col = None
    last_definitive_year = None
    
    for col in year_cols:
        if '+' in col:
            base_year = float(col.replace('+', ''))
            secular_variation_col = col
            last_definitive_year = base_year
        else:
            try:
                year_val = float(col)
                regular_years.append(year_val)
            except ValueError:
                continue
    
    regular_years.sort()
    
    if last_definitive_year is None and regular_years:
        last_definitive_year = regular_years[-1]
    
    if last_definitive_year is not None and target_year > last_definitive_year + 5.0:
        max_extrapolation_year = last_definitive_year + 5.0
        print(f"\n WARNING: IGRF extrapolation limit exceeded!")
        print(f"   Your current version of the IGRF model can only extrapolate to {last_definitive_year:.0f} + 5 years = {max_extrapolation_year:.0f}")
        print(f"   You entered a date for year {target_year:.2f}, which is beyond this limit.")
        print(f"   OTSO will use the furthest available extrapolation ({max_extrapolation_year:.0f}).")
        print(f"   We recommend updating your IGRF model if a newer version is available")
        print(f"   by running: OTSO.IGRFupdate in the terminal")
        print()
        
        target_year = max_extrapolation_year
    
    if target_year in regular_years:
        epoch_col = str(int(target_year))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    
    if target_year > last_definitive_year and secular_variation_col is not None:
        G, H = extrapolate_beyond_last_epoch(df, target_year, last_definitive_year, secular_variation_col)
        return G, H, target_year
    
    if target_year <= regular_years[0]:
        epoch_col = str(int(regular_years[0]))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    elif target_year >= regular_years[-1]:
        epoch_col = str(int(regular_years[-1]))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    else:
        for i in range(len(regular_years) - 1):
            if regular_years[i] <= target_year <= regular_years[i + 1]:
                year1, year2 = regular_years[i], regular_years[i + 1]
                col1, col2 = str(int(year1)), str(int(year2))
                
                G1, H1 = create_coefficient_arrays(df, col1)
                G2, H2 = create_coefficient_arrays(df, col2)
                
                t = (target_year - year1) / (year2 - year1)
                G_interp = G1 + t * (G2 - G1)
                H_interp = H1 + t * (H2 - H1)
                
                return G_interp, H_interp, target_year
    
    G, H = create_coefficient_arrays(df, year_cols[0])
    return G, H, target_year

def extrapolate_beyond_last_epoch(df: pd.DataFrame, target_year: float, last_year: float, 
                                 sv_col: str) -> Tuple[np.ndarray, np.ndarray]:

    dt = target_year - last_year
    
    last_year_col = str(int(last_year))
    G_base, H_base = create_coefficient_arrays(df, last_year_col)
    
    G_sv = np.zeros(137)
    H_sv = np.zeros(137)
    
    for _, row in df.iterrows():
        coeff_type = row['coeff']
        degree = int(row['SH_degree'])
        order = int(row['SH_order'])
        
        if sv_col not in df.columns:
            continue
            
        sv_value = row[sv_col]
        if pd.isna(sv_value):
            sv_value = 0.0
        else:
            sv_value = float(sv_value)
        
        index = geopack_index(degree, order)
        
    if 1 <= index <= 136:
            if coeff_type == 'g':
                G_sv[index] = sv_value
            elif coeff_type == 'h':
                H_sv[index] = sv_value
    
    G_extrap = G_base.copy()
    H_extrap = H_base.copy()
    for i in range(1, min(46, 137)):
        G_extrap[i] = G_base[i] + G_sv[i] * dt
        H_extrap[i] = H_base[i] + H_sv[i] * dt
    for i in range(46, 137):
        G_extrap[i] = G_base[i]
        H_extrap[i] = H_base[i]
    return G_extrap, H_extrap

def compute_gauss_coefficients(date_array):

    df = load_igrf_data('igrf_coefficients.csv')
    
    original_target_year = decimal_year_from_date_array(date_array)
    
    G, H, actual_year_used = interpolate_coefficients(df, original_target_year)

    G, H = schmidt_normalize(G, H, max_n=15)
    
    G_10 = abs(G[2])  # g(1,0) magnitude
    G_11 = G[3]       # g(1,1)
    H_11 = H[3]       # h(1,1)

    dipole_moment = np.sqrt(G[2]**2 + G[3]**2 + H[3]**2)

    dipole_latitude = np.degrees(np.arctan2(G[2], np.sqrt(G[3]**2 + H[3]**2)))
    dipole_longitude = np.degrees(np.arctan2(H[3], G[3]))

    return {
        'G_10': G_10,
        'G_coefficients': G[1:137],  # Return indices 1-136 (15th order)
        'H_coefficients': H[1:137],  # Return indices 1-136 (15th order)
        'dipole_latitude': dipole_latitude,
        'dipole_longitude': dipole_longitude,
        'dipole_moment': dipole_moment,
        'decimal_year': actual_year_used
    }


def schmidt_normalize(G, H, max_n=15):
    S = 1.0
    for N in range(2, max_n + 1):
        MN = N * (N - 1) // 2 + 1
        S *= (2 * N - 3) / (N - 1)
        G[MN] *= S
        H[MN] *= S
        P = S
        for M in range(2, N + 1):
            AA = 2.0 if M == 2 else 1.0
            P *= np.sqrt(AA * (N - M + 1) / (N + M - 2))
            MNN = MN + M - 1
            G[MNN] *= P
            H[MNN] *= P
    return G, H