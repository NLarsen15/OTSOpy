#!/usr/bin/env python3
"""
Simple IGRF coefficient interpolation for GEOPACK compatibility.
Loads IGRF coefficients from CSV and interpolates for a given date.
"""

import pandas as pd
import numpy as np
import os
from typing import Tuple

def decimal_year_from_date_array(date_array) -> float:
    """Convert DateArray [year, day, hour, minute, secs, seconds] to decimal year."""
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
    """Load IGRF coefficients from CSV file."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    subfolder_name = ''
    csv_file_name = csv_file
    
    csv_file_path = os.path.join(script_dir, subfolder_name, csv_file_name)
    return pd.read_csv(csv_file_path)

def get_coefficient_value(df: pd.DataFrame, coeff_type: str, degree: int, order: int, year_col: str) -> float:
    """Get a specific coefficient value from the dataframe."""
    mask = (df['coeff'] == coeff_type) & (df['SH_degree'] == degree) & (df['SH_order'] == order)
    rows = df[mask]
    
    if len(rows) == 0:
        return 0.0
    
    if year_col not in rows.columns:
        return 0.0
    
    value = rows[year_col].iloc[0]
    return float(value) if not pd.isna(value) else 0.0

def geopack_index(degree: int, order: int) -> int:
    """
    Calculate GEOPACK array index for coefficients.
    Based on sequential assignment: (1,0)→2, (1,1)→3, (2,0)→4, (2,1)→5, (2,2)→6, (3,0)→7, etc.
    """
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
    """Create G and H coefficient arrays for GEOPACK (indices 1-105)."""
    G = np.zeros(137)  # Index 0-136, but we use 1-136 for 15th order
    H = np.zeros(137)
    
    # Process all coefficients in the dataframe
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
        # Calculate GEOPACK index
        index = geopack_index(degree, order)
        if 1 <= index <= 136:
            if coeff_type == 'g':
                G[index] = value
            elif coeff_type == 'h':
                H[index] = value
    
    return G, H

def interpolate_coefficients(df: pd.DataFrame, target_year: float) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    Interpolate IGRF coefficients for a target year.
    Handles extrapolation beyond the last epoch using secular variation data.
    
    Returns:
        Tuple of (G_array, H_array, actual_year_used)
    """
    # Find available year columns (skip metadata columns)
    year_cols = [col for col in df.columns if col not in ['coeff', 'SH_degree', 'SH_order']]
    
    # Separate regular years from prediction columns and find secular variation column
    regular_years = []
    secular_variation_col = None
    last_definitive_year = None
    
    for col in year_cols:
        if '+' in col:
            # Found secular variation column - extract base year
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
    
    # If no secular variation column found, use the last regular year
    if last_definitive_year is None and regular_years:
        last_definitive_year = regular_years[-1]
    
    # Check for extrapolation limits (IGRF model valid for last_year + 5 years)
    if last_definitive_year is not None and target_year > last_definitive_year + 5.0:
        max_extrapolation_year = last_definitive_year + 5.0
        print(f"\n WARNING: IGRF extrapolation limit exceeded!")
        print(f"   Your current version of the IGRF model can only extrapolate to {last_definitive_year:.0f} + 5 years = {max_extrapolation_year:.0f}")
        print(f"   You entered a date for year {target_year:.2f}, which is beyond this limit.")
        print(f"   OTSO will use the furthest available extrapolation ({max_extrapolation_year:.0f}).")
        print(f"   We recommend updating your IGRF model if a newer version is available")
        print(f"   by running: OTSO.IGRFupdate in the terminal")
        print()
        
        # Cap the target year to the maximum safe extrapolation
        target_year = max_extrapolation_year
    
    # Check for exact match first
    if target_year in regular_years:
        epoch_col = str(int(target_year))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    
    # Handle extrapolation beyond the last year using secular variation
    if target_year > last_definitive_year and secular_variation_col is not None:
        G, H = extrapolate_beyond_last_epoch(df, target_year, last_definitive_year, secular_variation_col)
        return G, H, target_year
    
    # Find the two epochs to interpolate between
    if target_year <= regular_years[0]:
        # Use first epoch
        epoch_col = str(int(regular_years[0]))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    elif target_year >= regular_years[-1]:
        # Use last epoch if no extrapolation available
        epoch_col = str(int(regular_years[-1]))
        G, H = create_coefficient_arrays(df, epoch_col)
        return G, H, target_year
    else:
        # Interpolate between two epochs
        for i in range(len(regular_years) - 1):
            if regular_years[i] <= target_year <= regular_years[i + 1]:
                year1, year2 = regular_years[i], regular_years[i + 1]
                col1, col2 = str(int(year1)), str(int(year2))
                
                # Get coefficients for both epochs
                G1, H1 = create_coefficient_arrays(df, col1)
                G2, H2 = create_coefficient_arrays(df, col2)
                
                # Linear interpolation
                t = (target_year - year1) / (year2 - year1)
                G_interp = G1 + t * (G2 - G1)
                H_interp = H1 + t * (H2 - H1)
                
                return G_interp, H_interp, target_year
    
    # Fallback
    G, H = create_coefficient_arrays(df, year_cols[0])
    return G, H, target_year

def extrapolate_beyond_last_epoch(df: pd.DataFrame, target_year: float, last_year: float, 
                                 sv_col: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extrapolate coefficients beyond the last epoch using secular variation data.
    Similar to GEOPACK's extrapolation beyond 2025.
    
    Args:
        df: IGRF DataFrame
        target_year: Target year for extrapolation
        last_year: Last regular epoch year (e.g., 2025)
        sv_col: Secular variation column name (e.g., '2025+')
        
    Returns:
        Tuple of (G_array, H_array) extrapolated coefficients
    """
    # Calculate time difference in years
    dt = target_year - last_year
    
    # Get base coefficients from the last epoch
    last_year_col = str(int(last_year))
    G_base, H_base = create_coefficient_arrays(df, last_year_col)
    
    # Initialize arrays for secular variations
    G_sv = np.zeros(137)
    H_sv = np.zeros(137)
    
    # Extract secular variation coefficients
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
        
        # Calculate GEOPACK index
        index = geopack_index(degree, order)
        
    if 1 <= index <= 136:
            if coeff_type == 'g':
                G_sv[index] = sv_value
            elif coeff_type == 'h':
                H_sv[index] = sv_value
    
    # Apply extrapolation: coeff(target_year) = coeff(last_year) + sv * dt
    # Following GEOPACK, limit secular variation to first 45 coefficients
    G_extrap = G_base.copy()
    H_extrap = H_base.copy()
    # GEOPACK limits to first 45 coefficients for secular variation, keep that logic
    for i in range(1, min(46, 137)):
        G_extrap[i] = G_base[i] + G_sv[i] * dt
        H_extrap[i] = H_base[i] + H_sv[i] * dt
    # For coefficients beyond index 45, use base values without extrapolation
    for i in range(46, 137):
        G_extrap[i] = G_base[i]
        H_extrap[i] = H_base[i]
    return G_extrap, H_extrap

def compute_gauss_coefficients(date_array):
    """
    Compute IGRF Gauss coefficients for a given date array.
    
    Args:
        date_array: [year, day, hour, minute, secs, seconds]
        
    Returns:
        Dictionary with G and H coefficients and derived values
    """
    # Load IGRF data
    df = load_igrf_data('igrf_coefficients.csv')
    
    # Convert date array to decimal year
    original_target_year = decimal_year_from_date_array(date_array)
    
    # Get interpolated coefficients (returns actual year used)
    G, H, actual_year_used = interpolate_coefficients(df, original_target_year)
    
    # Calculate derived values
    G_10 = abs(G[2])  # g(1,0) magnitude
    G_11 = G[3]       # g(1,1)
    H_11 = H[3]       # h(1,1)
    # Calculate dipole parameters
    dipole_moment = np.sqrt(G[2]**2 + G[3]**2 + H[3]**2)
    # Dipole latitude and longitude (simplified calculation)
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
