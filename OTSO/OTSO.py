import sys, os, csv
import subprocess
import platform
import shutil
import pandas as pd
import urllib.request
from typing import List, Tuple

def clean():
    print("Cleaning OTSO...")
    Delete()
    print("OTSO cleaned")


def Delete():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    setupfile = os.path.join(script_dir, 'setup_complete.txt')

    if os.path.exists(setupfile):
        os.remove(setupfile)

    server_data_folder_path = os.path.join(script_dir, 'Parameters', 'functions', 'ServerData')
    if os.path.exists(server_data_folder_path):
        shutil.rmtree(server_data_folder_path)

def AddStation(Name, Latitude, Longitude):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    existing_row = df[df["Name"] == Name]

    if not existing_row.empty:
        print(f"Station '{Name}' already exists:")
        print(existing_row)

        user_input = input("Do you want to overwrite this entry? (y/n): ").strip().lower()
        
        if user_input != "y":
            print("No changes were made.")
            return

        df = df[df["Name"] != Name]

    new_station = pd.DataFrame([[Name, Latitude, Longitude]], columns=["Name", "Latitude", "Longitude"])
    df = pd.concat([df, new_station], ignore_index=True)

    df = df.sort_values(by="Name")

    df.to_csv(stationfile, index=False)

    print(f"Station '{Name}' has been added/updated successfully.")


def RemoveStation(Name):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    existing_row = df[df["Name"] == Name]

    if existing_row.empty:
        print(f"Station '{Name}' not found. No changes were made.")
        return

    print(f"Station '{Name}' found:")
    print(existing_row)

    user_input = input("Do you want to delete this entry? (y/n): ").strip().lower()
    
    if user_input != "y":
        print("No changes were made.")
        return

    df = df[df["Name"] != Name]

    df.to_csv(stationfile, index=False)

    print(f"Station '{Name}' has been removed successfully.")


def ListStations():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    stationfile = os.path.join(script_dir, 'Parameters', 'functions', 'StationList.csv')

    df = pd.read_csv(stationfile)

    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)

    print(df.to_string(index=False))


def find_latest_igrf_url(base_url: str = "https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/") -> str:
    latest_valid_generation = None
    latest_valid_url = None
    
    # Start from IGRF-12 and count up until we find an invalid one
    generation = 12
    while True:
        test_url = f"{base_url}igrf{generation}coeffs.txt"
        try:
            # Try to open the URL and read just the first few bytes to check if it exists
            with urllib.request.urlopen(test_url) as response:
                # If we can read the first 100 characters without error, the file exists
                first_bytes = response.read(100)
                if b'Generation International Geomagnetic Reference Field' in first_bytes or b'IGRF' in first_bytes:
                    latest_valid_generation = generation
                    latest_valid_url = test_url
                    generation += 1  # Check the next generation
                else:
                    # File exists but doesn't seem to be an IGRF file
                    break
        except Exception:
            # File doesn't exist or can't be accessed, we've found the latest
            break
    
    if latest_valid_url:
        return latest_valid_url, latest_valid_generation
    else:
        # Fallback to IGRF-14 if we can't find any valid version
        fallback_url = f"{base_url}igrf14coeffs.txt"
        return fallback_url, 14


def download_igrf_data(url: str = None) -> List[str]:
    if url is None:
        url, _ = find_latest_igrf_url()
    
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')
        return data.strip().split('\n')
    except Exception as e:
        print(f"Error downloading IGRF data from {url}: {e}")
        sys.exit(1)


def parse_igrf_header(lines: List[str]) -> Tuple[List[str], int]:
    header_line = None
    data_start_line = 0
    
    for i, line in enumerate(lines):
        # Skip comment lines that start with #
        if line.strip().startswith('#'):
            continue
        # Look for header line with g/h n m (with flexible spacing)
        stripped_line = line.strip()
        if stripped_line.startswith('g/h') and 'n' in stripped_line and 'm' in stripped_line:
            # Check if this looks like a header line (has numbers after g/h n m)
            parts = stripped_line.split()
            if len(parts) > 3:
                header_line = stripped_line
                data_start_line = i + 1
                break
    
    if header_line is None:
        raise ValueError("Could not find header line in IGRF data")
    
    # Extract years from header
    parts = header_line.split()
    years = []
    # Find where the year data starts (after g/h, n, m)
    start_index = 3
    for j, part in enumerate(parts):
        if part in ['g/h', 'n', 'm']:
            continue
        if j >= start_index:
            try:
                # Handle regular years that might have decimal places (like 2025.0)
                year_str = part.strip()
                if '.' in year_str:
                    year = int(float(year_str))
                    years.append(str(year))
                elif '-' in year_str:
                    # Handle prediction columns like "2025-30"
                    # Use the first year and add a '+' suffix
                    start_year = year_str.split('-')[0]
                    year = int(start_year)
                    years.append(f"{year}+")
                else:
                    year = int(year_str)
                    years.append(str(year))
            except ValueError:
                # Skip non-numeric parts (like 'SV' for secular variation)
                continue
    
    return years, data_start_line


def parse_igrf_data_all_years(lines: List[str], start_line: int, available_years: List[str]) -> List[Tuple[str, int, int, List[float]]]:
    coefficients = []
    zero_values = [0.0] * len(available_years)
    
    # Keep track of which (degree, order) combinations we've seen for each coefficient type
    seen_coeffs = set()
    
    # Parse the IGRF data first
    for line in lines[start_line:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
            
        parts = line.split()
        if len(parts) < 4:
            continue
            
        coeff_type = parts[0]  # 'g' or 'h'
        degree = int(parts[1])
        order = int(parts[2])
        
        # Extract values for all years
        values = []
        for year_index in range(len(available_years)):
            value_index = 3 + year_index
            
            if value_index < len(parts):
                try:
                    value = float(parts[value_index])
                    values.append(value)
                except ValueError:
                    values.append(0.0)  # Default to 0.0 for missing values
            else:
                values.append(0.0)  # Default to 0.0 for missing values
        
        coefficients.append((coeff_type, degree, order, values))
        seen_coeffs.add((coeff_type, degree, order))
    
    # Now add missing coefficients that should be zero
    # Find the maximum degree in the data
    max_degree = max(degree for _, degree, _, _ in coefficients) if coefficients else 0
    
    missing_coeffs = []
    
    # Add g,0,0 and h,0,0 (always zero)
    if ('g', 0, 0) not in seen_coeffs:
        missing_coeffs.append(('g', 0, 0, zero_values.copy()))
    if ('h', 0, 0) not in seen_coeffs:
        missing_coeffs.append(('h', 0, 0, zero_values.copy()))
    
    # Add missing h coefficients with order 0 for all degrees (these are always zero by spherical harmonics convention)
    for degree in range(1, max_degree + 1):
        if ('h', degree, 0) not in seen_coeffs:
            missing_coeffs.append(('h', degree, 0, zero_values.copy()))
    
    # Add the missing coefficients to the main list
    coefficients.extend(missing_coeffs)
    
    return coefficients


def write_csv_all_years(coefficients: List[Tuple[str, int, int, List[float]]], output_file: str, years: List[str]):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        
        # Write header
        header = ['coeff', 'SH_degree', 'SH_order'] + years
        writer.writerow(header)
        
        # Sort coefficients by degree, then order, then coefficient type (g before h)
        coefficients.sort(key=lambda x: (x[1], x[2], x[0]))
        
        # Write data rows
        for coeff_type, degree, order, values in coefficients:
            row = [coeff_type, degree, order] + values
            writer.writerow(row)


def IGRFupdate():
    # Get the path to the functions folder
    script_dir = os.path.dirname(os.path.abspath(__file__))
    functions_dir = os.path.join(script_dir, 'Parameters', 'functions')
    output_file = os.path.join(functions_dir, 'igrf_coefficients.csv')

    try:
        # Check for version argument
        if len(sys.argv) > 1:
            try:
                version = int(sys.argv[1])
                if version < 1 or version > 99:
                    raise ValueError
            except ValueError:
                print("Error: IGRF version must be an integer between 1 and 99 (e.g., 12 for IGRF12)")
                sys.exit(1)
            url = f"https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf{version}coeffs.txt"
            generation = version
            print(f"Downloading specified IGRF-{generation} coefficients...")
        else:
            # Download IGRF data (auto-detects latest version)
            url, generation = find_latest_igrf_url()
            print(f"Found and downloading IGRF-{generation}")

        lines = download_igrf_data(url)

        # Parse header to get available years
        available_years, data_start_line = parse_igrf_header(lines)

        # Process all years into a single CSV
        coefficients = parse_igrf_data_all_years(lines, data_start_line, available_years)

        write_csv_all_years(coefficients, output_file, available_years)

        print(f"Successfully updated IGRF-{generation} coefficients")

    except Exception as e:
        print(f"Error updating IGRF coefficients: {e}")
        sys.exit(1)


def addstation():
    if len(sys.argv) != 4:
        print("Usage: OTSO.addstation <Name> <Latitude> <Longitude>")
        print("Example: OTSO.addstation 'My Station' 45.0 -122.0")
        sys.exit(1)
    
    name = sys.argv[1]
    try:
        latitude = float(sys.argv[2])
        longitude = float(sys.argv[3])
    except ValueError:
        print("Error: Latitude and Longitude must be numbers")
        sys.exit(1)
    
    AddStation(name, latitude, longitude)


def removestation():
    if len(sys.argv) != 2:
        print("Usage: OTSO.removestation <Name>")
        print("Example: OTSO.removestation 'My Station'")
        sys.exit(1)
    
    name = sys.argv[1]
    RemoveStation(name)


def liststations():
    ListStations()




