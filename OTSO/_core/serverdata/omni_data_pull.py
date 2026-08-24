import pandas as pd
import numpy as np
import csv
from datetime import datetime, timedelta
import shutil
import os
import requests

from ..utils.tsy_params_utils import Pdyn_comp
from . import soho_data_pull
from .scratch_dir import SCRATCH_DIR
from .prior_year import load_prior_year_tail

def download_omni_data(year):
    url = f"https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/omni_5min{year}.asc"

    save_path = os.path.join(SCRATCH_DIR, f'omni_5min_{year}.lst')

    #print(f"Downloading OMNI 5-minute data for year {year}...")
    #print(f"URL: {url}")
    #print(f"Saving to: {save_path}")

    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        with open(save_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        #print(f"Downloaded successfully: {save_path}")
        parse_and_convert_to_csv_high_res(save_path, f'omni_{year}_high_res.csv', year=year)

    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return


def download_omni_low_res_data(year):
    base_url = "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/"
    file_name = f'omni2_{year}.dat'
    file_name2 = f'omni2_{year+1}.dat'

    endfile_name = f'omni2_{year}.dat'
    endfile_name2 = f'omni2_{year+1}.dat'

    file_path = os.path.join(SCRATCH_DIR, endfile_name)
    file_path2 = os.path.join(SCRATCH_DIR, endfile_name2)


    # Download primary year file
    url = base_url + file_name
    #print(f"Downloading {file_name}...")
    
    try:
        response = requests.get(url, stream=True, timeout=30)
        response.raise_for_status()

        with open(file_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        parse_and_convert_to_csv_low_res(file_path, f'omni_{year}_low_res.csv')

    except requests.exceptions.RequestException as e:
        print(f"Error downloading {file_name}: {e}")
        return

    # Only download next year's file if not processing the current year
    current_year = datetime.now().year
    if year < current_year:
        url2 = base_url + file_name2
        
        try:
            response2 = requests.get(url2, stream=True, timeout=30)
            response2.raise_for_status()

            with open(file_path2, 'wb') as f:
                for chunk in response2.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)

            #print(f"Downloaded successfully: {file_path2}")
            parse_and_convert_to_csv_low_res(file_path2, f'omni_{year+1}_low_res.csv')

        except requests.exceptions.RequestException as e:
            print(f"Could not download {file_name2}: {e} (this is normal if the file doesn't exist yet)")



def convert_to_datetime(year, decimal_day, hour, minute):
    """ Convert year, decimal day, and hour into a datetime string. """
    year = int(year)
    day_of_year = int(float(decimal_day))
    hour = int(hour)
    minute = int(minute)


    initial_date = datetime(year, 1, 1) + timedelta(days=day_of_year - 1)
    datetime_obj = initial_date + timedelta(hours=hour,minutes=minute)

    return datetime_obj.strftime('%Y-%m-%d %H:%M:%S')

def parse_and_convert_to_csv_high_res(input_file, output_file, year=None):
    """ Parse high-res OMNI data file and convert to CSV with specific columns. """
    with open(input_file, 'r') as datfile:
        lines = datfile.readlines()

    FILLER_VALUES = {"99999", "9999999", "-999.9", "9999.99", "9999.9", "999.9",
                      "99999.0", "-1.00e+05", "999999", "999.99", "99.99", "9999.0"}


    rows = [line.split() for line in lines]

    output_headers = ["Date", "Bx", "By", "Bz", "V", "Density", "Pdyn"]

    output_path = os.path.join(SCRATCH_DIR, output_file)

    with open(output_path, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(output_headers)

        for row in rows:
            try:
                datetime_str = convert_to_datetime(row[0], row[1], row[2],row[3])  # Year, DOY, Hour, Minute
                Bx_value = row[14]
                By_value = row[17]
                Bz_value = row[18]
                V_value = row[21]
                Density_value = row[25]
                Pdyn_value = row[27]

                values_to_check = {Bx_value, By_value, Bz_value, V_value, Density_value, Pdyn_value}

                # Skip the row if any value is a filler
                if any(val in FILLER_VALUES for val in values_to_check):
                    continue

                csv_writer.writerow([datetime_str, Bx_value, By_value, Bz_value, V_value, Density_value, Pdyn_value])
            except (IndexError, ValueError):
                continue

    add_rolling_averages(output_path, year=year)

def add_rolling_averages(csv_file_path, year=None):
    """
    Add 30-minute rolling averages for By and Bz columns to the CSV file.
    OMNI 5-minute data: 30 minutes = 6 data points for rolling average.

    Warm-starts the rolling window from the tail of the previous year's cached data
    (if available) so the first ~30 minutes of the year aren't computed from a
    truncated window - see prior_year.load_prior_year_tail.
    """
    try:
        df = pd.read_csv(csv_file_path)

        df['Date'] = pd.to_datetime(df['Date'])

        df = df.sort_values('Date')

        df['By'] = pd.to_numeric(df['By'], errors='coerce')
        df['Bz'] = pd.to_numeric(df['Bz'], errors='coerce')

        prior_tail = load_prior_year_tail(year, hours=1) if year is not None else None
        if prior_tail is not None and not prior_tail.empty:
            warmup = prior_tail[['Date', 'By', 'Bz']].copy()
            warmup['By'] = pd.to_numeric(warmup['By'], errors='coerce')
            warmup['Bz'] = pd.to_numeric(warmup['Bz'], errors='coerce')
            combined_By = pd.concat([warmup['By'], df['By']], ignore_index=True)
            combined_Bz = pd.concat([warmup['Bz'], df['Bz']], ignore_index=True)
            df['By_avg'] = combined_By.rolling(window=6, min_periods=1, center=False).mean().iloc[len(warmup):].values
            df['Bz_avg'] = combined_Bz.rolling(window=6, min_periods=1, center=False).mean().iloc[len(warmup):].values
        else:
            df['By_avg'] = df['By'].rolling(window=6, min_periods=1, center=False).mean()
            df['Bz_avg'] = df['Bz'].rolling(window=6, min_periods=1, center=False).mean()

        df['By_avg'] = df['By_avg'].round(2)
        df['Bz_avg'] = df['Bz_avg'].round(2)

        df.to_csv(csv_file_path, index=False)

    except Exception as e:
        print(f"Error adding rolling averages: {e}")
        print("Continuing without rolling averages...")

def parse_and_convert_to_csv_low_res(input_file, output_file):
    """ Parse the data file and convert it to a CSV with only specific columns. """
    with open(input_file, 'r') as datfile:
        lines = datfile.readlines()
        
        rows = [line.split() for line in lines]
        
    output_headers = ["Date", "Kp", "Kp_raw", "Dst", "Bx", "By", "Bz", "V", "Density", "Pdyn"]

    output_path = os.path.join(SCRATCH_DIR, output_file)
    
    with open(output_path, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(output_headers)

        for row in rows:
            try:
                datetime_str = convert_to_datetime(row[0], row[1], row[2],0)
                kp_value = process_kp_value(row[38])  # Word 39: Kp
                kp_raw_value = get_raw_kp_value(row[38])  # Word 39: Kp raw
                dst_value = row[40]  # Word 41: DST Index
                Bx_value = row[12]  # Word 13: Bx GSE
                By_value = row[15]  # Word 16: By GSM
                Bz_value = row[16]  # Word 17: Bz GSM
                V_value = row[24]   # Word 25: Plasma speed
                Density_value = row[23]  # Word 24: Proton Density
                Pdyn_value = row[28]  # Word 29: Flow Pressure

                csv_writer.writerow([datetime_str, kp_value, kp_raw_value, dst_value, Bx_value, By_value, Bz_value, V_value, Density_value, Pdyn_value])
            
            except (IndexError, ValueError):
                continue

def process_kp_value(kp_value):
    try:
        kp_value = int(kp_value)
        rounded_kp = round(kp_value, -1) 
        return int(rounded_kp / 10) 
    except ValueError:
        return kp_value

def get_raw_kp_value(kp_value):
    """Get the unrounded Kp value divided by 10."""
    try:
        kp_value = int(kp_value)
        return kp_value / 10.0
    except ValueError:
        return kp_value

def PullOMNI(year):
    download_omni_data(year)
    download_omni_low_res_data(year)

def PullOMNILowRes(year):
    download_omni_low_res_data(year)

def TSY01(File):
    """
    Computes rolling G1/G2/G3 (TSY01 storm-time indices) for each row from a trailing,
    time-based window: up to 12 rows (~1h at 5-min resolution) further restricted to
    Date > current_row_time - 1h. G1 and G2 are the window MEAN of a per-row term; G3
    is the window SUM (not mean) of a different per-row term divided by 2000 - this
    mirrors the original per-row implementation exactly (verified by fuzz-testing
    against it across thousands of synthetic rows, including gaps, filler V==9999.0,
    negative V, and the By==Bz==0 edge case).

    Vectorized because the original ran a full Python loop with a fresh pandas window
    slice per row - about 105k rows/year of 5-minute data, which made this the slowest
    step in downloading a new year's server data.

    One deliberate behavior change: if Bz is NaN, the original crashes with
    UnboundLocalError (its if/elif on Bz<0 / Bz>=0 leaves `Bs` unassigned for NaN).
    Here it propagates as NaN instead, poisoning just that row's window - NaN Bz is not
    expected in practice (OMNI files always populate the field, filler value or not),
    but crashing on it serves nothing.
    """
    filepath = os.path.join(SCRATCH_DIR, File)
    data = pd.read_csv(filepath)

    data['Date'] = pd.to_datetime(data['Date'])

    By = data['By'].to_numpy(dtype=float)
    Bz = data['Bz'].to_numpy(dtype=float)
    V = data['V'].to_numpy(dtype=float)
    N = data['Density'].to_numpy(dtype=float)

    V = np.where(V == 9999.0, np.nan, np.abs(V))

    B = np.sqrt(By ** 2 + Bz ** 2)
    h = (B / 40) ** 2 / (1 + B / 40)

    phi = np.arctan2(By, Bz)
    phi = np.where(phi <= 0, phi + 2 * np.pi, phi)
    phi = np.where((By == 0) & (Bz == 0), 0.0, phi)

    Bs = np.where(Bz < 0, np.abs(Bz), np.where(Bz >= 0, 0.0, np.nan))

    term1 = V * h * np.sin(phi / 2) ** 3  # feeds G1 (window mean)
    term2 = V * Bs                        # feeds G2 (window mean)
    term3 = N * V * Bs                    # feeds G3 (window sum)

    dt_index = pd.DatetimeIndex(data['Date'])
    row_count = pd.Series(1.0, index=dt_index).rolling('1h', min_periods=1).sum()

    def rolling_sum_poisoned(term):
        s = pd.Series(term, index=dt_index)
        window_sum = s.rolling('1h', min_periods=1).sum()
        any_nan_in_window = s.isna().astype(float).rolling('1h', min_periods=1).max() > 0
        return window_sum.where(~any_nan_in_window, np.nan)

    sum1 = rolling_sum_poisoned(term1)
    sum2 = rolling_sum_poisoned(term2)
    sum3 = rolling_sum_poisoned(term3)

    G1_values = np.round((sum1 / row_count).to_numpy(), 2)
    G2_values = np.round((0.005 * sum2 / row_count).to_numpy(), 2)
    G3_values = np.round((sum3 / 2000).to_numpy(), 2)

    pdyn_index = data.columns.get_loc('Pdyn') + 1

    data.insert(pdyn_index, 'G1', G1_values)
    data.insert(pdyn_index + 1, 'G2', G2_values)
    data.insert(pdyn_index + 2, 'G3', G3_values)

    data.to_csv(filepath, index=False)

def CombineLowRes(input_file, year):
    df = pd.read_csv(input_file)
    
    new_headers = ['Date', 'By', 'Bz', 'V', 'Density', 'Pdyn', 'Dst', 'Kp', 'Kp_raw',
                   'G1', 'G2', 'G3']
    

    new_df = pd.DataFrame(columns=new_headers)
    
    new_df['Date'] = df['Date']
    new_df['Kp'] = df['Kp']
    new_df['Kp_raw'] = df['Kp_raw'] if 'Kp_raw' in df.columns else df['Kp']
    new_df['Dst'] = df['Dst']
    
    for col in new_headers:
        if col not in ['Date', 'Kp', 'Kp_raw', 'Dst']:
            new_df[col] = 0

    new_df.to_csv(os.path.join(SCRATCH_DIR, f'{year}_TSY_Inputs.csv'), index=False)
    source_file = os.path.join(SCRATCH_DIR, f'{year}_TSY_Inputs.csv')
    destination_folder = os.path.join(os.path.dirname(__file__), "ServerData")
    os.makedirs(destination_folder, exist_ok=True)
    destination_file = os.path.join(destination_folder, os.path.basename(source_file))
    shutil.move(source_file, destination_file)

def Combine(high_res_file, low_res_file, TSY15file, year):
    base_dir = SCRATCH_DIR
    high_res_file = os.path.join(base_dir, high_res_file)
    low_res_file = os.path.join(base_dir, low_res_file)
    TSY15file = os.path.join(base_dir, TSY15file)
    futurefile = os.path.join(base_dir, f'omni_{year+1}_low_res.csv')

    high_res_df = pd.read_csv(high_res_file, parse_dates=['Date']).set_index('Date')
    low_res_df = pd.read_csv(low_res_file, parse_dates=['Date']).set_index('Date')
    
    tsy15_df = None
    if os.path.exists(TSY15file):
        try:
            tsy15_df = pd.read_csv(TSY15file, parse_dates=['DateTime']).set_index('DateTime')
            tsy15_df.index.name = 'Date'
            
            available_columns = []
            if 'N_index_normalised' in tsy15_df.columns:
                available_columns.append('N_index_normalised')
            elif 'N_index' in tsy15_df.columns:
                available_columns.append('N_index')
                
            if 'B_index' in tsy15_df.columns:
                available_columns.append('B_index')
                
            if 'SYM_H' in tsy15_df.columns:
                available_columns.append('SYM_H')
            elif 'SYM-H' in tsy15_df.columns:
                available_columns.append('SYM-H')
                
            if available_columns:
                tsy15_df = tsy15_df[available_columns]
                
                if 'N_index_normalised' in tsy15_df.columns:
                    tsy15_df = tsy15_df.rename(columns={'N_index_normalised': 'N_index'})
                
        except Exception as e:
            print(f"Error loading TSY15 file: {e}")
            tsy15_df = None
    else:
        print(f"TSY15 file not found: {TSY15file}")
    

    low_res_5min = low_res_df.resample('5min').ffill()
    

    if 'Dst' in low_res_df.columns:
        low_res_df.index = low_res_df.index + pd.Timedelta(hours=1)
        low_res_5min['Dst'] = low_res_df['Dst'].resample('5min').ffill().astype(int)

    # Copy By and Bz to By_avg and Bz_avg for low resolution data
    if 'By' in low_res_5min.columns:
        low_res_5min['By_avg'] = low_res_5min['By']
    if 'Bz' in low_res_5min.columns:
        low_res_5min['Bz_avg'] = low_res_5min['Bz']

    start_time = high_res_df.index.min()
    end_time = high_res_df.index.max()
    full_index = pd.date_range(start=start_time, end=end_time, freq='5min')

    # Reindex and combine with priority: high_res > partial > low_res > tsy15
    combined_df = pd.DataFrame(index=full_index)
    combined_df = combined_df.combine_first(high_res_df)
    combined_df = combined_df.combine_first(low_res_5min)
    if tsy15_df is not None:
        combined_df = combined_df.combine_first(tsy15_df)

    combined_df = combined_df.reset_index().rename(columns={'index': 'Date'})
    
    if 'V' in combined_df.columns:
        combined_df['V'] = combined_df['V'].abs()

    # Fill final Dec 31 values from following year if available
    # Only attempt if not the current year (future file won't exist yet)
    current_year = datetime.now().year
    if year < current_year and os.path.exists(futurefile):
        future_df = pd.read_csv(futurefile)
        dst_value = future_df.loc[0, 'Dst'] if 'Dst' in future_df.columns else None
        kp_value = future_df.loc[0, 'Kp'] if 'Kp' in future_df.columns else None
        kp_raw_value = future_df.loc[0, 'Kp_raw'] if 'Kp_raw' in future_df.columns else None
        if dst_value is not None and kp_value is not None:
            mask = (combined_df['Date'].dt.year == year) & \
                   (combined_df['Date'].dt.month == 12) & \
                   (combined_df['Date'].dt.day == 31) & \
                   (combined_df['Date'].dt.hour == 23) & \
                   (combined_df['Date'].dt.minute >= 5)
            combined_df.loc[mask, 'Dst'] = dst_value
            combined_df.loc[mask, 'Kp'] = kp_value
            if kp_raw_value is not None:
                combined_df.loc[mask, 'Kp_raw'] = kp_raw_value

    # Ensure column completeness and order
    desired_order = ['Date', 'Bx', 'By', 'Bz', 'By_avg', 'Bz_avg', 'V', 'Density', 'Pdyn', 'Dst', 'Kp', 'Kp_raw',
                     'G1', 'G2', 'G3',
                     'N_index', 'B_index', 'SYM_H']
    

    existing_columns = ['Date']  
    for col in desired_order[1:]:  
        if col in combined_df.columns:
            existing_columns.append(col)
    
    for col in existing_columns[1:]:
        if col not in combined_df.columns:
            combined_df[col] = pd.NA
    
    combined_df = combined_df[existing_columns]
    combined_df = combined_df.drop_duplicates(subset='Date', keep='first')

    output_file = os.path.join(base_dir, f'{year}_TSY_Inputs.csv')
    combined_df.to_csv(output_file, index=False)
    
    if soho_data_pull.celias_url_exists(year):
        if year >= 1996:
            SohoCombine(f'{year}_TSY_Inputs.csv', f'{year}_CELIAS_Proton_Monitor_5min.csv', year)

    tsy_inputs_path = os.path.join(base_dir, f'{year}_TSY_Inputs.csv')
    df = pd.read_csv(tsy_inputs_path, parse_dates=['Date'])
    df = df.set_index('Date')
    FILLER_VALUES = {"99999", "9999999", "-999.9", "9999.99", "9999.0",
                     "9999.9", "999.9", "99999.0", "-1.00e+05", "999999", "999.99", "99.99"}

    FILLER_VALUES_NUMERIC = set()
    for val in FILLER_VALUES:
        try:
            FILLER_VALUES_NUMERIC.add(float(val))
        except Exception:
            pass
    # Replace filler values (string and numeric) with NaN for all columns
    for col in df.columns:
        # Replace string versions
        df[col] = df[col].replace(FILLER_VALUES, pd.NA)
        # Replace numeric versions
        df[col] = df[col].replace(list(FILLER_VALUES_NUMERIC), pd.NA)

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    for col in df.columns:
        if df[col].isnull().any():
            isnan = df[col].isnull()
            mask = isnan.copy()
            start = None
            for i, val in enumerate(isnan):
                if val and start is None:
                    start = i
                if not val and start is not None:
                    end = i
                    if end - start <= 36:
                        mask.iloc[start:end] = True
                    else:
                        mask.iloc[start:end] = False
                    start = None

            if start is not None:
                end = len(isnan)
                if end - start <= 36:
                    mask.iloc[start:end] = True
                else:
                    mask.iloc[start:end] = False

            df.loc[mask, col] = df[col].interpolate(limit=36, limit_direction='both')[mask]
    df = df.reset_index()
    df.to_csv(tsy_inputs_path, index=False)

def CombineTSY04(TSY04file, year):
    base_dir = SCRATCH_DIR
    TSY04file = os.path.join(base_dir, TSY04file)
    mainfile = os.path.join(base_dir, f'{year}_TSY_Inputs.csv')


    main_df = pd.read_csv(mainfile, parse_dates=['Date'])
    tsy04_df = pd.read_csv(TSY04file, parse_dates=['Date'])


    merged = pd.merge(main_df, tsy04_df, on='Date', how='left', suffixes=('', '_TSY04'))

    cols = list(merged.columns)
    try:
        g3_idx = cols.index('G3')
    except ValueError:

        g3_idx = len(cols) - 1

    for w in ['W1','W2','W3','W4','W5','W6']:
        if w in cols:
            cols.remove(w)

    for i, w in enumerate(['W1','W2','W3','W4','W5','W6']):
        if w in merged.columns:
            cols.insert(g3_idx + 1 + i, w)

    merged = merged[cols]

    merged.to_csv(mainfile, index=False)

def move_to_server_data(file_name):
    destination_folder = os.path.join(os.path.dirname(__file__), "ServerData")
    output_file = os.path.join(SCRATCH_DIR, file_name)
    os.makedirs(destination_folder, exist_ok=True)
    shutil.move(output_file, os.path.join(destination_folder, os.path.basename(output_file)))

def SohoCombine(TSYfile, celias_file, year):
    base_dir = SCRATCH_DIR
    TSYfile = os.path.join(base_dir, TSYfile)
    celias_file = os.path.join(base_dir, celias_file)


    tsy_df = pd.read_csv(TSYfile, parse_dates=['Date']).set_index('Date')
    celias_df = pd.read_csv(celias_file, parse_dates=['Date']).set_index('Date')



    for idx, row in tsy_df.iterrows():
        if 'V' in tsy_df.columns and idx in celias_df.index:
            if row['V'] == 9999.0:
                new_v = celias_df.at[idx, 'V']
                tsy_df.at[idx, 'V'] = new_v
        if 'Density' in tsy_df.columns and idx in celias_df.index:
            if row['Density'] == 999.9:
                new_density = celias_df.at[idx, 'Density']
                tsy_df.at[idx, 'Density'] = new_density
        if 'Pdyn' in tsy_df.columns:
            if row['Pdyn'] == 99.99:
                new_v = tsy_df.at[idx, 'V']
                new_density = tsy_df.at[idx, 'Density']
                if pd.notna(new_v) and pd.notna(new_density):
                    try:
                        new_pdyn = Pdyn_comp(float(new_density), float(new_v))
                        tsy_df.at[idx, 'Pdyn'] = new_pdyn
                    except Exception:
                        pass

    combined_df = tsy_df.reset_index()

    output_file = os.path.join(base_dir, f'{year}_TSY_Inputs.csv')
    combined_df.to_csv(output_file, index=False)



def Omnidelete(OMNIYEAR):

    current_directory = SCRATCH_DIR

    filelist = [
        f'omni_{OMNIYEAR}_low_res.csv',
        f'omni_5min_{OMNIYEAR}.lst',
        f'omni2_{OMNIYEAR}.dat',
        f'{OMNIYEAR}_IMF_&_SW_gaps_le_3hrs_filled.txt',
        f'{OMNIYEAR}_IMF_gaps_le_3hrs_filled.txt',
        f'{OMNIYEAR}_Interval_list.txt',
        f'{OMNIYEAR}_OMNI_5m_with_TA15_drivers.dat',
        f'omni_{OMNIYEAR+1}_low_res.csv',
        f'omni2_{OMNIYEAR+1}.dat',
        f'omni_{OMNIYEAR}_high_res.csv',
        f'tempfile.csv',
        f'TSY15_{OMNIYEAR}.csv',
        f'omni_5min{OMNIYEAR}.asc',
        f'{OMNIYEAR}_CELIAS_Proton_Monitor_5min.csv',
        f'{OMNIYEAR}_CELIAS_Proton_Monitor_5min.zip',
        f'TSY04_W_parameters_{OMNIYEAR}.csv',
        f'{OMNIYEAR}_TSY_Inputs.csv',
    ]

    for i in filelist:
        file_path = os.path.join(current_directory, i)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Error: {e}")



def OmnideleteLowRes(OMNIYEAR):

    current_directory = SCRATCH_DIR

    file2 = f'omni_{OMNIYEAR}_low_res.csv'
    file3 = f'omni_{OMNIYEAR+1}_low_res.csv'
    file4 = f'omni2_{OMNIYEAR+1}.dat'
    file5 = f'omni2_{OMNIYEAR}.dat'
    
    filelist = [file2,file3,file4,file5]

    for i in filelist:
        file_path = os.path.join(current_directory, i)
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
            except Exception as e:
                print(f"Error: {e}")
