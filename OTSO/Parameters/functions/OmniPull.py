import pandas as pd
import csv, math
from datetime import datetime, timedelta
import shutil
import os
import numpy as np
import requests

def download_omni_data(year):
    url = f"https://spdf.gsfc.nasa.gov/pub/data/omni/high_res_omni/omni_5min{year}.asc"
    
    # Get the directory of the current script
    script_dir = os.path.join(os.path.dirname(__file__), "")
    save_path = os.path.join(script_dir, f'omni_5min_{year}.lst')

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
        parse_and_convert_to_csv_high_res(save_path, f'omni_{year}_high_res.csv')

    except requests.exceptions.RequestException as e:
        print(f"Error downloading file: {e}")
        return


def download_omni_low_res_data(year):
    base_url = "https://spdf.gsfc.nasa.gov/pub/data/omni/low_res_omni/"
    file_name = f'omni2_{year}.dat'
    file_name2 = f'omni2_{year+1}.dat'

    endfile_name = f'functions/omni2_{year}.dat'
    endfile_name2 = f'functions/omni2_{year+1}.dat'

    csvfile_name = f'omni_{year}_low_res.csv'
    csvfile_name2 = f'omni_{year+1}_low_res.csv'

    script_dir = os.path.dirname(os.path.dirname(__file__))
    file_path = os.path.join(script_dir, endfile_name)
    file_path2 = os.path.join(script_dir, endfile_name2)

    csvfile_path = os.path.join(script_dir, csvfile_name)
    csvfile_path2 = os.path.join(script_dir, csvfile_name2)
    
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

def parse_and_convert_to_csv_high_res(input_file, output_file):
    """ Parse high-res OMNI data file and convert to CSV with specific columns. """
    with open(input_file, 'r') as datfile:
        lines = datfile.readlines()

    FILLER_VALUES = {"99999", "9999999", "-999.9", "9999.99", "9999.9", "999.9", "99999.0", "-1.00e+05", "999999", "999.99"}


    rows = [line.split() for line in lines]

    output_headers = ["Date", "Bx", "By", "Bz", "V", "Density", "Pdyn"]

    script_dir = os.path.dirname(os.path.realpath(__file__))
    output_path = os.path.join(script_dir, output_file)

    with open(output_path, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(output_headers)

        for row in rows:
            try:
                datetime_str = convert_to_datetime(row[0], row[1], row[2],row[3])  # Year, DOY, Hour, Minute
                Bx_value = row[16]
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

    add_rolling_averages(output_path)

def add_rolling_averages(csv_file_path):
    """
    Add 30-minute rolling averages for By and Bz columns to the CSV file.
    OMNI 5-minute data: 30 minutes = 6 data points for rolling average.
    """
    try:
        df = pd.read_csv(csv_file_path)
        
        df['Date'] = pd.to_datetime(df['Date'])
        
        df = df.sort_values('Date')

        df['By'] = pd.to_numeric(df['By'], errors='coerce')
        df['Bz'] = pd.to_numeric(df['Bz'], errors='coerce')
        
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
        
    output_headers = ["Date", "Kp", "Dst", "Bx", "By", "Bz", "V", "Density", "Pdyn"]

    script_dir = os.path.dirname(os.path.realpath(__file__))
    output_path = os.path.join(script_dir, output_file)
    
    with open(output_path, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(output_headers)

        for row in rows:
            try:
                datetime_str = convert_to_datetime(row[0], row[1], row[2],0)
                kp_value = process_kp_value(row[38])
                dst_value = row[40]
                Bx_value = row[14]
                By_value = row[15]
                Bz_value = row[16]
                V_value = row[24]
                Density_value = row[23]
                Pdyn_value = row[28]

                csv_writer.writerow([datetime_str, kp_value, dst_value, Bx_value, By_value, Bz_value, V_value, Density_value, Pdyn_value])
            
            except (IndexError, ValueError):
                continue

def parse_and_convert_to_csv(input_file, output_file):
    """ Parse the data file and convert it to a CSV with only specific columns. """
    with open(input_file, 'r') as datfile:
        lines = datfile.readlines()
        
        rows = [line.split() for line in lines]
        
    output_headers = ["Date", "Kp", "Dst", "Bx", "By", "Bz", "V", "Density", "Pdyn"]

    script_dir = os.path.dirname(os.path.realpath(__file__))
    output_path = os.path.join(script_dir, output_file)
    
    with open(output_path, mode='w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(output_headers)

        for row in rows:
            try:
                datetime_str = convert_to_datetime(row[0], row[1], row[2],0)
                kp_value = process_kp_value(row[38])
                dst_value = row[40]
                Bx_value = row[14]
                By_value = row[15]
                Bz_value = row[16]
                V_value = row[24]
                Density_value = row[23]
                Pdyn_value = row[28]

                csv_writer.writerow([datetime_str, kp_value, dst_value, Bx_value, By_value, Bz_value, V_value, Density_value, Pdyn_value])
            
            except (IndexError, ValueError):
                continue

def process_kp_value(kp_value):
    try:
        kp_value = int(kp_value)
        rounded_kp = round(kp_value, -1) 
        return int(rounded_kp / 10) 
    except ValueError:
        return kp_value

def extract_row_from_csv(output_file_path, target_datetime):
    target_datetime_str = target_datetime.strftime('%d/%m/%y %H:00:00')

    with open(output_file_path, 'r') as csvfile:
        csv_reader = csv.DictReader(csvfile)

        for row in csv_reader:
            if row['Datetime'] == target_datetime_str:
                validate_row(row)
                return row

    print(f"No matching row found for the given datetime: {target_datetime_str}")
    return None

def validate_row(row):
    for key in row:
        if row[key] == 'NaN' or row[key] == '9999.' or row[key] == '999.9':
            print("ERROR: No valid OMNI data found for given date. Please check other sources. \nOTSO will now terminate.")
            exit()

def PullOMNI(year):
    download_omni_data(year)
    download_omni_low_res_data(year)

def PullOMNILowRes(year):
    download_omni_low_res_data(year)

def OMNI_to_csv(year):
    headers = [
        'BXGSM', 'By', 'Bz', 'V', 
        'VYGSE', 'VZGSE', 'Density', 'TEMP', 'IMFFLAG', 'ISWFLAG', 'TILT', 
        'Pdyn', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6', 'Date'
    ]

    file = os.path.join(os.path.dirname(__file__), f'{year}_OMNI_5m_with_TS05_variables.dat')
    data = pd.read_csv(file, sep=r'\s+', header=None)

    IYEAR = data[0]
    IDAY = data[1]
    IHOUR = data[2]
    MIN = data[3]

    def split_column(value):
       if '*' in str(value):
           return value.split('*')[0], 0.0
       return float(value), None

    split_values = data[13].apply(split_column)
    data[13] = split_values.apply(lambda x: x[0])

    if split_values.apply(lambda x: x[1]).notna().any():
        data.insert(14, 'new_column', split_values.apply(lambda x: x[1]).fillna(0.0))

    data[13] = data[13].astype(float)
    
    data['Date'] = [
        datetime(int(year), 1, 1) + pd.to_timedelta(int(day)-1 if int(day) > 0 else 0, unit='d') + \
        pd.to_timedelta(int(hour), unit='h') + pd.to_timedelta(int(minute), unit='m')
        for year, day, hour, minute in zip(IYEAR, IDAY, IHOUR, MIN)
    ]

    data = data.drop(columns=[0, 1, 2, 3])
    data.columns = headers

    columns_to_remove = ['BXGSM', 'VYGSE', 'VZGSE', 'TEMP', 'IMFFLAG', 'ISWFLAG', 'TILT']
    data = data.drop(columns=columns_to_remove, errors='ignore')

    remaining_columns = [col for col in headers if col != 'Date' and col not in columns_to_remove]
    data = data[['Date'] + remaining_columns] 
    
    csv_filename = os.path.join(os.path.dirname(__file__), f'{year}_TSY_Data.csv')
    data.to_csv(csv_filename, index=False)


def TSY01(File):
    filepath = os.path.join(os.path.dirname(__file__), File)
    data = pd.read_csv(filepath)

    data['Date'] = pd.to_datetime(data['Date'])

    results = []
    G1_values = []
    G2_values = []
    G3_values = []

    for i in range(len(data)):
        current_time = data['Date'].iloc[i]
        start_index = max(0, i - 11)
        relevant_rows = data.iloc[start_index:i + 1]
        relevant_rows = relevant_rows[relevant_rows['Date'] > (current_time - pd.Timedelta(hours=1))]
    

        G1,G2,G3 = TSY01_Constants(relevant_rows)
        G1_values.append(G1)
        G2_values.append(G2)
        G3_values.append(G3)

    pdyn_index = data.columns.get_loc('Pdyn') + 1

    data.insert(pdyn_index, 'G1', G1_values)
    data.insert(pdyn_index + 1, 'G2', G2_values)
    data.insert(pdyn_index + 2, 'G3', G3_values)

    data.to_csv(filepath, index=False)
    

def TSY01_Constants(data):

    IMFy = data["By"]
    IMFz = data["Bz"]
    Speed = data["V"]
    Density = data["Density"]
    G1 = 0
    G2 = 0
    G3 = 0

    for (By, Bz, V, N) in zip(IMFy, IMFz, Speed, Density):
        By = By
        Bz = Bz
        if V < 0:
            V = -1*V
        W = 1/len(IMFy)
    
    
        B = (By*By + Bz*Bz)**(0.5)
        h = (((B/40)**(2))/(1 + B/40))
    
        if(By == 0 and Bz == 0):
            phi = 0
        else:
            phi = math.atan2(By,Bz)
            if(phi <= 0):
                phi = phi + 2*math.pi
    
        if(Bz < 0):
            Bs = abs(Bz)
        elif Bz >= 0:
            Bs = 0
        
        G1 = G1 + (W)*V*h*(math.sin(phi/2)**3)
        G2 = G2 + (0.005)*(W)*(V)*Bs
        G3 = G3 + (N*V*Bs)/2000

    return round(G1, 2), round(G2, 2), round(G3, 2)

def CombineLowRes(input_file, year):
    df = pd.read_csv(input_file)
    
    new_headers = ['Date', 'By', 'Bz', 'V', 'Density', 'Pdyn', 'Dst', 'Kp', 
                   'G1', 'G2', 'G3', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    

    new_df = pd.DataFrame(columns=new_headers)
    
    new_df['Date'] = df['Date']
    new_df['Kp'] = df['Kp']
    new_df['Dst'] = df['Dst']
    
    for col in new_headers:
        if col not in ['Date', 'Kp', 'Dst']:
            new_df[col] = 0

    new_df.to_csv(f'{year}_TSY_Inputs.csv', index=False)
    source_file = f'{year}_TSY_Inputs.csv'
    destination_folder = os.path.join(os.path.dirname(__file__), "ServerData")
    os.makedirs(destination_folder, exist_ok=True)
    destination_file = os.path.join(destination_folder, os.path.basename(source_file))
    shutil.move(source_file, destination_file)
    
import os
import pandas as pd
import shutil

def Combine(TSYfile, high_res_file, low_res_file, TSY15file, year):
    base_dir = os.path.dirname(__file__)
    TSYfile = os.path.join(base_dir, TSYfile)
    high_res_file = os.path.join(base_dir, high_res_file)
    low_res_file = os.path.join(base_dir, low_res_file)
    TSY15file = os.path.join(base_dir, TSY15file)
    futurefile = os.path.join(base_dir, f'omni_{year+1}_low_res.csv')

    partial_df = pd.read_csv(TSYfile, parse_dates=['Date']).set_index('Date')
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
                
                #print(f"Loaded TSY15 data with columns: {available_columns}")
                #if 'N_index_normalised' in available_columns:
                    #print("Renamed N_index_normalised to N_index")
            #else:
                #print("Warning: No N_index, B_index, or SYM_H columns found in TSY15 file")
                #tsy15_df = None
        except Exception as e:
            print(f"Error loading TSY15 file: {e}")
            tsy15_df = None
    else:
        print(f"TSY15 file not found: {TSY15file}")
    

    low_res_5min = low_res_df.resample('5min').ffill()
    

    if 'Dst' in low_res_df.columns:
        low_res_df.index = low_res_df.index + pd.Timedelta(hours=1)
        low_res_5min['Dst'] = low_res_df['Dst'].resample('5min').bfill().astype(int)


    start_time = partial_df.index.min()
    end_time = partial_df.index.max()
    full_index = pd.date_range(start=start_time, end=end_time, freq='5min')

    # Reindex and combine with priority: high_res > partial > low_res > tsy15
    combined_df = pd.DataFrame(index=full_index)
    combined_df = combined_df.combine_first(high_res_df)
    combined_df = combined_df.combine_first(partial_df)
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
        if dst_value is not None and kp_value is not None:
            mask = (combined_df['Date'].dt.year == year) & \
                   (combined_df['Date'].dt.month == 12) & \
                   (combined_df['Date'].dt.day == 31) & \
                   (combined_df['Date'].dt.hour == 23) & \
                   (combined_df['Date'].dt.minute >= 5)
            combined_df.loc[mask, 'Dst'] = dst_value
            combined_df.loc[mask, 'Kp'] = kp_value

    # Ensure column completeness and order
    desired_order = ['Date', 'Bx', 'By', 'Bz', 'By_avg', 'Bz_avg', 'V', 'Density', 'Pdyn', 'Dst', 'Kp',
                     'G1', 'G2', 'G3', 'W1', 'W2', 'W3', 'W4', 'W5', 'W6',
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

    # Write output CSV
    output_file = os.path.join(base_dir, f'{year}_TSY_Inputs.csv')
    combined_df.to_csv(output_file, index=False)

    # Move to ServerData folder
    destination_folder = os.path.join(base_dir, "ServerData")
    os.makedirs(destination_folder, exist_ok=True)
    shutil.move(output_file, os.path.join(destination_folder, os.path.basename(output_file)))




def Omnidelete(OMNIYEAR):

    current_directory = os.path.join(os.path.dirname(__file__), "")

    file1 = f'{OMNIYEAR}_TSY_Data.csv'
    file2 = f'omni_{OMNIYEAR}_low_res.csv'
    file3 = f'omni_5min_{OMNIYEAR}.lst'
    file4 = f'omni2_{OMNIYEAR}.dat'
    file5 = f'{OMNIYEAR}_IMF_&_SW_gaps_le_3hrs_filled.txt'
    file6 = f'{OMNIYEAR}_IMF_gaps_le_3hrs_filled.txt'
    file7 = f'{OMNIYEAR}_Interval_list.txt'
    file8 = f'{OMNIYEAR}_OMNI_5m_with_TS05_variables.dat'
    file9 = f'omni_{OMNIYEAR+1}_low_res.csv'
    file10 = f'omni2_{OMNIYEAR+1}.dat'
    file11 = f'omni_{OMNIYEAR}_high_res.csv'
    file12 = f'tempfile.csv'
    file13 = f'TSY15_{OMNIYEAR}.csv'
    file14 = f'omni_5min{OMNIYEAR}.asc'
    
    if OMNIYEAR < datetime.now().year:
        filelist = [file1,file2,file3,file4,file5,file6,file7,file8,file9,file10,file11,file12,file13,file14]
    else:
        filelist = [file1,file2,file3,file4,file5,file6,file7,file8]

    for i in filelist:
        try:
            os.remove(os.path.join(current_directory, i))
        except Exception as e:
            print(f"Error: {e}")

    if os.path.exists(file9):
        os.remove(os.path.join(current_directory, file9))

def OmnideleteLowRes(OMNIYEAR):

    current_directory = os.path.dirname(__file__)

    file2 = f'omni_{OMNIYEAR}_low_res.csv'
    file3 = f'omni_{OMNIYEAR+1}_low_res.csv'
    file4 = f'omni2_{OMNIYEAR+1}.dat'
    file5 = f'omni2_{OMNIYEAR}.dat'
    
    if OMNIYEAR < datetime.now().year:
        filelist = [file2,file3,file4,file5]
    else:
        filelist = [file2,file3,file5]

    for i in filelist:
        try:
            os.remove(os.path.join(current_directory, i))
        except Exception as e:
            print(f"Error: {e}")

    if os.path.exists(file3):
        os.remove(os.path.join(current_directory, file3))
