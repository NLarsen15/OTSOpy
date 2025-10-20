
import pandas as pd
from datetime import datetime, timedelta
import re
import os
from OTSO.Parameters.functions.N_B_variables import N_index_normalized, B_index

output_directory = os.path.dirname(os.path.realpath(__file__))




def extract_30min_average(current_time):
    """
    Returns 30-minute average values for both solar wind and magnetic parameters prior to current_time,
    as well as N_normalized and B_normalized variables.
    Returns: (speed_avg, density_avg, By_avg, Bz_avg, Bx_avg, N_norm, B_norm)
    """
    # Solar wind
    sw_file = os.path.join(output_directory, 'space_data.csv')
    sw_df = pd.read_csv(sw_file)
    sw_df = sw_df.dropna(how='any')
    sw_df['time_tag'] = pd.to_datetime(sw_df['time_tag'])
    sw_df['fetch_timestamp'] = pd.to_datetime(sw_df['fetch_timestamp'])
    sw_avg = average_30min_prior(sw_df, current_time)
    density_avg = sw_avg['density'] if 'density' in sw_avg else None
    speed_avg = sw_avg['speed'] if 'speed' in sw_avg else None
    if speed_avg is not None:
        speed_avg = speed_avg * -1

    # Magnetic
    mag_file = os.path.join(output_directory, 'Magnetic_data.csv')
    mag_df = pd.read_csv(mag_file)
    mag_df = mag_df.dropna(how='any')
    mag_df['time_tag'] = pd.to_datetime(mag_df['time_tag'])
    mag_df['fetch_timestamp'] = pd.to_datetime(mag_df['fetch_timestamp'])
    mag_avg = average_30min_prior(mag_df, current_time)
    By_avg = mag_avg['by_gsm'] if 'by_gsm' in mag_avg else None
    Bz_avg = mag_avg['bz_gsm'] if 'bz_gsm' in mag_avg else None
    # Try to get Bx if available
    Bx_avg = mag_avg['bx_gsm'] if 'bx_gsm' in mag_avg else (mag_avg['bx'] if 'bx' in mag_avg else None)

    # Compute Bt and clock angle (thetac)
    Bt = None
    thetac = None
    if By_avg is not None and Bz_avg is not None:
        Bt = (By_avg**2 + Bz_avg**2)**0.5
        import math
        thetac = math.atan2(By_avg, Bz_avg)
        if thetac < 0:
             thetac += 2 * math.pi
 
    # Compute N_normalized and B_normalized if possible
    N_norm = None
    B_norm = None
    if speed_avg is not None and Bt is not None and thetac is not None:
        try:
            N_norm = N_index_normalized(abs(speed_avg), Bt, thetac)
        except Exception:
            N_norm = None
    if density_avg is not None and speed_avg is not None and Bt is not None and thetac is not None:
        try:
            B_norm = B_index(density_avg, abs(speed_avg), Bt, thetac)
        except Exception:
            B_norm = None

    # Ensure all returned values are real numbers (not complex)
    def to_real_rounded(val):
        if isinstance(val, complex):
            val = val.real
        if isinstance(val, float):
            return round(val, 3)
        return val

    return (
        to_real_rounded(speed_avg),
        to_real_rounded(density_avg),
        to_real_rounded(By_avg),
        to_real_rounded(Bz_avg),
        to_real_rounded(Bx_avg),
        to_real_rounded(N_norm),
        to_real_rounded(B_norm)
    )

def extract_solar_wind(current_time):
    """
    Returns 1-hour average values for solar wind parameters prior to current_time.
    Returns: (speed_avg, density_avg)
    """
    current_time = current_time - timedelta(hours=1)
    file_path = os.path.join(output_directory, 'space_data.csv')
    df = pd.read_csv(file_path)
    df = df.dropna(how='any')
    df['time_tag'] = pd.to_datetime(df['time_tag'])
    df['fetch_timestamp'] = pd.to_datetime(df['fetch_timestamp'])
    hourly_avg = hourly_average(df, current_time)
    density_avg = hourly_avg['density'] if 'density' in hourly_avg else None
    speed_avg = hourly_avg['speed'] if 'speed' in hourly_avg else None
    if speed_avg is not None:
        speed_avg = speed_avg * -1
    return speed_avg, density_avg



def extract_magnetic(current_time):


    current_time = current_time - timedelta(hours=1)
    file_path = os.path.join(output_directory, f'Magnetic_data.csv')
    df = pd.read_csv(file_path)
    df = df.dropna(how='any')
    df['time_tag'] = pd.to_datetime(df['time_tag'])
    df['fetch_timestamp'] = pd.to_datetime(df['fetch_timestamp'])
    hourly_avg = hourly_average(df, current_time)
    By_avg = hourly_avg['by_gsm'] if 'by_gsm' in hourly_avg else None
    Bz_avg = hourly_avg['bz_gsm'] if 'bz_gsm' in hourly_avg else None
    return By_avg, Bz_avg

def extract_magnetic_30min(current_time):
    """
    Returns 30-minute average values for magnetic parameters prior to current_time.
    """
    file_path = os.path.join(output_directory, f'Magnetic_data.csv')
    df = pd.read_csv(file_path)
    df = df.dropna(how='any')
    df['time_tag'] = pd.to_datetime(df['time_tag'])
    df['fetch_timestamp'] = pd.to_datetime(df['fetch_timestamp'])
    avg = average_30min_prior(df, current_time)
    By_avg = avg['by_gsm'] if 'by_gsm' in avg else None
    Bz_avg = avg['bz_gsm'] if 'bz_gsm' in avg else None
    return By_avg, Bz_avg


def hourly_average(df, lookup_time):
    lookup_time_ceil = pd.Timestamp(lookup_time).ceil('h')
    
    df_filtered = df[df['time_tag'].dt.floor('h') == lookup_time_ceil]
    
    hourly_avg = df_filtered.mean(numeric_only=True)

    return hourly_avg

def average_30min_prior(df, lookup_time):
    """
    Returns the average of all numeric columns for the 30 minutes prior to lookup_time.
    """
    end_time = pd.Timestamp(lookup_time)
    start_time = end_time - pd.Timedelta(minutes=30)
    mask = (df['time_tag'] > start_time) & (df['time_tag'] <= end_time)
    df_filtered = df.loc[mask]
    avg = df_filtered.mean(numeric_only=True)
    # Ensure all values are real numbers (not complex)
    avg = avg.apply(lambda x: x.real if isinstance(x, complex) else x)
    return avg

def extract_dst_value(file_path, current_time):
    current_day = current_time.day
    current_hour = current_time.hour

    with open(file_path, 'r') as file:
        lines = file.readlines()

    dst_values = {}
    daily_averages = {}

    # Initialize variables to store results in case no data is found
    prev_dst = None
    daily_avg = None

    for line in lines:
        if line.startswith('DST'):
            # Extract year_month and day from the line
            year_month = line[3:7]
            day = int(line[8:10])

            # Check if this line corresponds to the current day
            if day == current_day:
                # Use a regex to extract all numbers, including potential negative and multi-digit values
                dst_data = re.findall(r'-?\d+', line[11:])
                dst_data = dst_data[2:]
                
                # Ignore invalid values (e.g., "999") when extracting hourly data
                hourly_dst_data = [int(val) for val in dst_data[:-1] if int(val) != 999]
                
                # Check if the data for the current hour exists
                if len(hourly_dst_data) > current_hour:
                    prev_dst = hourly_dst_data[current_hour]
                    dst_values[(year_month, day)] = prev_dst

                # Extract the daily average, which is the last value
                try:
                    daily_average = int(dst_data[-1])
                    if daily_average != 999:  # Ignore invalid daily averages
                        daily_averages[(year_month, day)] = daily_average
                except (ValueError, IndexError):
                    pass

    # Retrieve the daily average for the current day if available
    daily_avg = daily_averages.get((year_month, current_day), None)

    return prev_dst, daily_avg

def extract_kp_index(current_time):
    
    kp_index_before = None
    kp_index_after = None
    date_time_before = None
    date_time_after = None

    file_path = os.path.join(output_directory, f'Kp_data.txt')

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            
            components = line.split()
            if len(components) < 9:
                continue
            
            date_str = f"{components[0]} {components[1]} {components[2]}"
            time_str = components[3]
            datetime_str = f"{date_str} {time_str}"
            current_datetime = datetime.strptime(datetime_str, '%Y %m %d %H.%M')

            kp_index = float(components[7])

            if current_datetime <= current_time:
                if kp_index_before is None or current_datetime > date_time_before:
                    kp_index_before = kp_index
                    date_time_before = current_datetime
            if current_datetime > current_time:
                if kp_index_after is None or current_datetime < date_time_after:
                    kp_index_after = kp_index
                    date_time_after = current_datetime

    if kp_index_before is not None and kp_index_after is None:
        kp_index_before = round(kp_index_before)
        if kp_index_before >= 6:
            IOPT = 7
        else:
            IOPT = kp_index_before + 1
            IOPT = round(IOPT)

        return kp_index_before, IOPT
    elif kp_index_before is not None and kp_index_after is not None:
        kp_index_before = round(kp_index_before)
        if kp_index_before >= 6:
            IOPT = 7
        else:
            IOPT = kp_index_before + 1
            IOPT = round(IOPT)
        return kp_index_before,IOPT
    else:
        return None
