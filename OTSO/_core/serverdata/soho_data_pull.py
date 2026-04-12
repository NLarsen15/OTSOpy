import requests
import zipfile
import io
import pandas as pd
import re
import os

def celias_url_exists(year: int) -> bool:
    data_url = f"https://l1.umd.edu/data/{year}_CELIAS_Proton_Monitor_5min.zip"
    try:
        response = requests.head(data_url, allow_redirects=True, timeout=10)
        return response.status_code == 200
    except Exception:
        return False

def download_and_unpack_celias(YEAR: int) -> None:
        
    data_url = f"https://l1.umd.edu/data/{YEAR}_CELIAS_Proton_Monitor_5min.zip"
    response = requests.get(data_url)
    response.raise_for_status()
    script_dir = os.path.join(os.path.dirname(__file__))

    # Save ZIP file to script directory
    zip_filename = os.path.basename(data_url)
    zip_path = os.path.join(script_dir, zip_filename)
    with open(zip_path, 'wb') as zip_file:
        zip_file.write(response.content)

    with zipfile.ZipFile(zip_path) as z:
        txt_files = [f for f in z.namelist() if f.lower().endswith('.txt')]
        if not txt_files:
            raise ValueError("No .txt file found in the ZIP archive.")
        txt_filename = txt_files[0]
        with z.open(txt_filename) as f:
            lines = f.read().decode('utf-8', errors='replace').splitlines()

    data_start = 0
    for i, line in enumerate(lines):
        if re.match(r"^\d{2} \w{3} ", line):
            data_start = i
            break

    data_lines = lines[data_start:]

    columns = [
        "YY", "MON", "DY", "DOY:HH:MM:SS", "SPEED", "Np", "Vth", "N/S", "V_He",
        "GSE_X", "GSE_Y", "GSE_Z", "RANGE", "HGLAT", "HGLONG", "CRN(E)"
    ]

    def parse_line(line):
        m = re.match(r"(\d{2}) (\w{3}) (\d{2}) (\d{3}:\d{2}:\d{2}:\d{2}) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+) +([\d\-.]+)", line)
        if not m:
            return None
        return [m.group(i) for i in range(1, 17)]

    parsed = [parse_line(l) for l in data_lines]
    parsed = [row for row in parsed if row]

    df = pd.DataFrame(parsed, columns=columns)

    month_map = {m: i for i, m in enumerate(['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'], 1)}
    def parse_datetime(row):
        try:
            year = int(row['YY'])
            year += 1900 if year > 50 else 2000
            month = month_map.get(row['MON'].title(), 1)
            day = int(row['DY'])
            doy, hour, minute, second = map(int, re.split('[:]', row['DOY:HH:MM:SS']))
            from datetime import datetime, timedelta
            dt = datetime(year, 1, 1) + timedelta(days=doy-1, hours=hour, minutes=minute, seconds=second)
            return dt
        except Exception:
            return pd.NaT
    df['datetime'] = df.apply(parse_datetime, axis=1)
    df = df.drop(['YY', 'MON', 'DY', 'DOY:HH:MM:SS'], axis=1)
    cols = ['datetime'] + [c for c in df.columns if c != 'datetime']
    df = df[cols]

    df = df[['datetime', 'SPEED', 'Np']]
    df = df.rename(columns={'SPEED': 'V', 'Np': 'Density'})
    df = df.set_index('datetime')
    df['V'] = pd.to_numeric(df['V'], errors='coerce')
    df['Density'] = pd.to_numeric(df['Density'], errors='coerce')

    df_resampled = df.resample('5min', label='right', closed='right').mean()
    df_rolling = df_resampled.rolling('1h', min_periods=1).mean()
    df_rolling = df_rolling.shift(freq='1h')
    start = df_rolling.index.min().replace(second=0, microsecond=0)
    if start.minute % 5 != 0:
        from datetime import timedelta
        start += timedelta(minutes=5 - start.minute % 5)
    end = df_rolling.index.max().replace(second=0, microsecond=0)
    new_index = pd.date_range(start, end, freq='5min')
    df_rolling = df_rolling.reindex(new_index)
    df_rolling.index.name = 'Date'

    df_rolling = df_rolling.round({'V': 3, 'Density': 3})
    output_filename = f"{YEAR}_CELIAS_Proton_Monitor_5min.csv"
    output_path = os.path.join(script_dir, output_filename)
    df_rolling.to_csv(output_path)
