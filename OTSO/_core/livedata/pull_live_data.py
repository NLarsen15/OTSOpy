import requests
import pandas as pd
import os
import json
import time
from datetime import datetime, timedelta

from .extract import *
from ..utils.tsy_params_utils import TSY01_Constants

output_directory = LIVEDATA_DIR

REQUEST_TIMEOUT = 15  # seconds, per HTTP request
MAX_ATTEMPTS = 3      # attempts per source before giving up
RETRY_BACKOFF = 2     # seconds, multiplied by attempt number
MAX_CACHE_AGE = timedelta(hours=1)  # reuse an on-disk file instead of redownloading if younger than this

def _fetch_with_retry(fetch_fn, url, label):
    """Retries a transient failure a few times; only prints/gives up after exhausting attempts."""
    last_error = None
    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            return fetch_fn(url)
        except requests.exceptions.RequestException as e:
            last_error = e
            if attempt < MAX_ATTEMPTS:
                time.sleep(RETRY_BACKOFF * attempt)
    print(f"Warning: failed to fetch live {label} data after {MAX_ATTEMPTS} attempts ({last_error})")
    return None

def fetch_data_Space(urlSpace: str) -> dict:
    response = requests.get(urlSpace, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json()

def fetch_data_Mag(urlMag: str) -> dict:
    response = requests.get(urlMag, timeout=REQUEST_TIMEOUT)
    response.raise_for_status()
    return response.json()

def fetch_data_Dst(urlDst: str) -> str:
    responseDst = requests.get(urlDst, timeout=REQUEST_TIMEOUT)
    responseDst.raise_for_status()
    return responseDst.text

def fetch_data_Kp(urlKp: str) -> str:
    responseKp = requests.get(urlKp, timeout=REQUEST_TIMEOUT)
    responseKp.raise_for_status()
    return responseKp.text

def save_to_json_space(data: dict) -> None:
    json_file = os.path.join(output_directory, f'space_data.json')
    with open(json_file, 'w') as file:
        json.dump(data, file, indent=4)

def save_to_json_magnetic(data: dict) -> None:
    json_file = os.path.join(output_directory, f'Magnetic_data.json')
    with open(json_file, 'w') as file:
        json.dump(data, file, indent=4)

def save_to_txt_Dst(data: str, dst_file: str) -> None:
    with open(dst_file, "w") as fileDst:
            fileDst.write(data)

def save_to_txt_Kp(data: str) -> None:
    Kp_file = os.path.join(output_directory, f'Kp_data.txt')
    with open(Kp_file, "w") as fileKp:
            fileKp.write(data)

def _active_rows(data: list) -> list:
    active = [row for row in data if row.get('active')]
    return active if active else data

def save_to_csv_Space(data: dict) -> None:
    csv_file = os.path.join(output_directory, f'space_data.csv')
    try:
        if isinstance(data, list) and all(isinstance(row, dict) for row in data):
            df = pd.DataFrame(_active_rows(data))
            df = df.rename(columns={'proton_speed': 'speed', 'proton_density': 'density'})
            df = df[['time_tag', 'speed', 'density']]
            df['fetch_timestamp'] = pd.Timestamp.now()
            df.to_csv(csv_file, index=False)

    except Exception as e:
        print(f"Error converting JSON to CSV: {e}")

def save_to_csv_Mag(data: dict) -> None:
    csv_file = os.path.join(output_directory, f'Magnetic_data.csv')
    try:
        if isinstance(data, list) and all(isinstance(row, dict) for row in data):
            df = pd.DataFrame(_active_rows(data))
            df = df[['time_tag', 'bx_gsm', 'by_gsm', 'bz_gsm']]
            df['fetch_timestamp'] = pd.Timestamp.now()
            df.to_csv(csv_file, index=False)

    except Exception as e:
        print(f"Error converting JSON to CSV: {e}")

def _is_fresh(path: str) -> bool:
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    age = time.time() - os.path.getmtime(path)
    return age < MAX_CACHE_AGE.total_seconds()

def Get_Data(current_time: datetime) -> tuple:
    script_dir = LIVEDATA_DIR
    os.makedirs(script_dir, exist_ok=True)

    space_csv = os.path.join(script_dir, 'space_data.csv')
    mag_csv = os.path.join(script_dir, 'Magnetic_data.csv')
    kp_txt = os.path.join(script_dir, 'Kp_data.txt')

    current_month = f"{current_time.month:02d}"
    current_year = current_time.year
    curretnt_year_two_digits = current_year % 100

    dststring1 = str(current_year) + str(current_month)
    dststring2 = str(curretnt_year_two_digits) + str(current_month)

    dst_txt = os.path.join(script_dir, f'Dst_data_{dststring1}.txt')

    urlSpace = 'https://services.swpc.noaa.gov/json/rtsw/rtsw_wind_1m.json'
    urlMag = 'https://services.swpc.noaa.gov/json/rtsw/rtsw_mag_1m.json'
    urlDst = 'https://wdc.kugi.kyoto-u.ac.jp/dst_realtime/' + dststring1 + '/dst' + dststring2 + '.for.request'
    urlKp = 'https://kp.gfz-potsdam.de/app/files/Kp_ap_nowcast.txt'

    if not (_is_fresh(space_csv) and _is_fresh(mag_csv)):
        data_Space = _fetch_with_retry(fetch_data_Space, urlSpace, "solar wind")
        if data_Space is not None:
            save_to_json_space(data_Space)
            save_to_csv_Space(data_Space)
        elif not os.path.exists(space_csv):
            raise RuntimeError("Live solar wind data could not be downloaded and no cached copy is available.")
        else:
            print("Warning: reusing stale cached solar wind data; live fetch failed.")

        data_Mag = _fetch_with_retry(fetch_data_Mag, urlMag, "magnetic field")
        if data_Mag is not None:
            save_to_json_magnetic(data_Mag)
            save_to_csv_Mag(data_Mag)
        elif not os.path.exists(mag_csv):
            raise RuntimeError("Live magnetic field data could not be downloaded and no cached copy is available.")
        else:
            print("Warning: reusing stale cached magnetic field data; live fetch failed.")

    if not _is_fresh(dst_txt):
        data_Dst = _fetch_with_retry(fetch_data_Dst, urlDst, "Dst index")
        if data_Dst is not None:
            save_to_txt_Dst(data_Dst, dst_txt)
        elif not os.path.exists(dst_txt):
            raise RuntimeError("Live Dst index data could not be downloaded and no cached copy is available.")
        else:
            print("Warning: reusing stale cached Dst index data; live fetch failed.")

    if not _is_fresh(kp_txt):
        data_Kp = _fetch_with_retry(fetch_data_Kp, urlKp, "Kp index")
        if data_Kp is not None:
            save_to_txt_Kp(data_Kp)
        elif not os.path.exists(kp_txt):
            raise RuntimeError("Live Kp index data could not be downloaded and no cached copy is available.")
        else:
            print("Warning: reusing stale cached Kp index data; live fetch failed.")

    Dst, daily_dst = extract_dst_value(dst_txt, current_time)
    Speed, Density = extract_solar_wind(current_time)
    speed_avg, density_avg, By_avg, Bz_avg, Bx_avg, N_norm, B_norm = extract_30min_average(current_time)
    By, Bz = extract_magnetic(current_time)
    Kp, IOPT = extract_kp_index(current_time)
    G1, G2, G3 = TSY01_Constants(By, Bz, -1 * Speed, Density)
    Speed = round(Speed, 3)
    Density = round(Density, 3)
    By = round(By, 3)
    Bz = round(Bz, 3)
    G1 = round(G1, 3)
    G2 = round(G2, 3)
    G3 = round(G3, 3)

    return Dst, Speed, Density, By, Bz, IOPT, G1, G2, G3, Kp, Bx_avg, By_avg, Bz_avg, N_norm, B_norm
