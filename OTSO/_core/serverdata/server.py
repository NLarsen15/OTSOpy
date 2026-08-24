import os
from datetime import timedelta
import pandas as pd
import numpy as np

from . import omni_data_pull
from . import tsy15_download
from . import soho_data_pull
from . import tsy04_param_generator
from .scratch_dir import SCRATCH_DIR

PACKAGE_DIR = os.path.dirname(os.path.realpath(__file__))
SERVERDATA_DIR = os.path.join(PACKAGE_DIR, "ServerData")

INVALID_OMNI_VALUES = [9999.0, 99999.0, 999.9, 9999.9, 9999.99, 99.99, 999.99, -999.9]

def is_invalid_value(value):
    """Check if a value is invalid (NaN or OMNI filler value)"""
    if pd.isna(value) or np.isnan(value):
        return True
    try:
        float_value = float(value)
        return any(abs(float_value - invalid) < 0.01 for invalid in INVALID_OMNI_VALUES)
    except (ValueError, TypeError):
        return True

def round_to_nearest_five_minutes(date):
    # If year is below 1981, round to nearest full hour
    if date.year < 1981:
        # Calculate minutes past the hour
        minutes = date.minute + date.second / 60 + date.microsecond / 60000000
        # If minutes < 30, round down; else round up
        if minutes < 30:
            return date.replace(minute=0, second=0, microsecond=0)
        else:
            rounded = date.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1)
            return rounded
    # Otherwise, round to nearest 5 minutes
    minutes = date.minute + date.second / 60 + date.microsecond / 60000000
    nearest = int(5 * round(minutes / 5))
    rounded = date.replace(minute=0, second=0, microsecond=0) + timedelta(minutes=nearest)
    return rounded

def GetServerData(Date, External, AdaptiveExternalModel=False):
    RoundedDate = round_to_nearest_five_minutes(Date)
    (Date, Bx, By, Bz, V, Density, Pdyn, Kp, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6,
     ByAvg, BzAvg, NIndex, BIndex, SymHCorrected, External, Kp_raw) = ExtractServerData(
        RoundedDate, External, AdaptiveExternalModel)

    return Bx, By, Bz, V, Density, Pdyn, Kp, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, ByAvg, BzAvg, NIndex, BIndex, SymHCorrected, External


def _tsy_file_path(year: int) -> str:
    return os.path.join(SERVERDATA_DIR, f'{year}_TSY_Inputs.csv')

def _load_year_dataframe(year: int) -> pd.DataFrame:
    """Loads and validates one year's TSY input CSV. Raises a clear, specific error if
    the file is missing/unreadable/malformed instead of returning None, which used to
    make callers crash several frames away with an opaque "cannot unpack NoneType" error."""
    TSY_file = _tsy_file_path(year)
    if not os.path.exists(TSY_file):
        raise FileNotFoundError(
            f"Server data for {year} not found at {TSY_file}. "
            f"Call DownloadServerFile({year}, g, h) first."
        )
    try:
        df = pd.read_csv(TSY_file)
    except Exception as e:
        raise RuntimeError(f"Error reading {TSY_file}: {e}") from e
    if 'Date' not in df.columns:
        raise ValueError(f"'Date' column not found in {TSY_file}")
    df['Date'] = pd.to_datetime(df['Date'], errors='coerce')
    return df

def _match_row(df: pd.DataFrame, RoundedDate) -> pd.Series:
    matching_row = df[df['Date'] == RoundedDate]
    if matching_row.empty:
        raise ValueError(f"No matching server data row found for datetime: {RoundedDate}")
    return matching_row.iloc[0]

# Magnetospheric model checks, newest to oldest. Each check tells whether the row has
# everything that model needs; used both for adaptive fallback (try each in turn) and
# for validating a single explicitly-requested model.
_MODEL_CHECKS = {
    11: ("TA16_RBF", lambda f: not (is_invalid_value(f['SymHCorrected']) or is_invalid_value(f['ByAvg']) or is_invalid_value(f['BzAvg']) or is_invalid_value(f['NIndex']))),
    10: ("TSY15B", lambda f: not (is_invalid_value(f['BIndex']) or is_invalid_value(f['ByAvg']) or is_invalid_value(f['BzAvg']))),
    9: ("TSY15N", lambda f: not (is_invalid_value(f['NIndex']) or is_invalid_value(f['ByAvg']) or is_invalid_value(f['BzAvg']))),
    7: ("TSY04", lambda f: not (is_invalid_value(f['W1']) or is_invalid_value(f['V']) or is_invalid_value(f['Bz']))),
    6: ("TSY01S", lambda f: not (is_invalid_value(f['G3']) or is_invalid_value(f['V']) or is_invalid_value(f['Bz']))),
    5: ("TSY01", lambda f: not (is_invalid_value(f['G1']) or is_invalid_value(f['G2']) or is_invalid_value(f['V']) or is_invalid_value(f['Bz']))),
    4: ("TSY96", lambda f: not (is_invalid_value(f['V']) or is_invalid_value(f['Bz']))),
    3: ("TSY89", lambda f: not is_invalid_value(f['Kp'])),
}

_STANDARD_VALIDATION_NAMES = {
    7: "TSY04",
    6: "TSY01S",
    5: "TSY01",
    4: "TSY96",
    9: "TSY15N",
    10: "TSY15B",
    11: "TA16_RBF",
}

def _process_row(row_data: pd.Series, External, AdaptiveExternalModel=False):
    fields = {
        'Bx': row_data['Bx'] if 'Bx' in row_data else None,
        'By': row_data['By'] if 'By' in row_data else None,
        'Bz': row_data['Bz'] if 'Bz' in row_data else None,
        'V': row_data['V'] if 'V' in row_data else None,
        'Density': row_data['Density'] if 'Density' in row_data else None,
        'Pdyn': row_data['Pdyn'] if 'Pdyn' in row_data else None,
        'Kp': row_data['Kp'] if 'Kp' in row_data else None,
        'Dst': row_data['Dst'] if 'Dst' in row_data else None,
        'G1': row_data['G1'] if 'G1' in row_data else None,
        'G2': row_data['G2'] if 'G2' in row_data else None,
        'G3': row_data['G3'] if 'G3' in row_data else None,
        'W1': row_data['W1'] if 'W1' in row_data else None,
        'W2': row_data['W2'] if 'W2' in row_data else None,
        'W3': row_data['W3'] if 'W3' in row_data else None,
        'W4': row_data['W4'] if 'W4' in row_data else None,
        'W5': row_data['W5'] if 'W5' in row_data else None,
        'W6': row_data['W6'] if 'W6' in row_data else None,
        'ByAvg': row_data['By_avg'] if 'By_avg' in row_data else None,
        'BzAvg': row_data['Bz_avg'] if 'Bz_avg' in row_data else None,
        'NIndex': row_data['N_index'] if 'N_index' in row_data else None,
        'BIndex': row_data['B_index'] if 'B_index' in row_data else None,
        'SymHCorrected': row_data['SYM_H'] if 'SYM_H' in row_data else None,
    }
    fields['Kp_raw'] = row_data['Kp_raw'] if 'Kp_raw' in row_data else fields['Kp']

    if AdaptiveExternalModel:
        fallback_order = sorted((k for k in _MODEL_CHECKS if k <= External), reverse=True)
        for model_num in fallback_order:
            model_name, check_func = _MODEL_CHECKS[model_num]
            if check_func(fields):
                if model_num != External:
                    print(f"WARNING: Parameters not available for requested model. Falling back to {model_name}.")
                External = model_num
                break
        else:
            raise ValueError("ERROR: No valid magnetospheric model parameters found for given time.")
    elif External in _STANDARD_VALIDATION_NAMES:
        _, check_func = _MODEL_CHECKS[External]
        if not check_func(fields):
            raise ValueError(f"ERROR: No {_STANDARD_VALIDATION_NAMES[External]} parameters found for given time.")

    if External in (9, 11) and not is_invalid_value(fields['NIndex']) and fields['NIndex'] > 2:
        print(f"WARNING: N-INDEX ({fields['NIndex']}) BEYOND TYPICAL RANGE (0-2). RESULTS MAY BE UNRELIABLE.")

    if External == 10 and not is_invalid_value(fields['BIndex']) and fields['BIndex'] > 2:
        print(f"WARNING: B-INDEX ({fields['BIndex']}) OUT OF ALLOWED RANGE (0-2). SETTING TO MAX ALLOWED VALUE OF 2")
        fields['BIndex'] = 2

    if not is_invalid_value(fields['V']) and fields['V'] > 0:
        fields['V'] = -1 * fields['V']

    defaults = {
        'Bx': 0, 'By': 5, 'Bz': 5, 'V': -500, 'Density': 1, 'Pdyn': 0,
        'Kp': 0, 'Kp_raw': 0, 'Dst': 0, 'G1': 0, 'G2': 0, 'G3': 0,
        'W1': 0, 'W2': 0, 'W3': 0, 'W4': 0, 'W5': 0, 'W6': 0,
        'ByAvg': 0, 'BzAvg': 0, 'NIndex': 0, 'BIndex': 0, 'SymHCorrected': 0,
    }
    for key, default in defaults.items():
        if is_invalid_value(fields[key]):
            fields[key] = default

    if External == 100:
        fields['Kp'] = fields['Kp_raw']

    return fields, External

def ExtractServerData(RoundedDate, External, AdaptiveExternalModel=False):
    df = _load_year_dataframe(RoundedDate.year)
    row_data = _match_row(df, RoundedDate)
    fields, External = _process_row(row_data, External, AdaptiveExternalModel)
    Date = row_data['Date'] if 'Date' in row_data else None

    return (Date, fields['Bx'], fields['By'], fields['Bz'], fields['V'], fields['Density'],
            fields['Pdyn'], fields['Kp'], fields['Dst'], fields['G1'], fields['G2'], fields['G3'],
            fields['W1'], fields['W2'], fields['W3'], fields['W4'], fields['W5'], fields['W6'],
            fields['ByAvg'], fields['BzAvg'], fields['NIndex'], fields['BIndex'],
            fields['SymHCorrected'], External, fields['Kp_raw'])


def DownloadServerFile(OMNIYEAR, g, h):
    TSY_file = _tsy_file_path(OMNIYEAR)

    if os.path.exists(TSY_file):
        return

    print(f"Data for {OMNIYEAR} does not exist in OTSO files.")
    print(f'Attempting to download data for {OMNIYEAR}')
    tsy15_download.process_year(OMNIYEAR, g, h)
    if soho_data_pull.celias_url_exists(OMNIYEAR) and OMNIYEAR >= 1996:
        soho_data_pull.download_and_unpack_celias(OMNIYEAR)
    omni_data_pull.PullOMNI(OMNIYEAR)
    omni_data_pull.TSY01(f'omni_{OMNIYEAR}_high_res.csv')
    omni_data_pull.TSY01(f'omni_{OMNIYEAR}_low_res.csv')
    omni_data_pull.Combine(f'omni_{OMNIYEAR}_high_res.csv', f'omni_{OMNIYEAR}_low_res.csv',
                     f'TSY15_{OMNIYEAR}.csv', OMNIYEAR)
    tsy04_param_generator.TSY04_param_generator(OMNIYEAR)
    omni_data_pull.CombineTSY04(f"TSY04_W_parameters_{OMNIYEAR}.csv", OMNIYEAR)
    omni_data_pull.move_to_server_data(f'{OMNIYEAR}_TSY_Inputs.csv')
    omni_data_pull.Omnidelete(OMNIYEAR)
    print(f'Finished downloading data for {OMNIYEAR}.')


def DownloadServerFileLowRes(OMNIYEAR):
    TSY_file = _tsy_file_path(OMNIYEAR)

    if os.path.exists(TSY_file):
        return

    FilePath = os.path.join(SCRATCH_DIR, f'omni_{OMNIYEAR}_low_res.csv')
    print(f"Data for {OMNIYEAR} does not exist in OTSO files.")
    print(f'Attempting to download data for {OMNIYEAR}')
    omni_data_pull.PullOMNILowRes(OMNIYEAR)
    omni_data_pull.CombineLowRes(FilePath, OMNIYEAR)
    omni_data_pull.OmnideleteLowRes(OMNIYEAR)
    print(f'Finished downloading data for {OMNIYEAR}.')


def GetServerDataFlight(Dates, External, AdaptiveExternalModel=False):
    unique_years = list(set(dt.year for dt in Dates))
    datadf = load_all_csvs(unique_years)

    serverdata = []
    for x in Dates:
        RoundedDate = round_to_nearest_five_minutes(x)
        data = ExtractServerFlightData(RoundedDate, External, datadf, AdaptiveExternalModel)
        serverdata.append(data)

    return serverdata


def load_all_csvs(unique_years):
    yeardfs = [_load_year_dataframe(year) for year in unique_years]
    return pd.concat(yeardfs, ignore_index=True)

def ExtractServerFlightData(RoundedDate, External, datadf, AdaptiveExternalModel=False):
    row_data = _match_row(datadf, RoundedDate)
    fields, External = _process_row(row_data, External, AdaptiveExternalModel)

    return [fields['Bx'], fields['By'], fields['Bz'], fields['V'], fields['Density'], fields['Pdyn'],
            fields['Kp'], fields['Dst'], fields['G1'], fields['G2'], fields['G3'],
            fields['W1'], fields['W2'], fields['W3'], fields['W4'], fields['W5'], fields['W6'],
            fields['ByAvg'], fields['BzAvg'], fields['NIndex'], fields['BIndex'],
            fields['SymHCorrected'], External]
