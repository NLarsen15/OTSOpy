from . import MiddleMan as OTSOLib
from . import OmniPull
from . import TSY15_OTSO_Download_v2 as TSY15
import os
from datetime import datetime, timedelta
import pandas as pd
import numpy as np
import sys
import platform

def is_invalid_value(value):
    """Check if a value is invalid (NaN or OMNI filler value)"""
    if pd.isna(value) or np.isnan(value):
        return True
    # Check for common OMNI filler values
    invalid_values = [9999.0, 99999.0, 999.9, 9999.9, 9999.99, 99.99, 999.99, -999.9]
    try:
        float_value = float(value)
        return any(abs(float_value - invalid) < 0.01 for invalid in invalid_values)
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
    OMNIYEAR = int(Date.year)
    RoundedDate = round_to_nearest_five_minutes(Date)
    Date, Bx, By, Bz, V, Density, Pdyn, Kp, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, ByAvg, BzAvg, NIndex, BIndex, SymHCorrected, External = ExtractServerData(RoundedDate, External, AdaptiveExternalModel)

    return Bx, By, Bz, V, Density, Pdyn, Kp, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, ByAvg, BzAvg, NIndex, BIndex, SymHCorrected, External


def ExtractServerData(RoundedDate, External, AdaptiveExternalModel=False):
    year = RoundedDate.year
    source_file = f'{year}_TSY_Inputs.csv'
    TSY_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "ServerData")
    TSY_file = os.path.join(TSY_folder, os.path.basename(source_file))

    if not os.path.exists(TSY_file):
        print(f"File for year {year} not found in the directory.")
        return

    try:
        df = pd.read_csv(TSY_file)
    except Exception as e:
        print(f"Error reading {source_file}: {e}")
        return

    if 'Date' not in df.columns:
        print(f"'Date' column not found in {source_file}")
        return

    df['Date'] = pd.to_datetime(df['Date'], errors='coerce')

    matching_row = df[df['Date'] == RoundedDate]

    #print(RoundedDate)

    if matching_row.empty:
        print(f"No matching row found for datetime: {RoundedDate}")
        return None
    else:
        row_data = matching_row.iloc[0]  

        Date = row_data['Date'] if 'Date' in row_data else None
        Bx = row_data['Bx'] if 'Bx' in row_data else None
        By = row_data['By'] if 'By' in row_data else None
        Bz = row_data['Bz'] if 'Bz' in row_data else None
        V = row_data['V'] if 'V' in row_data else None
        Density = row_data['Density'] if 'Density' in row_data else None
        Pdyn = row_data['Pdyn'] if 'Pdyn' in row_data else None
        Kp = row_data['Kp'] if 'Kp' in row_data else None
        Dst = row_data['Dst'] if 'Dst' in row_data else None
        G1 = row_data['G1'] if 'G1' in row_data else None
        G2 = row_data['G2'] if 'G2' in row_data else None
        G3 = row_data['G3'] if 'G3' in row_data else None
        W1 = row_data['W1'] if 'W1' in row_data else None
        W2 = row_data['W2'] if 'W2' in row_data else None
        W3 = row_data['W3'] if 'W3' in row_data else None
        W4 = row_data['W4'] if 'W4' in row_data else None
        W5 = row_data['W5'] if 'W5' in row_data else None
        W6 = row_data['W6'] if 'W6' in row_data else None
        ByAvg = row_data['By_avg'] if 'By_avg' in row_data else None
        BzAvg = row_data['Bz_avg'] if 'Bz_avg' in row_data else None
        NIndex = row_data['N_index'] if 'N_index' in row_data else None
        BIndex = row_data['B_index'] if 'B_index' in row_data else None
        SymHCorrected = row_data['SYM_H'] if 'SYM_H' in row_data else None

        # Apply adaptive fallback logic if AdaptiveExternalModel is True
        if AdaptiveExternalModel:
            # Define all model checks in order from newest to oldest
            model_checks = {
                11: ("TA16_RBF", lambda: not (is_invalid_value(SymHCorrected) or is_invalid_value(ByAvg) or is_invalid_value(BzAvg) or is_invalid_value(NIndex))),
                10: ("TSY15B", lambda: not (is_invalid_value(BIndex) or is_invalid_value(ByAvg) or is_invalid_value(BzAvg))),
                9: ("TSY15N", lambda: not (is_invalid_value(NIndex) or is_invalid_value(ByAvg) or is_invalid_value(BzAvg))),
                7: ("TSY04", lambda: not (is_invalid_value(W1) or is_invalid_value(V) or is_invalid_value(Bz))),
                6: ("TSY01S", lambda: not (is_invalid_value(G3) or is_invalid_value(V) or is_invalid_value(Bz))),
                5: ("TSY01", lambda: not (is_invalid_value(G1) or is_invalid_value(G2) or is_invalid_value(V) or is_invalid_value(Bz))),
                4: ("TSY96", lambda: not (is_invalid_value(V) or is_invalid_value(Bz))),
                3: ("TSY89", lambda: not (is_invalid_value(Kp)))
            }
            
            # Build fallback order starting from the requested model and going backwards
            fallback_order = []
            for model_num in sorted([k for k in model_checks.keys() if k <= External], reverse=True):
                model_name, check_func = model_checks[model_num]
                fallback_order.append((model_num, model_name, check_func))
            
            # Try each model in the fallback order
            original_external = External
            for model_num, model_name, check_func in fallback_order:
                if check_func():
                    if model_num != original_external:
                        print(f"WARNING: Parameters not available for requested model. Falling back to {model_name}.")
                    External = model_num
                    break
            else:
                raise ValueError("ERROR: No valid magnetospheric model parameters found for given time.")
        
        # Standard validation when AdaptiveExternalModel is False
        else:
            if External == 7 and ((is_invalid_value(W1)) or (is_invalid_value(V)) or (is_invalid_value(Bz))):
                raise ValueError("ERROR: No TSY04 parameters found for given time.")
            elif External == 6 and ((is_invalid_value(G3)) or (is_invalid_value(V)) or (is_invalid_value(Bz))):
                raise ValueError("ERROR: No TSY01S parameters found for given time.")
            elif External == 5 and ((is_invalid_value(G1)) or (is_invalid_value(G2)) or (is_invalid_value(V)) or (is_invalid_value(Bz))):
                raise ValueError("ERROR: No TSY01 parameters found for given time.")
            elif External == 4 and ((is_invalid_value(V)) or (is_invalid_value(Bz))):
                raise ValueError("ERROR: No TSY96 parameters found for given time.")
            elif External == 9 and ((is_invalid_value(NIndex)) or (is_invalid_value(ByAvg)) or (is_invalid_value(BzAvg))):
                raise ValueError("ERROR: No TSY15N parameters found for given time.")
            elif External == 10 and ((is_invalid_value(BIndex)) or (is_invalid_value(ByAvg)) or (is_invalid_value(BzAvg))):
                raise ValueError("ERROR: No TSY15B parameters found for given time.")
            elif External == 11 and ((is_invalid_value(SymHCorrected)) or (is_invalid_value(ByAvg)) or (is_invalid_value(BzAvg)) or (is_invalid_value(NIndex))):
                raise ValueError("ERROR: No TA16 parameters found for given time.")
        
        if External == 9 or External == 11:
             if not is_invalid_value(NIndex) and NIndex > 2:
                 print(f'WARNING: N-INDEX ({NIndex}) OUT OF ALLOWED RANGE (0-2). SETTING TO MAX ALLOWED VALUE OF 2')
                 NIndex = 2

        if External == 10:
            if not is_invalid_value(BIndex) and BIndex > 2:
                print(f'WARNING: B-INDEX ({BIndex}) OUT OF ALLOWED RANGE (0-2). SETTING TO MAX ALLOWED VALUE OF 2')
                BIndex = 2

        if not is_invalid_value(V) and V > 0:
            V = -1*V

        # Convert invalid values to default values after model selection
        Bx = 0 if is_invalid_value(Bx) else Bx
        By = 5 if is_invalid_value(By) else By
        Bz = 5 if is_invalid_value(Bz) else Bz
        V = -500 if is_invalid_value(V) else V
        Density = 1 if is_invalid_value(Density) else Density
        Pdyn = 0 if is_invalid_value(Pdyn) else Pdyn
        Kp = 0 if is_invalid_value(Kp) else Kp
        Dst = 0 if is_invalid_value(Dst) else Dst
        G1 = 0 if is_invalid_value(G1) else G1
        G2 = 0 if is_invalid_value(G2) else G2
        G3 = 0 if is_invalid_value(G3) else G3
        W1 = 0 if is_invalid_value(W1) else W1
        W2 = 0 if is_invalid_value(W2) else W2
        W3 = 0 if is_invalid_value(W3) else W3
        W4 = 0 if is_invalid_value(W4) else W4
        W5 = 0 if is_invalid_value(W5) else W5
        W6 = 0 if is_invalid_value(W6) else W6
        ByAvg = 0 if is_invalid_value(ByAvg) else ByAvg
        BzAvg = 0 if is_invalid_value(BzAvg) else BzAvg
        NIndex = 0 if is_invalid_value(NIndex) else NIndex
        BIndex = 0 if is_invalid_value(BIndex) else BIndex
        SymHCorrected = 0 if is_invalid_value(SymHCorrected) else SymHCorrected

        return Date, Bx, By, Bz, V, Density, Pdyn, Kp, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, ByAvg, BzAvg, NIndex, BIndex, SymHCorrected, External


def DownloadServerFile(OMNIYEAR, g, h):
    source_file = f'{OMNIYEAR}_TSY_Inputs.csv'
    TSY_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "ServerData")
    TSY_file = os.path.join(TSY_folder, source_file)

    if not os.path.exists(TSY_file):
        print(f"Data for {OMNIYEAR} does not exist in OTSO files.")
        print(f'Attempting to download data for {OMNIYEAR}')
        TSY15.process_year(OMNIYEAR, g, h)
        OmniPull.PullOMNI(OMNIYEAR)
        DIRECTORY = os.path.join(os.path.dirname(os.path.dirname(__file__)),"functions")
        length = len(DIRECTORY)
        if platform.system() == "Windows":
            OTSOLib.gettsy04datawindows(OMNIYEAR, length)
            OmniPull.OMNI_to_csv(OMNIYEAR)
            OmniPull.TSY01(f'{OMNIYEAR}_TSY_Data.csv')
            OmniPull.TSY01(f'omni_{OMNIYEAR}_high_res.csv')
            OmniPull.TSY01(f'omni_{OMNIYEAR}_low_res.csv')
            OmniPull.Combine(f'{OMNIYEAR}_TSY_Data.csv', f'omni_{OMNIYEAR}_high_res.csv', f'omni_{OMNIYEAR}_low_res.csv',
                             f'TSY15_{OMNIYEAR}.csv', OMNIYEAR)
            OmniPull.Omnidelete(OMNIYEAR)
            print(f'Finished downloading data for {OMNIYEAR}.')
        elif platform.system() == "Linux"  or "Darwin":
            OTSOLib.gettsy04datalinux(OMNIYEAR, DIRECTORY, length)
            OmniPull.OMNI_to_csv(OMNIYEAR)
            OmniPull.TSY01(f'{OMNIYEAR}_TSY_Data.csv')
            OmniPull.TSY01(f'omni_{OMNIYEAR}_high_res.csv')
            OmniPull.TSY01(f'omni_{OMNIYEAR}_low_res.csv')
            OmniPull.Combine(f'{OMNIYEAR}_TSY_Data.csv', f'omni_{OMNIYEAR}_high_res.csv', f'omni_{OMNIYEAR}_low_res.csv',
                              f'TSY15_{OMNIYEAR}.csv', OMNIYEAR)
            OmniPull.Omnidelete(OMNIYEAR)
            print(f'Finished downloading data for {OMNIYEAR}.')

    return


def DownloadServerFileLowRes(OMNIYEAR):
    source_file = f'{OMNIYEAR}_TSY_Inputs.csv'
    TSY_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "ServerData")
    TSY_file = os.path.join(TSY_folder, os.path.basename(source_file))
    FilePath = os.path.join(os.path.dirname(os.path.realpath(__file__)), f'omni_{OMNIYEAR}_low_res.csv')

    if not os.path.exists(TSY_file):
        print(f"Data for {OMNIYEAR} does not exist in OTSO files.")
        print(f'Attempting to download data for {OMNIYEAR}')
        OmniPull.PullOMNILowRes(OMNIYEAR)
        OmniPull.CombineLowRes(FilePath,OMNIYEAR)
        OmniPull.OmnideleteLowRes(OMNIYEAR)
        print(f'Finished downloading data for {OMNIYEAR}.')


    