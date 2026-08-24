import pandas as pd
import numpy as np
from tqdm import tqdm
import os
from numba import njit

from .scratch_dir import SCRATCH_DIR
from .prior_year import load_prior_year_tail

def TSY04_param_generator(YEAR: int) -> None:

    param_path = os.path.join(os.path.dirname(__file__), 'parameters.csv')
    params_df = pd.read_csv(param_path, delimiter='\t')
    params_df.columns = [c.strip() for c in params_df.columns]
    params = {str(row['Parameter']).strip(): row['Value'] for _, row in params_df.iterrows()}

    DT = [params[f'DT{i+1}']/60 for i in range(6)]
    Bs = [params[f'Bs{i+1}'] for i in range(6)]
    FAC_lam = [params[f'FAC{i+1}_gam'] for i in range(6)]
    FAC_beta = [params[f'FAC{i+1}_beta'] for i in range(6)]

    tsy_inputs_path = os.path.join(SCRATCH_DIR, f'{YEAR}_TSY_Inputs.csv')
    df = pd.read_csv(tsy_inputs_path)
    df['SYM_H'] = pd.to_numeric(df['SYM_H'], errors='coerce')

    raw_dt_per_hour = [params[f'DT{i+1}'] for i in range(6)]
    warmup_hours = min(30 * 24, (10 / min(raw_dt_per_hour)) * 1.1)
    prior_tail = load_prior_year_tail(YEAR, hours=warmup_hours)
    n_warmup = 0
    if prior_tail is not None and not prior_tail.empty:
        warmup_df = pd.DataFrame({
            'Date': prior_tail['Date'].dt.strftime('%Y-%m-%d %H:%M:%S'),
            'V': prior_tail['V'],
            'Density': prior_tail['Density'],
            'Bz': prior_tail['Bz'],
            'Pdyn': prior_tail['Pdyn'],
            'SYM_H': pd.to_numeric(prior_tail['SYM_H'], errors='coerce'),
        })
        n_warmup = len(warmup_df)
        df = pd.concat([warmup_df, df], ignore_index=True)

    window_size = 24
    valid_intervals = []
    in_interval = False
    interval_start = None
    for idx in range(window_size, len(df)):
        window_symh = df['SYM_H'].iloc[idx-window_size:idx]
        window_bz = df['Bz'].iloc[idx-window_size:idx]
        if window_symh.isnull().any() or window_bz.isnull().any():
            continue
        if not (window_bz >= 0).all():
            continue
        avg_symh = window_symh.mean()
        v_invalid = 'V' in df.columns and (df['V'].iloc[idx] == 9999.0 or np.isnan(df['V'].iloc[idx]))
        density_invalid = 'Density' in df.columns and (df['Density'].iloc[idx] == 999.9 or np.isnan(df['Density'].iloc[idx]))
        pdyn_invalid = 'Pdyn' in df.columns and (df['Pdyn'].iloc[idx] == 99.99 or np.isnan(df['Pdyn'].iloc[idx]))
        if in_interval and (v_invalid or density_invalid or pdyn_invalid):
            valid_intervals.append((interval_start, idx-1))
            in_interval = False
            interval_start = None
            continue
        if not in_interval and (-20 <= avg_symh <= 5):
            if not (v_invalid or density_invalid or pdyn_invalid):
                in_interval = True
                interval_start = idx
        if in_interval:
            if pd.isnull(df['SYM_H'].iloc[idx]):
                valid_intervals.append((interval_start, idx-1))
                in_interval = False
                interval_start = None
    if in_interval:
        valid_intervals.append((interval_start, len(df)-1))
    
    V = df['V'].values
    Density = df['Density'].values
    Bz = df['Bz'].values
    Date = df['Date'].values

    @njit
    def compute_Ws(V, Density, Bz, DT, Bs, FAC_lam, FAC_beta, start_idx, end_idx):
        n = end_idx - start_idx + 1
        Wmat = np.zeros((n, 6))
        for ind_pos in range(n):
            W = np.zeros(6)
            KEY = np.ones(6)
            for kk in range(ind_pos, -1, -1):
                vnorm = V[start_idx+kk] / 400.0
                dennorm = Density[start_idx+kk] * 1.16 / 5.0
                bsnorm = -Bz[start_idx+kk] / 5.0
                if bsnorm <= 0:
                    bs = np.zeros(6)
                else:
                    bs = np.array([bsnorm ** Bs[j] for j in range(6)])
                fac = np.array([
                    dennorm ** FAC_lam[j] * vnorm ** FAC_beta[j] * bs[j] for j in range(6)
                ])
                taumt = (ind_pos - kk) * 5
                arg = np.array([-taumt * DT[j] for j in range(6)])
                for n in range(6):
                    if arg[n] > -10 and KEY[n] == 1:
                        W[n] += fac[n] * np.exp(arg[n])
                    else:
                        KEY[n] = 0
                if np.all(KEY == 0):
                    break
            for n in range(6):
                W[n] = W[n] * DT[n] * 5
            Wmat[ind_pos, :] = W
        return Wmat
    
    results = []
    for start_idx, end_idx in valid_intervals:
        for i, ind in enumerate(range(start_idx, end_idx+1)):
            if i == 0:
                Wmat = compute_Ws(V, Density, Bz, np.array(DT), np.array(Bs), np.array(FAC_lam), np.array(FAC_beta), start_idx, end_idx)
            if ind < n_warmup:
                continue
            results.append({
                'Date': Date[ind],
                'W1': round(Wmat[i,0], 2), 'W2': round(Wmat[i,1], 2), 'W3': round(Wmat[i,2], 2),
                'W4': round(Wmat[i,3], 2), 'W5': round(Wmat[i,4], 2), 'W6': round(Wmat[i,5], 2)
            })
            
    output_filename = os.path.join(SCRATCH_DIR, f'TSY04_W_parameters_{YEAR}.csv')
    pd.DataFrame(results).to_csv(output_filename, index=False)
