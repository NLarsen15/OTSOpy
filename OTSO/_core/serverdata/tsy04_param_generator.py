import pandas as pd
import numpy as np
import os
from numba import njit

from .scratch_dir import SCRATCH_DIR
from .prior_year import load_prior_year_tail


@njit(cache=True)
def compute_Ws(V, Density, Bz, DT, Bs, FAC_lam, FAC_beta, start_idx, end_idx):
    """TSY04 (TS05) W1-W6 over one interval: Tsyganenko's truncated sums

        W_k(i) = 5 DT_k * sum_m fac_k(i - m) * exp(-5 m DT_k)

    over m = 0, 1, ... while 5 m DT_k < 10 and i - m stays inside the interval; a NaN input
    in that window makes W_k NaN. Same definition as the former per-row re-summation, but as
    a sliding window (decay the running sum, add the new term, drop the one leaving the
    window): O(n) instead of O(n x 3762) for W3.
    """
    n = end_idx - start_idx + 1
    fac = np.empty((n, 6))
    for i in range(n):
        vnorm = V[start_idx + i] / 400.0
        dennorm = Density[start_idx + i] * 1.16 / 5.0
        bsnorm = -Bz[start_idx + i] / 5.0
        for j in range(6):
            if bsnorm <= 0:
                bs = 0.0
            else:
                bs = bsnorm ** Bs[j]
            fac[i, j] = dennorm ** FAC_lam[j] * vnorm ** FAC_beta[j] * bs

    Wmat = np.empty((n, 6))
    for j in range(6):
        # number of terms kept: the same test as the per-row sum, -(m*5)*DT > -10
        M = 0
        while -((M + 1) * 5) * DT[j] > -10:
            M += 1
        decay = np.exp(-5 * DT[j])
        drop_weight = np.exp(-((M + 1) * 5) * DT[j])
        running = 0.0
        n_nan = 0
        for i in range(n):
            new = fac[i, j]
            if np.isnan(new):
                n_nan += 1
                new = 0.0
            running = running * decay + new
            leaving = i - (M + 1)
            if leaving >= 0:
                old = fac[leaving, j]
                if np.isnan(old):
                    n_nan -= 1
                else:
                    running -= old * drop_weight
            if n_nan > 0:
                Wmat[i, j] = np.nan
            else:
                Wmat[i, j] = running * DT[j] * 5
    return Wmat


def find_intervals(df: pd.DataFrame, window_size: int = 24) -> list:
    """Integration intervals (start, end) by OTSOpy's rules. The 2-hour window checks use
    cumulative counts over numpy arrays rather than pandas slices in the loop; the state
    machine itself is unchanged.

    For each row idx, using rows idx-24 .. idx-1: if the window has a missing SYM-H or Bz, or
    any Bz < 0, the row is skipped. Otherwise an open interval ends before idx if V, Density
    or Pdyn is invalid at idx; a new one starts at idx if the window's mean SYM-H is within
    [-20, 5] and idx is valid; an open interval also ends before idx if SYM-H is missing there.
    """
    n = len(df)
    symh = pd.to_numeric(df['SYM_H'], errors='coerce').to_numpy(dtype=float)
    bz = pd.to_numeric(df['Bz'], errors='coerce').to_numpy(dtype=float)

    def invalid(col, fill):
        if col not in df.columns:
            return np.zeros(n, dtype=bool)
        values = pd.to_numeric(df[col], errors='coerce').to_numpy(dtype=float)
        return (values == fill) | np.isnan(values)

    row_invalid = invalid('V', 9999.0) | invalid('Density', 999.9) | invalid('Pdyn', 99.99)
    good = (~np.isnan(symh) & ~np.isnan(bz) & (bz >= 0)).astype(np.int64)
    good_count = np.concatenate([[0], np.cumsum(good)])
    symh_missing = np.isnan(symh)

    valid_intervals = []
    in_interval = False
    interval_start = None
    for idx in range(window_size, n):
        if good_count[idx] - good_count[idx - window_size] != window_size:
            continue
        if in_interval and row_invalid[idx]:
            valid_intervals.append((interval_start, idx - 1))
            in_interval = False
            interval_start = None
            continue
        if not in_interval and not row_invalid[idx]:
            avg_symh = symh[idx - window_size:idx].mean()
            if -20 <= avg_symh <= 5:
                in_interval = True
                interval_start = idx
        if in_interval and symh_missing[idx]:
            valid_intervals.append((interval_start, idx - 1))
            in_interval = False
            interval_start = None
    if in_interval:
        valid_intervals.append((interval_start, n - 1))
    return valid_intervals


def TSY04_param_generator(YEAR: int) -> None:

    param_path = os.path.join(os.path.dirname(__file__), 'parameters.csv')
    params_df = pd.read_csv(param_path, delimiter='\t')
    params_df.columns = [c.strip() for c in params_df.columns]
    params = dict(zip(params_df['Parameter'].astype(str).str.strip(), params_df['Value']))

    DT = np.array([params[f'DT{i+1}']/60 for i in range(6)])
    Bs = np.array([params[f'Bs{i+1}'] for i in range(6)])
    FAC_lam = np.array([params[f'FAC{i+1}_gam'] for i in range(6)])
    FAC_beta = np.array([params[f'FAC{i+1}_beta'] for i in range(6)])

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

    valid_intervals = find_intervals(df)

    V = df['V'].to_numpy(dtype=float)
    Density = df['Density'].to_numpy(dtype=float)
    Bz = df['Bz'].to_numpy(dtype=float)
    Date = df['Date'].values

    dates, blocks = [], []
    for start_idx, end_idx in valid_intervals:
        Wmat = compute_Ws(V, Density, Bz, DT, Bs, FAC_lam, FAC_beta, start_idx, end_idx)
        keep = np.arange(start_idx, end_idx + 1) >= n_warmup
        dates.append(Date[start_idx:end_idx + 1][keep])
        blocks.append(Wmat[keep])

    columns = ['W1', 'W2', 'W3', 'W4', 'W5', 'W6']
    if blocks:
        # round() on numpy floats, as before, is numpy's rounding
        out = pd.DataFrame(np.round(np.concatenate(blocks), 2), columns=columns)
        out.insert(0, 'Date', np.concatenate(dates))
    else:
        out = pd.DataFrame(columns=['Date'] + columns)

    output_filename = os.path.join(SCRATCH_DIR, f'TSY04_W_parameters_{YEAR}.csv')
    out.to_csv(output_filename, index=False)
