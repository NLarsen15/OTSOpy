import os
import pandas as pd

SERVERDATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "ServerData")


def load_prior_year_tail(year: int, hours: float) -> pd.DataFrame | None:
    """
    Loads the last `hours` of the previous year's persisted {year-1}_TSY_Inputs.csv, if
    a cached copy exists.

    Several of the trailing-window computations that build a year's TSY input file need
    real history from before the year starts. Processing each year in
    isolation leaves the first stretch of every year computed from a truncated window.

    Returns None if no prior-year cache exists and reverts to original behavour.
    """
    prior_path = os.path.join(SERVERDATA_DIR, f'{year - 1}_TSY_Inputs.csv')
    if not os.path.exists(prior_path):
        return None
    df = pd.read_csv(prior_path, parse_dates=['Date'])
    if df.empty:
        return None
    cutoff = df['Date'].max() - pd.Timedelta(hours=hours)
    return df[df['Date'] > cutoff].reset_index(drop=True)
