import numpy as np
import pandas as pd

def Lshell(Trace: pd.DataFrame, CoordSystem: str) -> tuple[float, float]:
        
    Bx = Trace["Bx_GSM [nT]"]
    By = Trace["By_GSM [nT]"]
    Bz = Trace["Bz_GSM [nT]"]

    Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)
    idx_min = Bmag.idxmin()


    # Check for alt [km] column
    if "alt [km]" in Trace.columns:
        # Use altitude at min B, convert to Earth radii
        alt_km = Trace.loc[idx_min, "alt [km]"]
        L = (alt_km+6371.0) / 6371.0
    elif "radius [Re]" in Trace.columns:
        L = Trace.loc[idx_min, "radius [Re]"]
    else:
        Xmin = Trace.loc[idx_min, f"X_{CoordSystem} [Re]"]
        Ymin = Trace.loc[idx_min, f"Y_{CoordSystem} [Re]"]
        Zmin = Trace.loc[idx_min, f"Z_{CoordSystem} [Re]"]
        L = np.sqrt(Xmin**2 + Ymin**2 + Zmin**2)

    invlat = invariant_latitude_from_L(L)

    return L, invlat

def invariant_latitude_from_L(L: float) -> float:
    if L <= 1.0:
        return np.nan
    return np.degrees(np.arccos(1.0 / np.sqrt(L)))