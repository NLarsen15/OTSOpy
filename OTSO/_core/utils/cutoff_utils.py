import pandas as pd
import numpy as np
import os

from ..libs.MiddleMan import Middleman as OTSOLib


def find_cutoff_range(df, rigidity_step):
    # Ensure sorted by rigidity
    df = df.sort_values("R [GV]").reset_index(drop=True)
    R_step = rigidity_step
    decimals = len(str(R_step).split(".")[1])


    # Lowest rigidity where A == 1
    R_low = df.loc[df["Filter"] == 1, "R [GV]"].min()

    # Highest forbidden rigidity (A < 0)
    R_last_forbidden = df.loc[df["Filter"] < 0, "R [GV]"].max()

    # First rigidity above the last forbidden value
    R_high = df.loc[df["R [GV]"] > R_last_forbidden, "R [GV]"].min()

    # Penumbra region (exclude R_low and R_high)
    penumbra = df[(df["R [GV]"] >= R_low) & (df["R [GV]"] < R_high)]

    # Number of allowed values in the penumbra
    count_allowed = (penumbra["Filter"] == 1).sum()

    # Total number of rigidity values in the penumbra
    count_total = len(penumbra)

    # Penumbral transparency
    transparency = count_allowed / count_total if count_total > 0 else float("nan")
    if np.isnan(transparency):
        transparency = 0.0

    # Effective cutoff
    R_eff = R_high - R_step * count_allowed

    R_eff = round(R_eff, decimals)
    R_low = round(R_low, decimals)
    R_high = round(R_high, decimals)
    transparency = round(transparency, decimals)


    return (
        R_low,
        R_high,
        transparency,
        R_eff
    )

def compute_number_of_rigidities(start_rigidity, end_rigidity, rigidity_step):

    #print(f"Computing number of rigidities with start: {start_rigidity}, end: {end_rigidity}, step: {rigidity_step}")
    
    n = int((start_rigidity - end_rigidity) / rigidity_step)

    return n

def Rigidity_scan(FortranData, g, h, CutoffDataInstance, margin=0.75):
    smallscan = False
    FortranData.rigiditystep = 0.25

    #print("scan")
    #print(
    #    "PID:", os.getpid(),
    #    "CPUs:", os.sched_getaffinity(0)
    #)

    n = compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    FortranData.gyropercent = CutoffDataInstance.maxsteppercent

    Rigidities = np.zeros(n, dtype=np.float64)
    Allowed = np.zeros(n, dtype=np.int32)
    
    # Call new Fortran interface

    OTSOLib.cutoff(
        FortranData,
        g,
        h,
        Rigidities,
        Allowed
    )

    df = pd.DataFrame({
    "R [GV]": Rigidities,
    "Filter": Allowed})

    if np.all(Allowed == 1):
        smallscan = smallest_step_check(FortranData, g, h, CutoffDataInstance)
        if smallscan:
           return FortranData, smallscan

  # Sort by rigidity
    df = df.sort_values(by="R [GV]").reset_index(drop=True)

    R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

    
    # Create result DataFrame with rigidities

    if np.isnan(R_low):
        R_low = 0
    if np.isnan(R_high):
        R_high = R_low

    R_low = R_low - margin
    if (R_low < 0):
       R_low = 0

    R_high = R_high + margin
    if (R_high > FortranData.startrigidity):
       R_high = FortranData.startrigidity

    FortranData.rigiditystep = CutoffDataInstance.rigiditystep
    FortranData.startrigidity = R_high
    FortranData.endrigidity = R_low

    n = compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n

    return FortranData, smallscan

def cutoff_range_arrays(R, A, rigidity_step):
    """find_cutoff_range on arrays (same rules and rounding, no sorting or DataFrame needed)."""
    decimals = len(str(rigidity_step).split(".")[1])
    allowed = R[A == 1]
    R_low = allowed.min() if allowed.size else np.nan
    forbidden = R[A < 0]
    R_last_forbidden = forbidden.max() if forbidden.size else np.nan
    above = R[R > R_last_forbidden]
    R_high = above.min() if above.size else np.nan
    penumbra = (R >= R_low) & (R < R_high)
    count_allowed = int((A[penumbra] == 1).sum())
    count_total = int(penumbra.sum())
    # np.float64 so round() matches find_cutoff_range (numpy rounding, not Python float rounding)
    transparency = np.float64(count_allowed) / count_total if count_total > 0 else 0.0
    R_eff = R_high - rigidity_step * count_allowed
    return (round(R_low, decimals), round(R_high, decimals), round(transparency, decimals),
            round(R_eff, decimals))

def scan_range(FortranData, g, h, Data, end_rigidity, start_rigidity, step=None):
    """Trace end_rigidity + step, ..., start_rigidity (default step: the rigidity step).
    Returns the rigidities and allowed flags."""
    step = Data.rigiditystep if step is None else step
    FortranData.rigiditystep = step
    FortranData.endrigidity = end_rigidity
    FortranData.startrigidity = start_rigidity
    FortranData.n = int(round((start_rigidity - end_rigidity) / step))
    FortranData.gyropercent = Data.maxsteppercent
    Rigidities = np.zeros(FortranData.n, dtype=np.float64)
    Allowed = np.zeros(FortranData.n, dtype=np.int32)
    OTSOLib.cutoff(FortranData, g, h, Rigidities, Allowed)
    return Rigidities, Allowed

def adaptive_cutoff(FortranData, g, h, Data, margin=0.25, guard_below=1.25, guard_above=0.75, chunk=0.25,
                    guard_scale=0.25, guard_min=0.5):
    """
    Penumbra scan with adaptive guard bands (rigidityscan="ADAPTIVE").

    The coarse scan's bracket, padded by `margin`, is scanned at the rigidity step. The window
    is then extended by `chunk` (also scanned at the rigidity step) on each side until it holds
    at least the lower guard of only forbidden rigidities below the lowest allowed one and
    `guard_above` GV of only allowed rigidities above the highest forbidden one, or reaches the
    rigidity limits. Allowed bands narrower than the 0.25 GV coarse step can lie well below the
    coarse bracket (up to ~2 GV in tests), so a fixed margin either misses them or wastes time
    on cells that have none.

    The lower guard scales with the lowest allowed rigidity found so far, Rl:
    min(guard_below, max(guard_min, guard_scale * Rl)), i.e. the full 1.25 GV from Rl = 5 GV up
    and 0.5 GV below 2 GV. Penumbral structure narrows with the cutoff, and trajectories below
    ~1 GV are the most expensive, so a fixed 1.25 GV spent most of the run time on low-cutoff
    cells. On 17 regular 5° maps of the May 2024 storm (TSY89_refit + Boberg) this is 16% faster
    with no cutoff changed; a 0.25 GV minimum was 22% faster but missed allowed bands 0.3-0.5 GV
    below Rl in 2 of 45,917 cells.
    """
    start, end = Data.startrigidity, Data.endrigidity

    # Coarse scan to bracket the penumbra
    R, A = scan_range(FortranData, g, h, Data, end, start, step=0.25)
    if np.all(A == 1):
        FortranData.startrigidity, FortranData.endrigidity = start, end
        if smallest_step_check(FortranData, g, h, Data):
            return [0.0, 0.0, 0.0, 0.0]
    R_low, R_high, _, _ = cutoff_range_arrays(R, A, 0.25)
    if np.isnan(R_low):
        R_low = 0
    if np.isnan(R_high):
        R_high = R_low
    lo, hi = max(R_low - margin, 0), min(R_high + margin, start)

    # Fine scan, extended until both guard bands are clean
    R, A = scan_range(FortranData, g, h, Data, lo, hi)
    Rs, As = [R], [A]
    while True:
        allowed, forbidden = R[A == 1], R[A != 1]
        guard = guard_below
        if allowed.size:
            guard = min(guard_below, max(guard_min, guard_scale * allowed.min()))
        grow_down = lo > end and allowed.size and allowed.min() - lo < guard - 1e-9
        grow_up = hi < start and forbidden.size and hi - forbidden.max() < guard_above - 1e-9
        if not (grow_down or grow_up):
            break
        if grow_down:
            new_lo = max(lo - chunk, end)
            r, a = scan_range(FortranData, g, h, Data, new_lo, lo)
            Rs.append(r), As.append(a)
            lo = new_lo
        if grow_up:
            new_hi = min(hi + chunk, start)
            r, a = scan_range(FortranData, g, h, Data, hi, new_hi)
            Rs.append(r), As.append(a)
            hi = new_hi
        R, A = np.concatenate(Rs), np.concatenate(As)

    R_low, R_high, transparency, R_eff = cutoff_range_arrays(R, A, Data.rigiditystep)
    return [R_high, R_eff, R_low, transparency]

def standard_cutoff(FortranData, g, h, Data):
    smallscan = False

    if Data.rigidityscan == "ADAPTIVE":
        return adaptive_cutoff(FortranData, g, h, Data)

    if Data.rigidityscan == "ON":
        FortranData, smallscan = Rigidity_scan(FortranData, g, h, Data)

    if smallscan:
            returnvalues = [0.0, 0.0, 0.0, 0.0]
            R_low, R_high, transparency, R_eff = returnvalues
            
    else:
            Rigidities = np.zeros(FortranData.n, dtype=np.float64)
            Allowed = np.zeros(FortranData.n, dtype=np.int32)
            FortranData.gyropercent = Data.maxsteppercent
            
            # Call new Fortran interface
            OTSOLib.cutoff(
                FortranData,
                g,
                h,
                Rigidities,
                Allowed
            )

            df = pd.DataFrame({
            "R [GV]": Rigidities,
            "Filter": Allowed})

            # Sort by rigidity
            df = df.sort_values(by="R [GV]").reset_index(drop=True)

            R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

    returnvalues = [R_high, R_eff, R_low, transparency]

    return returnvalues

def Apparent_cutoff(Positioninfo, FortranData, g, h, Data):

    Zenith_and_Azimuths = [
        [0, 0],        # vertical
        [30, 0],
        [30, 45],
        [30, 90],
        [30, 135],
        [30, 180],
        [30, 225],
        [30, 270],
        [30, 315]
    ]

    # Store results from each direction
    results = []

    for za in Zenith_and_Azimuths:
        smallscan = False
        Position = [Positioninfo[0], Positioninfo[1], Positioninfo[2], za[0], za[1]]
        FortranData.positionin = Position
        FortranData.startrigidity = Data.startrigidity
        FortranData.endrigidity = Data.endrigidity

        if Data.rigidityscan == "ADAPTIVE":
            R_high, R_eff, R_low, transparency = adaptive_cutoff(FortranData, g, h, Data)
            results.append({"Direction": tuple(za), "Ru": R_high, "Rc": R_eff, "Rl": R_low, "TF": transparency})
            continue

        if Data.rigidityscan == "ON":
            FortranData, smallscan = Rigidity_scan(FortranData, g, h,  Data)

        if smallscan:
            R_high, R_eff, R_low, transparency = [0.0, 0.0, 0.0, 0.0]
        else:
            n = compute_number_of_rigidities(
                FortranData.startrigidity,
                FortranData.endrigidity,
                FortranData.rigiditystep
            )

            FortranData.n = n
            FortranData.gyropercent = Data.maxsteppercent

            Rigidities = np.zeros(n, dtype=np.float64)
            Allowed = np.zeros(n, dtype=np.int32)

            OTSOLib.cutoff(
                FortranData,
                g,
                h,
                Rigidities,
                Allowed
            )

            df = pd.DataFrame({
                "R [GV]": Rigidities,
                "Filter": Allowed
            }).sort_values(by="R [GV]").reset_index(drop=True)

            R_low, R_high, transparency, R_eff = find_cutoff_range(df, FortranData.rigiditystep)

        results.append({
            "Direction": tuple(za),
            "Ru": R_high,
            "Rc": R_eff,
            "Rl": R_low,
            "TF": transparency
            })

    rigidity_dataframe = pd.DataFrame(results).set_index("Direction")

    # Weights: vertical = 1/2, each of the other 8 = 1/16
    weights = np.array([0.5] + [1/16] * 8)

    R_step = FortranData.rigiditystep
    decimals = len(str(R_step).split(".")[1])

    weighted_mean = (
        rigidity_dataframe[["Ru", "Rc", "Rl", "TF"]]
        .multiply(weights, axis=0)
        .sum(axis=0)
        .round(decimals)
    )

    returnvalues = [weighted_mean["Ru"], weighted_mean["Rc"], 
                    weighted_mean["Rl"], weighted_mean["TF"]]

    return returnvalues


def smallest_step_check(FortranData, g, h, CutoffDataInstance):
    n = 1

    FortranData.gyropercent = CutoffDataInstance.maxsteppercent

    FortranData.n = n
    FortranData.rigiditystep = CutoffDataInstance.rigidityarray[2] + CutoffDataInstance.rigidityarray[1]
    #FortranData.rigiditystep = 0.05

    Rigidities = np.zeros(n, dtype=np.float64)
    Allowed = np.zeros(n, dtype=np.int32)
    
    # Call new Fortran interface
    OTSOLib.cutoff(
        FortranData,
        g,
        h,
        Rigidities,
        Allowed
    )

    df = pd.DataFrame({
    "R [GV]": Rigidities,
    "Filter": Allowed})

    if Allowed[0] == 1:
       result = True # The lowest tested rigidity is allowed
       #print("Smallest Step Escaped")
    else:
       result = False # The lowest tested rigidity is not allowed
       #print("Smallest Step Did Not Escape")

    return result