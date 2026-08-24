import numpy as np
import chaosmagpy as cp
from datetime import datetime, timezone
import re
import os
import sys
from pathlib import Path


# ------------------------------------------------------------
# helper: datetime -> MJD2000
# ------------------------------------------------------------
def datetime_to_mjd2000(dt):
    """
    Convert datetime to MJD2000 (days since 2000-01-01 12:00:00 UTC)
    """
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)

    mjd2000_epoch = datetime(2000, 1, 1, 12, 0, 0, tzinfo=timezone.utc)

    delta = dt - mjd2000_epoch
    return delta.total_seconds() / 86400.0


def get_latest_chaos_model(data_path='./data'):
    """
    Finds the CHAOS model with the highest version number
    in the specified data directory.

    For example:
        CHAOS-8.5.mat
        CHAOS-8.6.mat
        CHAOS-8.10.mat

    Returns:
        Path to the highest-version CHAOS model.
    """

    pattern = re.compile(r'^CHAOS-(\d+(?:\.\d+)?)\.mat$')

    models = []

    for filename in os.listdir(data_path):
        match = pattern.match(filename)

        if match:
            version = tuple(map(int, match.group(1).split('.')))
            models.append((version, filename))

    if not models:
        raise FileNotFoundError(
            f"No CHAOS model files found in '{data_path}'."
        )

    # Sort by version number and take the largest
    _, latest_filename = max(models, key=lambda x: x[0])

    return os.path.join(data_path, latest_filename)



def get_chaos_g_h(dt, model_path=None, nmax=15):
    """
    Returns CHAOS MF g and h coefficients as flat Python lists.

    Ordering:
        g = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...]
        h = [(0,0), (1,0), (1,1), (2,0), (2,1), (2,2), ...]

    g(0,0) and all h(n,0) terms are set to 0.0.
    """

    data_path = Path(__file__).resolve().parent / "data"

    if model_path == None:
        model_path = get_latest_chaos_model(data_path)

    #print(f"Using CHAOS model: {model_path}")

    # load model
    model = cp.load_CHAOS_matfile(model_path)

    # convert datetime -> MJD2000
    time_mjd = np.array([datetime_to_mjd2000(dt)])

    # compute MF coefficients
    coeffs = model.synth_coeffs_tdep(
        time_mjd,
        nmax=nmax,
        deriv=0
    )[0]

    # Output lists
    g = [0.0]  # g(0,0)
    h = [0.0]  # h(0,0)

    idx = 0

    for n in range(1, nmax + 1):

        for m in range(0, n + 1):

            # g(n,m)
            g.append(float(coeffs[idx]))
            idx += 1

            # h(n,0) is zero
            if m == 0:
                h.append(0.0)

            # h(n,m), m >= 1
            else:
                h.append(float(coeffs[idx]))
                idx += 1

    return g, h


def check_chaos_exists():
    data_path = Path(__file__).resolve().parent / "data"
    if not data_path.exists():
        print("CHAOS .mat files are not installed. \n Please install using the OTSO.chaosdownload CLI command.")
        sys.exit(0)