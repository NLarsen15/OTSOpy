"""Regression tests for bugs fixed in OTSO 1.3.15.

Each test checks one fix, so a regression shows which fix broke:
  - Kp rounding: OMNI Kp values such as 2.3 (meaning 2+, i.e. 7/3) must land in the 2+ bin, not the 2 bin.
  - Boberg DST cap: Dst above the mean Dst for Kp 0 is capped at that mean.
  - Sun position before 1900: GEOPACK/IRBEM SUN no longer returns zeros outside 1900-2099.
"""
import datetime

import numpy as np
import pytest

from OTSO import coordtrans, magfield
from OTSO._core.utils.tsy_params_utils import IOPTprocess_refit

DSTMEAN_KP0 = -2.57  # nT, DSTMEAN(1) in TSY89_refit.f95

LOCATIONS = [[6.0, 0.0, 0.0], [-8.0, 2.0, 1.0], [0.0, 5.0, 3.0]]


def _refit_field(kp, dst, bobergtype="CONTINUOUS", boberg=True):
    df = magfield(
        Locations=LOCATIONS,
        datetime_params={"year": 2024, "month": 5, "day": 1, "hour": 0},
        magfield_params={"internalmag": "IGRF", "externalmag": "TSY89_refit",
                         "boberg": boberg, "bobergtype": bobergtype},
        geomagnetic_params={"Dst": dst, "kp": kp, "n_index": 0, "b_index": 0, "sym_h_corrected": 0},
        coordinate_params={"inputcoord": "GSM", "coordsystem": "GSM", "coordout": "GSM"},
        computation_params={"corenum": 1, "Verbose": False},
    )[0]
    return df.iloc[:, -3:].to_numpy(dtype=float)


# ---------------------------------------------------------------- Kp rounding

@pytest.mark.parametrize("kp, expected", [
    (0.0, 1), (0.3, 2), (0.333, 2), (0.7, 3), (0.667, 3), (1.0, 4),
    (2.0, 7), (2.3, 8), (2.333, 8), (2.7, 9), (2.667, 9),
    (5.3, 17), (6.7, 21), (7.0, 22), (9.0, 22),
])
def test_iopt_refit_rounds_to_nearest_third(kp, expected):
    assert IOPTprocess_refit(kp) == expected


def test_refit_kp_plus_uses_plus_bin():
    """Kp 2.3 and 2.333 are the same bin (2+) and must differ from Kp 2.0 (the old truncation bug)."""
    b_20 = _refit_field(2.0, -20, boberg=False)
    b_23 = _refit_field(2.3, -20, boberg=False)
    b_2333 = _refit_field(7 / 3, -20, boberg=False)
    np.testing.assert_allclose(b_23, b_2333, rtol=0, atol=1e-9)
    assert not np.allclose(b_23, b_20, rtol=0, atol=1e-6)


# ---------------------------------------------------------------- Boberg DST cap

def test_boberg_dst_capped_at_kp0_mean():
    """Any Dst above DSTMEAN(Kp 0) must give the same field as Dst = DSTMEAN(Kp 0)."""
    at_cap = _refit_field(1.0, DSTMEAN_KP0)
    for dst in (0.0, 5.0, 30.0):
        np.testing.assert_allclose(_refit_field(1.0, dst), at_cap, rtol=0, atol=1e-9, err_msg=f"Dst={dst}")


def test_boberg_below_cap_still_depends_on_dst():
    assert not np.allclose(_refit_field(1.0, -10.0), _refit_field(1.0, -40.0), rtol=0, atol=1e-6)


# ---------------------------------------------------------------- Sun before 1900

def _geo_to_gse(date):
    df = coordtrans(Locations=[[1, 0, 0]], dates=[date], CoordIN="GEO", CoordOUT="GSE", corenum=1, Verbose=False)[0]
    return df.iloc[0, -3:].to_numpy(dtype=float)


@pytest.mark.parametrize("year", [1590, 1850, 1899, 1900, 2000, 2099, 2150])
def test_geo_to_gse_is_a_rotation(year):
    v = _geo_to_gse(datetime.datetime(year, 6, 21, 12))
    # The GEOPACK Sun vector is an approximation, so |v| drifts by ~1e-4 over 1600-2150.
    assert np.isclose(np.linalg.norm(v), 1.0, atol=1e-3), v
    # Before the fix, years outside 1900-2099 left the Sun direction at zero and GSE equalled GEO.
    assert not np.allclose(v, [1.0, 0.0, 0.0], atol=1e-3), v


def test_geo_to_gse_continuous_across_1900():
    """One day apart either side of 1900-01-01 the transform should change by well under a degree."""
    a = _geo_to_gse(datetime.datetime(1899, 12, 31, 12))
    b = _geo_to_gse(datetime.datetime(1900, 1, 1, 12))
    assert np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))) < 2.0
