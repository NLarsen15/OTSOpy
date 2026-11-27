"""Regression tests for the 1.3.15 fixes, one test per fix."""
import datetime

import numpy as np

from OTSO import coordtrans, magfield
from OTSO._core.utils.tsy_params_utils import IOPTprocess_refit


def _refit_field(kp, dst, boberg=True):
    df = magfield(
        Locations=[[6.0, 0.0, 0.0], [-8.0, 2.0, 1.0], [0.0, 5.0, 3.0]],
        datetime_params={"year": 2024, "month": 5, "day": 1, "hour": 0},
        magfield_params={"internalmag": "IGRF", "externalmag": "TSY89_refit",
                         "boberg": boberg, "bobergtype": "CONTINUOUS"},
        geomagnetic_params={"Dst": dst, "kp": kp, "n_index": 0, "b_index": 0, "sym_h_corrected": 0},
        coordinate_params={"inputcoord": "GSM", "coordsystem": "GSM", "coordout": "GSM"},
        computation_params={"corenum": 1, "Verbose": False},
    )[0]
    return df.iloc[:, -3:].to_numpy(dtype=float)


def _geo_to_gse(date):
    df = coordtrans(Locations=[[1, 0, 0]], dates=[date], CoordIN="GEO", CoordOUT="GSE", corenum=1, Verbose=False)[0]
    return df.iloc[0, -3:].to_numpy(dtype=float)


def test_kp_rounds_to_nearest_third():
    """OMNI-style Kp (2.3 = 2+) must use the 2+ parameter set, not be truncated to 2."""
    cases = {0.0: 1, 0.3: 2, 0.7: 3, 2.0: 7, 2.3: 8, 2.333: 8, 2.7: 9, 5.3: 17, 6.7: 21, 7.0: 22, 9.0: 22}
    assert {kp: IOPTprocess_refit(kp) for kp in cases} == cases
    np.testing.assert_allclose(_refit_field(2.3, -20, boberg=False), _refit_field(7 / 3, -20, boberg=False), atol=1e-9)
    assert not np.allclose(_refit_field(2.3, -20, boberg=False), _refit_field(2.0, -20, boberg=False), atol=1e-6)


def test_boberg_dst_cap():
    """Dst above the Kp 0 mean (-2.57 nT) is capped at that mean; below it Dst still matters."""
    at_cap = _refit_field(1.0, -2.57)
    for dst in (0.0, 30.0):
        np.testing.assert_allclose(_refit_field(1.0, dst), at_cap, atol=1e-9, err_msg=f"Dst={dst}")
    assert not np.allclose(_refit_field(1.0, -10.0), _refit_field(1.0, -40.0), atol=1e-6)


def test_sun_position_outside_1900_2099():
    """Sun-dependent systems must work outside 1900-2099 (previously GSE came back equal to GEO)."""
    for year in (1590, 1850, 2150):
        v = _geo_to_gse(datetime.datetime(year, 6, 21, 12))
        assert np.isclose(np.linalg.norm(v), 1.0, atol=1e-3), (year, v)  # GEOPACK Sun vector drifts ~1e-4
        assert not np.allclose(v, [1.0, 0.0, 0.0], atol=1e-3), (year, v)
    a, b = _geo_to_gse(datetime.datetime(1899, 12, 31, 12)), _geo_to_gse(datetime.datetime(1900, 1, 1, 12))
    assert np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))) < 2.0
