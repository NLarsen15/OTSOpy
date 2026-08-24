import pytest

DATE_INPUTS = {"year": 2020, "month": 1, "day": 1, "hour": 1, "minute": 0, "second": 0}

SOLAR_WIND_INPUTS = {
    "vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
    "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0,
}
GEOMAGNETIC_INPUTS = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0, "sym_h_corrected": 0}
TSYGANENKO_INPUTS = {
    "G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
    "W4": 0, "W5": 0, "W6": 0,
}
CUSTOM_FIELD_INPUTS = {"g": None, "h": None, "max_degree": 13, "MHDfile": None, "MHDcoordsys": None}

DATA_RETRIEVAL_INPUTS = {"serverdata": "OFF", "livedata": "OFF"}


@pytest.fixture
def date_inputs():
    return dict(DATE_INPUTS)


@pytest.fixture
def solar_wind_inputs():
    return dict(SOLAR_WIND_INPUTS)


@pytest.fixture
def geomagnetic_inputs():
    return dict(GEOMAGNETIC_INPUTS)


@pytest.fixture
def tsyganenko_inputs():
    return dict(TSYGANENKO_INPUTS)


@pytest.fixture
def custom_field_inputs():
    return dict(CUSTOM_FIELD_INPUTS)


@pytest.fixture
def data_retrieval_inputs():
    return dict(DATA_RETRIEVAL_INPUTS)


@pytest.fixture
def magfield_inputs():
    return {
        "internalmag": "IGRF", "externalmag": "TSY89c",
        "boberg": False, "bobergtype": "EXTENSION",
        "magnetopause": "Sphere", "spheresize": 10,
        "AdaptiveExternalModel": False,
    }
