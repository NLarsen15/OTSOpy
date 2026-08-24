import datetime

import pandas as pd
import pytest

from OTSO import (
    cone, coordtrans, cutoff, flight, magfield,
    planet, skymap, trace, trajectory, transmission,
)

COMPUTATION_FAST = {"corenum": 2, "threadnum": 1, "Verbose": False, "delim": ";"}

def assert_nonempty_df(df, name="dataframe"):
    assert isinstance(df, pd.DataFrame), f"{name} should be a DataFrame, got {type(df)}"
    assert not df.empty, f"{name} should not be empty"


def test_coordtrans():
    result = coordtrans(
        Locations=[[10, 10, 10]],
        dates=[datetime.datetime(2000, 10, 12, 8)],
        CoordIN="GEO",
        CoordOUT="GSM",
        corenum=1,
    )
    assert_nonempty_df(result[0], "coordtrans dataframe")
    assert isinstance(result[-1], str) and result[-1]


def test_magfield(date_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs):
    result = magfield(
        Locations=[[10, 10, 10]],
        coordinate_params={"inputcoord": "GDZ", "coordout": "GSM"},
        computation_params={"corenum": 1, "Verbose": False},
        datetime_params=date_inputs,
        magfield_params={"internalmag": "IGRF", "externalmag": "TSY89c", "boberg": False, "bobergtype": "EXTENSION"},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
    )
    assert_nonempty_df(result[0], "magfield dataframe")
    assert isinstance(result[-1], str) and result[-1]


def test_trace(date_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs):
    result = trace(
        coordsys="GEO",
        datetime_params=date_inputs,
        computation_params={"corenum": 2, "Verbose": False},
        magfield_params={
            "internalmag": "IGRF", "externalmag": "TSY89c",
            "boberg": False, "bobergtype": "EXTENSION",
            "magnetopause": "Sphere", "spheresize": 10,
        },
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        grid_params={"latstep": -30, "longstep": 180, "maxlat": 30, "minlat": -30, "maxlong": 360, "minlong": 0},
    )
    field_lines = result[0]
    assert isinstance(field_lines, dict) and field_lines
    for key, entry in field_lines.items():
        assert isinstance(entry, dict) and "Trace" in entry, f"trace field line {key} missing 'Trace'"
        assert_nonempty_df(entry["Trace"], f"trace field line {key}")
    assert isinstance(result[-1], str) and result[-1]


def test_trajectory(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                     custom_field_inputs, data_retrieval_inputs):
    result = trajectory(
        Stations=["OULU"],
        particle_params={"rigidity": 1, "Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params={"corenum": 1, "Verbose": False},
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        custom_field_params=custom_field_inputs,
        coordinate_params={"coordsystem": "GEO", "inputcoord": "GDZ"},
        data_retrieval_params=data_retrieval_inputs,
    )
    trajectories = result[0]
    assert isinstance(trajectories, list) and trajectories
    for entry in trajectories:
        assert isinstance(entry, dict) and "trajectory" in entry, "trajectory entry missing 'trajectory' key"
        assert_nonempty_df(entry["trajectory"], f"trajectory for {entry.get('station')}")
    assert isinstance(result[-1], str) and result[-1]


def test_cutoff(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                 custom_field_inputs, data_retrieval_inputs):
    result = cutoff(
        Stations=["OULU"],
        cutoff_comp="Vertical",
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 2, "endrigidity": 0, "rigiditystep": 0.1, "rigidityscan": "ON"},
        asymptotic_params={"asymptotic": "NO", "asymlevels": [1, 5, 10], "unit": "GeV"},
        transmission_params={"transmission": False, "transmissionRstep": 0.001, "transmissionsamples": 20},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"coordsystem": "GEO", "inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )
    assert_nonempty_df(result[0], "cutoff dataframe")
    assert isinstance(result[-1], str) and result[-1]


def test_cone(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
              custom_field_inputs, data_retrieval_inputs):
    result = cone(
        Stations=["OULU"],
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 1, "endrigidity": 0, "rigiditystep": 0.1},
        transmission_params={"transmission": False, "transmissionRstep": 0.001, "transmissionsamples": 20},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"coordsystem": "GSE", "inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )
    assert_nonempty_df(result[0], "cone asymptotic dataframe")
    assert_nonempty_df(result[1], "cone Ru/Rc/Rl dataframe")
    assert isinstance(result[-1], str) and result[-1]


def test_skymap(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                 custom_field_inputs, data_retrieval_inputs):
    result = skymap(
        Stations=["OULU"],
        skymap_params={"zenithstep": 30, "azimuthstep": 180, "maxzenith": 30, "minzenith": 0,
                        "maxazimuth": 360, "minazimuth": 0},
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 2, "endrigidity": 0, "rigiditystep": 0.1, "rigidityscan": "ON"},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES"},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )
    result_dict = result[0]
    assert isinstance(result_dict, dict) and result_dict
    for key, df in result_dict.items():
        assert_nonempty_df(df, f"skymap {key}")


def test_transmission(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                       custom_field_inputs, data_retrieval_inputs):
    result = transmission(
        Stations=["OULU"],
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 2, "endrigidity": 0, "rigiditystep": 0.1},
        transmission_params={"transmissionRstep": 0.01, "transmissionsamples": 5},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )
    assert_nonempty_df(result[0], "transmission dataframe")


def test_flight(magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                 custom_field_inputs, data_retrieval_inputs):
    dates = [datetime.datetime(2000, 10, 12, 8), datetime.datetime(2000, 10, 12, 9)]
    result = flight(
        latitudes=[65, 65],
        longitudes=[25, 25],
        dates=dates,
        altitudes=[30, 30],
        cutoff_comp="Vertical",
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 2, "endrigidity": 0, "rigiditystep": 0.1, "rigidityscan": "ON"},
        asymptotic_params={"asymptotic": "NO", "asymlevels": [1, 5, 10], "unit": "GeV"},
        transmission_params={"transmission": False, "transmissionRstep": 0.001, "transmissionsamples": 20},
        solar_wind_params={k: [v, v] for k, v in solar_wind_inputs.items()},
        geomagnetic_params={k: [v, v] for k, v in geomagnetic_inputs.items()},
        tsyganenko_params={k: [v, v] for k, v in tsyganenko_inputs.items()},
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"coordsystem": "GEO", "inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )
    assert_nonempty_df(result[0], "flight dataframe")
    assert isinstance(result[1], str) and result[1]


def test_planet(date_inputs, magfield_inputs, solar_wind_inputs, geomagnetic_inputs, tsyganenko_inputs,
                 custom_field_inputs, data_retrieval_inputs):
    result = planet(
        cutoff_comp="Vertical",
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params={"startrigidity": 2, "endrigidity": 0, "rigiditystep": 0.1, "rigidityscan": "ON"},
        asymptotic_params={"asymptotic": "NO", "asymlevels": [1, 5, 10], "unit": "GeV"},
        transmission_params={"transmission": False, "transmissionRstep": 0.001, "transmissionsamples": 20},
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params={"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                             "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                             "startaltitude": 20, "betaerror": 0.001,
                             "totalbetacheck": True, "adaptivestep": True, "maxsteps": 100000},
        particle_params={"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0},
        computation_params=COMPUTATION_FAST,
        coordinate_params={"coordsystem": "GEO", "inputcoord": "GDZ"},
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs,
        grid_params={"latstep": -30, "longstep": 120, "maxlat": 30, "minlat": -30, "maxlong": 240, "minlong": 0},
    )
    assert_nonempty_df(result[0], "planet cutoff grid dataframe")
    assert isinstance(result[-1], str) and result[-1]
