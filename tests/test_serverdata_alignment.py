import numpy as np
import pandas as pd
import pytest

from OTSO._core.serverdata import omni_data_pull, soho_data_pull


def _write_inputs(directory, year):
    """Two hours of 5-minute OMNI (with a one-hour gap) and three hourly OMNI records."""
    high_res = pd.DataFrame({"Date": pd.date_range(f"{year}-12-31 21:00", f"{year}-12-31 23:55", freq="5min")})
    high_res = high_res[high_res["Date"].dt.hour != 22]  # 22 UT falls back to hourly
    for col, value in {"Bx": 1.0, "By": 2.0, "Bz": -3.0, "V": 400.0, "Density": 5.0, "Pdyn": 1.5,
                       "By_avg": 2.0, "Bz_avg": -3.0, "G1": 1.0, "G2": 1.0, "G3": 1.0}.items():
        high_res[col] = value
    high_res.to_csv(directory / f"omni_{year}_high_res.csv", index=False)

    low_res = pd.DataFrame({
        "Date": pd.to_datetime([f"{year}-12-31 21:00", f"{year}-12-31 22:00", f"{year}-12-31 23:00"]),
        "Kp": [1, 2, 3], "Kp_raw": [0.7, 2.0, 3.3], "Dst": [-10, -20, -30],
        "Bx": 0.5, "By": 1.5, "Bz": -2.5, "V": 410.0, "Density": 4.0, "Pdyn": 1.2,
        "G1": 2.0, "G2": 2.0, "G3": 2.0,
    })
    low_res.to_csv(directory / f"omni_{year}_low_res.csv", index=False)
    return low_res.set_index("Date")


def _combine(tmp_path, monkeypatch, year=2024):
    monkeypatch.setattr(omni_data_pull, "SCRATCH_DIR", str(tmp_path))
    monkeypatch.setattr(omni_data_pull, "load_prior_year_tail", lambda year, hours: None)
    monkeypatch.setattr(soho_data_pull, "celias_url_exists", lambda year: False)
    hourly = _write_inputs(tmp_path, year)
    omni_data_pull.Combine(f"omni_{year}_high_res.csv", f"omni_{year}_low_res.csv",
                           f"TSY15_{year}.csv", year)
    out = pd.read_csv(tmp_path / f"{year}_TSY_Inputs.csv", parse_dates=["Date"])
    return out, hourly


def test_hourly_dst_and_kp_cover_their_own_hour(tmp_path, monkeypatch):
    """OMNI's hourly record stamped HH:00 covers HH:00-HH:59, so every 5-minute row must
    carry the Dst/Kp of its own hour - not the previous hour - including the final hour."""
    year = 2024
    out, hourly = _combine(tmp_path, monkeypatch, year)
    assert out["Date"].iloc[-1] == pd.Timestamp(f"{year}-12-31 23:55")
    expected = out["Date"].dt.floor("h").map(hourly["Dst"])
    assert (out["Dst"] == expected).all(), out.loc[out["Dst"] != expected, ["Date", "Dst"]]
    expected_kp = out["Date"].dt.floor("h").map(hourly["Kp"])
    assert (out["Kp"] == expected_kp).all()

    # the hourly fallback hour takes its solar wind from the same hourly record
    gap = out[out["Date"].dt.hour == 22]
    assert (gap["Source"] == 2).all() and (gap["V"] == 410.0).all()


def test_g_parameters_are_trailing_hour_means(tmp_path, monkeypatch):
    """G2 = <V Bs>/200 and G3 = <N V Bs>/2000 are means over the preceding hour
    (Tsyganenko 2002b; Tsyganenko, Singer & Kasper 2003), on 5-minute and hourly rows alike."""
    out, _ = _combine(tmp_path, monkeypatch)
    out = out.set_index("Date")
    G1_5min = 400 * ((13 ** 0.5 / 40) ** 2 / (1 + 13 ** 0.5 / 40)) * np.sin(np.arctan2(2, -3) / 2) ** 3

    # 21:55 - a full hour of 5-minute rows (By 2, Bz -3, V 400, N 5)
    full = out.loc["2024-12-31 21:55"]
    assert full["G1"] == pytest.approx(G1_5min, abs=0.01)
    assert full["G2"] == pytest.approx(400 * 3 / 200, abs=0.01)
    assert full["G3"] == pytest.approx(5 * 400 * 3 / 2000, abs=0.01)

    # 22:55 - a full hour of hourly fallback rows (Bz -2.5, V 410, N 4)
    hourly = out.loc["2024-12-31 22:55"]
    assert hourly["G2"] == pytest.approx(410 * 2.5 / 200, abs=0.01)
    assert hourly["G3"] == pytest.approx(4 * 410 * 2.5 / 2000, abs=0.01)

    # no jump at the 5-minute -> hourly switch: G moves gradually over the window
    switch = out.loc["2024-12-31 21:55":"2024-12-31 22:55", "G2"]
    assert switch.is_monotonic_decreasing
    assert switch.diff().abs().max() < 0.1
