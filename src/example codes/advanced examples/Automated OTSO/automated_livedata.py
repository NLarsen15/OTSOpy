import csv
import os
import time
import requests

from datetime import datetime, timedelta, timezone

from OTSO import cutoff


# ============================================================
# CONFIGURATION
# ============================================================
# Live data (OMNI/NOAA real-time feeds) lags behind the current time before
# it is published, so each computation targets the current UTC hour minus
# this offset rather than "now".
LIVE_DATA_LAG_HOURS = 2

FALLBACK_STATION = "OULU"

CSV_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "current_location_cutoff_log.csv")


magfield_inputs = {
    "internalmag": "IGRF",
    "externalmag": "TSY01S",   # Must not be TSY04 or TA16_RBF: unsupported with live data
    "boberg": False,
    "bobergtype": "EXTENSION"
}

rigidity_inputs = {
    "startrigidity": 20,
    "endrigidity": 0,
    "rigiditystep": 0.01,
    "rigidityscan": "ON"
}

integration_inputs = {
    "intmodel": "Boris-Buneman",
    "gyropercent": 15,
    "minaltitude": 20,
    "maxdistance": 100,
    "maxtime": 0,
    "mintrapdist": 6.6,
    "startaltitude": 20,
    "betaerror": 0.001,
    "totalbetacheck": True,
    "adaptivestep": False,
    "maxsteps": 100000
}

particle_inputs = {
    "Anum": 1,
    "anti": "YES",
    "zenith": 0,
    "azimuth": 0
}

computation_inputs = {
    "corenum": 1,
    "threadnum": 4,
    "Verbose": False,
    "delim": ";"
}

coord_inputs = {
    "coordsystem": "GEO",
    "inputcoord": "GDZ"
}

# Solar wind / geomagnetic / Tsyganenko inputs are ignored once livedata is
# "ON" (OTSO pulls real-time values instead), so the defaults are left as-is.
data_retrieval_inputs = {
    "serverdata": "OFF",
    "livedata": "ON"
}


# ============================================================
# LOCATION DETECTION
# ============================================================

def _sanitise_name(name):
    cleaned = "".join(c if c.isalnum() else "_" for c in name).strip("_")
    return cleaned or "MyLocation"


def _probe_ipinfo():
    response = requests.get("https://ipinfo.io/json", timeout=5)
    response.raise_for_status()
    data = response.json()

    lat_str, lon_str = data["loc"].split(",")
    city = data.get("city") or "MyLocation"

    return _sanitise_name(city), float(lat_str), float(lon_str)


def _probe_ip_api_com():
    response = requests.get("http://ip-api.com/json/", timeout=5)
    response.raise_for_status()
    data = response.json()

    if data.get("status") != "success":
        raise ValueError(data.get("message", "unknown error"))

    city = data.get("city") or "MyLocation"

    return _sanitise_name(city), float(data["lat"]), float(data["lon"])


def detect_location():
    """Detects an approximate location from the current IP address.

    IP-based geolocation only gives the location of the network's exit
    point (e.g. the nearest ISP hub), not a precise physical position.
    Returns (name, lat, lon), or None if every lookup fails.
    """
    for probe in (_probe_ipinfo, _probe_ip_api_com):
        try:
            return probe()
        except Exception as e:
            print(f"Location probe '{probe.__name__}' failed: {e}")

    return None


# ============================================================
# CUTOFF COMPUTATION
# ============================================================

def compute_cutoff(target_time_utc, station_name, custom_location):
    date_inputs = {
        "year": target_time_utc.year,
        "month": target_time_utc.month,
        "day": target_time_utc.day,
        "hour": target_time_utc.hour,
        "minute": 0,
        "second": 0
    }

    stations_list = [] if custom_location else [station_name]
    custom_location_list = [custom_location] if custom_location else None

    results = cutoff(
        Stations=stations_list,
        customlocations=custom_location_list,
        cutoff_comp="Vertical",
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params=rigidity_inputs,
        integration_params=integration_inputs,
        particle_params=particle_inputs,
        computation_params=computation_inputs,
        coordinate_params=coord_inputs,
        data_retrieval_params=data_retrieval_inputs
    )

    rigidity_df = results[0]

    Rl = rigidity_df.loc["Rl", station_name]
    Rc = rigidity_df.loc["Rc", station_name]
    Ru = rigidity_df.loc["Ru", station_name]
    PTF = rigidity_df.loc["PTF", station_name]

    return Rl, Rc, Ru, PTF


def append_result(date_str, Rl, Rc, Ru, PTF):
    write_header = not os.path.exists(CSV_FILE)

    with open(CSV_FILE, "a", newline="") as f:
        writer = csv.writer(f)
        if write_header:
            writer.writerow(["Date", "Rl", "Rc", "Ru", "PTF"])
        writer.writerow([date_str, Rl, Rc, Ru, PTF])


# ============================================================
# SCHEDULING
# ============================================================

def next_hour_mark(dt):
    return (dt.replace(minute=0, second=0, microsecond=0) + timedelta(hours=1))


def run_once(station_name, custom_location):
    now_utc = datetime.now(timezone.utc)
    target_time = now_utc.replace(minute=0, second=0, microsecond=0) - timedelta(hours=LIVE_DATA_LAG_HOURS)

    Rl, Rc, Ru, PTF = compute_cutoff(target_time, station_name, custom_location)

    date_str = target_time.strftime("%Y-%m-%d %H:%M:%S")
    append_result(date_str, Rl, Rc, Ru, PTF)

    print(f"[{date_str}] {station_name}: Rl={Rl:.4f} GV, Rc={Rc:.4f} GV, Ru={Ru:.4f} GV, PTF={PTF:.4f}")


def main():
    location = detect_location()

    if location is None:
        print(f"Could not detect current location; defaulting to the {FALLBACK_STATION} station.")
        station_name = FALLBACK_STATION
        custom_location = None
    else:
        station_name, lat, lon = location
        custom_location = [station_name, lat, lon]
        print(f"Detected location: {station_name} ({lat:.4f}, {lon:.4f})")

    print(f"Logging hourly cutoff rigidities to {CSV_FILE}")
    print(f"Each computation uses live data timestamped {LIVE_DATA_LAG_HOURS}h behind the current UTC hour.\n")

    while True:
        try:
            run_once(station_name, custom_location)
        except (Exception, SystemExit) as e:
            print(f"Cutoff computation failed for this hour, skipping: {e}")

        now_utc = datetime.now(timezone.utc)
        sleep_seconds = (next_hour_mark(now_utc) - now_utc).total_seconds()
        time.sleep(max(sleep_seconds, 1))


if __name__ == "__main__":
    main()
