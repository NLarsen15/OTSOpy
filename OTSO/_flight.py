from __future__ import annotations
from typing import TYPE_CHECKING, Sequence, Optional, Union

if TYPE_CHECKING:
    from ._core.otso_functions.otso_flight import OTSO_flight

def flight(**kwargs):
    import psutil    
    from dataclasses import fields
    from ._core.otso_functions import otso_flight
    from ._core.data_classes.flight_data import FlightData

    allowed = {f.name for f in fields(FlightData)}

    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(
            f"Unexpected keyword arguments: {', '.join(sorted(unknown))}"
    )

    # Set corenum if not provided in kwargs
    if 'corenum' not in kwargs or kwargs['corenum'] is None:
        kwargs['corenum'] = psutil.cpu_count(logical=True) - 2
        if kwargs['corenum'] <= 0:
            kwargs['corenum'] = 1
    
    # Handle None values for list parameters
    for key in ['g', 'h', 'MHDfile', 'MHDcoordsys']:
        if key in kwargs and kwargs[key] is None:
            kwargs[key] = []

    # Handle special flight-specific parameter defaults for arrays
    if kwargs.get("serverdata", "OFF") != "ON" and kwargs.get("livedata", "OFF") != "ON":
        default_values = {
            "vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 5.0, "bz": 5.0, "density": 1, "pdyn": 0, "Dst": 0,
            "G1": 0.00, "G2": 0.00, "G3": 0.00, "W1": 0, "W2": 0, "W3": 0, "W4": 0, "W5": 0, "W6": 0, "kp": 0,
            "by_avg": 0, "bz_avg": 0, "n_index": 0, "b_index": 0, "sym_h_corrected": 0
        }
        
        for var_name, default_value in default_values.items():
            if kwargs.get(var_name) is None or not kwargs.get(var_name):  # If None or empty list
                kwargs[var_name] = [default_value] * len(kwargs["latitudes"])
    else:
        # For serverdata="ON" or livedata="ON", ensure parameters exist but can be None
        default_values = {
            "vx": None, "vy": None, "vz": None, "bx": None, "by": None, "bz": None, "density": None, "pdyn": None, "Dst": None,
            "G1": None, "G2": None, "G3": None, "W1": None, "W2": None, "W3": None, "W4": None, "W5": None, "W6": None, "kp": None,
            "by_avg": None, "bz_avg": None, "n_index": None, "b_index": None, "sym_h_corrected": None
        }
        
        for var_name, default_value in default_values.items():
            if var_name not in kwargs:
                kwargs[var_name] = default_value

    FlightDataInstance = FlightData(
        latitudes = kwargs["latitudes"], 
        longitudes = kwargs["longitudes"], 
        dates = kwargs["dates"], 
        altitudes = kwargs["altitudes"],
        cutoff_comp = kwargs["cutoff_comp"], 
        minaltitude = kwargs["minaltitude"], 
        maxdistance = kwargs["maxdistance"], 
        maxtime = kwargs["maxtime"], 
        serverdata = kwargs["serverdata"], 
        livedata = kwargs["livedata"], 
        vx = kwargs["vx"], 
        vy = kwargs["vy"], 
        vz = kwargs["vz"], 
        bx = kwargs["bx"], 
        by = kwargs["by"], 
        bz = kwargs["bz"], 
        density = kwargs["density"], 
        pdyn = kwargs["pdyn"], 
        Dst = kwargs["Dst"], 
        G1 = kwargs["G1"], 
        G2 = kwargs["G2"], 
        G3 = kwargs["G3"], 
        W1 = kwargs["W1"], 
        W2 = kwargs["W2"], 
        W3 = kwargs["W3"], 
        W4 = kwargs["W4"], 
        W5 = kwargs["W5"], 
        W6 = kwargs["W6"], 
        kp = kwargs["kp"], 
        by_avg = kwargs["by_avg"], 
        bz_avg = kwargs["bz_avg"], 
        n_index = kwargs["n_index"], 
        b_index = kwargs["b_index"], 
        sym_h_corrected = kwargs["sym_h_corrected"], 
        Anum = kwargs["Anum"], 
        anti = kwargs["anti"], 
        internalmag = kwargs["internalmag"], 
        externalmag = kwargs["externalmag"],
        boberg = kwargs["boberg"], 
        bobergtype = kwargs["bobergtype"], 
        intmodel = kwargs["intmodel"], 
        startrigidity = kwargs["startrigidity"], 
        endrigidity = kwargs["endrigidity"], 
        rigiditystep = kwargs["rigiditystep"], 
        rigidityscan = kwargs["rigidityscan"], 
        coordsystem = kwargs["coordsystem"],
        gyropercent = kwargs["gyropercent"],
        fixedstep = kwargs["fixedstep"],
        magnetopause = kwargs["magnetopause"],
        corenum = kwargs["corenum"], 
        azimuth = kwargs["azimuth"], 
        zenith = kwargs["zenith"], 
        g = kwargs["g"], 
        h = kwargs["h"], 
        asymptotic = kwargs["asymptotic"], 
        asymlevels = kwargs["asymlevels"], 
        unit = kwargs["unit"],
        MHDfile = kwargs["MHDfile"], 
        MHDcoordsys = kwargs["MHDcoordsys"], 
        spheresize = kwargs["spheresize"], 
        inputcoord = kwargs["inputcoord"], 
        Verbose = kwargs["Verbose"], 
        AdaptiveExternalModel = kwargs["AdaptiveExternalModel"], 
        mintrapdist = kwargs["mintrapdist"],
        delim = kwargs["delim"],
        adaptivestep = kwargs['adaptivestep'], 
        betaerror = kwargs['betaerror'], 
        totalbetacheck = kwargs['totalbetacheck'],
        maxsteps = kwargs['maxsteps'],
        threadnum = kwargs['threadnum'],
        transmission = kwargs['transmission'],
        transmissionsamples = kwargs['transmissionsamples'],
        transmissionRstep = kwargs['transmissionRstep'],
        max_degree = kwargs["max_degree"]
    )
    
    flight = otso_flight.OTSO_flight(FlightDataInstance)

    
    return flight