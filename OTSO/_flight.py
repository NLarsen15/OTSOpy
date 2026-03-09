def flight(**kwargs):
    import psutil    
    from ._core.otso_functions import otso_flight
    from ._core.data_classes.flight_data import FlightData 

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
        kwargs["latitudes"], kwargs["longitudes"], kwargs["dates"], kwargs["altitudes"],
        kwargs["cutoff_comp"], kwargs["minaltitude"], kwargs["maxdistance"], kwargs["maxtime"], 
        kwargs["serverdata"], kwargs["livedata"], 
        kwargs["vx"], kwargs["vy"], kwargs["vz"], 
        kwargs["bx"], kwargs["by"], kwargs["bz"], 
        kwargs["density"], kwargs["pdyn"], kwargs["Dst"], 
        kwargs["G1"], kwargs["G2"], kwargs["G3"], 
        kwargs["W1"], kwargs["W2"], kwargs["W3"], kwargs["W4"], 
        kwargs["W5"], kwargs["W6"], kwargs["kp"], 
        kwargs["by_avg"], kwargs["bz_avg"], kwargs["n_index"], 
        kwargs["b_index"], kwargs["sym_h_corrected"], 
        kwargs["Anum"], kwargs["anti"], 
        kwargs["internalmag"], kwargs["externalmag"],
        kwargs["boberg"], kwargs["bobergtype"], 
        kwargs["intmodel"], kwargs["startrigidity"], kwargs["endrigidity"], kwargs["rigiditystep"], 
        kwargs["rigidityscan"], kwargs["coordsystem"], kwargs["gyropercent"], kwargs["magnetopause"], 
        kwargs["corenum"], kwargs["azimuth"], kwargs["zenith"], 
        kwargs["g"], kwargs["h"], kwargs["asymptotic"], kwargs["asymlevels"], kwargs["unit"],
        kwargs["MHDfile"], kwargs["MHDcoordsys"], kwargs["spheresize"], 
        kwargs["inputcoord"], kwargs["Verbose"], kwargs["AdaptiveExternalModel"], kwargs["mintrapdist"],
        kwargs["delim"],kwargs['adaptivestep'], kwargs['betaerror'], kwargs['totalbetacheck'],
        kwargs['maxsteps']
    )
    
    flight = otso_flight.OTSO_flight(FlightDataInstance)

    
    return flight