from OTSO import trace

date_inputs           = {"year": 2020, "month": 1, "day": 1, "hour":1,
                        "minute":0, "second":0}

magfield_inputs       = {"internalmag": "IGRF", "externalmag": "TSY89c",
                        "boberg": False, "bobergtype": "EXTENSION",
                        "magnetopause": "Sphere", "spheresize": 10}

solar_wind_inputs     = {"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
                         "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0}

geomagnetic_inputs    = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0,
                         "sym_h_corrected": 0}

tsyganenko_inputs     = {"G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
                         "W4": 0, "W5": 0, "W6": 0}

grid_inputs           = {"latstep": -5, "longstep": 30, "maxlat": 85, "minlat":-85,
                         "maxlong": 360, "minlong": 0}

integration_inputs    = {"minaltitude": 20, "maxdistance": 100, 
                         "maxtime": 0, "startaltitude": 20, 
                         "maxsteps": 100000}

computation_inputs    = {"corenum": 4, "Verbose": True}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDfile": None,
                         "MHDcoordsys": None}

data_retrieval_inputs = {"serverdata": "OFF", "livedata": "OFF"}

if __name__ == '__main__':

    coordinatesystem = "GEO"

    # Example using grouped parameters
    trace_results = trace(
        coordsys=coordinatesystem,
        datetime_params=date_inputs,
        computation_params=computation_inputs,
        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        grid_params=grid_inputs
    )

    print(trace_results[0]) # dictionary output containing positional information magnetic field lines generated over
                    # the globe
    print(trace_results[1]) # dataframe output containing L-shell and invariant latitude for each traced location
    print(trace_results[-1]) # text output of input variable information

