from OTSO import trajectory

date_inputs           = {"year": 2020, "month": 1, "day": 1, "hour":1,
                        "minute":0, "second":0}

magfield_inputs       = {"internalmag": "IGRF", "externalmag": "TSY89c",
                        "boberg": False, "bobergtype": "EXTENSION",
                        "magnetopause": "Sphere", "spheresize": 10,
                        "AdaptiveExternalModel": False}

solar_wind_inputs     = {"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
                         "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0}

geomagnetic_inputs    = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0,
                         "sym_h_corrected": 0}

tsyganenko_inputs     = {"G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
                         "W4": 0, "W5": 0, "W6": 0}

integration_inputs    = {"intmodel": "Boris-Buneman", "gyropercent": 5, "minaltitude": 20,
                         "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                         "startaltitude": 20, "betaerror": 0.001,
                         "totalbetacheck": True, "adaptivestep": False,
                         "maxsteps": 100000}

particle_inputs       = {"rigidity": 1, "Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0}

computation_inputs    = {"corenum": 1, "Verbose": True}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDfile": None,
                         "MHDcoordsys": None}

coord_inputs          = {"coordsystem": "GEO", "inputcoord": "GDZ"}
# Note: coordsystem will also control the coordinate system of asymptotic viewing direction results

data_retrieval_inputs = {"serverdata": "OFF", "livedata": "OFF"}

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)
    custom_location_list = [["Custom", 25, 25]] # list of lists ["name", latitude, longitude]

    # Example using grouped parameters
    trajectory_results = trajectory(
        Stations=stations_list,
        customlocations=custom_location_list,
        particle_params=particle_inputs,
        computation_params=computation_inputs,
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params=integration_inputs,
        custom_field_params=custom_field_inputs,
        coordinate_params=coord_inputs,
        data_retrieval_params=data_retrieval_inputs,
    )

    print(trajectory_results[0]) # dictionary output containing positional information for all trajectories generated starting
                         # from input stations
    print(trajectory_results[-1]) # text output of input variable information


