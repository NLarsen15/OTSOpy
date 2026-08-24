from OTSO import cutoff

date_inputs           = {"year": 2020, "month": 1, "day": 1, "hour":1,
                        "minute":0, "second":0}

magfield_inputs       = {"internalmag": "IGRF", "externalmag": "TSY89c",
                        "boberg": False, "bobergtype": "EXTENSION",
                        "magnetopause": "Sphere", "spheresize": 10,
                        "AdaptiveExternalModel": False}

rigidity_inputs       = {"startrigidity": 20, "endrigidity": 0,
                        "rigiditystep": 0.01, "rigidityscan": "ON"}

asymptotic_inputs     = {"asymptotic": "NO", "asymlevels": [0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000],
                        "unit": "GeV"}

transmission_inputs   = {"transmission": False, "transmissionRstep": 0.001,
                         "transmissionsamples": 20}

solar_wind_inputs     = {"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
                         "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0}

geomagnetic_inputs    = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0,
                         "sym_h_corrected": 0}

tsyganenko_inputs     = {"G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
                         "W4": 0, "W5": 0, "W6": 0}

integration_inputs    = {"intmodel": "Boris-Buneman", "gyropercent": 15, "minaltitude": 20,
                         "maxdistance": 100, "maxtime": 0, "mintrapdist": 6.6,
                         "startaltitude": 20, "betaerror": 0.001,
                         "totalbetacheck": True, "adaptivestep": True}

particle_inputs       = {"Anum": 1, "anti": "YES", "zenith": 0, "azimuth": 0}

computation_inputs    = {"corenum": 1, "threadnum": 1, "Verbose": True,
                         "delim": ";"}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDfile": None,
                         "MHDcoordsys": None}

coord_inputs          = {"coordsystem": "GEO", "inputcoord": "GDZ"}
# Note: coordsystem will also control the coordinate system of asymptotic viewing direction results

data_retrieval_inputs = {"serverdata": "OFF", "livedata": "OFF"}


if __name__ == '__main__':

    stations_list = ["OULU"]
    #cutoff_comp can be set as "Vertical, Apparent, and Custom"

    # Example using grouped parameters
    cutoff_results = cutoff(
        Stations=stations_list,
        cutoff_comp="Vertical",
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        rigidity_params=rigidity_inputs,
        asymptotic_params=asymptotic_inputs,
        transmission_params=transmission_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs,
        integration_params=integration_inputs,
        particle_params=particle_inputs,
        computation_params=computation_inputs,
        coordinate_params=coord_inputs,
        custom_field_params=custom_field_inputs,
        data_retrieval_params=data_retrieval_inputs
    )

    print(cutoff_results[0]) # dataframe output containing Ru, Rc, Rl for all input locations
    #print(cutoff_results[1]) # dataframe output containing asymptotic viewing direction results for all input locations
    #print(cutoff_results[2]) # dataframe output containing transmission functions for all input locations
    #print(cutoff_results[-1]) # text output of input variable information



