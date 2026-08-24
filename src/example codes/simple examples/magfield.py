from OTSO import magfield

date_inputs           = {"year": 2020, "month": 1, "day": 1, "hour":1,
                        "minute":0, "second":0}

magfield_inputs       = {"internalmag": "IGRF", "externalmag": "TSY89c",
                        "boberg": False, "bobergtype": "EXTENSION"}

solar_wind_inputs     = {"vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 1, "bz": 1,
                         "by_avg": 0, "bz_avg": 0, "density": 1, "pdyn": 0}

geomagnetic_inputs    = {"Dst": 0, "kp": 2, "n_index": 0, "b_index": 0,
                         "sym_h_corrected": 0}

tsyganenko_inputs     = {"G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0,
                         "W4": 0, "W5": 0, "W6": 0}

computation_inputs    = {"corenum": 1, "Verbose": True}

custom_field_inputs   = {"g": None, "h": None, "max_degree": 13, "MHDfile": None,
                         "MHDcoordsys": None}

coord_inputs          = {"inputcoord": "GDZ", "coordout": "GSM"}
# Note: coordsystem will also control the coordinate system of asymptotic viewing direction results

data_retrieval_inputs = {"serverdata": "OFF", "livedata": "OFF"}

if __name__ == '__main__':

    location_list = [[10,10,10]] # [[X,Y,Z]] Earth radii Geocentric coordinates in this instance

    # Example using grouped parameters
    magfield_results = magfield(
        Locations=location_list,
        coordinate_params=coord_inputs,
        computation_params=computation_inputs,
        datetime_params=date_inputs,
        magfield_params=magfield_inputs,
        solar_wind_params=solar_wind_inputs,
        geomagnetic_params=geomagnetic_inputs,
        tsyganenko_params=tsyganenko_inputs
    )

    print(magfield_results[0]) # dataframe of returned magnetic field vectors at inputted locations
    print(magfield_results[-1]) # text output of input variable information


