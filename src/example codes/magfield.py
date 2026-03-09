from OTSO import magfield

if __name__ == '__main__':

    location_list = [[10,10,10]] # [[X,Y,Z]] Earth radii Geocentric coordinates in this instance

    # Example using grouped parameters
    magfield_results = magfield(
        Locations=location_list,
        coordinate_params={"coordsystem": "GEO"},
        computation_params={"corenum": 1}
    )

    print(magfield_results[0]) # dataframe of returned magnetic field vectors at inputted locations
    print(magfield_results[1]) # text output of input variable information


