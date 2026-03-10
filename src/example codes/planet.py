from OTSO import planet

if __name__ == '__main__':
    
    # Example using grouped parameters
    planet_results = planet(
        cutoff_comp="Vertical",
        grid_params={"latstep": -10, "lonstep": 15},
        computation_params={"corenum": 6},
        datetime_params={"year": 2000},
        rigidity_params={"rigiditystep": 0.1},
        integration_params={"totalbetacheck": False}
    )

    print(planet_results[0]) # dataframe containing cutoff results for planet grid
    print(planet_results[1]) # text output of input variable information

