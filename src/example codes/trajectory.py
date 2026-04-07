from OTSO import trajectory

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    # Example using grouped parameters
    trajectory_results = trajectory(
        Stations=stations_list,
        particle_params={"rigidity": 5},
        computation_params={"corenum": 1}
    )

    print(trajectory_results[0]) # dictionary output containing positional information for all trajectories generated starting
                         # from input stations
    print(trajectory_results[1]) # text output of input variable information


