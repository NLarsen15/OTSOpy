from OTSO import cutoff

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    # Example using grouped parameters
    cutoff_results = cutoff(
        Stations=stations_list,
        computation_params={"corenum": 1},
        datetime_params={"year": 2000, "month": 1, "day": 1, "hour": 0}
    )

    print(cutoff_results[0]) # dataframe output containing Ru, Rc, Rl for all input locations
    print(cutoff_results[1]) # text output of input variable information



