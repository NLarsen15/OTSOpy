from OTSO import cone

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    cone_results = cone(Stations=stations_list,computation_params={"corenum":1},
                        datetime_params={"year":2000,"month":1,"day":1,"hour":0})

    print(cone_results[0]) # dataframe output containing asymptotic cones for all input locations
    print(cone_results[1]) # dataframe output containing Ru, Rc, Rl for all inputted locations
    print(cone_results[2]) # text output of input variable information



