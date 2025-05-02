import OTSO

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    cone = OTSO.cone(Stations=stations_list,corenum=1,year=2000,month=1,day=1,hour=0)

    print(cone[0]) # dataframe output containing asymptotic cones for all input locations
    print(cone[1]) # dataframe output containing Ru, Rc, Rl for all inputted locations
    print(cone[2]) # text output of input variable information




