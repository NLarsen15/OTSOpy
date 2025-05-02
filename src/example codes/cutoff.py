import OTSO

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    cutoff = OTSO.cutoff(Stations=stations_list,corenum=1,year=2000,month=1,day=1,hour=0)


    print(cutoff[0]) # dataframe output contating Ru, Rc, Rl for all input locations
    print(cutoff[1]) # text output of input variable information




