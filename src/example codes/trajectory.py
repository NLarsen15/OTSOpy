import OTSO

if __name__ == '__main__':

    stations_list = ["OULU","ROME","ATHN","CALG"] # list of neutron monitor stations (using their abbreviations)

    trajectory = OTSO.trajectory(Stations=stations_list,rigidity=5,corenum=1)

    print(trajectory[0]) # dictionary output containing positional information for all trajectories generated starting
                         # from input stations
    print(trajectory[1]) # text output of input variable information



