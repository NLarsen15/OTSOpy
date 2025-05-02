import OTSO

if __name__ == '__main__':

    location_list = [[10,10,10]] # [[X,Y,Z]] Earth radii Geocentric coordinates in this instance

    magfield = OTSO.magfield(Locations=location_list,coordsystem="GEO",corenum=1)

    print(magfield[0]) # dataframe of returned magnetic field vectors at inputted locations
    print(magfield[1]) # text output of input variable information



