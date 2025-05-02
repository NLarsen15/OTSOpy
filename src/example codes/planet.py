import OTSO

if __name__ == '__main__':
    

    planet = OTSO.planet(corenum=1, cutoff_comp="Vertical", year=2000, rigiditystep=0.1)

    print(planet[0]) # dataframe containing cutoff results for planet grid
    print(planet[1]) # text output of input variable information


