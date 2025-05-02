import OTSO
import datetime

if __name__ == '__main__':

    lat_lon_alt_list = [[10,10,10]] # [[Latitude,Longitude,Altitude]]
    date_list = [datetime.datetime(2000,10,12,8)] # [dates]
    
    Coords = OTSO.coordtrans(Locations=lat_lon_alt_list,dates=date_list,CoordIN="GEO",CoordOUT="GSM",corenum=1)

    print(Coords[0]) # dataframe output of converted coordinates
    print(Coords[1]) # text output detailing the initial and final conversion coordinate system
    
