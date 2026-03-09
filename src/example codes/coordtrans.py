from OTSO import coordtrans
import datetime

if __name__ == '__main__':

    lat_lon_alt_list = [[10,10,10]] # [[Latitude,Longitude,Altitude]]
    date_list = [datetime.datetime(2000,10,12,8)] # [dates]
    
    # Example using grouped parameters (coordtrans has fewer parameter groups)
    Coords = coordtrans(
        Locations=lat_lon_alt_list,
        dates=date_list,
        CoordIN="GEO",
        CoordOUT="GSM",
        corenum=1  # coordtrans uses individual parameters, not grouped ones
    )

    print(Coords[0]) # dataframe output of converted coordinates
    print(Coords[1]) # text output detailing the initial and final conversion coordinate system
    
