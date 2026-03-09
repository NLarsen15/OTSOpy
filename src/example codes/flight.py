from OTSO import flight
import datetime

if __name__ == '__main__':

    latitude_list = [10,15,20,25,30] # [Latitudes]
    longitude_list = [10,15,20,25,30] # [Longitudes]
    altitude_list = [30,40,50,60,80] # [Altitudes]
    date_list = [datetime.datetime(2000,10,12,8),datetime.datetime(2000,10,12,9),datetime.datetime(2000,10,12,10),
                 datetime.datetime(2000,10,12,11),datetime.datetime(2000,10,12,12)] # [dates]

    # Example using grouped parameters
    flight_results = flight(
        latitudes=latitude_list, 
        longitudes=longitude_list,
        dates=date_list,
        altitudes=altitude_list,
        cutoff_comp="Vertical",
        computation_params={"corenum": 1}
    )
    
    print(flight_results[0]) # dataframe output containing Ru, Rc, Rl along flightpath
    print(flight_results[1]) # text output of input variable information
    print(flight_results[2])  # dataframe output of input variables
