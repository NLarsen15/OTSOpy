import pandas as pd
import multiprocessing as mp
import numpy as np
import json

from ..custom_classes import date
from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..utils import cutoff_utils as cu
from ..utils import transmission_utils as tu
from ..utils import cpus_utils as cpu_util
from ..data_classes.flight_data import FlightData

def FortranFlight(Data: list, DateArray: list, IOPT: list, WindArray: list, GArray: list, HArray: list,
                  JsonFile: str, FData: FlightData, queue: mp.Queue, cpus) -> None:

    cpu_util.set_process_affinity(cpus)

    if FData.model[1] == 99:
      mhd_utils.MHDinitialise(FData.MHDfile)


    for x,y,z,I,G,H in zip(Data,DateArray,WindArray,IOPT,GArray,HArray):
        
        Position = [x[3],x[1],x[2],x[4],x[5]]
        Station = x[0]

        FortranData = fortran_data_utils.prepare_fortran_flight(x,  
                                                                FData,
                                                                y,z,I,G,H)

        if FData.cutoff_comp == "Apparent":
            returnvalues = cu.Apparent_cutoff(Position, FortranData, G, H, FData)
        else:
            returnvalues = cu.standard_cutoff(FortranData, G, H, FData)
        
        datetimeobj = date.convert_to_datetime(y)

        R_high, R_eff, R_low, transparency = returnvalues

        if FData.transmission:
            Transmissiondf = tu.transmission(FortranData, R_high, R_low, Station, FData)
            # Save transmission dataframe to CSV
            #print(Transmissiondf)

        if FData.asymptotic == "YES":
          Energy_List = FData.asymlevels.copy()
          P_List = []
          if FData.unit == "GeV":   
              E_0 = 0.938  # Rest mass energy of proton in GeV
              for i in Energy_List:
                  R = (i**2 + 2*i*E_0)**(0.5)
                  P_List.append(R)
          else:
              if FData.unit == "GV": 
                  P_List = Energy_List.copy()

          AsyAllowed = np.zeros(len(P_List), dtype=np.int32)
          AsyLat = np.zeros(len(P_List), dtype=np.float64)
          AsyLong = np.zeros(len(P_List), dtype=np.float64)
          Plist_array = np.array(P_List, dtype=np.float64)
          OTSOLib.trajectory(FortranData, G, H, Plist_array, 
                             len(P_List), AsyAllowed, AsyLat, AsyLong)

        if FData.inputcoord == "GDZ":
            altunit = "km"
        else:
            altunit = "Re"
            
        record = {
            "Date": datetimeobj.strftime("%Y-%m-%d %H:%M:%S"),
            "Latitude": x[1],
            "Longitude": x[2],
            f"Altitude [{altunit}]": x[3],
            "Ru [GV]": R_high,
            "Rc [GV]": R_eff,
            "Rl [GV]": R_low,
            "PTF": transparency
          }

        # Optional asymptotic
        if FData.asymptotic == "YES":
            asymlevels_with_units = [
                      f"{level} [{FData.unit}]"
                      for level in FData.asymlevels
                  ]
        
            formatted_asymptotic_data = [
                      f"{AsyAllowed[i]}{FData.delim}"
                      f"{round(AsyLat[i], 4)}{FData.delim}"
                      f"{round(AsyLong[i], 4)}"
                      for i in range(len(AsyAllowed))
                  ]
        
            asymptotic_df = pd.DataFrame(
                      [[Station] + formatted_asymptotic_data],
                      columns=["Station"] + asymlevels_with_units,
                  )
        
            record["asymptotic"] = asymptotic_df.to_dict(orient="records")
        
              # Optional transmission
        if FData.transmission:
            transmission = {}
            tf_column = Transmissiondf.columns[1]
            for _, row in Transmissiondf.iterrows():
                transmission[f'{row["R [GV]"]} [GV]'] = row[tf_column]
        
            record["transmission"] = transmission

        with open(JsonFile, "a", encoding="utf-8") as f:
              json.dump(record, f)
              f.write("\n")

        queue.put(1)
    return