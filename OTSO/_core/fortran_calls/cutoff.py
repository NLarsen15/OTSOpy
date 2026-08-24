import json
import pandas as pd
import multiprocessing as mp
import numpy as np

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.cutoff_data import CutoffData
from ..utils import cutoff_utils as cu
from ..utils import transmission_utils as tu
from ..utils import cpus_utils as cpu_util


def FortranCutoff(Data: list, CutoffDataInstance: CutoffData, queue: mp.Queue, cpus, JsonFile, lock) -> None:

    cpu_util.set_process_affinity(cpus)

    if CutoffDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(CutoffDataInstance.MHDfile)

    for x in Data:

      Station = x[0]
      Position = [x[3],x[1],x[2],x[4],x[5]]

      FortranData = fortran_data_utils.prepare_fortran_cutoff(x,  CutoffDataInstance)

      if CutoffDataInstance.cutoff_comp == "Apparent":
            returnvalues = cu.Apparent_cutoff(Position, FortranData, CutoffDataInstance.g,
                                               CutoffDataInstance.h, CutoffDataInstance)
      else:
            returnvalues = cu.standard_cutoff(FortranData, CutoffDataInstance.g, CutoffDataInstance.h,
                                               CutoffDataInstance)

      R_high, R_eff, R_low, transparency = returnvalues

      if CutoffDataInstance.transmission:
            Transmissiondf = tu.transmission(FortranData, R_high, R_low, Station, CutoffDataInstance)
            # Save transmission dataframe to CSV
            #print(Transmissiondf)

      if CutoffDataInstance.asymptotic == "YES":
          Energy_List = CutoffDataInstance.asymlevels.copy()
          P_List = []
          if CutoffDataInstance.unit == "GeV":   
              E_0 = 0.938  # Rest mass energy of proton in GeV
              for i in Energy_List:
                  R = (i**2 + 2*i*E_0)**(0.5)
                  P_List.append(R)
          else:
              if CutoffDataInstance.unit == "GV": 
                  P_List = Energy_List.copy()

          AsyAllowed = np.zeros(len(P_List), dtype=np.int32)
          AsyLat = np.zeros(len(P_List), dtype=np.float64)
          AsyLong = np.zeros(len(P_List), dtype=np.float64)
          Plist_array = np.array(P_List, dtype=np.float64)
          OTSOLib.trajectory(FortranData, CutoffDataInstance.g, CutoffDataInstance.h, 
                             Plist_array, len(P_List), AsyAllowed, AsyLat, AsyLong)


      
      Rigiditydataframe = pd.DataFrame({Station: returnvalues}, index=['Ru', 'Rc', 'Rl', 'PTF'])
      
      record = {
          "station": Station,
          "rigidities": Rigiditydataframe.to_dict(orient="split"),
      }

      # Optional transmission
      if CutoffDataInstance.transmission:
          record["transmission"] = Transmissiondf.to_dict(orient="split")
          # or orient="records" if that's more convenient

      # Optional asymptotic
      if CutoffDataInstance.asymptotic == "YES":
          asymlevels_with_units = [
              f"{level} [{CutoffDataInstance.unit}]"
              for level in CutoffDataInstance.asymlevels
          ]

          formatted_asymptotic_data = [
              f"{AsyAllowed[i]}{CutoffDataInstance.delim}"
              f"{round(AsyLat[i], 4)}{CutoffDataInstance.delim}"
              f"{round(AsyLong[i], 4)}"
              for i in range(len(AsyAllowed))
          ]

          asymptotic_df = pd.DataFrame(
              [[Station] + formatted_asymptotic_data],
              columns=["Station"] + asymlevels_with_units,
          )

          record["asymptotic"] = asymptotic_df.to_dict(orient="records")

      with lock:
          with open(JsonFile, "a", encoding="utf-8") as f:
              json.dump(record, f)
              f.write("\n")

      queue.put(1)

    return


