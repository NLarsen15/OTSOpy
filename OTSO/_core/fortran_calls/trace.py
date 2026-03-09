import pandas as pd
import os
import multiprocessing as mp
import time

from ..utils.maglines_utils import Lshell
from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.trace_data import TraceData

def FortranTrace(Data: list, TraceDataInstance: TraceData, queue: mp.Queue):
    
    if TraceDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(TraceDataInstance.MHDfile)
    
    for x in Data:
      
      Position = [x[3],x[1],x[2],x[4],x[5]]

      AtomicNum = TraceDataInstance.particlearray[0]
      AntiCheck = TraceDataInstance.particlearray[1]

      Rigidity = TraceDataInstance.rigidity
      CoordinateSystem = TraceDataInstance.coordsys
      MaxStepPercent = TraceDataInstance.maxsteppercent
      EndParams = TraceDataInstance.endparams
      Magnetopause = TraceDataInstance.magnetopause
      model = TraceDataInstance.model
      IntModel = TraceDataInstance.integrationmodel
      IOPT = TraceDataInstance.IOPT
      WindArray = TraceDataInstance.windarray
      g = TraceDataInstance.g
      h = TraceDataInstance.h
      MHDCoordSys = TraceDataInstance.MHDcoordsys
      spheresize = TraceDataInstance.spheresize
      inputcoord = TraceDataInstance.inputcoord
      DateArray = TraceDataInstance.datearray
      trapdist = TraceDataInstance.mindist

      Filename = f"{x[1]}_{x[2]}.csv"
      name = f"{x[1]}_{x[2]}"

      OTSOLib.fieldtrace(Position, Rigidity, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, MaxStepPercent,
                          EndParams, Filename, g, h, MHDCoordSys,spheresize, inputcoord)
      max_retries = 5
      retry_delay = 0.2  # seconds
      attempt = 0
      while attempt < max_retries:
        if not os.path.exists(Filename):
          #print(f"[fortrancallTrace] File not found: {Filename} (attempt {attempt+1}/{max_retries})")
          time.sleep(retry_delay)
          attempt += 1
          continue
        try:
          Trace = pd.read_csv(Filename)
          break
        except Exception as e:
          #print(f"[fortrancallTrace] Error reading {Filename} (attempt {attempt+1}/{max_retries}): {e}")
          time.sleep(retry_delay)
          attempt += 1
      else:
        #print(f"[fortrancallTrace] Failed to read {Filename} after {max_retries} attempts. Skipping this trace.")
        continue
      coordsystem2 = "GSM"
      columns = Trace.columns

      new_columns = [f"{col}_{CoordinateSystem} [Re]" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]
      if CoordinateSystem == "GDZ":
        column_names = ["alt [km]", "latitude", "longitude"]
        new_columns = [f"{column_names[i]}" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]
      elif CoordinateSystem == "SPH":
        column_names = ["radius [Re]", "latitude", "longitude"]
        new_columns = [f"{column_names[i]}" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]
  
      Trace.columns = new_columns

      Trace.columns = new_columns

      L, invlat = Lshell(Trace, CoordinateSystem)
      # Round to 4 decimal places and convert to float
      L = round(float(L), 4)
      invlat = round(float(invlat), 4)
      dataframe_dict = {
        name: {
          "Trace": Trace,
          "L_shell": L,
          "Invariant_Latitude": invlat
        }}  
      queue.put(dataframe_dict)
      os.remove(Filename)
    
    return