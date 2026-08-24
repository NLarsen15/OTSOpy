import pandas as pd
import os
import multiprocessing as mp
import time
import tempfile
import csv
import json

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils.maglines_utils import Lshell
from ..utils import fortran_data_utils
from ..utils import mhd_utils
from ..data_classes.trace_data import TraceData

def FortranTrace(Data: list, TraceDataInstance: TraceData, queue: mp.Queue, JsonFile):
    
    if TraceDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(TraceDataInstance.MHDfile)
    
    for x in Data:

      File = tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name
      with open(File, mode='a', newline='', encoding='utf-8') as file:  # Open in write ('w') mode
                  writer = csv.writer(file)
                  coordsystem = TraceDataInstance.coordsys
                  headers = [f"X_{coordsystem} [Re]", f"Y_{coordsystem} [Re]", f"Z_{coordsystem} [Re]", 
                             f"Bx_GSM [nT]", f"By_GSM [nT]", f"Bz_GSM [nT]"]
                  writer.writerow(headers)
    
      CoordinateSystem = TraceDataInstance.coordsys

      FortranData = fortran_data_utils.prepare_fortran_trace(x,  TraceDataInstance)


      FileNameLen = len(File)
      name = f"{x[1]}_{x[2]}"

      OTSOLib.fieldtrace(FortranData, File, FileNameLen, TraceDataInstance.g, TraceDataInstance.h)

      max_retries = 5
      retry_delay = 0.2
      attempt = 0

      while attempt < max_retries:
        if not os.path.exists(File):
          time.sleep(retry_delay)
          attempt += 1
          continue

        try:
          Trace = pd.read_csv(File)
          break

        except Exception:
          time.sleep(retry_delay)
          attempt += 1

      else:
        if os.path.exists(File):
          os.remove(File)
        continue


      # Rename columns based on coordinate system
      coordsystem2 = "GSM"
      columns = Trace.columns
      new_columns = None

      if CoordinateSystem == "GDZ":
          column_names = ["alt [km]", "latitude", "longitude"]
          new_columns = [f"{column_names[i]}" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]

      elif CoordinateSystem == "SPH":
          column_names = ["radius [Re]", "latitude", "longitude"]
          new_columns = [f"{column_names[i]}" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]

      if new_columns:
        Trace.columns = new_columns


      # Calculate L-shell and invariant latitude
      L, invlat = Lshell(Trace, CoordinateSystem)

      L = round(float(L), 4)
      invlat = round(float(invlat), 4)

      if TraceDataInstance.inputcoord == "GDZ":
        altunit = "km"
      else:
        altunit = "Re"


      # Convert dataframe into JSON-compatible format
      record = {
          "location": name,
          f"altitude [{altunit}]": TraceDataInstance.startaltitude,
          "Trace": {
              "columns": Trace.columns.tolist(),
              "data": Trace.values.tolist()
          },
          "L_shell": L,
          "Invariant_Latitude": invlat
      }

      with open(JsonFile, "a", encoding="utf-8") as f:
          json.dump(record, f)
          f.write("\n")

      os.remove(File)
      queue.put(1)
    
    return