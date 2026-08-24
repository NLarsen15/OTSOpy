import pandas as pd
import multiprocessing as mp
import numpy as np
import os
import tempfile
import csv
import json

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.trajectory_data import TrajectoryData


def FortranTrajectory(Data: list, TrajectoryDataInstance: TrajectoryData, queue: mp.Queue, JsonFile: str, lock) -> None:
    
    if TrajectoryDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(TrajectoryDataInstance.MHDfile)

    for x in Data:

      File = tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name

      Station = x[0]

      FortranData = fortran_data_utils.prepare_fortran_trajectory(x,  TrajectoryDataInstance)

      filter = np.array([0], dtype=np.int32)
      alat = np.array([0.0], dtype=np.float64)
      along = np.array([0.0], dtype=np.float64)

      filter, alat, along = OTSOLib.trajectory_full(FortranData,
                TrajectoryDataInstance.g,
                TrajectoryDataInstance.h,
                TrajectoryDataInstance.rigidity,
                File,
                int(len(File))
            )

      coordsystem = TrajectoryDataInstance.coordsystem
      headers = [
          f"{coordsystem}_X [Re]",
          f"{coordsystem}_Y [Re]",
          f"{coordsystem}_Z [Re]",
          f"{coordsystem}_Vx [km/s]",
          f"{coordsystem}_Vy [km/s]",
          f"{coordsystem}_Vz [km/s]",
      ]

      with open(File, "r", newline="", encoding="utf-8") as f:
          rows = list(csv.reader(f))

      with open(File, "w", newline="", encoding="utf-8") as f:
          writer = csv.writer(f)
          writer.writerow(headers)
          writer.writerows(rows)

      trajectory_df = pd.read_csv(File)

      #print(f"Trajectory computation for station {Station} completed.")
      #print(trajectory_df)

      record = {
        "station": Station,
        "rigidity": TrajectoryDataInstance.rigidity,
        "Filter": filter,
        "Alat": alat,
        "Along": along,
        "trajectory": trajectory_df.to_dict(orient="records")
      }

      with lock:        
          with open(JsonFile, "a", encoding="utf-8") as f:
              json.dump(record, f)
              f.write("\n")


      os.remove(File)

      queue.put(1)
    exit(1)
    return


