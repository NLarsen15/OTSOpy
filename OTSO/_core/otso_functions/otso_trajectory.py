import os
import time
from datetime import datetime
import multiprocessing as mp
import sys
import queue
from tqdm import tqdm
import tempfile
import csv
import json
import pandas as pd

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import trajectory_inputs
from ..fortran_calls import trajectory
from ..readme_generators import trajectory_readme
from ..data_classes.trajectory_data import TrajectoryData

def OTSO_trajectory(TrajectoryDataInstance: 'TrajectoryData') -> list:

    trajectory_inputs.TrajectoryInputs(TrajectoryDataInstance)

    lock = mp.Lock()

    ChildProcesses = []

    UsedCores = cores.Cores(TrajectoryDataInstance.station_array, TrajectoryDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    File = tempfile.NamedTemporaryFile(delete=False, suffix=".json").name

    start = time.time()
    if TrajectoryDataInstance.Verbose:
        print("OTSO Trajectory Computation Started")

    total_stations = len(TrajectoryDataInstance.station_array)
    results = []

    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass
    ProcessQueue = mp.Manager().Queue()
    for Data in Positionlists:
        Child = mp.Process(target=trajectory.FortranTrajectory,  args=(Data, TrajectoryDataInstance, ProcessQueue, File, lock))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    processed = 0

    progress_bar = None
    if TrajectoryDataInstance.Verbose:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit="trajectory")

    while processed < total_stations:
      try:
        result_df = ProcessQueue.get(timeout=0.001)
        results.append(result_df)
        processed += 1
  
        if TrajectoryDataInstance.Verbose:
            if progress_bar is not None:
                progress_bar.update(1)
            else:
                percent_complete = (processed / total_stations) * 100
                sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                sys.stdout.flush()

      except queue.Empty:
        pass
      
      time.sleep(0.0001)

    if progress_bar is not None:
        progress_bar.close()

    for b in ChildProcesses:
        b.join()

    result_list = []


    with open(File, "r", encoding="utf-8") as f:
        for line in f:
            result = json.loads(line)
            result["trajectory"] = pd.DataFrame(result["trajectory"])
            result_list.append(result)

    stop = time.time()
    Printtime = round((stop-start),3)
    
    if TrajectoryDataInstance.Verbose:
        print("\nOTSO Trajectory Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")

    EventDate = datetime(TrajectoryDataInstance.year,TrajectoryDataInstance.month,TrajectoryDataInstance.day,TrajectoryDataInstance.hour,
                         TrajectoryDataInstance.minute,TrajectoryDataInstance.second)

    readme = trajectory_readme.READMETrajectory(TrajectoryDataInstance, EventDate, Printtime)

    if os.path.exists(File):
        os.remove(File)

    if TrajectoryDataInstance.livedata == "ON" or TrajectoryDataInstance.livedata == 1:
        file_clean.remove_files()

    return [result_list, readme]