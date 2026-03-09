import time
from datetime import datetime
import multiprocessing as mp
import sys
import queue
from tqdm import tqdm

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import trajectory_inputs
from ..fortran_calls import trajectory
from ..readme_generators import trajectory_readme
from ..data_classes.trajectory_data import TrajectoryData

def OTSO_trajectory(TrajectoryDataInstance: 'TrajectoryData') -> list:

    trajectory_inputs.TrajectoryInputs(TrajectoryDataInstance)

    ChildProcesses = []

    UsedCores = cores.Cores(TrajectoryDataInstance.station_array, TrajectoryDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()
    if TrajectoryDataInstance.Verbose:
        print("OTSO Trajectory Computation Started")

    total_stations = len(TrajectoryDataInstance.station_array)
    results = []

    if TrajectoryDataInstance.corenum == 1:
        num_batches = len(Positionlists)
        
        progress_bar = None
        if TrajectoryDataInstance.Verbose:
            progress_bar = tqdm(total=num_batches, desc="OTSO Running", unit=" batch")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        for Data in Positionlists:
            trajectory.FortranTrajectory(Data, TrajectoryDataInstance, simple_queue)
            
            processed += 1
            if TrajectoryDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(1)
                else:
                    percent_complete = (processed / num_batches) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()
        
        results = simple_queue.get_all()
        
        # Close progress bar if it was created
        if progress_bar is not None:
            progress_bar.close()

    else:
        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass
        ProcessQueue = mp.Manager().Queue()
        for Data in Positionlists:
            Child = mp.Process(target=trajectory.FortranTrajectory,  args=(Data, TrajectoryDataInstance, ProcessQueue))
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

    result_list = results

    stop = time.time()
    Printtime = round((stop-start),3)
    
    if TrajectoryDataInstance.Verbose:
        print("\nOTSO Trajectory Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")

    EventDate = datetime(TrajectoryDataInstance.year,TrajectoryDataInstance.month,TrajectoryDataInstance.day,TrajectoryDataInstance.hour,
                         TrajectoryDataInstance.minute,TrajectoryDataInstance.second)

    readme = trajectory_readme.READMETrajectory(TrajectoryDataInstance, EventDate, Printtime)

    if TrajectoryDataInstance.livedata == "ON" or TrajectoryDataInstance.livedata == 1:
        file_clean.remove_files()

    return [result_list, readme]