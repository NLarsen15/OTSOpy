import time
from datetime import datetime
import multiprocessing as mp
import pandas as pd
import sys
import queue
from tqdm import tqdm

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import cone_inputs
from ..fortran_calls import cone
from ..readme_generators import cone_readme
from ..data_classes.cone_data import ConeData


def OTSO_cone(CoreDataInstance: ConeData) -> list:

    cone_inputs.ConeInputs(CoreDataInstance)

    ChildProcesses = []
    results = []

    UsedCores = cores.Cores(CoreDataInstance.station_array, CoreDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()

    if CoreDataInstance.Verbose:
        print("OTSO Cone Computation Started")

    total_stations = len(CoreDataInstance.station_array)
    results = []

    if CoreDataInstance.corenum == 1:
        num_batches = len(Positionlists)
        
        progress_bar = None
        if CoreDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=num_batches, desc="OTSO Running", unit=" batch")
        elif CoreDataInstance.Verbose:
            print(f"Processing {total_stations} stations...")

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
            cone.FortranCone(Data, CoreDataInstance, simple_queue)
            
            processed += 1
            if CoreDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(1)
                else:
                    percent_complete = (processed / num_batches) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()
        
        results = simple_queue.get_all()
        
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
            Child = mp.Process(target=cone.FortranCone,  args=(Data, CoreDataInstance, ProcessQueue))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if CoreDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cone")
        elif CoreDataInstance.Verbose:
            print(f"Processing {total_stations} stations...")

        while processed < total_stations:
          try:
            result_df = ProcessQueue.get(timeout=0.001)
            results.append(result_df)
            processed += 1

            if CoreDataInstance.Verbose:
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

    conedf_list = []
    Rigiditylist = []

    for x in results:
        conedf_list.append(x[0])
        Rigiditylist.append(x[1])

    merged_cone_df = conedf_list[0]
    for df in conedf_list[1:]:
        merged_cone_df = pd.merge(merged_cone_df, df, on='R [GV]')
    cols = ['R [GV]'] + [col for col in merged_cone_df.columns if col != 'R [GV]']
    merged_cone_df = merged_cone_df[cols]

    merged_cone_df = merged_cone_df[['R [GV]'] + sorted(merged_cone_df.columns.drop('R [GV]'))]

    merged_R_df = pd.concat(Rigiditylist, axis=1)
    merged_R_df = merged_R_df.sort_index(axis=1)

    stop = time.time()
    Printtime = round((stop-start),3)

    if CoreDataInstance.Verbose:
        print("\nOTSO Cone Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(CoreDataInstance.year,CoreDataInstance.month,CoreDataInstance.day,
                         CoreDataInstance.hour,CoreDataInstance.minute,CoreDataInstance.second)
    README = cone_readme.READMECone(CoreDataInstance, EventDate, Printtime)

    if CoreDataInstance.livedata == "ON" or CoreDataInstance.livedata == 1:
        file_clean.remove_files()
    
    return [merged_cone_df, merged_R_df, README]