import time
from datetime import datetime
import multiprocessing as mp
import pandas as pd
import sys
import queue
import numpy as np
from tqdm import tqdm

from ..livedata import file_clean
from ..inputs import magfield_inputs
from ..fortran_calls import magfield
from ..readme_generators import magfield_readme
from ..data_classes.magfield_data import MagfieldData

def OTSO_magfield(MagfieldDataInstance: MagfieldData) -> list:

    magfield_inputs.MagFieldInputs(MagfieldDataInstance)

    ChildProcesses = []
    results = []

    LocationsList = np.array_split(MagfieldDataInstance.locations,
                                   MagfieldDataInstance.corenum)

    start = time.time()

    if MagfieldDataInstance.Verbose:
        print("OTSO Magfield Computation Started")

    total_stations = len(MagfieldDataInstance.locations)
    results = []

    if MagfieldDataInstance.corenum == 1:
        
        progress_bar = None
        if MagfieldDataInstance.Verbose:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" field calculation")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        for Data in LocationsList:
            magfield.FortranMagfield(Data, MagfieldDataInstance, simple_queue)
            
            processed += len(Data)
            if MagfieldDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(len(Data))
                else:
                    percent_complete = (processed / total_stations) * 100
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
        for Data in LocationsList:
            Child = mp.Process(target=magfield.FortranMagfield,  args=(Data, MagfieldDataInstance, ProcessQueue))
            
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if MagfieldDataInstance.Verbose:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" field calculation")

        while processed < total_stations:
          try:
            result_df = ProcessQueue.get(timeout=0.001)
            results.append(result_df)
            processed += 1

            if MagfieldDataInstance.Verbose:
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

    combined_df = pd.concat(results, ignore_index=True)
    sorted_df = combined_df.sort_values(by=combined_df.columns[:3].tolist())

    stop = time.time()
    Printtime = round((stop-start),3)

    if MagfieldDataInstance.Verbose:
      print("\nOTSO Magfield Computation Complete")
      print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(MagfieldDataInstance.year, MagfieldDataInstance.month, 
                         MagfieldDataInstance.day, MagfieldDataInstance.hour, 
                         MagfieldDataInstance.minute, MagfieldDataInstance.second)
    
    README = magfield_readme.READMEMagfield(EventDate, MagfieldDataInstance, Printtime)

    if MagfieldDataInstance.livedata == "ON" or MagfieldDataInstance.livedata == 1:
        file_clean.remove_files()
    
    return [sorted_df,README]