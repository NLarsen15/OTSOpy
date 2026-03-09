import time
from datetime import datetime
import multiprocessing as mp
import pandas as pd
import sys
import queue
from tqdm import tqdm

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import cutoff_inputs
from ..fortran_calls import cutoff
from ..readme_generators import cutoff_readme
from ..data_classes.cutoff_data import CutoffData

def OTSO_cutoff(CutoffDataInstance: 'CutoffData') -> object:

    cutoff_inputs.CutoffInputs(CutoffDataInstance)

    ChildProcesses = []

    UsedCores = cores.Cores(CutoffDataInstance.station_array, CutoffDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()
    if CutoffDataInstance.Verbose:
        print("OTSO Cutoff Computation Started")

    total_stations = len(CutoffDataInstance.station_array)
    results = []
    asymptotic_results = []

    if CutoffDataInstance.corenum == 1:
        num_batches = len(Positionlists)
        
        progress_bar = None
        if CutoffDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=num_batches, desc="OTSO Running", unit=" batch")
        elif CutoffDataInstance.Verbose:
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
            result_list = cutoff.FortranCutoff(Data, CutoffDataInstance, simple_queue)
            # Since we're not using multiprocessing, we need to handle the return directly
            if result_list:  # FortranCutoff returns None, but puts results in queue
                pass
            
            processed += 1
            if CutoffDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(1)
                else:
                    percent_complete = (processed / num_batches) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()
        
        # Process all results from simple_queue
        all_results = simple_queue.get_all()
        for result_list in all_results:
            cutoff_df, asymptotic_df = result_list
            results.append(cutoff_df)
            if asymptotic_df is not None:
                asymptotic_results.append(asymptotic_df)
        
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
            Child = mp.Process(target=cutoff.FortranCutoff,  args=(Data, CutoffDataInstance, ProcessQueue))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if CutoffDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cutoff")
        elif CutoffDataInstance.Verbose:
            print(f"Processing {total_stations} stations...")

        while processed < total_stations:
          try:
            result_list = ProcessQueue.get(timeout=0.001)
            # result_list contains [cutoff_df, asymptotic_df] where asymptotic_df can be None
            cutoff_df, asymptotic_df = result_list
            results.append(cutoff_df)
            
            # Store asymptotic results separately if they exist
            if asymptotic_df is not None:
                asymptotic_results.append(asymptotic_df)
            
            processed += 1

            if CutoffDataInstance.Verbose:
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
    merged_df = pd.concat(results, axis=1)
    merged_df = merged_df.sort_index(axis=1)
    
    # Merge asymptotic results if any exist
    merged_asymptotic_df = None
    if asymptotic_results:
        merged_asymptotic_df = pd.concat(asymptotic_results, ignore_index=True)

    stop = time.time()
    Printtime = round((stop-start),3)

    if CutoffDataInstance.Verbose:
        print("\nOTSO Cutoff Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(CutoffDataInstance.year,CutoffDataInstance.month,CutoffDataInstance.day,
                         CutoffDataInstance.hour,CutoffDataInstance.minute,CutoffDataInstance.second)
    README = cutoff_readme.READMECutoff(CutoffDataInstance, EventDate, Printtime)
    
    if CutoffDataInstance.livedata == "ON" or CutoffDataInstance.livedata == 1:
        file_clean.remove_files()
        
    if CutoffDataInstance.asymptotic == "YES":
        return [merged_df, merged_asymptotic_df, README]
    else:
        return [merged_df, README]