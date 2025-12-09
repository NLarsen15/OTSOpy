import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators, misc, date
import pandas as pd
import sys
import queue
import numpy as np
from tqdm import tqdm
from .igrf_process import compute_gauss_coefficients

def OTSO_coordtrans(Locations,Dates,CoordIN,CoordOUT,corenum,Verbose):
    
    if CoordIN not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordIN: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    if CoordOUT not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordOUT: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

    ChildProcesses = []
    results = []
    DateArrayList = []
    hlist = []
    glist = []

    for x in Dates:
          DateCreate = date.Date(x)
          DateArray = DateCreate.GetDate()
          DateArrayList.append(DateArray)
          gausscoefs = compute_gauss_coefficients(DateArray)
          hlist.append(gausscoefs['H_coefficients'])
          glist.append(gausscoefs['G_coefficients'])

    LocationsList = np.array_split(Locations, corenum)
    DateArrayList = np.array_split(DateArrayList, corenum)


    start = time.time()

    if Verbose:
        print("OTSO Coordtrans Computation Started")

    total_stations = len(Locations)
    results = []

    if corenum == 1:
        num_batches = len(LocationsList)
        
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" transformation")
        elif Verbose:
            print(f"Processing {total_stations} coordinates...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        for Data, Date, g, h in zip(LocationsList, DateArrayList, glist, hlist):
            fortran_calls.fortrancallCoordtrans(Data, Date, CoordIN, CoordOUT, simple_queue, g, h)
            
            processed += len(Data)
            if Verbose:
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
        for Data,Date,g,h in zip(LocationsList,DateArrayList,glist,hlist):
            Child = mp.Process(target=fortran_calls.fortrancallCoordtrans,  args=(Data, Date, CoordIN, CoordOUT, ProcessQueue, g, h))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" transformation")
        elif Verbose:
            print(f"Processing {total_stations} coordinates...")

        while processed < total_stations:
          try:
            result_df = ProcessQueue.get(timeout=0.001)
            results.append(result_df)
            processed += 1

            if Verbose:
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
    sorted_df = combined_df.sort_values(by=combined_df.columns[:4].tolist())

    stop = time.time()
    Printtime = round((stop-start),3)
    if Verbose:
        print("\nOTSO Coordtrans Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    README = readme_generators.READMECoordtrans(CoordIN,CoordOUT,Printtime)
    
    return [sorted_df, README]