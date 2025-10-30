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


    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:

        pass
# Create a shared message queue for the processes to produce/consume data
    ProcessQueue = mp.Manager().Queue()
    for Data,Date,g,h in zip(LocationsList,DateArrayList,glist,hlist):
        Child = mp.Process(target=fortran_calls.fortrancallCoordtrans,  args=(Data, Date, CoordIN, CoordOUT, ProcessQueue, g, h))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

# Wait for child processes to complete

    results = []
    total_stations = len(Locations)
    processed = 0

    # Initialize progress bar if tqdm is available and Verbose is True
    progress_bar = None
    if Verbose and tqdm is not None:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" transformation")
    elif Verbose:
        # Fallback to simple counter if tqdm is not available
        print(f"Processing {total_stations} coordinates...")

    while processed < total_stations:
      try:
        # Check if the ProcessQueue has any new results
        result_df = ProcessQueue.get(timeout=0.001)  # Use timeout to avoid blocking forever
        results.append(result_df)
        processed += 1

        # Update progress
        if Verbose:
            if progress_bar is not None:
                progress_bar.update(1)
            else:
                # Fallback to percentage if tqdm is not available
                percent_complete = (processed / total_stations) * 100
                sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                sys.stdout.flush()

      except queue.Empty:
        # Queue is empty, but processes are still running, so we continue checking
        pass
      
      time.sleep(0.0001)

    # Close progress bar if it was created
    if progress_bar is not None:
        progress_bar.close()

    # Ensure that all processes have completed
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