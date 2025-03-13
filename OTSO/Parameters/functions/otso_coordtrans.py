import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators, misc, date
import pandas as pd
import sys
import queue
import numpy as np

def OTSO_coordtrans(Locations,Dates,CoordIN,CoordOUT,corenum):
    
    if CoordIN not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordIN: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    if CoordOUT not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid CoordOUT: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

    ChildProcesses = []
    results = []
    DateArrayList = []

    for x in Dates:
          DateCreate = date.Date(x)
          DateArray = DateCreate.GetDate()
          DateArrayList.append(DateArray)

    LocationsList = np.array_split(Locations, corenum)
    DateArrayList = np.array_split(DateArrayList, corenum)


    start = time.time()
    print("OTSO Coordtrans Computation Started")
    sys.stdout.write(f"\r{0:.2f}% complete")


    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:

        pass
# Create a shared message queue for the processes to produce/consume data
    ProcessQueue = mp.Manager().Queue()
    for Data,Date in zip(LocationsList,DateArrayList):
        Child = mp.Process(target=fortran_calls.fortrancallCoordtrans,  args=(Data, Date, CoordIN, CoordOUT, ProcessQueue))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

# Wait for child processes to complete

    results = []
    total_stations = len(Locations)
    processed = 0

    while processed < total_stations:
      try:
        # Check if the ProcessQueue has any new results
        result_df = ProcessQueue.get(timeout=0.001)  # Use timeout to avoid blocking forever
        results.append(result_df)
        processed += 1

        # Calculate and print the progress
        percent_complete = (processed / total_stations) * 100
        sys.stdout.write(f"\r{percent_complete:.2f}% complete")
        sys.stdout.flush()

      except queue.Empty:
        # Queue is empty, but processes are still running, so we continue checking
        pass
      
      time.sleep(0.0001)

    # Ensure that all processes have completed
    for b in ChildProcesses:
        b.join()

    combined_df = pd.concat(results, ignore_index=True)
    sorted_df = combined_df.sort_values(by=combined_df.columns[:4].tolist())

    print("\nOTSO Coordtrans Computation Complete")
    stop = time.time()
    Printtime = round((stop-start),3)
    print("Whole Program Took: " + str(Printtime) + " seconds")
    
    README = readme_generators.READMECoordtrans(CoordIN,CoordOUT,Printtime)
    
    return [sorted_df, README]