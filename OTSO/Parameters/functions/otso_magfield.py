import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, magfield_inputs
import pandas as pd
import sys
import queue
import numpy as np

def OTSO_magfield(Locations,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,year,
           month,day,hour,minute,second,internalmag,externalmag,
           coordsystemIN,g,h,corenum,MHDfile,MHDcoordsys):

    magfieldInputArray = magfield_inputs.MagFieldInputs(Locations,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,year,
           month,day,hour,minute,second,internalmag,externalmag,
           coordsystemIN,g,h,corenum,MHDfile,MHDcoordsys)

    Locations = magfieldInputArray[0]
    DateArray = magfieldInputArray[1]
    Model = magfieldInputArray[2]
    IOPT = magfieldInputArray[3]
    WindArray = magfieldInputArray[4]
    CoordinateSystem = magfieldInputArray[5]
    Kp = magfieldInputArray[6]
    CoreNum = magfieldInputArray[7]
    LiveData = magfieldInputArray[8]
    serverdata = magfieldInputArray[9]
    g = magfieldInputArray[10]
    h = magfieldInputArray[11]

    ChildProcesses = []
    results = []

    LocationsList = np.array_split(Locations, CoreNum)

    start = time.time()
    print("OTSO Magfield Computation Started")
    sys.stdout.write(f"\r{0:.2f}% complete")


    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:

        pass
# Create a shared message queue for the processes to produce/consume data
    ProcessQueue = mp.Manager().Queue()
    for Data in LocationsList:
        Child = mp.Process(target=fortran_calls.fortrancallMagfield,  args=(Data, DateArray, Model, IOPT, WindArray, CoordinateSystem, ProcessQueue,g,h,
                                                                            MHDfile, MHDcoordsys))
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
    sorted_df = combined_df.sort_values(by=combined_df.columns[:3].tolist())

    print("\nOTSO Magfield Computation Complete")
    stop = time.time()
    Printtime = round((stop-start),3)
    print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    README = readme_generators.READMEMagfield(EventDate, Model, IOPT, WindArray,
                                          CoordinateSystem, Printtime, LiveData, serverdata, kp)

    if LiveData == 1:
        misc.remove_files()
    
    return [sorted_df,README]