import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, trajectory_inputs
import pandas as pd
import sys
import queue
import random
import numpy as np
from tqdm import tqdm

def OTSO_trajectory(Stations,rigidity, customlocations,startaltitude,
           minaltitude,zenith,azimuth,maxdistance,maxtime,serverdata,livedata,vx,vy,vz,bx,by,bz,density,
           pdyn,Dst,G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,intmodel,
           coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose):

    TrajectoryInputArray = trajectory_inputs.TrajectoryInputs(Stations,rigidity, customlocations,startaltitude,
           minaltitude,zenith,azimuth,maxdistance,maxtime,serverdata,livedata,vx,vy,vz,bx,by,bz,density,
           pdyn,Dst,G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,intmodel,
           coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile,MHDcoordsys,inputcoord)

    Rigidity = TrajectoryInputArray[0]
    DateArray = TrajectoryInputArray[1]
    Model = TrajectoryInputArray[2]
    IntModel = TrajectoryInputArray[3]
    ParticleArray = TrajectoryInputArray[4]
    IOPT = TrajectoryInputArray[5]
    WindArray = TrajectoryInputArray[6]
    Magnetopause = TrajectoryInputArray[7]
    CoordinateSystem = TrajectoryInputArray[8]
    MaxStepPercent = TrajectoryInputArray[9]/100
    EndParams = TrajectoryInputArray[10]
    Station_Array = TrajectoryInputArray[11]
    InputtedStations = TrajectoryInputArray[12]
    Kp = TrajectoryInputArray[13]
    CoreNum = TrajectoryInputArray[14]
    LiveData = TrajectoryInputArray[15]
    ServerData = TrajectoryInputArray[16]
    g = TrajectoryInputArray[17]
    h = TrajectoryInputArray[18]


    ChildProcesses = []

    UsedCores = cores.Cores(Station_Array, CoreNum)
    CoreList = UsedCores.getCoreList()
    Positionlists = UsedCores.getPositions()

    start = time.time()
    InputtedStations.find_non_matching_stations()
    if Verbose:
        print("OTSO Trajectory Computation Started")

# Set the process creation method to 'forkserver'
    try:
        # Check if the start method is already set
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        # If the start method is already set, a RuntimeError will be raised
        # You can log or handle this as needed
        pass
# Create a shared message queue for the processes to produce/consume data
    ProcessQueue = mp.Manager().Queue()
    for Data,Core in zip(Positionlists,CoreList):
        Child = mp.Process(target=fortran_calls.fortrancallTrajectory,  args=(Data, Core, Rigidity, DateArray, Model, IntModel, ParticleArray, IOPT, 
                                                                              WindArray, Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, 
                                                                              ProcessQueue,g,h, MHDfile, MHDcoordsys,spheresize,inputcoord))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    total_stations = len(Station_Array)
    processed = 0
    results = []

    # Initialize progress bar if tqdm is available and Verbose is True
    progress_bar = None
    if Verbose:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit="trajectory")

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

    for b in ChildProcesses:
        b.join()

# Wait for child processes to complete
    combined_dict = {}
    for d in results:
        combined_dict.update(d)

    stop = time.time()
    Printtime = round((stop-start),3)
    
    if Verbose:
        print("\nOTSO Trajectory Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")

    EventDate = datetime(year,month,day,hour,minute,second)
    AntiCheck = ParticleArray[0]
    Anum = ParticleArray[1]
    readme = readme_generators.READMETrajectory(Station_Array, Rigidity, EventDate, Model, IntModel, Anum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, Printtime, MaxStepPercent*100, EndParams, LiveData, Kp, serverdata)

    if livedata == "ON" or livedata == 1:
        misc.remove_files()

    return [combined_dict, readme]