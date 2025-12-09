import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, cutoff_inputs
import pandas as pd
import sys
import queue
from tqdm import tqdm

def OTSO_cutoff(Stations,customlocations,startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h, MHDfile, MHDcoordsys,spheresize,inputcoord,Verbose,
           AdaptiveExternalModel):

    CutoffInputArray = cutoff_inputs.CutoffInputs(Stations,customlocations,startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h, MHDfile, MHDcoordsys,inputcoord,AdaptiveExternalModel)

    RigidityArray = CutoffInputArray[0]
    DateArray = CutoffInputArray[1]
    Model = CutoffInputArray[2]
    IntModel = CutoffInputArray[3]
    ParticleArray = CutoffInputArray[4]
    IOPT = CutoffInputArray[5]
    WindArray = CutoffInputArray[6]
    Magnetopause = CutoffInputArray[7]
    CoordinateSystem = CutoffInputArray[8]
    MaxStepPercent = CutoffInputArray[9]/100
    EndParams = CutoffInputArray[10]
    Station_Array = CutoffInputArray[11]
    InputtedStations = CutoffInputArray[12]
    Rcomp = CutoffInputArray[13]
    Rscan = CutoffInputArray[14]
    Kp = CutoffInputArray[15]
    corenum = CutoffInputArray[16]
    LiveData = CutoffInputArray[17]
    g = CutoffInputArray[18]
    h = CutoffInputArray[19]

    AntiCheck = ParticleArray[1]

    ChildProcesses = []

    UsedCores = cores.Cores(Station_Array, corenum)
    CoreList = UsedCores.getCoreList()
    Positionlists = UsedCores.getPositions()
    InputtedStations.find_non_matching_stations()

    start = time.time()
    if Verbose:
        print("OTSO Cutoff Computation Started")

    total_stations = len(Station_Array)
    results = []

    if corenum == 1:
        # Single core processing - avoid multiprocessing overhead
        # When corenum=1, all stations are in one batch
        num_batches = len(Positionlists)
        
        # Initialize progress bar if tqdm is available and Verbose is True
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=num_batches, desc="OTSO Running", unit=" batch")
        elif Verbose:
            # Fallback to simple counter if tqdm is not available
            print(f"Processing {total_stations} stations...")

        # Create a simple queue-like list for single-core processing
        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        # Process all stations directly without multiprocessing
        for Data, Core in zip(Positionlists, CoreList):
            fortran_calls.fortrancallCutoff(Data, Core, RigidityArray, DateArray, Model, IntModel, ParticleArray, IOPT,
                WindArray, Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, Kp, simple_queue, g, h, MHDfile, MHDcoordsys, spheresize, inputcoord)
            
            # Update progress after each station
            processed += 1
            if Verbose:
                if progress_bar is not None:
                    progress_bar.update(1)
                else:
                    # Fallback to percentage if tqdm is not available
                    percent_complete = (processed / num_batches) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                    sys.stdout.flush()
        
        # Get all results
        results = simple_queue.get_all()
        
        # Close progress bar if it was created
        if progress_bar is not None:
            progress_bar.close()

    else:
        # Multi-core processing
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
            Child = mp.Process(target=fortran_calls.fortrancallCutoff,  args=(Data, Core, RigidityArray, DateArray, Model, IntModel, ParticleArray, IOPT,
             WindArray, Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, Kp, ProcessQueue, g, h,  MHDfile, MHDcoordsys,spheresize,inputcoord))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        # Wait for child processes to complete
        processed = 0

        # Initialize progress bar if tqdm is available and Verbose is True
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cutoff")
        elif Verbose:
            # Fallback to simple counter if tqdm is not available
            print(f"Processing {total_stations} stations...")

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


    # Concatenate the results on the index
    merged_df = pd.concat(results, axis=1)
    merged_df = merged_df.sort_index(axis=1)

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
        print("\nOTSO Cutoff Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    README = readme_generators.READMECutoff(Station_Array, RigidityArray, EventDate, Model, IntModel, Anum, AntiCheck, IOPT, WindArray, Magnetopause, 
                                            CoordinateSystem, Printtime, MaxStepPercent*100, EndParams, cutoff_comp, Rscan, LiveData, serverdata, Kp)
    
    if livedata == "ON" or livedata == 1:
        misc.remove_files()
        
    return [merged_df, README]