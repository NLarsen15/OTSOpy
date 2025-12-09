import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, cone_inputs
import pandas as pd
import sys
import queue
from tqdm import tqdm


def OTSO_cone(Stations,customlocations,startaltitude,minaltitude,zenith,azimuth,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,
           coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose,
           AdaptiveExternalModel):

    ConeInputArray = cone_inputs.ConeInputs(Stations,customlocations,startaltitude,minaltitude,zenith,azimuth,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,
           coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile,MHDcoordsys,inputcoord,AdaptiveExternalModel)

    RigidityArray = ConeInputArray[0]
    DateArray = ConeInputArray[1]
    Model = ConeInputArray[2]
    IntModel = ConeInputArray[3]
    ParticleArray = ConeInputArray[4]
    IOPT = ConeInputArray[5]
    WindArray = ConeInputArray[6]
    Magnetopause = ConeInputArray[7]
    CoordinateSystem = ConeInputArray[8]
    MaxStepPercent = ConeInputArray[9]/100
    EndParams = ConeInputArray[10]
    Station_Array = ConeInputArray[11]
    InputtedStations = ConeInputArray[12]
    Kp = ConeInputArray[13]
    CoreNum = ConeInputArray[14]
    LiveData = ConeInputArray[15]
    g = ConeInputArray[16]
    h = ConeInputArray[17]

    AntiCheck = ParticleArray[1]

    ChildProcesses = []
    results = []

    UsedCores = cores.Cores(Station_Array, CoreNum)
    CoreList = UsedCores.getCoreList()
    Positionlists = UsedCores.getPositions()

    start = time.time()
    InputtedStations.find_non_matching_stations()

    if Verbose:
        print("OTSO Cone Computation Started")

    total_stations = len(Station_Array)
    results = []

    if CoreNum == 1:
        # Single core processing - avoid multiprocessing overhead
        # When CoreNum=1, all stations are in one batch
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
            fortran_calls.fortrancallCone(Data, Core, RigidityArray, DateArray, Model, IntModel, ParticleArray, IOPT, WindArray, 
                                        Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, simple_queue, g, h, MHDfile, MHDcoordsys,
                                        spheresize, inputcoord)
            
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
        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass
        
        # Create a shared message queue for the processes to produce/consume data
        ProcessQueue = mp.Manager().Queue()
        for Data,Core in zip(Positionlists,CoreList):
            Child = mp.Process(target=fortran_calls.fortrancallCone,  args=(Data, Core, RigidityArray, DateArray, Model, IntModel, ParticleArray, IOPT, WindArray, 
                                                                            Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, ProcessQueue,g,h,MHDfile, MHDcoordsys,
                                                                            spheresize, inputcoord))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        # Wait for child processes to complete
        processed = 0

        # Initialize progress bar if tqdm is available and Verbose is True
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cone")
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

    if Verbose:
        print("\nOTSO Cone Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    README = readme_generators.READMECone(Station_Array, RigidityArray, EventDate, Model, IntModel, Anum, AntiCheck, IOPT, WindArray, Magnetopause, 
                                          CoordinateSystem, Printtime, MaxStepPercent*100, EndParams, LiveData, serverdata, Kp)

    if livedata == "ON" or livedata == 1:
        misc.remove_files()
    
    return [merged_cone_df, merged_R_df, README]