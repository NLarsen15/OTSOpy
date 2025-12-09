import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, trace_inputs
import pandas as pd
import sys
import queue
import random
import numpy as np
from tqdm import tqdm


def OTSO_trace(startaltitude,Coordsys,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,year,
           month,day,hour,minute,second,internalmag,externalmag,
           gyropercent,magnetopause,corenum,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose):

    Anum = 1
    TraceInputArray = trace_inputs.TraceInputs(startaltitude,Coordsys,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,year,
           month,day,hour,minute,second,internalmag,externalmag,
           gyropercent,magnetopause,corenum,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,inputcoord)

    LongitudeList = TraceInputArray[0]
    LatitudeList = TraceInputArray[1]
    Rigidity = TraceInputArray[2]
    ParticleArray = TraceInputArray[3]
    DateArray = TraceInputArray[4]
    Model = TraceInputArray[5]
    IntModel = TraceInputArray[6]
    IOPT = TraceInputArray[7]
    WindArray = TraceInputArray[8]
    Magnetopause = TraceInputArray[9]
    MaxStepPercent = TraceInputArray[10]/100
    EndParams = TraceInputArray[11]
    Zenith = TraceInputArray[12]
    Azimuth = TraceInputArray[13]
    CoreNum = TraceInputArray[14]
    Alt = TraceInputArray[15]
    LiveData = TraceInputArray[16]
    AntiCheck = TraceInputArray[17]
    g = TraceInputArray[18]
    h = TraceInputArray[19]

    ChildProcesses = []

    totalprocesses = len(LongitudeList)*len(LatitudeList)

    combined_coordinates = [(lat, lon) for lat in LatitudeList for lon in LongitudeList]

    NewCoreNum = misc.CheckCoreNumPlanet(CoreNum)
    FileNamesPlanet = []

    combined_coordinates = [(lat, lon) for lat in LatitudeList for lon in LongitudeList]
    for list in combined_coordinates:
        FileNamesPlanet.append(str(list[0]) + "_" + str(list[1]))
    DataPlanet = []
    i = 1
    for point,name in zip(combined_coordinates, FileNamesPlanet):
        Core = "Core " + str(i)
        DataPlanet.append([name,point[0],point[1],Alt,Zenith,Azimuth,Core])
        i = i + 1

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, CoreNum)
    CoreList = np.arange(1, CoreNum + 1)
    start = time.time()

    if Verbose:
        print("OTSO Trace Computation Started")

    results = {}

    if CoreNum == 1:
        # Single core processing - avoid multiprocessing overhead
        # When CoreNum=1, all points processed in batches
        num_batches = len(DataLists)
        
        # Initialize progress bar if tqdm is available and Verbose is True
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" trace")
        elif Verbose:
            # Fallback to simple counter if tqdm is not available
            print(f"Processing {totalprocesses} grid points...")

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
        
        # Process all data directly without multiprocessing
        # When single-core, process each item individually to show progress
        for Data, Core in zip(DataLists, CoreList):
            for single_item in Data:
                # Process one coordinate at a time
                fortran_calls.fortrancallTrace([single_item], Rigidity, DateArray, Model, IntModel, 
                                              ParticleArray, IOPT, WindArray, 
                                              Magnetopause, Coordsys, MaxStepPercent, EndParams,
                                              simple_queue, g, h, MHDfile, MHDcoordsys,
                                              spheresize, inputcoord)
                
                # Collect results from the simple queue and update progress
                num_items = len(simple_queue.items)
                for x in simple_queue.items:
                    results.update(x)
                    processed += 1
                
                # Update progress after each coordinate
                if Verbose:
                    if progress_bar is not None:
                        progress_bar.update(num_items)
                    else:
                        # Fallback to percentage if tqdm is not available
                        percent_complete = (processed / totalprocesses) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                        sys.stdout.flush()
                
                # Clear the simple queue for next item
                simple_queue.items = []
        
        # Close progress bar if it was created
        if progress_bar is not None:
            progress_bar.close()

    else:
        # Multi-core processing
        # Set the process creation method to 'forkserver'
        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass

        # Create a shared message queue for the processes to produce/consume data
        ProcessQueue = mp.Manager().Queue()
        for Data,Core in zip(DataLists, CoreList):
                Child = mp.Process(target=fortran_calls.fortrancallTrace,  args=(Data, Rigidity, DateArray, Model, IntModel, 
                                                                                  ParticleArray, IOPT, WindArray, 
                                                                                  Magnetopause, Coordsys, MaxStepPercent, EndParams,
                                                                                  ProcessQueue,g,h, MHDfile, MHDcoordsys,
                                                                                  spheresize,inputcoord))
                ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        # Initialize progress bar if tqdm is available and Verbose is True
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" trace")
        elif Verbose:
            # Fallback to simple counter if tqdm is not available
            print(f"Processing {totalprocesses} grid points...")

        processed = 0
        while processed < totalprocesses:
            try:
                # Collect all results available in the queue at this moment
                result_collector = []
                while True:
                    try:
                        result_df = ProcessQueue.get_nowait()  # Non-blocking, does not wait
                        result_collector.append(result_df)
                        processed += 1
                    except queue.Empty:
                        break
        
                # Append all collected results to the main results list
                for x in result_collector:
                    results.update(x)
        
                # Update progress
                if Verbose:
                    if progress_bar is not None:
                        progress_bar.update(len(result_collector))
                    else:
                        # Fallback to percentage if tqdm is not available
                        percent_complete = (processed / totalprocesses) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                        sys.stdout.flush()
        
            except queue.Empty:
                # Queue is empty, but processes are still running, so we continue checking
                pass
            
            # Wait for 5 seconds before the next iteration
            time.sleep(0.0001)

        # Close progress bar if it was created
        if progress_bar is not None:
            progress_bar.close()

        for b in ChildProcesses:
            b.join()

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
        print("\nOTSO Trace Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")

    # Sorting the dictionary by its keys
    sorted_results = dict(sorted(results.items(), key=lambda item: parse_key(item[0])))
    
    EventDate = datetime(year,month,day,hour,minute,second)
    readme = readme_generators.READMETrace(Data, EventDate, Model, IntModel, 
                                             AntiCheck, IOPT, WindArray, Magnetopause, Printtime,
                                             maxlat,maxlong,minlat,minlong, latstep, longstep,
                                             MaxStepPercent*100, EndParams, 
                                             LiveData, serverdata, kp)

    if livedata == "ON" or livedata == 1:
        misc.remove_files()

    return [sorted_results,readme]


def parse_key(key):
    lat, lon = key.split('_')
    return (int(lat), int(lon))
