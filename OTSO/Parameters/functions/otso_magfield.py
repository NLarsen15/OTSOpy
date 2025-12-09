import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, magfield_inputs
import pandas as pd
import sys
import queue
import numpy as np
from tqdm import tqdm

def OTSO_magfield(Locations,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,year,
           month,day,hour,minute,second,internalmag,externalmag,
           coordsystemIN,g,h,corenum,MHDfile,MHDcoordsys,Verbose):

    magfieldInputArray = magfield_inputs.MagFieldInputs(Locations,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,year,
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

    if Verbose:
        print("OTSO Magfield Computation Started")

    total_stations = len(Locations)
    results = []

    if CoreNum == 1:
        num_batches = len(LocationsList)
        
        progress_bar = None
        if Verbose:
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
            fortran_calls.fortrancallMagfield(Data, DateArray, Model, IOPT, WindArray, CoordinateSystem, simple_queue, g, h,
                                            MHDfile, MHDcoordsys)
            
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
        for Data in LocationsList:
            Child = mp.Process(target=fortran_calls.fortrancallMagfield,  args=(Data, DateArray, Model, IOPT, WindArray, CoordinateSystem, ProcessQueue,g,h,
                                                                                MHDfile, MHDcoordsys))
            ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        processed = 0

        progress_bar = None
        if Verbose:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" field calculation")

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
    sorted_df = combined_df.sort_values(by=combined_df.columns[:3].tolist())

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
      print("\nOTSO Magfield Computation Complete")
      print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    README = readme_generators.READMEMagfield(EventDate, Model, IOPT, WindArray,
                                          CoordinateSystem, Printtime, LiveData, serverdata, kp)

    if livedata == "ON" or livedata == 1:
        misc.remove_files()
    
    return [sorted_df,README]