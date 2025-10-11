import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, flight_inputs
import pandas as pd
import sys
import queue
import numpy as np
import tempfile
from tqdm import tqdm


def OTSO_flight(latitudes,longitudes,dates,altitudes,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,Anum,anti,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h,asymptotic,asymlevels,unit,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose):

    FlightInputArray = flight_inputs.FlightInputs(latitudes,longitudes,dates,altitudes,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,Anum,anti,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h,asymptotic,asymlevels,unit,MHDfile,MHDcoordsys,inputcoord)

    RigidityArray = FlightInputArray[0]
    DateArray = FlightInputArray[1]
    Model = FlightInputArray[2]
    IntModel = FlightInputArray[3]
    ParticleArray = FlightInputArray[4]
    IOPT = FlightInputArray[5]
    WindArray = FlightInputArray[6]
    Magnetopause = FlightInputArray[7]
    CoordinateSystem = FlightInputArray[8]
    MaxStepPercent = FlightInputArray[9]/100
    EndParams = FlightInputArray[10]
    Station_Array = FlightInputArray[11]
    Rcomp = FlightInputArray[12]
    Rscan = FlightInputArray[13]
    KpList = FlightInputArray[14]
    corenum = FlightInputArray[15]
    LiveData = FlightInputArray[16]
    serverdata = FlightInputArray[17]
    g = FlightInputArray[18]
    h = FlightInputArray[19]

    AntiCheck = ParticleArray[1]

    ChildProcesses = []

    UsedCores = cores.Cores(Station_Array, corenum)
    CoreList = UsedCores.getCoreList()
    Positionlists = UsedCores.getPositions()
    WindLists = np.array_split(WindArray, corenum)
    IOPTLists = np.array_split(IOPT, corenum)
    DateArrayLists = np.array_split(DateArray, corenum)


    current_dir = os.path.dirname(os.path.realpath(__file__))
    flight_list = [tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name for _ in range(corenum)]

    start = time.time()
    if Verbose:
        print("OTSO Flight Computation Started")

    results = []
    resultsfinal = []
    processed = 0
    totalp = 0
    total_stations = len(Station_Array)

    # Initialize progress bar if tqdm is available and Verbose is True
    progress_bar = None
    if Verbose and tqdm is not None:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" locations", position=0)
    elif Verbose:
        # Fallback to simple counter if tqdm is not available
        print(f"Processing flight paths for {total_stations} stations...")

    try:
        if not mp.get_start_method(allow_none=True):
             mp.set_start_method('spawn')
    except RuntimeError:
         pass

    ProcessQueue = mp.Manager().Queue()
    for Data, Core, Date, I, Wind, flightFile in zip(Positionlists,CoreList,DateArrayLists, IOPTLists, WindLists, flight_list):
        Child = mp.Process(target=fortran_calls.fortrancallFlight,  args=(Data, RigidityArray, Date, Model, IntModel, 
                                                                              ParticleArray, I, Wind, 
                                                                              Magnetopause, MaxStepPercent, EndParams, 
                                                                              Rcomp, Rscan, asymptotic, asymlevels, unit,
                                                                              ProcessQueue,g,h,CoordinateSystem, flightFile,  MHDfile, MHDcoordsys,
                                                                              spheresize,inputcoord))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    while processed < total_stations:
        try:
            result_collector = []
            while True:
                try:
                    countint = ProcessQueue.get(timeout=0.001)
                    result_collector.append(countint)
                    processed += 1
                except queue.Empty:
                    break
    
            # Update totalp with the sum of items processed by cores
            if result_collector:
                totalp = totalp + sum(result_collector)
            
            # Update progress
            if Verbose:
                if progress_bar is not None:
                    # Update batch progress bar
                    progress_bar.update(len(result_collector) if result_collector else 0)
                    #progress_bar.set_description(f"Flight batches ({totalp} points total)")
                else:
                    # Fallback to percentage if tqdm is not available
                    percent_complete = (processed / total_stations) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% batches complete ({totalp} points processed)")
                    sys.stdout.flush()

    
        except queue.Empty:
            pass
      
        time.sleep(0.1)  # Update every 100ms for smoother progress display

    # Close progress bar if it was created
    if progress_bar is not None:
        progress_bar.close()

    for b in ChildProcesses:
        b.join()
        b.close()

    for x in flight_list:
        df = pd.read_csv(x)
        resultsfinal.append(df)
        os.remove(x)

    merged_df = pd.concat(resultsfinal, ignore_index=True)
    merged_df = merged_df.sort_index(axis=0)

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
        print("\nOTSO Flight Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")


    readme = readme_generators.READMEFlight(Data, RigidityArray, Model, IntModel,
                                            AntiCheck, IOPT, WindArray, Magnetopause, Printtime,
                                            MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                            LiveData, asymptotic, asymlevels, unit, serverdata, kp)

    
    datareadme = readme_generators.READMEFlightData(DateArray,WindArray,KpList)
    
    if livedata == 1:
        misc.remove_files()

    return [merged_df,readme,datareadme]