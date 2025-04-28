import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, flight_inputs
import pandas as pd
import sys
import queue
import numpy as np

def OTSO_flight(latitudes,longitudes,dates,altitudes,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,Anum,anti,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h,asymptotic,asymlevels,unit,MHDfile,MHDcoordsys):

    FlightInputArray = flight_inputs.FlightInputs(latitudes,longitudes,dates,altitudes,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,Anum,anti,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           coordsystem,gyropercent,magnetopause,corenum,azimuth,zenith,g,h,asymptotic,asymlevels,unit,MHDfile,MHDcoordsys)

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
    flight_list = [os.path.join(current_dir, f"Flight{i}.csv") for i in range(1, (corenum) + 1)]
    
    print("OTSO Flight Computation Started")
    start = time.time()
    sys.stdout.write(f"\r{0:.2f}% complete")

    results = []
    resultsfinal = []
    processed = 0
    totalp = 0
    total_stations = len(Station_Array)

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
                                                                              ProcessQueue,g,h,CoordinateSystem, flightFile,  MHDfile, MHDcoordsys))
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
                    totalp = totalp + sum(result_collector)
                    result_collector = []
                except queue.Empty:
                    break
    
            percent_complete = (totalp / total_stations) * 100
            sys.stdout.write(f"\r{percent_complete:.2f}% complete")
            sys.stdout.flush()

    
        except queue.Empty:
            pass
      
        time.sleep(0.01)

    for b in ChildProcesses:
        b.join()
        b.close()

    for x in flight_list:
        df = pd.read_csv(x)
        resultsfinal.append(df)
        os.remove(x)

    merged_df = pd.concat(resultsfinal, ignore_index=True)
    merged_df = merged_df.sort_index(axis=0)

    print("\nOTSO Flight Computation Complete")
    stop = time.time()
    Printtime = round((stop-start),3)
    print("Whole Program Took: " + str(Printtime) + " seconds")

    

    readme = readme_generators.READMEFlight(Data, RigidityArray, Model, IntModel, 
                                            AntiCheck, IOPT, WindArray, Magnetopause, Printtime,
                                            MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                            LiveData, asymptotic, asymlevels, unit, serverdata, kp)

    
    datareadme = readme_generators.READMEFlightData(DateArray,WindArray,KpList)
    
    if livedata == 1:
        misc.remove_files()

    return [merged_df,readme,datareadme]