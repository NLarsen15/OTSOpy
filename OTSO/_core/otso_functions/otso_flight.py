import time
import multiprocessing as mp
import os
import pandas as pd
import sys
import queue
import numpy as np
import tempfile
from tqdm import tqdm

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import flight_inputs
from ..fortran_calls import flight
from ..readme_generators import flight_readme
from ..data_classes.flight_data import FlightData


def OTSO_flight(FlightDataInstance: FlightData) -> list:

    flight_inputs.FlightInputs(FlightDataInstance)

    #RigidityArray = FlightInputArray[0]
    #DateArray = FlightInputArray[1]
    #Model = FlightInputArray[2]
    #IntModel = FlightInputArray[3]
    #ParticleArray = FlightInputArray[4]
    #IOPT = FlightInputArray[5]
    #WindArray = FlightInputArray[6]
    #Magnetopause = FlightInputArray[7]
    #CoordinateSystem = FlightInputArray[8]
    #MaxStepPercent = FlightInputArray[9]/100
    #EndParams = FlightInputArray[10]
    #Station_Array = FlightInputArray[11]
    #Rcomp = FlightInputArray[12]
    #Rscan = FlightInputArray[13]
    #KpList = FlightInputArray[14]
    #corenum = FlightInputArray[15]
    #LiveData = FlightInputArray[16]
    #serverdata = FlightInputArray[17]
    #g = FlightInputArray[18]
    #h = FlightInputArray[19]

    Glist = FlightDataInstance.glist
    Hlist = FlightDataInstance.hlist

    ChildProcesses = []

    UsedCores = cores.Cores(FlightDataInstance.station_array, 
                            FlightDataInstance.corenum)
    Positionlists = UsedCores.getPositions()
    WindLists = np.array_split(FlightDataInstance.windarraylist, FlightDataInstance.corenum)
    IOPTLists = np.array_split(FlightDataInstance.IOPTlist, FlightDataInstance.corenum)
    DateArrayLists = np.array_split(FlightDataInstance.datearraylist, FlightDataInstance.corenum)
    GArrayLists = np.array_split(Glist, FlightDataInstance.corenum)
    HArrayLists = np.array_split(Hlist, FlightDataInstance.corenum)

    flight_list = [tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name for _ in range(FlightDataInstance.corenum)]

    start = time.time()
    if FlightDataInstance.Verbose:
        print("OTSO Flight Computation Started")

    resultsfinal = []
    processed = 0
    totalp = 0
    total_stations = len(FlightDataInstance.station_array)

    if FlightDataInstance.corenum == 1:
        num_batches = 1
        
        progress_bar = None
        if FlightDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=num_batches, desc="OTSO Running", unit=" batch", position=0)
        elif FlightDataInstance.Verbose:
            print(f"Processing flight paths for {total_stations} stations...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        
        for Data, Date, I, Wind, flightFile, G, H in zip(Positionlists, DateArrayLists, IOPTLists, WindLists, flight_list, GArrayLists, HArrayLists):
            flight.FortranFlight(Data, Date, I, Wind, G, H, flightFile, FlightDataInstance, simple_queue)
            
            processed += 1
            if simple_queue.items:
                totalp = totalp + simple_queue.items[-1]
            
            if FlightDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(1)
                else:
                    percent_complete = (processed / num_batches) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% batches complete ({totalp} points processed)")
                    sys.stdout.flush()
        
        if progress_bar is not None:
            progress_bar.close()

    else:
        progress_bar = None
        if FlightDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" location", position=0)
        elif FlightDataInstance.Verbose:
            print(f"Processing flight paths for {total_stations} stations...")

        try:
            if not mp.get_start_method(allow_none=True):
                 mp.set_start_method('spawn')
        except RuntimeError:
             pass

        ProcessQueue = mp.Manager().Queue()
        for Data, Date, I, Wind, flightFile, G, H in zip(Positionlists,DateArrayLists, IOPTLists, WindLists, flight_list, GArrayLists, HArrayLists):
            Child = mp.Process(target=flight.FortranFlight,  args=(Data, Date, I, Wind, G, H, flightFile, FlightDataInstance, ProcessQueue))
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
        
                if result_collector:
                    totalp = totalp + sum(result_collector)
                
                if FlightDataInstance.Verbose:
                    if progress_bar is not None:
                        progress_bar.update(len(result_collector) if result_collector else 0)
                    else:
                        percent_complete = (processed / total_stations) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% batches complete ({totalp} points processed)")
                        sys.stdout.flush()

        
            except queue.Empty:
                pass
          
            time.sleep(0.0001)

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

    if FlightDataInstance.Verbose:
        print("\nOTSO Flight Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")


    readme = flight_readme.READMEFlight(FlightDataInstance, Printtime)

    datareadme = flight_readme.READMEFlightData(FlightDataInstance.datearraylist,FlightDataInstance.windarraylist,FlightDataInstance.Kplist)
    
    if FlightDataInstance.livedata == "ON" or FlightDataInstance.livedata == 1:
        file_clean.remove_files()

    return [merged_df,readme,datareadme]