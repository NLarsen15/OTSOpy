import time
from datetime import datetime
import multiprocessing as mp
import sys
import queue
import random
import numpy as np
from tqdm import tqdm

from ..livedata import file_clean
from ..inputs import trace_inputs
from ..fortran_calls import trace
from ..readme_generators import trace_readme
from ..data_classes.trace_data import TraceData

def OTSO_trace(TraceDataInstance: TraceData) -> list:

    trace_inputs.TraceInputs(TraceDataInstance)

    #LongitudeList = TraceInputArray[0]
    #LatitudeList = TraceInputArray[1]
    #Rigidity = TraceInputArray[2]
    #ParticleArray = TraceInputArray[3]
    #DateArray = TraceInputArray[4]
    #Model = TraceInputArray[5]
    #IntModel = TraceInputArray[6]
    #IOPT = TraceInputArray[7]
    #WindArray = TraceInputArray[8]
    #Magnetopause = TraceInputArray[9]
    #MaxStepPercent = TraceInputArray[10]/100
    #EndParams = TraceInputArray[11]
    #Zenith = TraceInputArray[12]
    #Azimuth = TraceInputArray[13]
    #CoreNum = TraceInputArray[14]
    #Alt = TraceInputArray[15]
    #LiveData = TraceInputArray[16]
    #AntiCheck = TraceInputArray[17]
    #g = TraceInputArray[18]
    #h = TraceInputArray[19]

    ChildProcesses = []

    totalprocesses = len(TraceDataInstance.longitudelist)*len(TraceDataInstance.latitudelist)

    combined_coordinates = [(lat, lon) for lat in TraceDataInstance.latitudelist for lon in TraceDataInstance.longitudelist]

    FileNamesPlanet = []

    combined_coordinates = [(lat, lon) for lat in TraceDataInstance.latitudelist for lon in TraceDataInstance.longitudelist]
    for list in combined_coordinates:
        FileNamesPlanet.append(str(list[0]) + "_" + str(list[1]))
    DataPlanet = []
    i = 1
    for point,name in zip(combined_coordinates, FileNamesPlanet):
        DataPlanet.append([name,point[0],point[1],TraceDataInstance.startaltitude,TraceDataInstance.zenith,TraceDataInstance.azimuth])
        i = i + 1

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, TraceDataInstance.corenum)
    start = time.time()

    if TraceDataInstance.Verbose:
        print("OTSO Trace Computation Started")

    results = {}

    if TraceDataInstance.corenum == 1:
        
        progress_bar = None
        if TraceDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" trace")
        elif TraceDataInstance.Verbose:
            print(f"Processing {totalprocesses} grid points...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        processed = 0
        
        for Data in DataLists:
            for single_item in Data:
                trace.FortranTrace([single_item], TraceDataInstance, simple_queue)
                
                num_items = len(simple_queue.items)
                for x in simple_queue.items:
                    results.update(x)
                    processed += 1
                
                if TraceDataInstance.Verbose:
                    if progress_bar is not None:
                        progress_bar.update(num_items)
                    else:
                        percent_complete = (processed / totalprocesses) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                        sys.stdout.flush()
                
                simple_queue.items = []
        
        if progress_bar is not None:
            progress_bar.close()

    else:
        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass

        ProcessQueue = mp.Manager().Queue()
        for Data in DataLists:
                Child = mp.Process(target=trace.FortranTrace,  args=(Data, TraceDataInstance, ProcessQueue))
                ChildProcesses.append(Child)

        for a in ChildProcesses:
            a.start()

        progress_bar = None
        if TraceDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" trace")
        elif TraceDataInstance.Verbose:
            print(f"Processing {totalprocesses} grid points...")

        processed = 0
        while processed < totalprocesses:
            try:
                result_collector = []
                while True:
                    try:
                        result_df = ProcessQueue.get_nowait()
                        result_collector.append(result_df)
                        processed += 1
                    except queue.Empty:
                        break
        
                for x in result_collector:
                    results.update(x)
        
                if TraceDataInstance.Verbose:
                    if progress_bar is not None:
                        progress_bar.update(len(result_collector))
                    else:
                        percent_complete = (processed / totalprocesses) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% complete")
                        sys.stdout.flush()
        
            except queue.Empty:
                pass
            
            time.sleep(0.0001)

        if progress_bar is not None:
            progress_bar.close()

        for b in ChildProcesses:
            b.join()

    stop = time.time()
    Printtime = round((stop-start),3)

    if TraceDataInstance.Verbose:
        print("\nOTSO Trace Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")

    # Sorting the dictionary by its keys
    sorted_results = dict(sorted(results.items(), key=lambda item: parse_key(item[0])))
    
    EventDate = datetime(TraceDataInstance.year,TraceDataInstance.month,TraceDataInstance.day,
                         TraceDataInstance.hour,TraceDataInstance.minute,
                         TraceDataInstance.second)
    
    readme = trace_readme.READMETrace(TraceDataInstance, EventDate, Printtime)

    if TraceDataInstance.livedata == "ON" or TraceDataInstance.livedata == 1:
        file_clean.remove_files()

    return [sorted_results,readme]


def parse_key(key):
    lat, lon = key.split('_')
    return (int(lat), int(lon))
