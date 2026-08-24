import time
from datetime import datetime
import multiprocessing as mp
import sys
import queue
import random
import numpy as np
import gc
import os
import tempfile
from tqdm import tqdm
import pandas as pd
import json

from ..livedata import file_clean
from ..inputs import planet_inputs
from ..inputs import trace_inputs
from ..fortran_calls import trace
from ..readme_generators import trace_readme
from ..data_classes.trace_data import TraceData

def OTSO_trace(TraceDataInstance: TraceData) -> list:

    trace_inputs.TraceInputs(TraceDataInstance)

    gc.enable()

    available_cpus = sorted(os.sched_getaffinity(0))
    num_available = len(available_cpus)

    next_cpu = 0

    combined_coordinates = [(lat, lon) for lat in TraceDataInstance.latitudelist for lon in TraceDataInstance.longitudelist]

    ChildProcesses = []

    totalprocesses = len(TraceDataInstance.longitudelist)*len(TraceDataInstance.latitudelist)

    NewCoreNum = planet_inputs.CheckCoreNumPlanet(TraceDataInstance.corenum)
    
    DataTrace = []
    i = 0
    for point in combined_coordinates:
        DataTrace.append([f"{point[0]}_{point[1]}", point[0], point[1], TraceDataInstance.startaltitude, TraceDataInstance.zenith, TraceDataInstance.azimuth]) 
        i = i + 1

    if totalprocesses == 0:
        print("\nWarning: No coordinate pairs provided or generated. Skipping planet computation.")
        EventDate = datetime(TraceDataInstance.year,TraceDataInstance.month,TraceDataInstance.day,TraceDataInstance.hour,TraceDataInstance.minute,TraceDataInstance.second)
        readme = trace_readme.READMETrace(TraceDataInstance, EventDate, None,
                                        custom_coords_provided=(TraceDataInstance.array_of_lats_and_longs is not None))
        return [pd.DataFrame(), readme]
        
    actual_cores_to_use = min(NewCoreNum, totalprocesses)

    worker_files = generate_worker_files(actual_cores_to_use)

    shuffled_list = DataTrace.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, actual_cores_to_use)

    start = time.time()

    processed = 0
    totalp = 0

    if TraceDataInstance.Verbose:
        print("OTSO Trace Computation Started")

    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass

    ProcessQueue = mp.Queue()

    ChildProcesses = []
    for worker_id, Data in enumerate(DataLists):

            Child = mp.Process(
                    target=trace.FortranTrace,
                    args=(Data, TraceDataInstance, ProcessQueue, worker_files[worker_id])
                )
            ChildProcesses.append(Child)
        
    for a in ChildProcesses:
        a.start()

    progress_bar = None
    if TraceDataInstance.Verbose and tqdm is not None:
        progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
    elif TraceDataInstance.Verbose:
        print(f"Processing {totalprocesses} grid points...")
 
    while processed < totalprocesses:
        if all(not p.is_alive() for p in ChildProcesses):
            break
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
            
            gc.collect()
            if TraceDataInstance.Verbose:
                if progress_bar is not None:
                    progress_bar.update(sum(result_collector) if result_collector else 0)
                    progress_bar.set_description(f"OTSO Running ({totalp}/{totalprocesses})")
                else:
                    percent_complete = (totalp / totalprocesses) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete ({totalp}/{totalprocesses} points)")
                    sys.stdout.flush()

    
        except queue.Empty:
            pass
        
        time.sleep(0.0001)

    if progress_bar is not None:
        progress_bar.close()

    for b in ChildProcesses:
        b.join()

    processed = 0
    
    ChildProcesses.clear()

    trace_dict = {}
    lshell_records = []

    if TraceDataInstance.inputcoord == "GDZ":
        altunit = "km"
    else:
        altunit = "Re"


    for file in worker_files:

        with open(file, "r", encoding="utf-8") as f:

            for line in f:

                record = json.loads(line)

                name = record["location"]
                altitude = record[f"altitude [{altunit}]"]

                # Rebuild trace dataframe
                Trace = pd.DataFrame(
                    record["Trace"]["data"],
                    columns=record["Trace"]["columns"]
                )

                # Store full trace dictionary
                trace_dict[name] = {
                    f"altitude [{altunit}]": altitude,
                    "Trace": Trace
                }

                # Store L-shell summary information
                lshell_records.append({
                    "location": name,
                    f"altitude [{altunit}]": altitude,
                    "L_shell": record["L_shell"],
                    "Invariant_Latitude": record["Invariant_Latitude"]
                })

        os.remove(file)


    # Convert summary records to dataframe, sorted by location (latitude then longitude,
    # matching the planet grid output convention) rather than by Invariant_Latitude.
    lshell_df = pd.DataFrame(lshell_records)
    sort_keys = lshell_df["location"].apply(parse_key)
    lshell_df = (
        lshell_df
        .assign(_lat_sort=sort_keys.map(lambda k: k[0]), _lon_sort=sort_keys.map(lambda k: k[1]))
        .sort_values(["_lat_sort", "_lon_sort"])
        .drop(columns=["_lat_sort", "_lon_sort"])
        .reset_index(drop=True)
    )

    stop = time.time()
    Printtime = round((stop-start),3)

    if TraceDataInstance.Verbose:
        print("\nOTSO Trace Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(TraceDataInstance.year,TraceDataInstance.month,TraceDataInstance.day,
                         TraceDataInstance.hour,TraceDataInstance.minute,
                         TraceDataInstance.second)
    
    readme = trace_readme.READMETrace(TraceDataInstance, EventDate, Printtime)

    if TraceDataInstance.livedata == "ON" or TraceDataInstance.livedata == 1:
        file_clean.remove_files()

    output = [trace_dict, lshell_df, readme]

    return output


def parse_key(key):
    lat, lon = key.split('_')
    return (int(lat), int(lon))



def generate_worker_files(num_files: int, suffix=".json") -> list:
    """
    Generate unique temporary files for multiprocessing workers.

    Parameters
    ----------
    num_files : int
        Number of worker files to create.
    suffix : str
        File extension.

    Returns
    -------
    list
        List of unique file paths.
    """

    files = []

    for i in range(num_files):
        tmp = tempfile.NamedTemporaryFile(
            delete=False,
            suffix=f"_{i}{suffix}"
        )

        files.append(tmp.name)
        tmp.close()

    return files