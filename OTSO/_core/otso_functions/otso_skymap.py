import time
from datetime import datetime
import multiprocessing as mp
import os
import pandas as pd
import sys
import queue
import random
import numpy as np
import gc
import tempfile
from tqdm import tqdm
import json

from ..livedata import file_clean
from ..inputs import skymap_inputs
from ..fortran_calls import skymap
from ..readme_generators import skymap_readme
from ..data_classes.skymap_data import SkymapData

def OTSO_skymap(SkymapData: SkymapData) -> list:

    gc.enable()

    available_cpus = sorted(os.sched_getaffinity(0))
    num_available = len(available_cpus)

    next_cpu = 0

    skymap_inputs.SkymapInputs(SkymapData)
    
    ChildProcesses = []

    totalprocesses = len(SkymapData.station_array)

    NewCoreNum = SkymapData.corenum
    
    i = 0

    DataSkymap = []

    for station in SkymapData.station_array:
        name, lat, lon, alt, _, _ = station

        for zenith, azimuth in SkymapData.ZA_pairs:
            DataSkymap.append([
                f"{name}",
                lat,
                lon,
                alt,
                zenith,
                azimuth
            ])
        
    actual_cores_to_use = min(NewCoreNum, totalprocesses)

    worker_files = generate_worker_files(actual_cores_to_use)

    station_chunks = np.array_split(
        SkymapData.station_array,
        actual_cores_to_use
    )  # type: ignore

    start = time.time()

    processed = 0
    totalp = 0

    if SkymapData.Verbose:
        print("OTSO Skymap Computation Started")

    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass

    ProcessQueue = mp.Queue()

    ChildProcesses = []

    for worker_id, station_chunk in enumerate(station_chunks):

        # Build this worker's list of station/ZA combinations
        Data = []

        for name, lat, lon, alt, _, _ in station_chunk:
            for zenith, azimuth in SkymapData.ZA_pairs:
                Data.append([
                    name,
                    lat,
                    lon,
                    alt,
                    zenith,
                    azimuth
                ])

        threads = SkymapData.threadnum

        cpus = {
            available_cpus[(next_cpu + i) % num_available]
            for i in range(threads)
        }

        next_cpu = (next_cpu + threads) % num_available

        Child = mp.Process(
            target=skymap.FortranSkymap,
            args=(
                Data,
                SkymapData,
                ProcessQueue,
                cpus,
                worker_files[worker_id]
            )
        )

        ChildProcesses.append(Child)

    for child in ChildProcesses:
        child.start()

    progress_bar = None
    if SkymapData.Verbose and tqdm is not None:
        progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
    elif SkymapData.Verbose:
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
            if SkymapData.Verbose:
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

    skymap_results = {}

    for file in worker_files:

        worker_results = read_skymap_json(file)

        # merge stations from this worker
        skymap_results.update(worker_results)

        os.remove(file)

    stop = time.time()
    Printtime = round((stop-start),3)

    if SkymapData.Verbose:
        print("\nOTSO Skymap Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(SkymapData.year, SkymapData.month, SkymapData.day, SkymapData.hour, SkymapData.minute, SkymapData.second)

    readme = skymap_readme.READMESkymap(SkymapData, EventDate, Printtime)

    if SkymapData.livedata == "ON" or SkymapData.livedata == 1:
        file_clean.remove_files()

    output = [skymap_results]

    output.append(readme)

    return output

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

def read_skymap_json(JsonFile: str) -> dict[str, pd.DataFrame]:
    """
    Read worker skymap JSON and convert each station into a DataFrame.

    Returns
    -------
    dict[str, pd.DataFrame]
        Dictionary keyed by station name.
    """

    with open(JsonFile, "r", encoding="utf-8") as f:
        json_data = json.load(f)

    station_dfs = {}

    for station, records in json_data.items():
        station_dfs[station] = (
        pd.DataFrame(records)
        .sort_values(by=["Zenith", "Azimuth"])
        .reset_index(drop=True))

    return station_dfs