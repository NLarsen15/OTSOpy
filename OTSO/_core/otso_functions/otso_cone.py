import time
from datetime import datetime
import multiprocessing as mp
import pandas as pd
import sys
import queue
from tqdm import tqdm
import os
import tempfile
import json

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import cone_inputs
from ..fortran_calls import cone
from ..readme_generators import cone_readme
from ..data_classes.cone_data import ConeData


def OTSO_cone(ConeDataInstance: ConeData) -> list:

    File = tempfile.NamedTemporaryFile(delete=False, suffix=".json").name
    lock = mp.Lock()

    cone_inputs.ConeInputs(ConeDataInstance)

    ChildProcesses = []
    results = []

    total_cpus = os.cpu_count()

    next_cpu = 0

    UsedCores = cores.Cores(ConeDataInstance.station_array, ConeDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()

    if ConeDataInstance.Verbose:
        print("OTSO Cone Computation Started")

    total_stations = len(ConeDataInstance.station_array)
    results = []


    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass
    
    ProcessQueue = mp.Manager().Queue()
    for Data in Positionlists:
        threads = ConeDataInstance.threadnum
        cpus = set(range(next_cpu, next_cpu + threads))
        next_cpu += threads
        Child = mp.Process(target=cone.FortranCone,  args=(Data, ConeDataInstance, ProcessQueue, cpus, File, lock))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    processed = 0

    progress_bar = None
    if ConeDataInstance.Verbose and tqdm is not None:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cone")
    elif ConeDataInstance.Verbose:
        print(f"Processing {total_stations} stations...")

    while processed < total_stations:
        try:
            ProcessQueue.get(timeout=0.001)
            processed += 1

            if ConeDataInstance.Verbose:
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

    conedf_list = []
    Rigiditylist = []

    with open(File, "r", encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)

            conerec = record["cone"]
            conedf_list.append(
                pd.DataFrame(
                    conerec["data"],
                    index=conerec["index"],
                    columns=conerec["columns"],
                )
            )

            rig = record["rigidities"]
            Rigiditylist.append(
                pd.DataFrame(
                    rig["data"],
                    index=rig["index"],
                    columns=rig["columns"],
                )
            )

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

    if ConeDataInstance.Verbose:
        print("\nOTSO Cone Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(ConeDataInstance.year,ConeDataInstance.month,ConeDataInstance.day,
                         ConeDataInstance.hour,ConeDataInstance.minute,ConeDataInstance.second)
    README = cone_readme.READMECone(ConeDataInstance, EventDate, Printtime)

    if ConeDataInstance.livedata == "ON" or ConeDataInstance.livedata == 1:
        file_clean.remove_files()

    if os.path.exists(File):
        os.remove(File)
    
    return [merged_cone_df, merged_R_df, README]