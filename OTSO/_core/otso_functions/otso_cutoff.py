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
from ..inputs import cutoff_inputs
from ..fortran_calls import cutoff
from ..readme_generators import cutoff_readme
from ..data_classes.cutoff_data import CutoffData

def OTSO_cutoff(CutoffDataInstance: 'CutoffData') -> object:

    cutoff_inputs.CutoffInputs(CutoffDataInstance)

    File = tempfile.NamedTemporaryFile(delete=False, suffix=".json").name
    lock = mp.Lock()

    ChildProcesses = []

    UsedCores = cores.Cores(CutoffDataInstance.station_array, CutoffDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()
    if CutoffDataInstance.Verbose:
        print("OTSO Cutoff Computation Started")

    total_stations = len(CutoffDataInstance.station_array)
    results = []
    asymptotic_results = []

    next_cpu = 0

    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass
    
    ProcessQueue = mp.Manager().Queue()
    for Data in Positionlists:
        threads = CutoffDataInstance.threadnum
        cpus = set(range(next_cpu, next_cpu + threads))
        next_cpu += threads
        Child = mp.Process(target=cutoff.FortranCutoff,  args=(Data, CutoffDataInstance, ProcessQueue, cpus, File, lock))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    processed = 0

    progress_bar = None
    if CutoffDataInstance.Verbose and tqdm is not None:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cutoff")
    elif CutoffDataInstance.Verbose:
        print(f"Processing {total_stations} stations...")

    while processed < total_stations:
      try:
        result_list = ProcessQueue.get(timeout=0.001)
        
        processed += 1

        if CutoffDataInstance.Verbose:
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


    results = []
    asymptotic_results = []
    transmission_results = []

    with open(File, "r", encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)

            # Recover rigidity dataframe
            rig = record["rigidities"]
            Rigiditydataframe = pd.DataFrame(
                rig["data"],
                index=rig["index"],
                columns=rig["columns"],
            )
            results.append(Rigiditydataframe)

            # Recover asymptotic dataframe if present
            if "asymptotic" in record:
                asymptotic_results.append(
                    pd.DataFrame(record["asymptotic"])
                )

            # Recover transmission dataframe if present
            if "transmission" in record:
                trans = record["transmission"]
                transmission_results.append(
                    pd.DataFrame(
                        trans["data"],
                        index=trans["index"],
                        columns=trans["columns"],
                    )
                )

        # ------------------------------------------------------------------
        # Merge rigidity results
        # ------------------------------------------------------------------

        merged_df = pd.concat(results, axis=1)
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]
        merged_df = merged_df.sort_index(axis=1)

        # ------------------------------------------------------------------
        # Merge asymptotic results
        # ------------------------------------------------------------------

        merged_asymptotic_df = None
        if asymptotic_results:
            merged_asymptotic_df = pd.concat(asymptotic_results, ignore_index=True)

        # ------------------------------------------------------------------
        # Merge transmission results
        # ------------------------------------------------------------------

        merged_transmission_df = None
        if transmission_results:
            merged_transmission_df = transmission_results[0]

            for df in transmission_results[1:]:
                merged_transmission_df = merged_transmission_df.merge(
                    df,
                    on="R [GV]",
                    how="outer",
                )

            merged_transmission_df = (
                merged_transmission_df
                .sort_values("R [GV]")
                .reset_index(drop=True)
            )

    stop = time.time()
    Printtime = round((stop-start),3)

    if CutoffDataInstance.Verbose:
        print("\nOTSO Cutoff Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(CutoffDataInstance.year,CutoffDataInstance.month,CutoffDataInstance.day,
                         CutoffDataInstance.hour,CutoffDataInstance.minute,CutoffDataInstance.second)
    README = cutoff_readme.READMECutoff(CutoffDataInstance, EventDate, Printtime)

    if os.path.exists(File):
        os.remove(File)
    
    if CutoffDataInstance.livedata == "ON" or CutoffDataInstance.livedata == 1:
        file_clean.remove_files()
        
    output = [merged_df]

    if merged_asymptotic_df is not None:
        output.append(merged_asymptotic_df)
    else:
        output.append("Asymptotic directions were not generated")

    if merged_transmission_df is not None:
        output.append(merged_transmission_df)
    else:
        output.append("Transmission functions were not generated")

    output.append(README)

    return output