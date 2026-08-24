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
from ..inputs import transmission_inputs
from ..fortran_calls import transmission
from ..readme_generators import transmission_readme
from ..data_classes.transmission_data import TransmissionData


def OTSO_transmission(TransmissionDataInstance: TransmissionData) -> list:

    File = tempfile.NamedTemporaryFile(delete=False, suffix=".json").name
    lock = mp.Lock()

    transmission_inputs.TransmissionInputs(TransmissionDataInstance)

    ChildProcesses = []

    next_cpu = 0

    UsedCores = cores.Cores(TransmissionDataInstance.station_array, TransmissionDataInstance.corenum)
    Positionlists = UsedCores.getPositions()

    start = time.time()

    if TransmissionDataInstance.Verbose:
        print("OTSO Transmission Computation Started")

    total_stations = len(TransmissionDataInstance.station_array)
    results = []


    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass
    
    ProcessQueue = mp.Manager().Queue()
    for Data in Positionlists:
        threads = TransmissionDataInstance.threadnum
        cpus = set(range(next_cpu, next_cpu + threads))
        next_cpu += threads
        Child = mp.Process(target=transmission.FortranTransmission,  args=(Data, TransmissionDataInstance, ProcessQueue, cpus, File, lock))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    processed = 0

    progress_bar = None
    if TransmissionDataInstance.Verbose and tqdm is not None:
        progress_bar = tqdm(total=total_stations, desc="OTSO Running", unit=" cone")
    elif TransmissionDataInstance.Verbose:
        print(f"Processing {total_stations} stations...")

    while processed < total_stations:
        try:
            ProcessQueue.get(timeout=0.001)
            processed += 1

            if TransmissionDataInstance.Verbose:
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

    Transmissionlist = []

    with open(File, "r", encoding="utf-8") as f:
        for line in f:
            record = json.loads(line)

            transmissionrecord = record["transmission"]

            df = pd.DataFrame(
                transmissionrecord["data"],
                index=transmissionrecord["index"],
                columns=transmissionrecord["columns"],
            )

            # Rename the TF column if needed (already contains station name)
            Transmissionlist.append(df)


    # Combine all station transmissions
    merged_Transmission_df = Transmissionlist[0]

    for df in Transmissionlist[1:]:
        merged_Transmission_df = pd.merge(
            merged_Transmission_df,
            df,
            on="R [GV]"
        )


    # Put rigidity first
    cols = ["R [GV]"] + [
        col for col in merged_Transmission_df.columns 
        if col != "R [GV]"
    ]

    merged_Transmission_df = merged_Transmission_df[cols]

    # Sort station columns alphabetically
    merged_Transmission_df = (
        merged_Transmission_df[
            ["R [GV]"] +
            sorted(merged_Transmission_df.columns.drop("R [GV]"))
        ]
    )

    stop = time.time()
    Printtime = round((stop-start),3)

    if TransmissionDataInstance.Verbose:
        print("\nOTSO Transmission Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(TransmissionDataInstance.year,TransmissionDataInstance.month,TransmissionDataInstance.day,
                         TransmissionDataInstance.hour,TransmissionDataInstance.minute,TransmissionDataInstance.second)
    README = transmission_readme.READMETransmission(TransmissionDataInstance, EventDate, Printtime)

    if TransmissionDataInstance.livedata == "ON" or TransmissionDataInstance.livedata == 1:
        file_clean.remove_files()

    if os.path.exists(File):
        os.remove(File)
    
    return [merged_Transmission_df, README]