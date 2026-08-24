import time
import multiprocessing as mp
import os
import pandas as pd
import sys
import queue
import numpy as np
import tempfile
from tqdm import tqdm
import json

from ..custom_classes import cores
from ..livedata import file_clean
from ..inputs import flight_inputs
from ..fortran_calls import flight
from ..readme_generators import flight_readme
from ..data_classes.flight_data import FlightData


def OTSO_flight(FlightDataInstance: FlightData) -> list:

    flight_inputs.FlightInputs(FlightDataInstance)

    available_cpus = sorted(os.sched_getaffinity(0))
    num_available = len(available_cpus)

    next_cpu = 0

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


    worker_files = generate_worker_files(FlightDataInstance.corenum)

    start = time.time()
    if FlightDataInstance.Verbose:
        print("OTSO Flight Computation Started")

    resultsfinal = []
    processed = 0
    totalp = 0
    total_stations = len(FlightDataInstance.station_array)

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

    ProcessQueue = mp.Queue()
    for Data, Date, I, Wind, worker_file, G, H in zip(Positionlists,DateArrayLists, IOPTLists, WindLists, worker_files, GArrayLists, HArrayLists):
        threads = FlightDataInstance.threadnum

        cpus = {
                available_cpus[(next_cpu + i) % num_available]
                for i in range(threads)
                }

        next_cpu = (next_cpu + threads) % num_available

        Child = mp.Process(target=flight.FortranFlight,  args=(Data, Date, I, Wind, G, H, worker_file, FlightDataInstance, ProcessQueue, cpus))
        ChildProcesses.append(Child)

    for a in ChildProcesses:
        a.start()

    while processed < total_stations:
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

    rigidity_df_list = []
    asymptotic_df_list = []
    transmission_df_list = []

    
    for file in worker_files:

        rigidity_df, asymptotic_df, transmission_df = read_and_parse_jsons(file, FlightDataInstance)

        rigidity_df_list.append(rigidity_df)

        if asymptotic_df is not None:
            asymptotic_df_list.append(asymptotic_df)

        if transmission_df is not None:
            transmission_df_list.append(transmission_df)

        os.remove(file)

        # Combine rigidity
    final_flight = pd.concat(
        rigidity_df_list,
        ignore_index=True
    )

    final_flight = (
        final_flight
        .sort_values(["Date"])
        .reset_index(drop=True)
    )


    # Combine asymptotic
    final_asymptotic = None

    if asymptotic_df_list:
        final_asymptotic = (
        pd.concat(asymptotic_df_list, ignore_index=True)
        .sort_values("Date")
        .reset_index(drop=True)
    )


    # Combine transmission
    final_transmission = None

    if transmission_df_list:
        final_transmission = (
        pd.concat(transmission_df_list, ignore_index=True)
        .sort_values("Date")
        .reset_index(drop=True)
    )

    stop = time.time()
    Printtime = round((stop-start),3)

    if FlightDataInstance.Verbose:
        print("\nOTSO Flight Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")


    readme = flight_readme.READMEFlight(FlightDataInstance, Printtime)

    datareadme = flight_readme.READMEFlightData(FlightDataInstance.datearraylist,FlightDataInstance.windarraylist,FlightDataInstance.Kplist)
    
    if FlightDataInstance.livedata == "ON" or FlightDataInstance.livedata == 1:
        file_clean.remove_files()

    output = [final_flight]

    if final_asymptotic is not None:
        output.append(final_asymptotic)
    else:
        output.append("Asymptotic directions were not generated")

    if final_transmission is not None:
        output.append(final_transmission)
    else:
        output.append("Transmission functions were not generated")

    output.append(readme)
    output.append(datareadme)

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

def read_and_parse_jsons(File, Data):
    rigidity_results = []
    asymptotic_results = []
    transmission_results = []

    with open(File, "r", encoding="utf-8") as f:
        rigidity_results = []
        asymptotic_results = []
        transmission_results = []
        for line in f:

            record = json.loads(line)

            lat = record["Latitude"]
            lon = record["Longitude"]

            # --------------------------------------------------
            # Rigidity dataframe
            # --------------------------------------------------

            if Data.inputcoord == "GDZ":
                altunit = "km"
            else:
                altunit = "Re"


            rigidity_results.append(
            {
                "Date": record["Date"],
                "Latitude": record["Latitude"],
                "Longitude": record["Longitude"],
                f"Altitude [{altunit}]": record[f"Altitude [{altunit}]"],
                "Ru [GV]": record["Ru [GV]"],
                "Rc [GV]": record["Rc [GV]"],
                "Rl [GV]": record["Rl [GV]"],
                "PTF": record["PTF"],
            }
            )


            # --------------------------------------------------
            # Asymptotic dataframe
            # --------------------------------------------------

            if "asymptotic" in record:

                asym_record = record["asymptotic"][0]

                row = {
                    "Date": record["Date"],
                    "Latitude": record["Latitude"],
                    "Longitude": record["Longitude"],
                    f"Altitude [{altunit}]": record[f"Altitude [{altunit}]"]
                }

                for level, value in asym_record.items():

                    if level == "Station":
                        continue

                    row[level] = value

                asymptotic_results.append(row)


            # --------------------------------------------------
            # Transmission dataframe
            # --------------------------------------------------

            if "transmission" in record:

                row = {
                    "Date": record["Date"],
                    "Latitude": lat,
                    "Longitude": lon,
                    f"Altitude [{altunit}]": record[f"Altitude [{altunit}]"]
                }

                row.update(record["transmission"])

                transmission_results.append(row)
    
    
        # --------------------------------------------------
        # Final dataframe construction
        # --------------------------------------------------
    
        rigidity_df = pd.DataFrame(rigidity_results)
    
    
        asymptotic_df = None
        if asymptotic_results:
            asymptotic_df = pd.DataFrame(asymptotic_results)
    
    
        transmission_df = None
        if transmission_results:
            transmission_df = pd.DataFrame(transmission_results)

    return rigidity_df, asymptotic_df, transmission_df