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
from ..inputs import planet_inputs
from ..fortran_calls import planet
from ..readme_generators import planet_readme
from ..data_classes.planet_data import PlanetData

def OTSO_planet(PlanetDataInstance: 'PlanetData') -> list:

    # ----------------------------------------------------------
    # Planet inputs
    # ----------------------------------------------------------

    planet_inputs.PlanetInputs(PlanetDataInstance)

    PlanetDataInstance.custom_coords_provided = (
        PlanetDataInstance.array_of_lats_and_longs is not None
    )

    # ----------------------------------------------------------
    # CPU information
    # ----------------------------------------------------------

    available_cpus = sorted(
        os.sched_getaffinity(0)
    )

    num_available = len(available_cpus)

    threads = PlanetDataInstance.threadnum

    if threads < 1:
        threads = 1

    # Maximum number of workers that can receive
    # non-overlapping CPU groups.
    max_workers = num_available // threads

    # ----------------------------------------------------------
    # Coordinates
    # ----------------------------------------------------------

    totalprocesses = len(
        PlanetDataInstance.coordinate_pairs
    )

    NewCoreNum = planet_inputs.CheckCoreNumPlanet(
        PlanetDataInstance.corenum
    )

    DataPlanet = []

    for point in PlanetDataInstance.coordinate_pairs:

        DataPlanet.append([
            f"{point[0]}_{point[1]}",
            point[0],
            point[1],
            PlanetDataInstance.startaltitude,
            PlanetDataInstance.zenith,
            PlanetDataInstance.azimuth
        ])

    random.shuffle(DataPlanet)

    # ----------------------------------------------------------
    # No coordinates
    # ----------------------------------------------------------

    if totalprocesses == 0:

        print(
            "\nWarning: No coordinate pairs provided or generated. "
            "Skipping planet computation."
        )

        EventDate = datetime(
            PlanetDataInstance.year,
            PlanetDataInstance.month,
            PlanetDataInstance.day,
            PlanetDataInstance.hour,
            PlanetDataInstance.minute,
            PlanetDataInstance.second
        )

        readme = planet_readme.READMEPlanet(
            PlanetDataInstance,
            EventDate,
            None,
            custom_coords_provided=(
                PlanetDataInstance.array_of_lats_and_longs
                is not None
            )
        )

        return [
            pd.DataFrame(),
            readme
        ]

    # ----------------------------------------------------------
    # Determine actual worker count
    # ----------------------------------------------------------

    actual_cores_to_use = min(
        NewCoreNum,
        totalprocesses,
        max_workers
    )

    if actual_cores_to_use < 1:

        raise RuntimeError(
            "No workers can be created. "
            f"Available CPUs={num_available}, "
            f"threads per worker={threads}"
        )

    # ----------------------------------------------------------
    # Worker JSON files
    # ----------------------------------------------------------

    worker_files = generate_worker_files(
        actual_cores_to_use
    )

    # ----------------------------------------------------------
    # Multiprocessing start method
    # ----------------------------------------------------------

    try:

        if not mp.get_start_method(
            allow_none=True
        ):

            mp.set_start_method(
                "spawn"
            )

    except RuntimeError:
        pass

    # ----------------------------------------------------------
    # Create queues
    # ----------------------------------------------------------

    # TaskQueue:
    #     Contains coordinates waiting to be processed.
    #
    # ProcessQueue:
    #     Contains progress messages from workers.

    TaskQueue = mp.Queue()
    ProcessQueue = mp.Queue()

    # ----------------------------------------------------------
    # Populate dynamic task queue
    # ----------------------------------------------------------

    #
    # IMPORTANT:
    #
    # We do NOT split DataPlanet into DataLists.
    #
    # Every location goes into one shared queue.
    #

    for point in DataPlanet:
        TaskQueue.put(point)

    # ----------------------------------------------------------
    # Add termination signals
    # ----------------------------------------------------------

    #
    # Each worker gets exactly one None.
    #
    # Once a worker receives None, it exits.
    #

    for _ in range(actual_cores_to_use):
        TaskQueue.put(None)

    # ----------------------------------------------------------
    # Create worker processes
    # ----------------------------------------------------------

    ChildProcesses = []

    for worker_id in range(actual_cores_to_use):

        # Assign a unique block of CPUs to this worker.
        #
        # Example:
        #
        # threads = 2
        #
        # worker 0 -> CPUs [0, 1]
        # worker 1 -> CPUs [2, 3]
        # worker 2 -> CPUs [4, 5]
        # worker 3 -> CPUs [6, 7]

        start_cpu = worker_id * threads

        cpus = set(
            available_cpus[
                start_cpu:start_cpu + threads
            ]
        )

        Child = mp.Process(
            target=planet.FortranPlanet,
            args=(
                TaskQueue,
                PlanetDataInstance,
                ProcessQueue,
                cpus,
                worker_files[worker_id]
            )
        )

        ChildProcesses.append(Child)

    # ----------------------------------------------------------
    # Start timing
    # ----------------------------------------------------------

    start = time.perf_counter()

    if PlanetDataInstance.Verbose:

        print(
            "OTSO Planet Computation Started"
        )

    # ----------------------------------------------------------
    # Start workers
    # ----------------------------------------------------------

    for process in ChildProcesses:
        process.start()

    # ----------------------------------------------------------
    # Progress bar
    # ----------------------------------------------------------

    processed = 0

    progress_bar = None

    if (
        PlanetDataInstance.Verbose
        and tqdm is not None
    ):

        progress_bar = tqdm(
            total=totalprocesses,
            desc="OTSO Running",
            unit=" location"
        )

    elif PlanetDataInstance.Verbose:

        print(
            f"Processing "
            f"{totalprocesses} grid points..."
        )

    # ----------------------------------------------------------
    # Collect progress
    # ----------------------------------------------------------

    while processed < totalprocesses:

        try:

            # A worker sends one message after completing
            # one coordinate.

            ProcessQueue.get(
                timeout=0.001
            )

            processed += 1

            if progress_bar is not None:

                progress_bar.update(1)

                progress_bar.set_description(
                    f"OTSO Running "
                    f"({processed}/{totalprocesses})"
                )

            elif PlanetDataInstance.Verbose:

                percent_complete = (
                    processed
                    / totalprocesses
                    * 100
                )

                sys.stdout.write(
                    f"\r{percent_complete:.2f}% complete "
                    f"({processed}/{totalprocesses} points)"
                )

                sys.stdout.flush()

        except queue.Empty:

            # Check whether all workers have died
            # unexpectedly.

            if all(
                not p.is_alive()
                for p in ChildProcesses
            ):

                # If all workers stopped but we haven't
                # processed everything, something failed.

                if processed < totalprocesses:

                    print(
                        "\nWarning: All worker processes "
                        "have stopped before all locations "
                        "were completed."
                    )

                break

    # ----------------------------------------------------------
    # Close progress bar
    # ----------------------------------------------------------

    if progress_bar is not None:
        progress_bar.close()

    # ----------------------------------------------------------
    # Wait for workers
    # ----------------------------------------------------------

    for process in ChildProcesses:
        process.join()

    # ----------------------------------------------------------
    # Check worker exit codes
    # ----------------------------------------------------------

    failed_workers = []

    for worker_id, process in enumerate(
        ChildProcesses
    ):

        if process.exitcode != 0:

            failed_workers.append(
                (
                    worker_id,
                    process.exitcode
                )
            )

    if failed_workers:

        print(
            "\nWarning: Worker processes failed:"
        )

        for worker_id, exitcode in failed_workers:

            print(
                f"  Worker {worker_id}: "
                f"exit code {exitcode}"
            )

    # ----------------------------------------------------------
    # Parse worker JSON files
    # ----------------------------------------------------------

    rigidity_df_list = []
    asymptotic_df_list = []
    transmission_df_list = []

    for file in worker_files:

        try:

            (
                rigidity_df,
                asymptotic_df,
                transmission_df
            ) = read_and_parse_jsons(
                file,
                PlanetDataInstance
            )

            rigidity_df_list.append(
                rigidity_df
            )

            if asymptotic_df is not None:

                asymptotic_df_list.append(
                    asymptotic_df
                )

            if transmission_df is not None:

                transmission_df_list.append(
                    transmission_df
                )

        finally:

            # Remove temporary worker file
            if os.path.exists(file):
                os.remove(file)

    # ----------------------------------------------------------
    # Combine rigidity
    # ----------------------------------------------------------

    if rigidity_df_list:

        final_planet = pd.concat(
            rigidity_df_list,
            ignore_index=True
        )

    else:

        final_planet = pd.DataFrame()

    final_planet = (
        final_planet
        .drop_duplicates(
            subset=[
                "Latitude",
                "Longitude"
            ]
        )
        .sort_values(
            [
                "Latitude",
                "Longitude"
            ]
        )
        .reset_index(
            drop=True
        )
    )

    # ----------------------------------------------------------
    # Combine asymptotic
    # ----------------------------------------------------------

    final_asymptotic = None

    if asymptotic_df_list:

        final_asymptotic = apply_polar_regions(
            pd.concat(
                asymptotic_df_list,
                ignore_index=True
            ),
            PlanetDataInstance.custom_coords_provided
        )

    # ----------------------------------------------------------
    # Combine transmission
    # ----------------------------------------------------------

    final_transmission = None

    if transmission_df_list:

        final_transmission = apply_polar_regions(
            pd.concat(
                transmission_df_list,
                ignore_index=True
            ),
            PlanetDataInstance.custom_coords_provided
        )

    # ----------------------------------------------------------
    # Apply polar regions to rigidity
    # ----------------------------------------------------------

    final_planet = apply_polar_regions(
        final_planet,
        PlanetDataInstance.custom_coords_provided
    )

    # ----------------------------------------------------------
    # Timing
    # ----------------------------------------------------------

    stop = time.perf_counter()

    Printtime = round(
        stop - start,
        3
    )

    if PlanetDataInstance.Verbose:

        print(
            "\nOTSO Planet Computation Complete"
        )

        print(
            "Whole Program Took: "
            + str(Printtime)
            + " seconds"
        )

    # ----------------------------------------------------------
    # README
    # ----------------------------------------------------------

    EventDate = datetime(
        PlanetDataInstance.year,
        PlanetDataInstance.month,
        PlanetDataInstance.day,
        PlanetDataInstance.hour,
        PlanetDataInstance.minute,
        PlanetDataInstance.second
    )

    readme = planet_readme.READMEPlanet(
        PlanetDataInstance,
        EventDate,
        Printtime
    )

    # ----------------------------------------------------------
    # Live data cleanup
    # ----------------------------------------------------------

    if (
        PlanetDataInstance.livedata == "ON"
        or PlanetDataInstance.livedata == 1
    ):

        file_clean.remove_files()

    # ----------------------------------------------------------
    # Return
    # ----------------------------------------------------------

    output = [
        final_planet
    ]

    if final_asymptotic is not None:

        output.append(
            final_asymptotic
        )

    else:

        output.append(
            "Asymptotic directions were not generated"
        )

    if final_transmission is not None:

        output.append(
            final_transmission
        )

    else:

        output.append(
            "Transmission functions were not generated"
        )

    output.append(
        readme
    )

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

            rigidity_results.append(
            {
                "Latitude": record["Latitude"],
                "Longitude": record["Longitude"],
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
                    "Latitude": record["Latitude"],
                    "Longitude": record["Longitude"],
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
                    "Latitude": lat,
                    "Longitude": lon,
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


def apply_polar_regions(
    dataframe,
    custom_coords_provided=False,
):
    """
    Restore polar regions after JSON parsing.

    For grid calculations:
        Duplicate pole rows across all longitudes.

    For custom coordinates:
        Keep only the explicitly supplied pole coordinates.
    """

    if dataframe is None or dataframe.empty:
        return dataframe

    dataframe = dataframe.copy()

    dataframe = dataframe.drop_duplicates(
        subset=["Latitude", "Longitude"]
    )

    dataframe["Latitude"] = pd.to_numeric(
        dataframe["Latitude"],
        errors="coerce",
    )

    dataframe["Longitude"] = pd.to_numeric(
        dataframe["Longitude"],
        errors="coerce",
    )

    dataframe = dataframe.dropna(
        subset=["Latitude", "Longitude"]
    )

    unique_longitudes = sorted(
        dataframe["Longitude"].unique()
    )

    north_pole = dataframe[
        dataframe["Latitude"] == 90.0
    ]

    south_pole = dataframe[
        dataframe["Latitude"] == -90.0
    ]

    # Remove poles temporarily
    dataframe = dataframe[
        ~dataframe["Latitude"].isin([90.0, -90.0])
    ]

    polar_rows = []

    if custom_coords_provided:

        # Keep only the user-supplied poles
        if not north_pole.empty:
            polar_rows.append(
                north_pole.iloc[0].copy()
            )

        if not south_pole.empty:
            polar_rows.append(
                south_pole.iloc[0].copy()
            )

    else:

        # Duplicate poles across every longitude
        for longitude in unique_longitudes:

            if not north_pole.empty:
                row = north_pole.iloc[0].copy()
                row["Longitude"] = longitude
                polar_rows.append(row)

            if not south_pole.empty:
                row = south_pole.iloc[0].copy()
                row["Longitude"] = longitude
                polar_rows.append(row)

    if polar_rows:

        dataframe = pd.concat(
            [
                dataframe,
                pd.DataFrame(polar_rows),
            ],
            ignore_index=True,
        )

    return (
        dataframe
        .drop_duplicates(subset=["Latitude", "Longitude"])
        .sort_values(
            by=["Latitude", "Longitude"],
            ascending=[True, True],
            ignore_index=True,
        )
    )