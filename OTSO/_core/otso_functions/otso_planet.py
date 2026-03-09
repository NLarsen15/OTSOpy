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
import csv
import tempfile
from tqdm import tqdm

from ..livedata import file_clean
from ..inputs import planet_inputs
from ..fortran_calls import planet
from ..readme_generators import planet_readme
from ..data_classes.planet_data import PlanetData

def OTSO_planet(PlanetDataInstance: 'PlanetData') -> list:

    gc.enable()

    planet_inputs.PlanetInputs(PlanetDataInstance)
    PlanetDataInstance.custom_coords_provided=(PlanetDataInstance.array_of_lats_and_longs is not None)
    
    ChildProcesses = []

    totalprocesses = len(PlanetDataInstance.coordinate_pairs)

    NewCoreNum = planet_inputs.CheckCoreNumPlanet(PlanetDataInstance.corenum)
    FileNamesPlanet = []

    for coord_pair in PlanetDataInstance.coordinate_pairs:
        FileNamesPlanet.append(f"{coord_pair[0]:.8f}_{coord_pair[1]:.8f}".replace(".", "dot"))
    
    DataPlanet = []
    i = 1
    for point, name in zip(PlanetDataInstance.coordinate_pairs, FileNamesPlanet):
        DataPlanet.append([name, point[0], point[1], PlanetDataInstance.startaltitude, PlanetDataInstance.zenith, PlanetDataInstance.azimuth]) 
        i = i + 1

    if totalprocesses == 0:
        print("\nWarning: No coordinate pairs provided or generated. Skipping planet computation.")
        EventDate = datetime(PlanetDataInstance.year,PlanetDataInstance.month,PlanetDataInstance.day,PlanetDataInstance.hour,PlanetDataInstance.minute,PlanetDataInstance.second)
        readme = planet_readme.READMEPlanet(PlanetDataInstance, EventDate, None,
                                        custom_coords_provided=(PlanetDataInstance.array_of_lats_and_longs is not None))
        return [pd.DataFrame(), readme]
        
    actual_cores_to_use = min(NewCoreNum, totalprocesses)

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, actual_cores_to_use)

    start = time.time()

    planet_list = [tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name for _ in range(PlanetDataInstance.corenum)]

    for planet_file in planet_list:
        with open(planet_file, mode='a', newline='', encoding='utf-8') as file:  # Open in write ('w') mode
            writer = csv.writer(file)
            
            if PlanetDataInstance.asymptotic == "YES":
                asymlevels_with_units = [f"{level} [{PlanetDataInstance.unit}]" for level in PlanetDataInstance.asymlevels]
                default_headers = ["Latitude", "Longitude", "Rc GV", "Rc Asym"]
                headers = default_headers + asymlevels_with_units
            else:
                headers = ["Latitude", "Longitude", "Ru", "Rc", "Rl"]
            writer.writerow(headers)

    processed = 0
    totalp = 0

    if PlanetDataInstance.Verbose:
        print("OTSO Planet Computation Started")

    if actual_cores_to_use == 1:       
        progress_bar = None
        if PlanetDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
        elif PlanetDataInstance.Verbose:
            print(f"Processing {totalprocesses} grid points...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        
        for Data, planetfile in zip(DataLists, planet_list):
            for single_item in Data:
                planet.FortranPlanet([single_item], PlanetDataInstance, planetfile, simple_queue)
                
                if simple_queue.items:
                    batch_count = sum(simple_queue.items)
                    totalp += batch_count
                    processed += batch_count
                    
                    if PlanetDataInstance.Verbose:
                        if progress_bar is not None:
                            progress_bar.update(batch_count)
                            progress_bar.set_description(f"OTSO Running ({totalp}/{totalprocesses})")
                        else:
                            percent_complete = (totalp / totalprocesses) * 100
                            sys.stdout.write(f"\r{percent_complete:.2f}% complete ({totalp}/{totalprocesses} points)")
                            sys.stdout.flush()
                    
                    simple_queue.items = []
            
            gc.collect()
        
        if progress_bar is not None:
            progress_bar.close()

    else:
        ProcessQueue = mp.Manager().Queue()

        try:
            if not mp.get_start_method(allow_none=True):
                mp.set_start_method('spawn')
        except RuntimeError:
            pass

        ChildProcesses = []
        for Data, planetfile in zip(DataLists, planet_list):
                Child = mp.Process(target=planet.FortranPlanet,  args=(Data, PlanetDataInstance, planetfile, ProcessQueue))
                ChildProcesses.append(Child)
            
        for a in ChildProcesses:
            a.start()

        progress_bar = None
        if PlanetDataInstance.Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
        elif PlanetDataInstance.Verbose:
            print(f"Processing {totalprocesses} grid points...")
     
        while processed < totalprocesses:
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
                if PlanetDataInstance.Verbose:
                    if progress_bar is not None:
                        progress_bar.update(sum(result_collector) if result_collector else 0)
                        progress_bar.set_description(f"OTSO Running ({totalp}/{totalprocesses})")
                    else:
                        percent_complete = (totalp / totalprocesses) * 100
                        sys.stdout.write(f"\r{percent_complete:.2f}% complete ({totalp}/{totalprocesses} points)")
                        sys.stdout.flush()

        
            except queue.Empty:
                pass
            
            time.sleep(0.5)

        if progress_bar is not None:
            progress_bar.close()

        for b in ChildProcesses:
            b.join()
            b.close()

    processed = 0
    
    ChildProcesses.clear()

    final_planet = combine_planet_files(planet_list)

    stop = time.time()
    Printtime = round((stop-start),3)

    if PlanetDataInstance.Verbose:
        print("\nOTSO Planet Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(PlanetDataInstance.year, PlanetDataInstance.month, PlanetDataInstance.day, PlanetDataInstance.hour, PlanetDataInstance.minute, PlanetDataInstance.second)
    
    if not DataPlanet:
         print("\nError: DataPlanet list is empty before final README generation.")
         readme_context = None 
    else:
        readme_context = DataPlanet[0]

    readme = planet_readme.READMEPlanet(PlanetDataInstance, EventDate, Printtime)

    if PlanetDataInstance.livedata == "ON" or PlanetDataInstance.livedata == 1:
        file_clean.remove_files()

    return [final_planet, readme]


def combine_planet_files(planet_list):
    resultsfinal = []
    for x in planet_list:
        df = pd.read_csv(x, header=0)
        for col in df.columns[:3]:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df = df.dropna(subset=df.columns[:3])
        df = df.drop_duplicates(subset=['Latitude', 'Longitude'])
        resultsfinal.append(df)
        os.remove(x)
    combined_planet = pd.concat(resultsfinal, ignore_index=True)
    combined_planet = combined_planet.drop_duplicates(subset=['Latitude', 'Longitude'])
    combined_planet['Longitude'] = pd.to_numeric(combined_planet['Longitude'], errors='coerce')
    combined_planet['Latitude'] = pd.to_numeric(combined_planet['Latitude'], errors='coerce')
    
    unique_longitudes = sorted(combined_planet['Longitude'].unique())
    polar_rows = []
    
    north_pole_data = combined_planet[combined_planet['Latitude'] == 90.0]
    south_pole_data = combined_planet[combined_planet['Latitude'] == -90.0]
    
    combined_planet = combined_planet[~((combined_planet['Latitude'] == 90.0) | (combined_planet['Latitude'] == -90.0))]
    
    if not north_pole_data.empty:
        for longitude in unique_longitudes:
            for _, row in north_pole_data.iterrows():
                new_row = row.copy()
                new_row['Longitude'] = longitude
                polar_rows.append(new_row)
    
    if not south_pole_data.empty:
        for longitude in unique_longitudes:
            for _, row in south_pole_data.iterrows():
                new_row = row.copy()
                new_row['Longitude'] = longitude
                polar_rows.append(new_row)
    
    if polar_rows:
        polar_df = pd.DataFrame(polar_rows)
        combined_planet = pd.concat([combined_planet, polar_df], ignore_index=True)
    
    final_planet = combined_planet.sort_values(by=["Latitude", "Longitude"], ignore_index=True)
    return final_planet

def get_unique_filename(filepath):
    base, ext = os.path.splitext(filepath)
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    return new_filepath