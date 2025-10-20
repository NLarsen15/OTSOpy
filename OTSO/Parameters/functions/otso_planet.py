import time
from datetime import datetime
import multiprocessing as mp
import os
from . import fortran_calls, readme_generators,cores, misc, planet_inputs
import pandas as pd
import sys
import queue
import random
import numpy as np
import gc
import csv
from collections import defaultdict
import tempfile
from tqdm import tqdm

def OTSO_planet(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum, azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose,
           array_of_lats_and_longs=None,
           grid_params_user_set=False,):

    gc.enable()

    Anum = 1
    PlanetInputArray = planet_inputs.PlanetInputs(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum, azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,inputcoord,
           array_of_lats_and_longs=array_of_lats_and_longs,
           grid_params_user_set=grid_params_user_set)

    combined_coordinates = PlanetInputArray[0]
    RigidityArray = PlanetInputArray[1]
    DateArray = PlanetInputArray[2]
    Model = PlanetInputArray[3]
    IntModel = PlanetInputArray[4]
    ParticleArray = PlanetInputArray[5]
    IOPT = PlanetInputArray[6]
    WindArray = PlanetInputArray[7]
    Magnetopause = PlanetInputArray[8]
    MaxStepPercent = PlanetInputArray[9]/100
    EndParams = PlanetInputArray[10]
    Rcomp = PlanetInputArray[11]
    Rscan = PlanetInputArray[12]
    Zenith = PlanetInputArray[13]
    Azimuth = PlanetInputArray[14]
    CoreNum = PlanetInputArray[15]
    asymptotic = PlanetInputArray[16]
    asymlevels = PlanetInputArray[17]
    Alt = PlanetInputArray[18]
    LiveData = PlanetInputArray[19]
    AntiCheck = PlanetInputArray[20]
    g = PlanetInputArray[21]
    h = PlanetInputArray[22]

    kp = WindArray[17]

    del(PlanetInputArray)

    ChildProcesses = []

    totalprocesses = len(combined_coordinates)

    NewCoreNum = misc.CheckCoreNumPlanet(CoreNum)
    FileNamesPlanet = []

    # Generate filenames based on the coordinate pairs
    for coord_pair in combined_coordinates:
        # Ensure formatting is consistent, handle potential floats
        FileNamesPlanet.append(f"{coord_pair[0]:.8f}_{coord_pair[1]:.8f}".replace(".", "dot")) # Using more precision and replacing dots
    
    DataPlanet = []
    i = 1
    # Populate DataPlanet directly from combined_coordinates
    for point, name in zip(combined_coordinates, FileNamesPlanet):
        Core = "Core " + str(i)
        # point[0] is latitude, point[1] is longitude
        DataPlanet.append([name, point[0], point[1], Alt, Zenith, Azimuth, Core]) 
        i = i + 1

    # --- Fix Start: Cap cores and handle empty input --- 
    if totalprocesses == 0:
        print("\nWarning: No coordinate pairs provided or generated. Skipping planet computation.")
        # Need to decide what to return. An empty DataFrame and a basic README?
        EventDate = datetime(year,month,day,hour,minute,second)
        readme = readme_generators.READMEPlanet(None, RigidityArray, EventDate, Model, IntModel, 
                                        AntiCheck, IOPT, WindArray, Magnetopause, 0,
                                        maxlat,maxlong,minlat,minlong, latstep, longstep,
                                        MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                        LiveData, asymptotic, asymlevels, unit, serverdata, kp,
                                        custom_coords_provided=(array_of_lats_and_longs is not None)) # Added flag back temporarily
        return [pd.DataFrame(), readme] # Return empty dataframe and readme
        
    # Ensure the number of processes doesn't exceed the number of points
    actual_cores_to_use = min(NewCoreNum, totalprocesses)
    # --- Fix End --- 

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    # Use actual_cores_to_use for splitting
    DataLists = np.array_split(shuffled_list, actual_cores_to_use) 

    # Adjust core list generation based on actual_cores_to_use
    CoreList = np.arange(1, actual_cores_to_use + 1) 
    start = time.time()

    current_dir = os.path.dirname(os.path.realpath(__file__))

    planet_list = [tempfile.NamedTemporaryFile(delete=False, suffix=".csv").name for _ in range(corenum)]

    for planet_file in planet_list:
        with open(planet_file, mode='a', newline='', encoding='utf-8') as file:  # Open in write ('w') mode
            writer = csv.writer(file)
            
            if asymptotic == "YES":
                asymlevels_with_units = [f"{level} [{unit}]" for level in asymlevels]
                default_headers = ["Latitude", "Longitude", "Rc GV", "Rc Asym"]
                headers = default_headers + asymlevels_with_units
            else:
                headers = ["Latitude", "Longitude", "Ru", "Rc", "Rl"]
            writer.writerow(headers)

    ProcessQueue = mp.Manager().Queue()

    results = []
    resultsfinal = []
    processed = 0
    totalp = 0

    if Verbose:
        print("OTSO Planet Computation Started")

    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:
        pass

    ChildProcesses = []
    for Data, Core, planetfile in zip(DataLists, CoreList, planet_list):
            Child = mp.Process(target=fortran_calls.fortrancallPlanet,  args=(Data, RigidityArray, DateArray, Model, IntModel, 
                                                                              ParticleArray, IOPT, WindArray, 
                                                                              Magnetopause, MaxStepPercent, EndParams, 
                                                                              Rcomp, Rscan, asymptotic, asymlevels, unit,
                                                                              ProcessQueue,g,h,planetfile, MHDfile, MHDcoordsys,
                                                                              spheresize,inputcoord))
            ChildProcesses.append(Child)
        
    for a in ChildProcesses:
        a.start()

    # Initialize progress bar if tqdm is available and Verbose is True
    progress_bar = None
    if Verbose and tqdm is not None:
        progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
    elif Verbose:
        # Fallback to simple counter if tqdm is not available
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
    
            # Update totalp with the sum of items processed by cores
            if result_collector:
                totalp = totalp + sum(result_collector)
            
            gc.collect()
            # Update progress
            if Verbose:
                if progress_bar is not None:
                    # Update progress bar with the actual number of items processed
                    progress_bar.update(sum(result_collector) if result_collector else 0)
                    # Update the description to show current progress
                    progress_bar.set_description(f"OTSO Running ({totalp}/{totalprocesses})")
                else:
                    # Fallback to percentage if tqdm is not available
                    percent_complete = (totalp / totalprocesses) * 100
                    sys.stdout.write(f"\r{percent_complete:.2f}% complete ({totalp}/{totalprocesses} points)")
                    sys.stdout.flush()

    
        except queue.Empty:
            pass
        
        time.sleep(0.5)

    # Close progress bar if it was created
    if progress_bar is not None:
        progress_bar.close()

    for b in ChildProcesses:
        b.join()
        b.close()

    processed = 0
    
    ChildProcesses.clear()

    planet = combine_planet_files(planet_list)

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
        print("\nOTSO Planet Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    
    # --- Fix: Use the first element of DataPlanet for README context --- 
    if not DataPlanet:
         # This case should ideally be caught earlier, but handle defensively
         print("\nError: DataPlanet list is empty before final README generation.")
         readme_context = None 
    else:
        readme_context = DataPlanet[0] # Use the first point's data
    # --- End Fix --- 

    # Pass the representative context to READMEPlanet
    readme = readme_generators.READMEPlanet(readme_context, RigidityArray, EventDate, Model, IntModel, 
                                            AntiCheck, IOPT, WindArray, Magnetopause, Printtime,
                                            maxlat,maxlong,minlat,minlong, latstep, longstep,
                                            MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                            LiveData, asymptotic, asymlevels, unit, serverdata, kp,
                                            custom_coords_provided=(array_of_lats_and_longs is not None))

    if livedata == "ON" or livedata == 1:
        misc.remove_files()

    return [planet, readme]


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
    
    # Duplicate polar points for all longitude values
    unique_longitudes = sorted(combined_planet['Longitude'].unique())
    polar_rows = []
    
    # Find existing polar points (90° and -90°)
    north_pole_data = combined_planet[combined_planet['Latitude'] == 90.0]
    south_pole_data = combined_planet[combined_planet['Latitude'] == -90.0]
    
    # Remove existing polar points from the dataframe
    combined_planet = combined_planet[~((combined_planet['Latitude'] == 90.0) | (combined_planet['Latitude'] == -90.0))]
    
    # Duplicate north pole data for each longitude
    if not north_pole_data.empty:
        for longitude in unique_longitudes:
            for _, row in north_pole_data.iterrows():
                new_row = row.copy()
                new_row['Longitude'] = longitude
                polar_rows.append(new_row)
    
    # Duplicate south pole data for each longitude
    if not south_pole_data.empty:
        for longitude in unique_longitudes:
            for _, row in south_pole_data.iterrows():
                new_row = row.copy()
                new_row['Longitude'] = longitude
                polar_rows.append(new_row)
    
    # Add the duplicated polar rows back to the dataframe
    if polar_rows:
        polar_df = pd.DataFrame(polar_rows)
        combined_planet = pd.concat([combined_planet, polar_df], ignore_index=True)
    
    planet = combined_planet.sort_values(by=["Latitude", "Longitude"], ignore_index=True)
    return planet

def get_unique_filename(filepath):
    """If the file exists, add a numerical suffix to make it unique."""
    base, ext = os.path.splitext(filepath)
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    return new_filepath