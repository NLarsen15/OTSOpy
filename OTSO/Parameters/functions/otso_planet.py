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
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,spheresize,inputcoord,Verbose,AdaptiveExternalModel,
           array_of_lats_and_longs=None,
           grid_params_user_set=False):

    gc.enable()

    Anum = 1
    PlanetInputArray = planet_inputs.PlanetInputs(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum, azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile,MHDcoordsys,inputcoord,AdaptiveExternalModel,
           array_of_lats_and_longs=array_of_lats_and_longs,
           grid_params_user_set=grid_params_user_set,)

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

    for coord_pair in combined_coordinates:
        FileNamesPlanet.append(f"{coord_pair[0]:.8f}_{coord_pair[1]:.8f}".replace(".", "dot"))
    
    DataPlanet = []
    i = 1
    for point, name in zip(combined_coordinates, FileNamesPlanet):
        Core = "Core " + str(i)
        DataPlanet.append([name, point[0], point[1], Alt, Zenith, Azimuth, Core]) 
        i = i + 1

    if totalprocesses == 0:
        print("\nWarning: No coordinate pairs provided or generated. Skipping planet computation.")
        EventDate = datetime(year,month,day,hour,minute,second)
        readme = readme_generators.READMEPlanet(None, RigidityArray, EventDate, Model, IntModel, 
                                        AntiCheck, IOPT, WindArray, Magnetopause, 0,
                                        maxlat,maxlong,minlat,minlong, latstep, longstep,
                                        MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                        LiveData, asymptotic, asymlevels, unit, serverdata, kp,
                                        custom_coords_provided=(array_of_lats_and_longs is not None))
        return [pd.DataFrame(), readme]
        
    actual_cores_to_use = min(NewCoreNum, totalprocesses)

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, actual_cores_to_use)

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

    results = []
    resultsfinal = []
    processed = 0
    totalp = 0

    if Verbose:
        print("OTSO Planet Computation Started")

    if actual_cores_to_use == 1:
        num_batches = len(DataLists)
        
        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
        elif Verbose:
            print(f"Processing {totalprocesses} grid points...")

        class SimpleQueue:
            def __init__(self):
                self.items = []
            def put(self, item):
                self.items.append(item)
            def get_all(self):
                return self.items

        simple_queue = SimpleQueue()
        
        for Data, Core, planetfile in zip(DataLists, CoreList, planet_list):
            for single_item in Data:
                fortran_calls.fortrancallPlanet([single_item], RigidityArray, DateArray, Model, IntModel, 
                                              ParticleArray, IOPT, WindArray, 
                                              Magnetopause, MaxStepPercent, EndParams, 
                                              Rcomp, Rscan, asymptotic, asymlevels, unit,
                                              simple_queue, g, h, planetfile, MHDfile, MHDcoordsys,
                                              spheresize, inputcoord)
                
                if simple_queue.items:
                    batch_count = sum(simple_queue.items)
                    totalp += batch_count
                    processed += batch_count
                    
                    if Verbose:
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

        progress_bar = None
        if Verbose and tqdm is not None:
            progress_bar = tqdm(total=totalprocesses, desc="OTSO Running", unit=" location")
        elif Verbose:
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
                if Verbose:
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

    planet = combine_planet_files(planet_list)

    stop = time.time()
    Printtime = round((stop-start),3)

    if Verbose:
        print("\nOTSO Planet Computation Complete")
        print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    
    if not DataPlanet:
         print("\nError: DataPlanet list is empty before final README generation.")
         readme_context = None 
    else:
        readme_context = DataPlanet[0]

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
    
    planet = combined_planet.sort_values(by=["Latitude", "Longitude"], ignore_index=True)
    return planet

def get_unique_filename(filepath):
    base, ext = os.path.splitext(filepath)
    counter = 1
    new_filepath = filepath
    while os.path.exists(new_filepath):
        new_filepath = f"{base}_{counter}{ext}"
        counter += 1
    return new_filepath