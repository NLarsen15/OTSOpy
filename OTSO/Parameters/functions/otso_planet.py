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

def OTSO_planet(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum, azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h):
    gc.enable()

    Anum = 1
    PlanetInputArray = planet_inputs.PlanetInputs(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum, azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h)

    LongitudeList = PlanetInputArray[0]
    LatitudeList = PlanetInputArray[1]
    RigidityArray = PlanetInputArray[2]
    DateArray = PlanetInputArray[3]
    Model = PlanetInputArray[4]
    IntModel = PlanetInputArray[5]
    ParticleArray = PlanetInputArray[6]
    IOPT = PlanetInputArray[7]
    WindArray = PlanetInputArray[8]
    Magnetopause = PlanetInputArray[9]
    MaxStepPercent = PlanetInputArray[10]/100
    EndParams = PlanetInputArray[11]
    Rcomp = PlanetInputArray[12]
    Rscan = PlanetInputArray[13]
    Zenith = PlanetInputArray[14]
    Azimuth = PlanetInputArray[15]
    CoreNum = PlanetInputArray[16]
    asymptotic = PlanetInputArray[17]
    asymlevels = PlanetInputArray[18]
    Alt = PlanetInputArray[19]
    LiveData = PlanetInputArray[20]
    AntiCheck = PlanetInputArray[21]
    g = PlanetInputArray[22]
    h = PlanetInputArray[23]

    del(PlanetInputArray)

    ChildProcesses = []

    totalprocesses = len(LongitudeList)*len(LatitudeList)

    NewCoreNum = misc.CheckCoreNumPlanet(CoreNum)
    FileNamesPlanet = []

    combined_coordinates = [(round(lat, 6), round(lon, 6)) for lat in LatitudeList for lon in LongitudeList]
    
    lat_grid, lon_grid = np.meshgrid(LatitudeList, LongitudeList, indexing='ij')  

    combined_coordinates = list(zip(lat_grid.ravel(), lon_grid.ravel()))
    combined_coordinates = list(set((round(lat, 6), round(lon, 6)) for lat, lon in combined_coordinates))

    coord_dict = defaultdict(list)

    for coordlist in combined_coordinates:
        FileNamesPlanet.append(str(coordlist[0]) + "_" + str(coordlist[1]))
    DataPlanet = []
    i = 1
    for point,name in zip(combined_coordinates, FileNamesPlanet):
        Core = "Core " + str(i)
        DataPlanet.append([name,point[0],point[1],Alt,Zenith,Azimuth,Core])
        i = i + 1

    shuffled_list = DataPlanet.copy()
    random.shuffle(shuffled_list)
    DataLists = np.array_split(shuffled_list, CoreNum)

    CoreList = np.arange(1, CoreNum + 1)
    start = time.time()

    current_dir = os.path.dirname(os.path.realpath(__file__))

    planet_list = [os.path.join(current_dir, f"Planet{i}.csv") for i in range(1, (CoreNum) + 1)]
    for planet_file in planet_list:
        with open(planet_file, mode='a', newline='', encoding='utf-8') as file:  # Open in write ('w') mode
            writer = csv.writer(file)
            
            if asymptotic == "YES":
                asymlevels_with_units = [f"{level} [{unit}]" for level in asymlevels]
                default_headers = ["Latitude", "Longitude", "Rc GV", "Rc Asym"]
                headers = default_headers + asymlevels_with_units
            else:
                headers = ["Latitude", "Longitude", "Rl", "Rc", "Ru"]
            writer.writerow(headers)
    
    
    ProcessQueue = mp.Manager().Queue()

    results = []
    resultsfinal = []
    processed = 0
    totalp = 0

    print("OTSO Planet Computation Started")
    sys.stdout.write(f"\r{0:.2f}% complete")
    try:
        if not mp.get_start_method(allow_none=True):
            mp.set_start_method('spawn')
    except RuntimeError:

        pass

    i = 0
    
    ChildProcesses = []
    for Data, Core, planetfile in zip(DataLists, CoreList, planet_list):
            Child = mp.Process(target=fortran_calls.fortrancallPlanet,  args=(Data, RigidityArray, DateArray, Model, IntModel, 
                                                                              ParticleArray, IOPT, WindArray, 
                                                                              Magnetopause, MaxStepPercent, EndParams, 
                                                                              Rcomp, Rscan, asymptotic, asymlevels, unit,
                                                                              ProcessQueue,g,h,planetfile))
            ChildProcesses.append(Child)
        
    for a in ChildProcesses:
        a.start()
 
    while processed < totalprocesses:
        try:
            result_collector = []
            while True:
                try:
                    countint = ProcessQueue.get(timeout=0.001)
                    result_collector.append(countint)
                    processed += 1
                    totalp = totalp + sum(result_collector)
                    result_collector = []
                except queue.Empty:
                    break
    
            
            gc.collect()
            percent_complete = (totalp / totalprocesses) * 100
            sys.stdout.write(f"\r{percent_complete:.2f}% complete")
            sys.stdout.flush()

    
        except queue.Empty:
            pass
        
        time.sleep(2)

    for b in ChildProcesses:
        b.join()
        b.close()

    processed = 0
    
    ChildProcesses.clear()

    planet = combine_planet_files(planet_list)
   
    print("\nOTSO Planet Computation Complete")
    stop = time.time()
    Printtime = round((stop-start),3)
    print("Whole Program Took: " + str(Printtime) + " seconds")
    
    EventDate = datetime(year,month,day,hour,minute,second)
    readme = readme_generators.READMEPlanet(Data, RigidityArray, EventDate, Model, IntModel, 
                                            AntiCheck, IOPT, WindArray, Magnetopause, Printtime,
                                            maxlat,maxlong,minlat,minlong, latstep, longstep,
                                            MaxStepPercent*100, EndParams, cutoff_comp, Rscan, 
                                            LiveData, asymptotic, asymlevels, unit, serverdata, kp)

    if LiveData == 1:
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
    planet = combined_planet.sort_values(by=["Latitude", "Longitude"], ignore_index=True)
    return planet