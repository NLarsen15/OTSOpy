from . import MiddleMan as OTSOLib
import time
import os
import shutil
from datetime import datetime
import pandas as pd
import glob
import csv
import numpy as np
import multiprocessing as mp
from . import date
import gc
import sys
import psutil
import tracemalloc
from . import MHDinit

process = psutil.Process()

def fortrancallCutoff(Data, Core, RigidityArray, DateArray, model, IntModel, ParticleArray, IOPT, 
WindArray, Magnetopause, CoordinateSystem, MaxStepPercent, 
EndParams, Rcomp, Rscan, Kp, queue, g, h, MHDfile, MHDCoordSys):
    
    if model[1] == 99:
      MHDinit.MHDinitialise(MHDfile)

    for x in Data:
      
      newstart = time.time()
      Position = [x[3],x[1],x[2],x[4],x[5]]


      Station = x[0]

      #print(FileDescriptors)

      NMname = Station

      StartRigidity = RigidityArray[0]
      EndRigidity = RigidityArray[1]
      RigidityStep = RigidityArray[2]
      AtomicNum = ParticleArray[0]
      AntiCheck = ParticleArray[1]

      FileName = NMname + ".csv"
      Rigidities = OTSOLib.cutoff(Position, StartRigidity, EndRigidity, RigidityStep, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause,
                                   FileName, CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, g, h, MHDCoordSys)

      
      Rigiditydataframe = pd.DataFrame({Station: Rigidities}, index=['Ru', 'Rc', 'Rl'])
      queue.put(Rigiditydataframe)

    return

def fortrancallCone(Data, Core, RigidityArray, DateArray, model, IntModel, ParticleArray, IOPT, WindArray, Magnetopause,
                     CoordinateSystem, MaxStepPercent, EndParams, queue, g, h, MHDfile, MHDCoordSys):
    
    if model[1] == 99:
      MHDinit.MHDinitialise(MHDfile)
    
    for x in Data:
      
      Position = [x[3],x[1],x[2],x[4],x[5]]
      Station = x[0]

      NMname = Station

      StartRigidity = RigidityArray[0]
      EndRigidity = RigidityArray[1]
      RigidityStep = RigidityArray[2]
      AtomicNum = ParticleArray[0]
      AntiCheck = ParticleArray[1]

      length = int(StartRigidity/RigidityStep)

      FileName = NMname + ".csv"
      Cone, Rigidities = OTSOLib.cone(Position, StartRigidity, EndRigidity, RigidityStep, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray,
                                       Magnetopause, FileName, CoordinateSystem, MaxStepPercent, EndParams, length, g, h, MHDCoordSys)
      decoded_lines = [line.decode('utf-8').strip() for sublist in Cone for line in sublist]

      data = [line.split() for line in decoded_lines]

      Conedf = pd.DataFrame(data, columns=['R [GV]', 'Filter', 'ALat', 'ALong'])
      Conedf[NMname] = Conedf[['Filter', 'ALat', 'ALong']].astype(str).agg(';'.join, axis=1)
      Conedf.drop(columns=['Filter', 'ALat', 'ALong'], inplace=True)

      Rigiditydataframe = pd.DataFrame({Station: Rigidities}, index=['Ru', 'Rc', 'Rl'])
 
      queue.put([Conedf, Rigiditydataframe])
    
    return


def fortrancallPlanet(Data, Rigidity, DateArray, model, IntModel, ParticleArray, IOPT, WindArray, Magnetopause, MaxStepPercent, 
                      EndParams, Rcomp, Rscan, asymptotic, asymlevels, unit, queue, g, h, PlanetFile,MHDfile, MHDCoordSys):
  with open(PlanetFile, mode='a', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)

    if model[1] == 99:
      MHDinit.MHDinitialise(MHDfile)

    for x in Data:
        
        Position = [x[3],x[1],x[2],x[4],x[5]]
        Station = x[0]

        AtomicNum = ParticleArray[0]
        AntiCheck = ParticleArray[1]
  
        NMname = Station
        Rigidities = [0,0,0]
  
        FileName = NMname + ".csv"
        Rigidities = OTSOLib.planet(Position, Rigidity, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, FileName,
                                     MaxStepPercent, EndParams, Rcomp, Rscan, g, h, MHDCoordSys)
        CoordinateSystem = "GEO"
  
        lat_long_pairs = []
        P_List = []
  
        if asymptotic == "YES":
          Energy_List = asymlevels.copy()
          if unit == "GeV":   
            E_0 = 0.938
            for i in Energy_List:
              R = (i**2 + 2*i*E_0)**(0.5)
              P_List.append(R)
            P_List.insert(0, Rigidities[1])
          else:
            if unit == "GV": 
              P_List = Energy_List.copy()
            P_List.insert(0, Rigidities[1])
          for P in P_List:
              if P == P_List[0]:  # Check if it's the first value in P_List
                  while True:
                      bool, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, FileName,
                                                            CoordinateSystem, MaxStepPercent, EndParams, g, h, MHDCoordSys)
                      P += Rigidity[2]
                      if bool == 1:
                          lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
                          break 
              else:
                  bool, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, FileName,
                                                        CoordinateSystem, MaxStepPercent, EndParams, g, h, MHDCoordSys)
                  lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
  
          formatted_list = [f"{bool};{lat};{long}" for bool, lat, long in lat_long_pairs]
  
          formatted_list.insert(0, Position[1])
          formatted_list.insert(1, Position[2])
          formatted_list.insert(2, Rigidities[1])
          writer.writerow(formatted_list)
          queue.put(1)
  
        else:
           data = [x[1],x[2],Rigidities[0],Rigidities[1], Rigidities[2]]
           writer.writerow(data)
           queue.put(1)
           
    file.close()

  return

  

def fortrancallTrajectory(Data, Core, Rigidity, DateArray, model, IntModel, ParticleArray, IOPT, WindArray, 
                          Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, queue, g, h, MHDfile, MHDCoordSys):
  
  if model[1] == 99:
    MHDinit.MHDinitialise(MHDfile)
  
  for x in Data:
    
    Position = [x[3],x[1],x[2],x[4],x[5]]
    Station = x[0]

    AtomicNum = ParticleArray[0]
    AntiCheck = ParticleArray[1]

    NMname = Station

    FileName = NMname + ".csv"
    OTSOLib.trajectory_full(Position, Rigidity, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, FileName, 
                            CoordinateSystem, MaxStepPercent, EndParams, g, h, MHDCoordSys)
    Trajectory = pd.read_csv(FileName)
    Trajectory.columns = [f"{col}_Re [{CoordinateSystem}]" for col in Trajectory.columns]
    
    dataframe_dict = {NMname: Trajectory}

    queue.put(dataframe_dict)
    os.remove(FileName)

  return

def fortrancallFlight(Data, Rigidity, DateArray, model, IntModel, ParticleArray, IOPT, WindArray, Magnetopause, 
                      MaxStepPercent, EndParams, Rcomp, Rscan, asymptotic, asymlevels, unit, queue, g, h, 
                      CoordinateSystem, FlightFile, MHDfile, MHDCoordSys):
  with open(FlightFile, mode='a', newline='', encoding='utf-8') as file:
    headers = ["Date","Latitude", "Longitude","Altitude", "Rl", "Rc", "Ru"]
    writer = csv.writer(file)
    writer.writerow(headers)

    if model[1] == 99:
      MHDinit.MHDinitialise(MHDfile)


    for x,y,z,I in zip(Data,DateArray,WindArray,IOPT):
        
        Position = [x[3],x[1],x[2],x[4],x[5]]
        Station = x[0]
        
        datetimeobj = date.convert_to_datetime(y)
        Wind = z
  
        StartRigidity = Rigidity[0]
        EndRigidity = Rigidity[1]
        RigidityStep = Rigidity[2]
        AtomicNum = ParticleArray[0]
        AntiCheck = ParticleArray[1]
        AtomicNum = ParticleArray[0]
        AntiCheck = ParticleArray[1]
  
        NMname = Station
        Rigidities = [0,0,0]
  
        FileName = NMname + ".csv"
        Rigidities = Rigidities = OTSOLib.cutoff(Position, StartRigidity, EndRigidity, RigidityStep, y, model, IntModel, AtomicNum, AntiCheck, I, Wind, Magnetopause, 
                                                 FileName, CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, g, h, MHDCoordSys)
  
        lat_long_pairs = []
        P_List = []
  
        if asymptotic == "YES":
          Energy_List = asymlevels.copy()
          if unit == "GeV":   
            E_0 = 0.938
            for i in Energy_List:
              R = (i**2 + 2*i*E_0)**(0.5)
              P_List.append(R)
            P_List.insert(0, Rigidities[1])
          else:
            if unit == "GV": 
              P_List = Energy_List.copy()
            P_List.insert(0, Rigidities[1])
          for P in P_List:
              if P == P_List[0]:
                  while True:
                      bool, Lat, Long = OTSOLib.trajectory(Position, P, y, model, IntModel, AtomicNum, AntiCheck, I, Wind, Magnetopause, FileName, 
                                                           CoordinateSystem, MaxStepPercent, EndParams, g, h, MHDCoordSys)
                      P += RigidityStep
                      if bool == 1:
                          lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
                          break 
              else:
                  bool, Lat, Long = OTSOLib.trajectory(Position, P, y, model, IntModel, AtomicNum, AntiCheck, I, Wind, Magnetopause, FileName, 
                                                       CoordinateSystem, MaxStepPercent, EndParams, g, h, MHDCoordSys)
                  lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
          
          asymlevels_with_units = [f"{level} [{unit}]" for level in asymlevels]
          defualt_headers = ["Date","Latitude","Longitude","Altitude","Rc GV","Rc Asym"]
          headers = defualt_headers + asymlevels_with_units
          df = pd.DataFrame(columns=headers)
  
          formatted_list = [f"{bool};{lat};{long}" for bool, lat, long in lat_long_pairs]
  
          formatted_list.insert(0, datetimeobj)
          formatted_list.insert(1, Position[1])
          formatted_list.insert(2, Position[2])
          formatted_list.insert(3, Position[0])
          formatted_list.insert(4, Rigidities[1])
          writer.writerow(formatted_list)
          queue.put(1)
  
        else:
           headers = ["Date","Latitude", "Longitude","Altitude", "Rl", "Rc", "Ru"]
           data = [datetimeobj,x[1],x[2],x[3],Rigidities[0],Rigidities[1], Rigidities[2]]
           writer.writerow(data)
           queue.put(1)
    file.close()
  return


def fortrancallMagfield(Data, DateArray, model, IOPT, WindArray, CoordinateSystem, queue,g,h,MHDfile, MHDCoordSys):
    
  if model[1] == 99:
    MHDinit.MHDinitialise(MHDfile)

  for x in Data:
    Position = x
    
    Bfield = OTSOLib.magstrength(Position, DateArray, model, IOPT, WindArray, CoordinateSystem,MHDCoordSys, g, h)
    Bfield = Bfield*10**9
    combined_array = np.concatenate((Position, Bfield))
    coord_suffix = f"_{CoordinateSystem}"  # Add the coordinate system suffix
    columns = [f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]", f"GSM_Bx [nT]", f"GSM_By [nT]", f"GSM_Bz [nT]"]
    df = pd.DataFrame(columns=columns)
    df.loc[len(df)] = combined_array
    queue.put(df)

  return

def fortrancallCoordtrans(Data, DateArray, CoordIN, CoordOUT, queue):
    for x,y in zip(Data,DateArray):
      Position = [1+(x[2]/6371.0),x[0],x[1]]

      datetimeobj = date.convert_to_datetime(y)

      year = y[0]
      day = y[1]
      hour = y[2]
      minute = y[3]
      secint = y[4]
      sectot = y[5]
      
      Coords = OTSOLib.coordtrans(Position,year,day,hour,minute,secint,sectot,CoordIN,CoordOUT)
      combined_array = np.concatenate(([datetimeobj], Position, Coords))

      coord_suffix = f"_{CoordIN}"
      coord_suffix2 = f"_{CoordOUT}"
      columns = ["Date",f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]", f"X{coord_suffix2} [Re]", f"Y{coord_suffix2} [Re]", f"Z{coord_suffix2} [Re]"]
      if CoordOUT == "GDZ" or CoordOUT == "SPH":
         columns = ["Date",f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]", f"altitude{coord_suffix2} [km]", f"latitude{coord_suffix2}", f"longitude{coord_suffix2}"]
      df = pd.DataFrame(columns=columns)
      df.loc[len(df)] = combined_array
      queue.put(df)

    return

def fortrancallTrace(Data, Rigidity, DateArray, model, IntModel, ParticleArray, IOPT, WindArray, Magnetopause, 
                     CoordinateSystem, MaxStepPercent, EndParams, queue, g, h, MHDfile, MHDCoordSys):
    
    if model[1] == 99:
      MHDinit.MHDinitialise(MHDfile)
    
    for x in Data:
      
      Position = [x[3],x[1],x[2],x[4],x[5]]

      AtomicNum = ParticleArray[0]
      AntiCheck = ParticleArray[1]

      Filename = f"{x[1]}_{x[2]}.csv"
      name = f"{x[1]}_{x[2]}"

      OTSOLib.fieldtrace(Position, Rigidity, DateArray, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, MaxStepPercent, EndParams, Filename, g, h, MHDCoordSys)
      Trace = pd.read_csv(Filename)
      coordsystem2 = "GSM"
      columns = Trace.columns

      new_columns = [f"{col}_{CoordinateSystem} [Re]" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]
      if CoordinateSystem == "GDZ" or CoordinateSystem == "SPH":
        column_names = ["alt [km]", "latitude", "longitude"]
        new_columns = [f"{col}_{column_names[i]}" if i < 3 else f"{col}_{coordsystem2} [nT]" for i, col in enumerate(columns)]

      Trace.columns = new_columns

      Trace.columns = new_columns
            
      dataframe_dict = {name: Trace}
  
      queue.put(dataframe_dict)
      os.remove(Filename)
    
    return