import csv
import pandas as pd
import multiprocessing as mp

from ..custom_classes import date
from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.flight_data import FlightData

def FortranFlight(Data: list, DateArray: list, IOPT: list, WindArray: list, GArray: list, HArray: list,
                  FlightFile: str, FlightDataInstance: FlightData, queue: mp.Queue) -> None:
  
  with open(FlightFile, mode='a', newline='', encoding='utf-8') as file:
    if FlightDataInstance.asymptotic == "YES":
        asymlevels_with_units = [f"{level} [{FlightDataInstance.unit}]" for level in FlightDataInstance.asymlevels]
        default_headers = ["Date","Latitude","Longitude","Altitude","Rc GV","Rc Asym"]
        headers = default_headers + asymlevels_with_units
    else:
        headers = ["Date","Latitude", "Longitude","Altitude", "Ru", "Rc", "Rl"]

    writer = csv.writer(file)
    writer.writerow(headers)

    if FlightDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(FlightDataInstance.MHDfile)


    for x,y,z,I,G,H in zip(Data,DateArray,WindArray,IOPT,GArray,HArray):
        
        Position = [x[3],x[1],x[2],x[4],x[5]]
        Station = x[0]
        
        datetimeobj = date.convert_to_datetime(y)
        Wind = z
  
        StartRigidity = FlightDataInstance.rigidityarray[0]
        EndRigidity = FlightDataInstance.rigidityarray[1]
        RigidityStep = FlightDataInstance.rigidityarray[2]
        AtomicNum = FlightDataInstance.particlearray[0]
        AntiCheck = FlightDataInstance.particlearray[1]
        AtomicNum = FlightDataInstance.particlearray[0]
        AntiCheck = FlightDataInstance.particlearray[1]

        CoordinateSystem = FlightDataInstance.coordsystem
        MaxStepPercent = FlightDataInstance.maxsteppercent
        EndParams = FlightDataInstance.endparams
        Magnetopause = FlightDataInstance.magnetopause
        Rcomp = FlightDataInstance.Rcomp
        Rscan = FlightDataInstance.Rscan
        model = FlightDataInstance.model
        IntModel = FlightDataInstance.integrationmodel
        g = G
        h = H
        MHDCoordSys = FlightDataInstance.MHDcoordsys
        spheresize = FlightDataInstance.spheresize
        inputcoord = FlightDataInstance.inputcoord
        trapdist = FlightDataInstance.mindist
        adapt = FlightDataInstance.adapt
        Berr = FlightDataInstance.Berr
        totalbetacheck = FlightDataInstance.totalbetacheck
  
        NMname = Station
        Rigidities = [0,0,0]
  
        FileName = NMname + ".csv"
        Rigidities = Rigidities = OTSOLib.cutoff(Position, StartRigidity, EndRigidity, RigidityStep, y, 
                                                 model, IntModel, AtomicNum, AntiCheck, I, Wind, Magnetopause, 
                                                 CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, g, h, 
                                                 MHDCoordSys,spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
  
        lat_long_pairs = []
        P_List = []
  
        if FlightDataInstance.asymptotic == "YES":
          Energy_List = FlightDataInstance.asymlevels.copy()
          if FlightDataInstance.unit == "GeV":   
            E_0 = 0.938
            for i in Energy_List:
              R = (i**2 + 2*i*E_0)**(0.5)
              P_List.append(R)
            P_List.insert(0, round(Rigidities[1], 5))
          else:
            if FlightDataInstance.unit == "GV": 
              P_List = Energy_List.copy()
            P_List.insert(0, Rigidities[1])
          for P in P_List:
              if P == P_List[0]:
                  while True:
                      bool, Lat, Long = OTSOLib.trajectory(Position, P, y, model, IntModel, AtomicNum, AntiCheck, 
                                                           I, Wind, Magnetopause, 
                                                           CoordinateSystem, MaxStepPercent, EndParams, g, h, 
                                                           MHDCoordSys,spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
                      P += RigidityStep
                      if bool == 1:
                          lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
                          break 
              else:
                  bool, Lat, Long = OTSOLib.trajectory(Position, P, y, model, IntModel, AtomicNum, AntiCheck, I, 
                                                       Wind, Magnetopause, 
                                                       CoordinateSystem, MaxStepPercent, EndParams, g, h, 
                                                       MHDCoordSys,spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
                  lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
          
          asymlevels_with_units = [f"{level} [{FlightDataInstance.unit}]" for level in FlightDataInstance.asymlevels]
          default_headers = ["Date","Latitude","Longitude","Altitude","Rc GV","Rc Asym"]
          headers = default_headers + asymlevels_with_units
  
          formatted_list = [f"{bool}{FlightDataInstance.delim}{lat}{FlightDataInstance.delim}{long}" for bool, lat, long in lat_long_pairs]
  
          formatted_list.insert(0, datetimeobj)
          formatted_list.insert(1, round(float(Position[1]),5))
          formatted_list.insert(2, round(float(Position[2]),5))
          formatted_list.insert(3, round(float(Position[0]),5))
          formatted_list.insert(4, round(Rigidities[1],5))
          writer.writerow(formatted_list)
          queue.put(1)
  
        else:
           headers = ["Date","Latitude", "Longitude","Altitude", "Ru", "Rc", "Rl"]
           data = [datetimeobj,x[1],x[2],x[3],Rigidities[0],Rigidities[1], Rigidities[2]]
           writer.writerow(data)
           queue.put(1)
    file.close()
  return