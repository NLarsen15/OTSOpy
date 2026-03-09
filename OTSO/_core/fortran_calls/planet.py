import csv
import multiprocessing as mp

from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.planet_data import PlanetData

def FortranPlanet(Data, PlanetDataInstance: 'PlanetData', PlanetFile: str, queue: mp.Queue) -> None:
  with open(PlanetFile, mode='a', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)

    if PlanetDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(PlanetDataInstance.MHDfile)

    for x in Data:
        
        Position = [x[3],x[1],x[2],x[4],x[5]]
        Station = x[0]

        AtomicNum = PlanetDataInstance.particlearray[0]
        AntiCheck = PlanetDataInstance.particlearray[1]
        Rigidity = PlanetDataInstance.rigidityarray
        AtomicNum = PlanetDataInstance.particlearray[0]
        AntiCheck = PlanetDataInstance.particlearray[1]
        DateArray = PlanetDataInstance.datearray
        model = PlanetDataInstance.model
        IntModel = PlanetDataInstance.integrationmodel
        IOPT = PlanetDataInstance.IOPT
        WindArray = PlanetDataInstance.windarray
        Magnetopause = PlanetDataInstance.magnetopause
        MaxStepPercent = PlanetDataInstance.maxsteppercent
        EndParams = PlanetDataInstance.endparams
        Rcomp = PlanetDataInstance.Rcomp
        Rscan = PlanetDataInstance.Rscan
        g = PlanetDataInstance.g
        h = PlanetDataInstance.h
        MHDCoordSys = PlanetDataInstance.MHDcoordsys
        spheresize = PlanetDataInstance.spheresize
        inputcoord = PlanetDataInstance.inputcoord
        trapdist = PlanetDataInstance.mindist
        adapt = PlanetDataInstance.adapt
        Berr = PlanetDataInstance.Berr
        totalbetacheck = PlanetDataInstance.totalbetacheck
  
        Rigidities = [0,0,0]

        Rigidities = OTSOLib.planet(Position, Rigidity, DateArray, model, IntModel, AtomicNum, 
                                    AntiCheck, IOPT, WindArray, Magnetopause,
                                     MaxStepPercent, EndParams, Rcomp, Rscan, g, h, 
                                     MHDCoordSys,spheresize, inputcoord,trapdist, adapt, Berr, totalbetacheck)
        #print(f"Rigidities returned from Fortran: {Rigidities}")
        CoordinateSystem = "GEO"
  
        lat_long_pairs = []
        P_List = []
  
        if PlanetDataInstance.asymptotic == "YES":
          Energy_List = PlanetDataInstance.asymlevels.copy()
          if PlanetDataInstance.unit == "GeV":   
            E_0 = 0.938
            for i in Energy_List:
              R = (i**2 + 2*i*E_0)**(0.5)
              P_List.append(R)
            P_List.insert(0, Rigidities[1])
          else:
            if PlanetDataInstance.unit == "GV": 
              P_List = Energy_List.copy()
            P_List.insert(0, Rigidities[1])
          for P in P_List:
              if P == P_List[0]:  # Check if it's the first value in P_List
                  while True:
                      bool, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum,
                                                            AntiCheck, IOPT, WindArray, Magnetopause,
                                                            CoordinateSystem, MaxStepPercent, EndParams, 
                                                            g, h, MHDCoordSys,spheresize, 
                                                            inputcoord, trapdist, adapt, Berr, totalbetacheck)
                      P += Rigidity[2]
                      if bool == 1:
                          lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
                          break 
              else:
                  bool, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum, 
                                                       AntiCheck, IOPT, WindArray, Magnetopause,
                                                        CoordinateSystem, MaxStepPercent, EndParams, 
                                                        g, h, MHDCoordSys,spheresize, 
                                                        inputcoord, trapdist, adapt, Berr, totalbetacheck)
                  lat_long_pairs.append([bool, round(Lat, 3), round(Long, 3)])
  
          formatted_list = [f"{bool}{PlanetDataInstance.delim}{lat}{PlanetDataInstance.delim}{long}" for bool, lat, long in lat_long_pairs]
  
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