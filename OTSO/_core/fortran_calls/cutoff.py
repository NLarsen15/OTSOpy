import pandas as pd
import multiprocessing as mp

from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.cutoff_data import CutoffData

def FortranCutoff(Data: list, CutoffDataInstance: CutoffData, queue: mp.Queue) -> None:
    
    if CutoffDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(CutoffDataInstance.MHDfile)

    for x in Data:
      
      Position = [x[3],x[1],x[2],x[4],x[5]]


      Station = x[0]

      #print(FileDescriptors)

      StartRigidity = CutoffDataInstance.rigidityarray[0]
      EndRigidity = CutoffDataInstance.rigidityarray[1]
      RigidityStep = CutoffDataInstance.rigidityarray[2]
      AtomicNum = CutoffDataInstance.particlearray[0]
      AntiCheck = CutoffDataInstance.particlearray[1]
      DateArray = CutoffDataInstance.datearray
      model = CutoffDataInstance.model
      IntModel = CutoffDataInstance.integrationmodel
      IOPT = CutoffDataInstance.IOPT
      WindArray = CutoffDataInstance.windarray
      Magnetopause = CutoffDataInstance.magnetopause
      CoordinateSystem = CutoffDataInstance.coordsystem
      MaxStepPercent = CutoffDataInstance.maxsteppercent
      EndParams = CutoffDataInstance.endparams
      Rcomp = CutoffDataInstance.Rcomp
      Rscan = CutoffDataInstance.Rscan
      g = CutoffDataInstance.g
      h = CutoffDataInstance.h
      MHDCoordSys = CutoffDataInstance.MHDcoordsys
      spheresize = CutoffDataInstance.spheresize
      inputcoord = CutoffDataInstance.inputcoord
      trapdist = CutoffDataInstance.mindist
      adapt = CutoffDataInstance.adapt
      Berr = CutoffDataInstance.Berr
      totalbetacheck = CutoffDataInstance.totalbetacheck

      Rigidities = OTSOLib.cutoff(Position, StartRigidity, EndRigidity, RigidityStep, DateArray, model, IntModel, 
                                  AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause,
                                  CoordinateSystem, MaxStepPercent, EndParams, Rcomp, Rscan, g, h, 
                                  MHDCoordSys,spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
      # Add asymptotic computation similar to planet.py implementation
      asymptotic_data = []
      
      if CutoffDataInstance.asymptotic == "YES":
          print("Starting Asymptotic Computation...")
          Energy_List = CutoffDataInstance.asymlevels.copy()
          P_List = []
          
          if CutoffDataInstance.unit == "GeV":   
              E_0 = 0.938  # Rest mass energy of proton in GeV
              for i in Energy_List:
                  R = (i**2 + 2*i*E_0)**(0.5)
                  P_List.append(R)
              P_List.insert(0, Rigidities[1])  # Insert cutoff rigidity (Rc) at beginning
          else:
              if CutoffDataInstance.unit == "GV": 
                  P_List = Energy_List.copy()
              P_List.insert(0, Rigidities[1])  # Insert cutoff rigidity (Rc) at beginning
          
          for P in P_List:
              if P == P_List[0]:  # First value is the cutoff rigidity
                  # Find the first rigidity above cutoff that allows escape
                  while True:
                      bool_result, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum, 
                                                                  AntiCheck, IOPT, WindArray, Magnetopause,
                                                                CoordinateSystem, MaxStepPercent, EndParams, g, h, 
                                                                MHDCoordSys, spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
                      P += RigidityStep
                      if bool_result == 1:
                          asymptotic_data.append([bool_result, round(Lat, 3), round(Long, 3), P - RigidityStep])
                          break 
              else:
                  bool_result, Lat, Long = OTSOLib.trajectory(Position, P, DateArray, model, IntModel, AtomicNum, 
                                                              AntiCheck, IOPT, WindArray, Magnetopause,
                                                            CoordinateSystem, MaxStepPercent, EndParams, g, h, 
                                                            MHDCoordSys, spheresize, inputcoord, trapdist, adapt, Berr, totalbetacheck)
                  asymptotic_data.append([bool_result, round(Lat, 3), round(Long, 3), P])
      
      # Create result DataFrame with rigidities
      Rigiditydataframe = pd.DataFrame({Station: Rigidities}, index=['Ru', 'Rc', 'Rl'])
      
      # If asymptotic computation was performed, create separate asymptotic DataFrame
      if CutoffDataInstance.asymptotic and asymptotic_data:
          # Create asymptotic level headers with units (add "Rc Asym" for cutoff rigidity + asymlevels)
          asymlevels_with_units = ["Rc Asym"] + [f"{level} [{CutoffDataInstance.unit}]" for level in CutoffDataInstance.asymlevels]
          
          # Format asymptotic data as strings "bool;lat;long"
          formatted_asymptotic_data = [f"{data[0]}{CutoffDataInstance.delim}{data[1]}{CutoffDataInstance.delim}{data[2]}" for data in asymptotic_data]
          
          # Create row data: Station name + asymptotic data for each level
          row_data = [Station] + formatted_asymptotic_data
          
          # Column headers: Station + asymptotic levels with units
          column_headers = ["Station"] + asymlevels_with_units
          
          # Create DataFrame with single row
          asymptotic_df = pd.DataFrame([row_data], columns=column_headers)
          
          # Return both dataframes
          queue.put([Rigiditydataframe, asymptotic_df])
      else:
          # Return cutoff dataframe and None for asymptotic
          queue.put([Rigiditydataframe, None])
    return