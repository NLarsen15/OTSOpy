import pandas as pd
import multiprocessing as mp

from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.cone_data import ConeData

def FortranCone(Data: list, ConeDataInstance: ConeData, queue: mp.Queue) -> None:
    
    if ConeDataInstance.model[1] == 99:
      mhd_utils.MHDinitialise(ConeDataInstance.MHDfile)
    
    for x in Data:
      
      Position = [x[3],x[1],x[2],x[4],x[5]]
      Station = x[0]

      NMname = Station

      StartRigidity = ConeDataInstance.rigidityarray[0]
      EndRigidity = ConeDataInstance.rigidityarray[1]
      RigidityStep = ConeDataInstance.rigidityarray[2]
      AtomicNum = ConeDataInstance.particlearray[0]
      AntiCheck = ConeDataInstance.particlearray[1]
      DateArray = ConeDataInstance.datearray
      model = ConeDataInstance.model
      IntModel = ConeDataInstance.integrationmodel
      IOPT = ConeDataInstance.IOPT
      WindArray = ConeDataInstance.windarray
      Magnetopause = ConeDataInstance.magnetopause
      CoordinateSystem = ConeDataInstance.coordsystem
      MaxStepPercent = ConeDataInstance.maxsteppercent
      EndParams = ConeDataInstance.endparams
      g = ConeDataInstance.g
      h = ConeDataInstance.h
      MHDCoordSys = ConeDataInstance.MHDcoordsys
      spheresize = ConeDataInstance.spheresize
      inputcoord = ConeDataInstance.inputcoord
      trapdist = ConeDataInstance.mindist
      adapt = ConeDataInstance.adapt
      Berr = ConeDataInstance.Berr
      totalbetacheck = ConeDataInstance.totalbetacheck


      length = int((StartRigidity-EndRigidity)/RigidityStep)

      Cone, Rigidities = OTSOLib.cone(Position, StartRigidity, EndRigidity, RigidityStep, DateArray, model, 
                                      IntModel, AtomicNum, AntiCheck, IOPT, WindArray,
                                       Magnetopause, CoordinateSystem, MaxStepPercent, 
                                       EndParams, length, g, h, MHDCoordSys,spheresize, 
                                       inputcoord, trapdist, adapt, Berr, totalbetacheck)
      
      decoded_lines = [line.decode('utf-8').strip() for sublist in Cone for line in sublist]

      data = [line.split() for line in decoded_lines]

      Conedf = pd.DataFrame(data, columns=['R [GV]', 'Filter', 'ALat', 'ALong'])
      Conedf[NMname] = Conedf[['Filter', 'ALat', 'ALong']].astype(str).agg(f'{ConeDataInstance.delim}'.join, axis=1)
      Conedf.drop(columns=['Filter', 'ALat', 'ALong'], inplace=True)

      Rigiditydataframe = pd.DataFrame({Station: Rigidities}, index=['Ru', 'Rc', 'Rl'])
 
      queue.put([Conedf, Rigiditydataframe])
    
    return