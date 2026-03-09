import pandas as pd
import os
import multiprocessing as mp

from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.trajectory_data import TrajectoryData

def FortranTrajectory(Data: list, TrajectoryDataInstance: 'TrajectoryData', queue: mp.Queue) -> None:
  
  if TrajectoryDataInstance.model[1] == 99:
    mhd_utils.MHDinitialise(TrajectoryDataInstance.MHDfile)
  
  for x in Data:
    
    Position = [x[3],x[1],x[2],x[4],x[5]]
    Station = x[0]

    AtomicNum = TrajectoryDataInstance.particlearray[0]
    AntiCheck = TrajectoryDataInstance.particlearray[1]
    Rigidity = TrajectoryDataInstance.rigidity
    AtomicNum = TrajectoryDataInstance.particlearray[0]
    AntiCheck = TrajectoryDataInstance.particlearray[1]
    DateArray = TrajectoryDataInstance.datearray
    model = TrajectoryDataInstance.model
    IntModel = TrajectoryDataInstance.integrationmodel
    IOPT = TrajectoryDataInstance.IOPT
    WindArray = TrajectoryDataInstance.windarray
    Magnetopause = TrajectoryDataInstance.magnetopause
    CoordinateSystem = TrajectoryDataInstance.coordsystem
    MaxStepPercent = TrajectoryDataInstance.maxsteppercent
    EndParams = TrajectoryDataInstance.endparams
    g = TrajectoryDataInstance.g
    h = TrajectoryDataInstance.h
    MHDCoordSys = TrajectoryDataInstance.MHDcoordsys
    spheresize = TrajectoryDataInstance.spheresize
    inputcoord = TrajectoryDataInstance.inputcoord
    trapdist = TrajectoryDataInstance.mindist
    adapt = TrajectoryDataInstance.adapt
    Berr = TrajectoryDataInstance.Berr
    totalbetacheck = TrajectoryDataInstance.totalbetacheck

    NMname = Station

    FileName = NMname + ".csv"
    bool_val, lat, long = OTSOLib.trajectory_full(Position, Rigidity, DateArray, 
                            model, IntModel, AtomicNum, 
                            AntiCheck, IOPT, WindArray, 
                            Magnetopause, FileName, 
                            CoordinateSystem, MaxStepPercent, 
                            EndParams, g, h, MHDCoordSys,spheresize, 
                            inputcoord,trapdist, adapt, Berr, totalbetacheck)
    
    Trajectory = pd.read_csv(FileName)
    Trajectory.columns = [f"{col}_Re [{CoordinateSystem}]" for col in Trajectory.columns]

    result = {
        "NMname": NMname,
        "trajectory": Trajectory,
        "Filter": bool_val,
        "AsymLat": lat,
        "AsymLong": long
    }

    queue.put(result)
    os.remove(FileName)

  return