import numpy as np

from OTSO._core.libs.MiddleMan import Middleman as Mid
from ..data_classes.cutoff_data import CutoffData
from ..data_classes.flight_data import FlightData
from ..data_classes.trace_data import TraceData
from ..data_classes.skymap_data import SkymapData
from ..data_classes.trajectory_data import TrajectoryData
from ..utils import cutoff_utils as cu

##########################################################################################################
def prepare_fortran_cutoff(Positioninfo: list, DataInstance: CutoffData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    StartRigidity = DataInstance.rigidityarray[0]
    EndRigidity = DataInstance.rigidityarray[1]
    RigidityStep = DataInstance.rigidityarray[2]
    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsystem
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    Rcomp = DataInstance.Rcomp
    Rscan = DataInstance.Rscan
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    FortranThreads = DataInstance.threadnum
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)

    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()

    # Copy Python dataclass values into the Fortran object
    FortranData.startrigidity = StartRigidity
    FortranData.endrigidity = EndRigidity
    FortranData.rigiditystep = RigidityStep

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams

    FortranData.intmode = IntModel

    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)

    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep

    FortranData.rcomputation = Rcomp
    FortranData.scanchoice = Rscan

    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord

    FortranData.trapdist = trapdist

    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False

    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck

    FortranData.fortranthreads = FortranThreads

    n = cu.compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    return FortranData

def prepare_fortran_cone(Positioninfo: list, DataInstance: CutoffData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    StartRigidity = DataInstance.rigidityarray[0]
    EndRigidity = DataInstance.rigidityarray[1]
    RigidityStep = DataInstance.rigidityarray[2]
    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsystem
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    FortranThreads = DataInstance.threadnum
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)

    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()

    # Copy Python dataclass values into the Fortran object
    FortranData.startrigidity = StartRigidity
    FortranData.endrigidity = EndRigidity
    FortranData.rigiditystep = RigidityStep

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams

    FortranData.intmode = IntModel

    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)

    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep

    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord

    FortranData.trapdist = trapdist

    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False

    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck

    FortranData.fortranthreads = FortranThreads

    n = cu.compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    return FortranData
#############################################################################################################################
#############################################################################################################################
def prepare_fortran_trajectory(Positioninfo: list, DataInstance: TrajectoryData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsystem
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)

    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams

    FortranData.intmode = IntModel

    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)

    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep

    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord

    FortranData.trapdist = trapdist

    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False

    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck

    return FortranData
#############################################################################################################################
#############################################################################################################################
def prepare_fortran_magfield(DataInstance: CutoffData):

    DateArray = DataInstance.datearray
    model = DataInstance.model
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)
    
    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()
    
    # Copy Python dataclass values into the Fortran object

    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.iopt = IOPT

    FortranData.mode = np.asarray(model, dtype=np.int32)
    
    FortranData.mhdcoordsys = MHDCoordSys

    return FortranData
#############################################################################################################################
#############################################################################################################################
def prepare_fortran_Coordtrans(Date, maxdegree, g, h):

    
    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()
    
    # Copy Python dataclass values into the Fortran object

    FortranData.date = Date
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = len(g)
    

    return FortranData
##########################################################################################################
def prepare_fortran_flight(Positioninfo: list, 
                           DataInstance: FlightData,
                           DateArrayinput,
                           WindArrayinput,
                           IOPTinput,
                           GArrayinput,
                           HArrayinput):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    StartRigidity = DataInstance.rigidityarray[0]
    EndRigidity = DataInstance.rigidityarray[1]
    RigidityStep = DataInstance.rigidityarray[2]
    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DateArrayinput
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = IOPTinput
    WindArray = WindArrayinput
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsystem
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    Rcomp = DataInstance.Rcomp
    Rscan = DataInstance.Rscan
    g = GArrayinput
    h = HArrayinput
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    FortranThreads = DataInstance.threadnum
    maxdegree = DataInstance.max_degree
    gaussianlength = len(GArrayinput)
    
    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()
    
    # Copy Python dataclass values into the Fortran object
    FortranData.startrigidity = StartRigidity
    FortranData.endrigidity = EndRigidity
    FortranData.rigiditystep = RigidityStep

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams
    
    FortranData.intmode = IntModel
    
    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)
    
    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep

    FortranData.rcomputation = Rcomp
    FortranData.scanchoice = Rscan

    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord

    FortranData.trapdist = trapdist

    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False

    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck

    FortranData.fortranthreads = FortranThreads

    n = cu.compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    return FortranData
################################################################################################################################
################################################################################################################################
def prepare_fortran_trace(Positioninfo: list, DataInstance: TraceData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsys
    MaxStepPercent = DataInstance.maxsteppercent
    EndParams = DataInstance.endparams
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)
    
    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams
    
    FortranData.intmode = IntModel
    
    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)
    
    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem
    
    FortranData.gyropercent = MaxStepPercent
    
    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord
    

    return FortranData
################################################################################################################################
def prepare_fortran_transmission(Positioninfo: list, DataInstance: CutoffData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    StartRigidity = DataInstance.rigidityarray[0]
    EndRigidity = DataInstance.rigidityarray[1]
    RigidityStep = DataInstance.rigidityarray[2]
    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    #CoordinateSystem = "GSM"
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    FortranThreads = DataInstance.threadnum
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)
    
    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()
    
    # Copy Python dataclass values into the Fortran object
    FortranData.startrigidity = StartRigidity
    FortranData.endrigidity = EndRigidity
    FortranData.rigiditystep = RigidityStep

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams
    
    FortranData.intmode = IntModel
    
    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)
    
    FortranData.pause = Magnetopauseinput

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep
    #FortranData.coordsystem = CoordinateSystem
    
    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord
    
    FortranData.trapdist = trapdist
    
    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False
    
    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck
    
    FortranData.fortranthreads = FortranThreads

    n = cu.compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    return FortranData
#############################################################################################################################
def prepare_fortran_skymap(Positioninfo: list, DataInstance: SkymapData):

    Position = [Positioninfo[3], Positioninfo[1], Positioninfo[2], Positioninfo[4], Positioninfo[5]]

    StartRigidity = DataInstance.rigidityarray[0]
    EndRigidity = DataInstance.rigidityarray[1]
    RigidityStep = DataInstance.rigidityarray[2]
    AtomicNum = DataInstance.particlearray[0]
    AntiCheck = DataInstance.particlearray[1]
    DateArray = DataInstance.datearray
    model = DataInstance.model
    IntModel = DataInstance.integrationmodel
    IOPT = DataInstance.IOPT
    WindArray = DataInstance.windarray
    Magnetopauseinput = DataInstance.magnetopauseinput
    CoordinateSystem = DataInstance.coordsystem
    MaxStepPercent = DataInstance.maxsteppercent
    FixedStep = DataInstance.fixedstep
    EndParams = DataInstance.endparams
    Rcomp = DataInstance.Rcomp
    Rscan = DataInstance.Rscan
    g = DataInstance.g
    h = DataInstance.h
    MHDCoordSys = DataInstance.MHDcoordsys
    spheresize = DataInstance.spheresize
    inputcoord = DataInstance.inputcoord
    trapdist = DataInstance.mintrapdist
    adapt = DataInstance.adaptivestep
    Berr = DataInstance.betaerror
    totalbetacheck = DataInstance.totalbetacheck
    FortranThreads = DataInstance.threadnum
    maxdegree = DataInstance.max_degree
    gaussianlength = len(DataInstance.g)

    # Create the f90wrap Fortran derived type object
    FortranData = Mid.FortranData()

    # Copy Python dataclass values into the Fortran object
    FortranData.startrigidity = StartRigidity
    FortranData.endrigidity = EndRigidity
    FortranData.rigiditystep = RigidityStep

    FortranData.positionin = Position
    FortranData.wind = WindArray
    FortranData.date = DateArray
    FortranData.g = g
    FortranData.h = h
    FortranData.maxdegree = maxdegree
    FortranData.gaussianlength = gaussianlength
    FortranData.end = EndParams

    FortranData.intmode = IntModel

    FortranData.atomicnumber = AtomicNum
    FortranData.anti = AntiCheck
    FortranData.iopt = IOPT
    FortranData.mode = np.asarray(model, dtype=np.int32)

    FortranData.pause = Magnetopauseinput
    FortranData.coordsystem = CoordinateSystem

    FortranData.gyropercent = MaxStepPercent
    FortranData.fixedstepsize = FixedStep

    FortranData.scanchoice = Rscan
    
    FortranData.mhdcoordsys = MHDCoordSys
    FortranData.sphere = spheresize
    FortranData.inputcoord = inputcoord
    
    FortranData.trapdist = trapdist
    
    if trapdist > 0:
       FortranData.trapdistcheck = True
    else:
       FortranData.trapdistcheck = False
    
    FortranData.adapt = adapt
    FortranData.berr = Berr
    FortranData.totalbetacheck = totalbetacheck
    
    FortranData.fortranthreads = FortranThreads

    n = cu.compute_number_of_rigidities(FortranData.startrigidity, FortranData.endrigidity, FortranData.rigiditystep)

    FortranData.n = n
    return FortranData