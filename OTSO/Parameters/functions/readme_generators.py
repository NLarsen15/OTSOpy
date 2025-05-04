from datetime import date
import os
import shutil
from . import misc
from . import date as d
import pandas as pd

def READMECone(UsedStationstemp, RigidityArray, EventDate, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, Printtime, MaxStepPercent, EndParams, LiveData, serverdata, kp):
    result = []

    pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

    if LiveData == 1:
       OnlineData = "Online Space Weather Data Used"
    elif serverdata == "ON":
       OnlineData = "Server Data Generated Using OMNI Database Used"
    else:
       OnlineData = "User Inputted Data Used"
    
    particle = "anti-particle" if AntiCheck == 1 else "Normal Particle"
    
    IntegrationMethods = ["4th Order Runge-Kutta Method", "Boris Method", "Vay Method", "Higuera-Cary Method"]
    IntegrationMethod = IntegrationMethods[IntModel-1] if 0 <= IntModel <= 4 else "Unknown Integration Method"

    InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
    Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"

    ExternalModels = [
        "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
        "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
    ]
    External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
    if model[1] == 99:
         External = "MHD file"

    PauseModels = [
        "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
        "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
        "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
    ]
    PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Output Coordinate System:\n{CoordinateSystem}\n\n")
    result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Max Time Step [% of gyrofrequency]: {MaxStepPercent}\n\n")
    result.append(f"Minimum Altitude: {EndParams[0]}km\n")
    result.append(f"Max Distance: {EndParams[1]}Re\n")
    result.append(f"Max Time: {EndParams[2]}s\n\n")
    result.append(f"Start Altitude = {UsedStationstemp[0][3]}km \n")
    result.append(f"Zenith = {UsedStationstemp[0][4]}\n")
    result.append(f"Azimuth = {UsedStationstemp[0][5]}\n\n")
    result.append(f"Kp = {kp}\n")
    result.append(f"IOPT = {IOPT}\n\n")
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(f"Vx = {abs(WindArray[0])}\n")
    result.append(f"Vy = {WindArray[1]}\n")
    result.append(f"Vz = {WindArray[2]}\n\n")
    result.append(f"IMF [nT]:\n")
    result.append(f"By = {WindArray[3]}\n")
    result.append(f"Bz = {WindArray[4]}\n\n")
    result.append(f"Density = {WindArray[5]} cm^-3\n")
    result.append(f"Pdyn = {pressure} nPa\n\n")
    result.append(f"Dst = {WindArray[7]} nT\n\n")
    result.append(f"G1 = {WindArray[8]}\n")
    result.append(f"G2 = {WindArray[9]}\n")
    result.append(f"G3 = {WindArray[10]}\n\n")
    result.append(f"W1 = {WindArray[11]}\n")
    result.append(f"W2 = {WindArray[12]}\n")
    result.append(f"W3 = {WindArray[13]}\n")
    result.append(f"W4 = {WindArray[14]}\n")
    result.append(f"W5 = {WindArray[15]}\n")
    result.append(f"W6 = {WindArray[16]}\n\n")
    result.append(f"Atomic Number = {AtomicNum}\n\n")
    result.append(f"Particle Type = {particle}\n\n")
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"External Model = {External}\n\n")
    result.append(f"Magnetopause Model = {PauseModel}\n\n")
    result.append(f"Rigidity:\n")
    result.append(f"Start = {RigidityArray[0]} [GV]\n")
    result.append(f"End = {RigidityArray[1]} [GV]\n")
    result.append(f"Step = {RigidityArray[2]} [GV]\n\n")
    result.append(f"Stations:\n")

    for station in UsedStationstemp:
        result.append(f"{station[0]}, Latitude: {station[1]}, Longitude: {station[2]}\n")

    return "".join(result)

def READMECutoff(UsedStationstemp, RigidityArray, EventDate, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, Printtime, MaxStepPercent, EndParams, Rcomp, Rscan, LiveData, serverdata, kp):
    
    result = []

    pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

    if LiveData == 1:
      OnlineData = "Online Space Weather Data Used"
    elif serverdata == "ON":
       OnlineData = "Server Data Generated Using OMNI Database Used"
    else:
       OnlineData = "User Inputted Data Used"
    
    if Rcomp == "Vertical":
      CutoffComp = "Vertical Cutoff Rigidity"
    elif Rcomp == "Apparent":
      CutoffComp = "Apparent Cutoff Rigidity"
    else:
      CutoffComp = "Custom Cutoff"

    particle = "anti-particle" if AntiCheck == 1 else "Normal Particle"
    RigidityScan = "Rigidity Scan Used" if Rscan != 0 else "No Rigidity Scan"
    
    IntegrationMethods = ["4th Order Runge-Kutta Method", "Boris Method", "Vay Method", "Higuera-Cary Method"]
    IntegrationMethod = IntegrationMethods[IntModel-1] if 0 <= IntModel <= 4 else "Unknown Integration Method"

    InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
    Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"

    ExternalModels = [
        "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
        "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
    ]
    External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
    if model[1] == 99:
         External = "MHD file"

    PauseModels = [
        "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
        "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
        "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
    ]
    PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Cutoff Computed: {CutoffComp}\n\n")
    result.append(f"Output Coordinate System:\n{CoordinateSystem}\n\n")
    result.append(f"Rigidity Scan:\n{RigidityScan}\n\n")
    result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Max Time Step [% of gyrofrequency]: {MaxStepPercent}\n\n")
    result.append(f"Minimum Altitude: {EndParams[0]}km\n")
    result.append(f"Max Distance: {EndParams[1]}Re\n")
    result.append(f"Max Time: {EndParams[2]}s\n\n")
    result.append(f"Start Altitude = {UsedStationstemp[0][3]}km \n")
    result.append(f"Zenith = {UsedStationstemp[0][4]}\n")
    result.append(f"Azimuth = {UsedStationstemp[0][5]}\n\n")
    result.append(f"Kp = {kp}\n")
    result.append(f"IOPT = {IOPT}\n\n")
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(f"Vx = {abs(WindArray[0])}\n")
    result.append(f"Vy = {WindArray[1]}\n")
    result.append(f"Vz = {WindArray[2]}\n\n")
    result.append(f"IMF [nT]:\n")
    result.append(f"By = {WindArray[3]}\n")
    result.append(f"Bz = {WindArray[4]}\n\n")
    result.append(f"Density = {WindArray[5]} cm^-3\n")
    result.append(f"Pdyn = {pressure} nPa\n\n")
    result.append(f"Dst = {WindArray[7]} nT\n\n")
    result.append(f"G1 = {WindArray[8]}\n")
    result.append(f"G2 = {WindArray[9]}\n")
    result.append(f"G3 = {WindArray[10]}\n\n")
    result.append(f"W1 = {WindArray[11]}\n")
    result.append(f"W2 = {WindArray[12]}\n")
    result.append(f"W3 = {WindArray[13]}\n")
    result.append(f"W4 = {WindArray[14]}\n")
    result.append(f"W5 = {WindArray[15]}\n")
    result.append(f"W6 = {WindArray[16]}\n\n")
    result.append(f"Atomic Number = {AtomicNum}\n\n")
    result.append(f"Particle Type = {particle}\n\n")
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"External Model = {External}\n\n")
    result.append(f"Magnetopause Model = {PauseModel}\n\n")
    result.append(f"Rigidity:\n")
    result.append(f"Start = {RigidityArray[0]} [GV]\n")
    result.append(f"End = {RigidityArray[1]} [GV]\n")
    result.append(f"Step = {RigidityArray[2]} [GV]\n\n")
    result.append(f"Stations:\n")

    for station in UsedStationstemp:
        result.append(f"{station[0]}, Latitude: {station[1]}, Longitude: {station[2]}\n")

    return "".join(result)

def READMEFlight(Data, Rigidity, model, IntModel, AntiCheck, IOPT, WindArray, Magnetopause, Printtime, MaxStepPercent, EndParams, Rcomp, Rscan, LiveData, asymptotic, asymlevels, unit, serverdata, kp):
      result = []

      if LiveData == 1:
         OnlineData = "Online Space Weather Data Used"
      elif serverdata == "ON":
         OnlineData = "Server Data Generated Using OMNI Database Used"
      else:
         OnlineData = "User Inputted Data Used"
      
      if Rcomp == "Vertical":
         CutoffComp = "Vertical Cutoff Rigidity"
      elif Rcomp == "Apparent":
         CutoffComp = "Apparent Cutoff Rigidity"
      else:
         CutoffComp = "Custom Cutoff"
      
      particle = "anti-particle" if AntiCheck == 1 else "Normal Particle"
      RigidityScan = "Rigidity Scan Used" if Rscan != 0 else "No Rigidity Scan"
      
      IntegrationMethods = ["4th Order Runge-Kutta Method", "Boris Method", "Vay Method", "Higuera-Cary Method"]
      IntegrationMethod = IntegrationMethods[IntModel-1] if 0 <= IntModel <= 4 else "Unknown Integration Method"
  
      InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
      Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"
  
      ExternalModels = [
          "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
          "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
      ]
      External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
      if model[1] == 99:
         External = "MHD file"
  
      PauseModels = [
          "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
          "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
          "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
      ]
      PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"
  
      today = date.today()
      result.append(f"\n")
      result.append(f"Date of OTSO computation: {today}\n")
      result.append(f"Total computation time: {Printtime} seconds\n\n")
      result.append(f"Cutoff Computed: {CutoffComp}\n\n")
      result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
      result.append(f"Rigidity Scan:\n{RigidityScan}\n\n")
      result.append(f"Input Variables:\n\n")
      result.append(f"Data Used: {OnlineData}\n\n")
      result.append(f"Max Time Step [% of gyrofrequency]: {MaxStepPercent}\n\n")
      result.append(f"Minimum Altitude: {EndParams[0]}km\n")
      result.append(f"Max Distance: {EndParams[1]}Re\n")
      result.append(f"Max Time: {EndParams[2]}s\n\n")
      result.append(f"Zenith = {Data[0][4]}\n")
      result.append(f"Azimuth = {Data[0][5]}\n\n")
      result.append(f"Particle Type = {particle}\n\n")
      result.append(f"Magnetic Field Models:\n")
      result.append(f"Internal Model = {Internal}\n")
      result.append(f"External Model = {External}\n\n")
      result.append(f"Magnetopause Model = {PauseModel}\n\n")
      result.append(f"Rigidity:\n")
      result.append(f"Start = {Rigidity[0]} [GV]\n")
      result.append(f"End = {Rigidity[1]} [GV]\n")
      result.append(f"Step = {Rigidity[2]} [GV]\n\n")
      if asymptotic == "YES":
         result.append("Asymptotic Directions: " + str(asymptotic) + " \n")
         result.append("Asymptotic Levels: " + str(asymlevels) + " \n")
         result.append("Asymptotic Levels Unit: " + str(unit) + " \n")
    
      return "".join(result)

def READMEFlightData(DateArray,WindArray,kparray):

   headers = ["Date","kp","Vx [km/s]", "Vy [km/s]", "Vz [km/s]", "By [nT]", "Bz [nT]", "Density [cm^-3]", "Pdyn [nPa]", "Dst [nT]", "G1", "G2",
              "G3", "W1", "W2", "W3", "W4", "W5", "W6"]
   
   rows = []

   for date,Wind,kp in zip(DateArray,WindArray,kparray):
      datetimeobj = d.convert_to_datetime(date)
      pressure = misc.Pdyn_comp(Wind[5],abs(Wind[0]))
      Vx = abs(Wind[0])
      Vy = Wind[1]
      Vz = Wind[2]
      By = Wind[3]
      Bz = Wind[4]
      Density = Wind[5]
      Pdyn = pressure
      Dst = Wind[7]
      G1 = Wind[8]
      G2 = Wind[9]
      G3 = Wind[10]
      W1 = Wind[11]
      W2 = Wind[12]
      W3 = Wind[13]
      W4 = Wind[14]
      W5 = Wind[15]
      W6 = Wind[16]
      data = {
      "Date": datetimeobj,
      "kp": kp,
      "Vx [km/s]": Vx,
      "Vy [km/s]": Vy,
      "Vz [km/s]": Vz,
      "By [nT]": By,
      "Bz [nT]": Bz,
      "Density [cm^-3]": Density,
      "Pdyn [nPa]": Pdyn,
      "Dst [nT]": Dst,
      "G1": G1,
      "G2": G2,
      "G3": G3,
      "W1": W1,
      "W2": W2,
      "W3": W3,
      "W4": W4,
      "W5": W5,
      "W6": W6,}
      
      rows.append(data)

   df = pd.DataFrame(rows, columns=headers)

   return df

def READMETrajectory(UsedStationstemp, Rigidity, EventDate, model, IntModel, AtomicNum, AntiCheck, IOPT, WindArray, Magnetopause, CoordinateSystem, Printtime, MaxStepPercent, EndParams, LiveData, kp, serverdata):
    result = []

    pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

    if LiveData == 1:
      OnlineData = "Online Space Weather Data Used"
    elif serverdata == "ON":
       OnlineData = "Server Data Generated Using OMNI Database Used"
    else:
       OnlineData = "User Inputted Data Used"

    particle = "anti-particle" if AntiCheck == 1 else "Normal Particle"
    
    IntegrationMethods = ["4th Order Runge-Kutta Method", "Boris Method", "Vay Method", "Higuera-Cary Method"]
    IntegrationMethod = IntegrationMethods[IntModel-1] if 0 <= IntModel <= 4 else "Unknown Integration Method"

    InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
    Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"

    ExternalModels = [
        "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
        "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
    ]
    External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
    if model[1] == 99:
         External = "MHD file"

    PauseModels = [
        "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
        "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
        "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
    ]
    PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Output Coordinate System:\n{CoordinateSystem}\n\n")
    result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Max Time Step [% of gyrofrequency]: {MaxStepPercent}\n\n")
    result.append(f"Minimum Altitude: {EndParams[0]}km\n")
    result.append(f"Max Distance: {EndParams[1]}Re\n")
    result.append(f"Max Time: {EndParams[2]}s\n\n")
    result.append(f"Start Altitude = {UsedStationstemp[0][3]}km \n")
    result.append(f"Zenith = {UsedStationstemp[0][4]}\n")
    result.append(f"Azimuth = {UsedStationstemp[0][5]}\n\n")
    result.append(f"Kp = {kp}\n")
    result.append(f"IOPT = {IOPT}\n\n")
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(f"Vx = {abs(WindArray[0])}\n")
    result.append(f"Vy = {WindArray[1]}\n")
    result.append(f"Vz = {WindArray[2]}\n\n")
    result.append(f"IMF [nT]:\n")
    result.append(f"By = {WindArray[3]}\n")
    result.append(f"Bz = {WindArray[4]}\n\n")
    result.append(f"Density = {WindArray[5]} cm^-3\n")
    result.append(f"Pdyn = {pressure} nPa\n\n")
    result.append(f"Dst = {WindArray[7]} nT\n\n")
    result.append(f"G1 = {WindArray[8]}\n")
    result.append(f"G2 = {WindArray[9]}\n")
    result.append(f"G3 = {WindArray[10]}\n\n")
    result.append(f"W1 = {WindArray[11]}\n")
    result.append(f"W2 = {WindArray[12]}\n")
    result.append(f"W3 = {WindArray[13]}\n")
    result.append(f"W4 = {WindArray[14]}\n")
    result.append(f"W5 = {WindArray[15]}\n")
    result.append(f"W6 = {WindArray[16]}\n\n")
    result.append(f"Atomic Number = {AtomicNum}\n\n")
    result.append(f"Particle Type = {particle}\n\n")
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"External Model = {External}\n\n")
    result.append(f"Magnetopause Model = {PauseModel}\n\n")
    result.append(f"Rigidity = {Rigidity}\n\n")
    result.append(f"Stations:\n")

    for station in UsedStationstemp:
        result.append(f"{station[0]}, Latitude: {station[1]}, Longitude: {station[2]}\n")

    return "".join(result)

def READMEPlanet(Data, Rigidity, EventDate, model, IntModel, AntiCheck, IOPT, WindArray, Magnetopause, Printtime,maxlat,maxlong,minlat,minlong,LatStep, LongStep, MaxStepPercent, EndParams, Rcomp, Rscan, LiveData, asymptotic, asymlevels, unit, serverdata, kp, custom_coords_provided=False):
      result = []

      pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

      if LiveData == 1:
         OnlineData = "Online Space Weather Data Used"
      elif serverdata == "ON":
         OnlineData = "Server Data Generated Using OMNI Database Used"
      else:
         OnlineData = "User Inputted Data Used"
      
      if Rcomp == "Vertical":
         CutoffComp = "Vertical Cutoff Rigidity"
      elif Rcomp == "Apparent":
         CutoffComp = "Apparent Cutoff Rigidity"
      else:
         CutoffComp = "Custom Cutoff"
      
      particle = "anti-particle" if AntiCheck == 1 else "Normal Particle"
      RigidityScan = "Rigidity Scan Used" if Rscan != 0 else "No Rigidity Scan"
      
      IntegrationMethods = ["4th Order Runge-Kutta Method", "Boris Method", "Vay Method", "Higuera-Cary Method"]
      IntegrationMethod = IntegrationMethods[IntModel-1] if 0 <= IntModel <= 4 else "Unknown Integration Method"
  
      InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
      Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"
  
      ExternalModels = [
          "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
          "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
      ]
      External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
      if model[1] == 99:
         External = "MHD file"
      
  
      PauseModels = [
          "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
          "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
          "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
      ]
      PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"
  
      today = date.today()
      result.append(f"\n")
      result.append(f"Date of OTSO computation: {today}\n")
      result.append(f"Total computation time: {Printtime} seconds\n\n")
      result.append(f"Cutoff Computed: {CutoffComp}\n\n")
      result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
      result.append(f"Rigidity Scan:\n{RigidityScan}\n\n")
      result.append(f"Input Variables:\n\n")
      result.append(f"Data Used: {OnlineData}\n\n")
      result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
      result.append(f"Max Time Step [% of gyrofrequency]: {MaxStepPercent}\n\n")
      result.append(f"Minimum Altitude: {EndParams[0]}km\n")
      result.append(f"Max Distance: {EndParams[1]}Re\n")
      result.append(f"Max Time: {EndParams[2]}s\n\n")
      
      if Data and len(Data[0]) > 5:
          result.append(f"Start Altitude = {Data[0][3]}km \n")
          result.append(f"Zenith = {Data[0][4]}\n")
          result.append(f"Azimuth = {Data[0][5]}\n\n")
      else:
          result.append(f"Start Altitude = Not Available\n")
          result.append(f"Zenith = Not Available\n")
          result.append(f"Azimuth = Not Available\n\n")

      result.append(f"Kp = {kp}\n")
      result.append(f"IOPT = {IOPT}\n\n")
      result.append(f"Solar Wind Speed [km/s]:\n")
      result.append(f"Vx = {abs(WindArray[0])}\n")
      result.append(f"Vy = {WindArray[1]}\n")
      result.append(f"Vz = {WindArray[2]}\n\n")
      result.append(f"IMF [nT]:\n")
      result.append(f"By = {WindArray[3]}\n")
      result.append(f"Bz = {WindArray[4]}\n\n")
      result.append(f"Density = {WindArray[5]} cm^-3\n")
      result.append(f"Pdyn = {pressure} nPa\n\n")
      result.append(f"Dst = {WindArray[7]} nT\n\n")
      result.append(f"G1 = {WindArray[8]}\n")
      result.append(f"G2 = {WindArray[9]}\n")
      result.append(f"G3 = {WindArray[10]}\n\n")
      result.append(f"W1 = {WindArray[11]}\n")
      result.append(f"W2 = {WindArray[12]}\n")
      result.append(f"W3 = {WindArray[13]}\n")
      result.append(f"W4 = {WindArray[14]}\n")
      result.append(f"W5 = {WindArray[15]}\n")
      result.append(f"W6 = {WindArray[16]}\n\n")
      result.append(f"Particle Type = {particle}\n\n")
      result.append(f"Magnetic Field Models:\n")
      result.append(f"Internal Model = {Internal}\n")
      result.append(f"External Model = {External}\n\n")
      result.append(f"Magnetopause Model = {PauseModel}\n\n")
      result.append(f"Rigidity:\n")
      result.append(f"Start = {Rigidity[0]} [GV]\n")
      result.append(f"End = {Rigidity[1]} [GV]\n")
      result.append(f"Step = {Rigidity[2]} [GV]\n\n")

      if custom_coords_provided:
          result.append("Coordinate Input: Custom list of (Latitude, Longitude) pairs provided.\n\n")
      else:
          result.append("Max and Min Latitude and Longitude:"+ "\n")
          result.append("Latitude: Max = " + str(maxlat) + " Min = " + str(minlat) + "\n")
          result.append("Longitude: Max = " + str(maxlong) + " Min = " + str(minlong) + "\n\n")
          result.append("Latitude and Longitude Steps:"+ "\n")
          lat_step_str = str(abs(LatStep)) if LatStep is not None else "N/A"
          long_step_str = str(abs(LongStep)) if LongStep is not None else "N/A"
          result.append("Latitude = " + lat_step_str + " degree steps" + "\n")
          result.append("Longitude = " + long_step_str + " degree steps" + "\n\n")
      
      if asymptotic == "YES":
         result.append("Asymptotic Directions: " + str(asymptotic) + " \n")
         result.append("Asymptotic Levels: " + str(asymlevels) + " \n")
         result.append("Asymptotic Levels Unit: " + str(unit) + " \n")
    
      return "".join(result)

def READMEMagfield(EventDate, model, IOPT, WindArray, CoordinateSystem, Printtime, LiveData, serverdata,kp):
    result = []

    pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

    if LiveData == 1:
      OnlineData = "Online Space Weather Data Used"
    elif serverdata == "ON":
       OnlineData = "Server Data Generated Using OMNI Database Used"
    else:
       OnlineData = "User Inputted Data Used"

    InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
    Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"

    ExternalModels = [
        "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
        "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
    ]
    External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
    if model[1] == 99:
         External = "MHD file"

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Input Coordinate System:\n{CoordinateSystem}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Kp = {kp}\n")
    result.append(f"IOPT = {IOPT}\n\n")
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(f"Vx = {abs(WindArray[0])}\n")
    result.append(f"Vy = {WindArray[1]}\n")
    result.append(f"Vz = {WindArray[2]}\n\n")
    result.append(f"IMF [nT]:\n")
    result.append(f"By = {WindArray[3]}\n")
    result.append(f"Bz = {WindArray[4]}\n\n")
    result.append(f"Density = {WindArray[5]} cm^-3\n")
    result.append(f"Pdyn = {pressure} nPa\n\n")
    result.append(f"Dst = {WindArray[7]} nT\n\n")
    result.append(f"G1 = {WindArray[8]}\n")
    result.append(f"G2 = {WindArray[9]}\n")
    result.append(f"G3 = {WindArray[10]}\n\n")
    result.append(f"W1 = {WindArray[11]}\n")
    result.append(f"W2 = {WindArray[12]}\n")
    result.append(f"W3 = {WindArray[13]}\n")
    result.append(f"W4 = {WindArray[14]}\n")
    result.append(f"W5 = {WindArray[15]}\n")
    result.append(f"W6 = {WindArray[16]}\n\n")
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"External Model = {External}\n\n")

    return "".join(result)

def READMECoordtrans(CoordIN, CoordOUT, Printtime):
    result = []

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Input Coordinate System:\n{CoordIN}\n\n")
    result.append(f"Output Coordinate System:\n{CoordOUT}\n\n")
    

    return "".join(result)

def READMETrace(Data, EventDate, model, IntModel, AntiCheck, IOPT, WindArray, Magnetopause, Printtime,maxlat,maxlong,minlat,minlong,LatStep, LongStep, MaxStepPercent, EndParams, LiveData, serverdata, kp):
      result = []

      pressure = misc.Pdyn_comp(WindArray[5],abs(WindArray[0]))

      if LiveData == 1:
         OnlineData = "Online Space Weather Data Used"
      elif serverdata == "ON":
         OnlineData = "Server Data Generated Using OMNI Database Used"
      else:
         OnlineData = "User Inputted Data Used"
  
      InternalModels = ["None", "IGRF", "Dipole", "Custom Gaussian Coefficients", "Custom Gaussian Coefficients"]
      Internal = InternalModels[model[0]] if 0 <= model[0] <= 4 else "Unknown Internal Model"
  
      ExternalModels = [
          "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89",
          "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89 Boberg Extension"
      ]
      External = ExternalModels[model[1]] if 0 <= model[1] <= 8 else "Unknown External Model"
      if model[1] == 99:
         External = "MHD file"
  
      PauseModels = [
          "25Re Sphere", "Aberrated Formisano Model", "Sibeck Model", "Kobel Model",
          "Tsyganenko 96 Magnetopause Model", "Tsyganenko 01 Magnetopause Model",
          "Tsyganenko 01 Storm Magnetopause Model", "Tsyganenko 04 Magnetopause Model"
      ]
      PauseModel = PauseModels[Magnetopause] if 0 <= Magnetopause <= 7 else "Unknown Magnetopause Model"
  
      today = date.today()
      result.append(f"\n")
      result.append(f"Date of OTSO computation: {today}\n")
      result.append(f"Total computation time: {Printtime} seconds\n\n")
      result.append(f"Input Variables:\n\n")
      result.append(f"Data Used: {OnlineData}\n\n")
      result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
      result.append(f"Kp = {kp}\n")
      result.append(f"IOPT = {IOPT}\n\n")
      result.append(f"Solar Wind Speed [km/s]:\n")
      result.append(f"Vx = {abs(WindArray[0])}\n")
      result.append(f"Vy = {WindArray[1]}\n")
      result.append(f"Vz = {WindArray[2]}\n\n")
      result.append(f"IMF [nT]:\n")
      result.append(f"By = {WindArray[3]}\n")
      result.append(f"Bz = {WindArray[4]}\n\n")
      result.append(f"Density = {WindArray[5]} cm^-3\n")
      result.append(f"Pdyn = {pressure} nPa\n\n")
      result.append(f"Dst = {WindArray[7]} nT\n\n")
      result.append(f"G1 = {WindArray[8]}\n")
      result.append(f"G2 = {WindArray[9]}\n")
      result.append(f"G3 = {WindArray[10]}\n\n")
      result.append(f"W1 = {WindArray[11]}\n")
      result.append(f"W2 = {WindArray[12]}\n")
      result.append(f"W3 = {WindArray[13]}\n")
      result.append(f"W4 = {WindArray[14]}\n")
      result.append(f"W5 = {WindArray[15]}\n")
      result.append(f"W6 = {WindArray[16]}\n\n")
      result.append(f"Magnetic Field Models:\n")
      result.append(f"Internal Model = {Internal}\n")
      result.append(f"External Model = {External}\n\n")
      result.append(f"Magnetopause Model = {PauseModel}\n\n")
      result.append("Max and Min Latitude and Longitude:"+ "\n")
      result.append("Latitude: Max = " + str(maxlat) + " Min = " + str(minlat) + "\n")
      result.append("Longitude: Max = " + str(maxlong) + " Min = " + str(minlong) + "\n\n")
      result.append("Latitude and Longitude Steps:"+ "\n")
      result.append("Latitude = " + str(abs(LatStep)) + " degree steps" + "\n")
      result.append("Longitude = " + str(abs(LongStep)) + " degree steps" + "\n\n")
    
      return "".join(result)