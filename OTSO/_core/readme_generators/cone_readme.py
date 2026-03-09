from datetime import date

from .readme_utils import *
from ..data_classes.cone_data import ConeData

def READMECone(Data: ConeData, EventDate, Printtime) -> str:
    
    result = []

    OnlineData = OnlineDataStatus(Data.livedata, Data.serverdata)
    
    particle = ParticleCheck(Data.particlearray[1])
    
    IntegrationMethod = IntegrationMethodCheck(Data.integrationmodel)

    Internal = InternalModelCheck(Data.model)

    External = ExternalModelCheck(Data.model)

    PauseModel = MagnetopauseModelCheck(Data.magnetopause)

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Output Coordinate System:\n{Data.coordsystem}\n\n")
    result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n")
    if OnlineData == "Online Space Weather Data Used":
      result.append("NOAA preliminary data used, these values are not final and may differ from the finalised OMNI database used in the ServerData function. The OMNI database results take precident over the preliminary NOAA values used in LiveData.\n\n")
    else:
      result.append("\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Max Time Step [% of gyrofrequency]: {Data.maxsteppercent*100}\n")
    result = beta_readme_section(result, Data.totalbetacheck, Data.Berr, Data.adapt)
    result = end_conditions_readme_section(result, Data.endparams)
    result.append(f"Start Altitude = {Data.station_array[0][3]}km \n")
    result.append(f"Zenith = {Data.station_array[0][4]}\n")
    result.append(f"Azimuth = {Data.station_array[0][5]}\n\n")
    result.append(f"Kp = {Data.Kp}\n")
    result.append(f"IOPT = {Data.IOPT}\n\n")
    result = solar_wind_readme_section(result, Data.windarray)
    result.append(f"Atomic Number = {Data.particlearray[0]}\n\n")
    result.append(f"Particle Type = {particle}\n\n")
    result = field_models_readme_section(result, Internal, External, PauseModel)
    result = boberg_readme_section(result, Data.boberg, Data.bobergtype)
    result.append(f"Rigidity:\n")
    result.append(f"Start = {Data.rigidityarray[0]} [GV]\n")
    result.append(f"End = {Data.rigidityarray[1]} [GV]\n")
    result.append(f"Step = {Data.rigidityarray[2]} [GV]\n\n")
    result.append(f"Stations:\n")

    for station in Data.station_array:
        result.append(f"{station[0]}, Latitude: {station[1]}, Longitude: {station[2]}\n")

    return "".join(result)