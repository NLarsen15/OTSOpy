from datetime import date

from .readme_utils import *
from ..data_classes.planet_data import PlanetData

def READMEPlanet(Data: 'PlanetData', EventDate, Printtime, custom_coords_provided=False):
    
    result = []

    OnlineData = OnlineDataStatus(Data.livedata, Data.serverdata)
    
    particle = ParticleCheck(Data.particlearray[1])

    IntegrationMethod = IntegrationMethodCheck(Data.integrationmodel)

    Internal = InternalModelCheck(Data.model)

    External = ExternalModelCheck(Data.model)
    
    PauseModel = MagnetopauseModelCheck(Data.magnetopause)

    CutoffComp = CutoffCompCheck(Data.Rcomp)

    RigidityScan = RigidityScanCheck(Data.Rscan)

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Cutoff Computed: {CutoffComp}\n\n")
    result.append(f"Integration Method:\n{IntegrationMethod}\n\n")
    result.append(f"Rigidity Scan:\n{RigidityScan}\n\n")
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
    
    result.append(f"Start Altitude = {Data.startaltitude}km \n")
    if Data.Rcomp == "Apparent":
        result.append(f"\n\n")
    else:
        result.append(f"Zenith = {Data.zenith}\n")
        result.append(f"Azimuth = {Data.azimuth}\n\n")

    result.append(f"Kp = {Data.Kp}\n")
    result.append(f"IOPT = {Data.IOPT}\n\n")
    result = solar_wind_readme_section(result, Data.windarray)
    result.append(f"Particle Type = {particle}\n\n")
    result = field_models_readme_section(result, Internal, External, PauseModel)
    result = boberg_readme_section(result, Data.boberg, Data.bobergtype)
    result.append(f"Rigidity:\n")
    result.append(f"Start = {Data.rigidityarray[0]} [GV]\n")
    result.append(f"End = {Data.rigidityarray[1]} [GV]\n")
    result.append(f"Step = {Data.rigidityarray[2]} [GV]\n\n")

    if custom_coords_provided:
        result.append("Coordinate Input: Custom list of (Latitude, Longitude) pairs provided.\n\n")
    else:
        result.append("Max and Min Latitude and Longitude:"+ "\n")
        result.append("Latitude: Max = " + str(Data.maxlat) + " Min = " + str(Data.minlat) + "\n")
        result.append("Longitude: Max = " + str(Data.maxlong) + " Min = " + str(Data.minlong) + "\n\n")
        result.append("Latitude and Longitude Steps:"+ "\n")
        lat_step_str = str(abs(Data.latstep)) if Data.latstep is not None else "N/A"
        long_step_str = str(abs(Data.longstep)) if Data.longstep is not None else "N/A"
        result.append("Latitude = " + lat_step_str + " degree steps" + "\n")
        result.append("Longitude = " + long_step_str + " degree steps" + "\n\n")
    
    if Data.asymptotic == "YES":
       result.append("Asymptotic Directions: " + str(Data.asymptotic) + " \n")
       result.append("Asymptotic Levels: " + str(Data.asymlevels) + " \n")
       result.append("Asymptotic Levels Unit: " + str(Data.unit) + " \n")
  
    return "".join(result)