from datetime import date

from .readme_utils import *
from ..data_classes.magfield_data import MagfieldData

def READMEMagfield(EventDate: date, MagfieldDataInstance: MagfieldData, Printtime: float) -> str:
    
    result = []

    OnlineData = OnlineDataStatus(MagfieldDataInstance.livedata, MagfieldDataInstance.serverdata)

    Internal = InternalModelCheck(MagfieldDataInstance.model)

    External = ExternalModelCheck(MagfieldDataInstance.model)

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Input Coordinate System:\n{MagfieldDataInstance.inputcoord}\n\n")
    result.append(f"Output Coordinate System:\n{MagfieldDataInstance.coordout}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n")
    if OnlineData == "Online Space Weather Data Used":
      result.append("NOAA preliminary data used, these values are not final and may differ from the finalised OMNI database used in the ServerData function. The OMNI database results take precident over the preliminary NOAA values used in LiveData.\n\n")
    else:
      result.append("\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Kp = {MagfieldDataInstance.Kp}\n")
    result.append(f"IOPT = {MagfieldDataInstance.IOPT}\n\n")
    result = solar_wind_readme_section(result, MagfieldDataInstance.windarray)
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"External Model = {External}\n\n")
    result = boberg_readme_section(result, MagfieldDataInstance.boberg, MagfieldDataInstance.bobergtype)

    return "".join(result)