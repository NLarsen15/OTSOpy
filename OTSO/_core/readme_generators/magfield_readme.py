from datetime import date

from .readme_utils import *
from ..data_classes.magfield_data import MagfieldData

def READMEMagfield(EventDate: date, Data: MagfieldData, Printtime: float) -> str:
    
    result = []

    OnlineData = OnlineDataStatus(Data.livedata, Data.serverdata)

    Internal = InternalModelCheck(Data.internalmag)

    External = ExternalModelCheck(Data.model)

    today = date.today()
    result.append(f"\n")
    result.append(f"OTSO Version: {OTSOVersion()}\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Input Coordinate System:\n{Data.inputcoord}\n\n")
    result.append(f"Output Coordinate System:\n{Data.coordout}\n\n")
    result.append(f"Input Variables:\n\n")
    result.append(f"Data Used: {OnlineData}\n")
    if OnlineData == "Online Space Weather Data Used":
      result.append("NOAA preliminary data used, these values are not final and may differ from the finalised OMNI database used in the ServerData function. The OMNI database results take precident over the preliminary NOAA values used in LiveData.\n\n")
    else:
      result.append("\n")
    result.append(f"Simulation Date: {EventDate.strftime('%d/%m/%Y, %H:%M:%S')}\n\n")
    result.append(f"Kp = {Data.Kp}\n")
    result.append(f"IOPT = {Data.IOPT}\n\n")
    result = solar_wind_readme_section(result, Data.windarray)
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"Max Spherical Harmonics Order/Degree = {Data.max_degree}\n")
    result.append(f"External Model = {External}\n\n")
    result = boberg_readme_section(result, Data.boberg, Data.bobergtype)

    result.append(Data.citationstring)
    result.append("\n")

    return "".join(result)