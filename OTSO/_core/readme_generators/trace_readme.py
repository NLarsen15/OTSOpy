from datetime import date

from .readme_utils import *
from ..data_classes.trace_data import TraceData

def READMETrace(Data: TraceData, EventDate: date, Printtime: float) -> str:
  
  result = []

  OnlineData = OnlineDataStatus(Data.livedata, Data.serverdata)

  Internal = InternalModelCheck(Data.model)

  External = ExternalModelCheck(Data.model)

  PauseModel = MagnetopauseModelCheck(Data.magnetopause)

  today = date.today()
  result.append(f"\n")
  result.append(f"Date of OTSO computation: {today}\n")
  result.append(f"Total computation time: {Printtime} seconds\n\n")
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
  result = field_models_readme_section(result, Internal, External, PauseModel)
  result = boberg_readme_section(result, Data.boberg, Data.bobergtype)
  result.append("Max and Min Latitude and Longitude:"+ "\n")
  result.append("Latitude: Max = " + str(Data.maxlat) + " Min = " + str(Data.minlat) + "\n")
  result.append("Longitude: Max = " + str(Data.maxlong) + " Min = " + str(Data.minlong) + "\n\n")
  result.append("Latitude and Longitude Steps:"+ "\n")
  result.append("Latitude = " + str(abs(Data.latstep)) + " degree steps" + "\n")
  result.append("Longitude = " + str(abs(Data.longstep)) + " degree steps" + "\n\n")

  return "".join(result)