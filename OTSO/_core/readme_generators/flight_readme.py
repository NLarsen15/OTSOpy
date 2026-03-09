from datetime import date
import pandas as pd

from .readme_utils import *
from ..custom_classes import date as d
from ..data_classes.flight_data import FlightData

def READMEFlight(Data: FlightData, Printtime: float) -> str:
   
   result = []

   OnlineData = OnlineDataStatus(Data.livedata, Data.serverdata)

   CutoffComp = CutoffCompCheck(Data.Rcomp)

   particle = ParticleCheck(Data.particlearray[1])
   
   RigidityScan = RigidityScanCheck(Data.Rscan)
   
   IntegrationMethod = IntegrationMethodCheck(Data.integrationmodel)

   Internal = InternalModelCheck(Data.model)

   External = ExternalModelCheck(Data.model)

   PauseModel = MagnetopauseModelCheck(Data.magnetopause)

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
   result.append(f"Max Time Step [% of gyrofrequency]: {Data.maxsteppercent*100}\n")
   result = beta_readme_section(result, Data.totalbetacheck, Data.Berr, Data.adapt)
   result = end_conditions_readme_section(result, Data.endparams)
   result.append(f"Zenith = {Data.zenith}\n")
   result.append(f"Azimuth = {Data.azimuth}\n\n")
   result.append(f"Particle Type = {particle}\n\n")
   result = field_models_readme_section(result, Internal, External, PauseModel)
   result = boberg_readme_section(result, Data.boberg, Data.bobergtype)
   result.append(f"Rigidity:\n")
   result.append(f"Start = {Data.rigidityarray[0]} [GV]\n")
   result.append(f"End = {Data.rigidityarray[1]} [GV]\n")
   result.append(f"Step = {Data.rigidityarray[2]} [GV]\n\n")
   if Data.asymptotic == "YES":
      result.append("Asymptotic Directions: " + str(Data.asymptotic) + " \n")
      result.append("Asymptotic Levels: " + str(Data.asymlevels) + " \n")
      result.append("Asymptotic Levels Unit: " + str(Data.unit) + " \n")
 
   return "".join(result)

def READMEFlightData(DateArray,WindArray,kparray):

   headers = ["Date","kp","Vx [km/s]", "Vy [km/s]", "Vz [km/s]", "Bx [nT]", "By [nT]", "Bz [nT]", "By Average (30mins) [nT]", "Bz Average (30mins)[nT]", "Density [cm^-3]", "Pdyn [nPa]", "Dst [nT]", "G1", "G2",
              "G3", "W1", "W2", "W3", "W4", "W5", "W6", "N", "B", "SYM-H"]
   
   rows = []

   for date,Wind,kp in zip(DateArray,WindArray,kparray):
      datetimeobj = d.convert_to_datetime(date)
      Vx = abs(Wind[0])
      Vy = Wind[1]
      Vz = Wind[2]
      Bx = Wind[3]
      By = Wind[4]
      Bz = Wind[5]
      By_av = Wind[20]
      Bz_av = Wind[21]
      Density = Wind[6]
      Pdyn = Wind[7]
      Dst = Wind[8]
      G1 = Wind[9]
      G2 = Wind[10]
      G3 = Wind[11]
      W1 = Wind[12]
      W2 = Wind[13]
      W3 = Wind[14]
      W4 = Wind[15]
      W5 = Wind[16]
      W6 = Wind[17]
      n_index = Wind[22]
      b_index = Wind[23]
      sym_h_corrected = Wind[24]
      data = {
      "Date": datetimeobj,
      "kp": kp,
      "Vx [km/s]": Vx,
      "Vy [km/s]": Vy,
      "Vz [km/s]": Vz,
      "Bx [nT]": Bx,
      "By [nT]": By,
      "Bz [nT]": Bz,
      "By Average (30mins) [nT]": By_av,
      "Bz Average (30mins) [nT]": Bz_av,
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
      "W6": W6,
      "N": n_index,
      "B": b_index,
      "SYM-H": sym_h_corrected,}
      
      rows.append(data)

   df = pd.DataFrame(rows, columns=headers)

   return df