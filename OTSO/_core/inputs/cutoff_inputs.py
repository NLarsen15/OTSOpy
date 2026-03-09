import numpy as np
from datetime import datetime

from ..utils import input_utils, tsy_params_utils
from ..custom_classes import date, solar_wind, stations
from ..livedata import pull_live_data
from ..serverdata import server
from ..data_classes.cutoff_data import CutoffData

def CutoffInputs(Data: 'CutoffData') -> None:
    
     EventDate = datetime(Data.year,Data.month,Data.day,Data.hour,Data.minute,Data.second)
     DateCreate = date.Date(EventDate)
     Data.datearray = DateCreate.GetDate()

     AntiCheck = input_utils.anti_check(Data.anti)
    
     Data.magnetopause = input_utils.magnetopause_check(Data.mpause)

     Data.integrationmodel = input_utils.intmodel_check(Data.intmodel)

     ServerData = input_utils.serverdata_check(Data.serverdata)

     LiveData = input_utils.livedata_check(Data.livedata)
    
     Internal, Data.g, Data.h = input_utils.internalmag_check(Data.internalmag, Data.datearray, Data.g, Data.h)
      
     External = input_utils.externalmag_check(Data.externalmag, Data.MHDfile)

     Bobon, bobtype = input_utils.BobergCheck(Data.boberg, Data.bobergtype)

     input_utils.coordsystem_check(Data.coordsystem, Data.inputcoord)
 
     input_utils.DataCheck(ServerData,LiveData,EventDate)

     input_utils.CustomLocationsCheck(Data.customlocations)

     Data.Rcomp, Data.zenith, Data.azimuth = input_utils.cutoff_comp_check(Data.cutoff_comp, Data.zenith, Data.azimuth)
     
     Data.Rscan = input_utils.rigidityscan_check(Data.rigidityscan)

     input_utils.DataCheck(ServerData,LiveData,EventDate)

     Data.endparams = [Data.minaltitude,Data.maxdistance,Data.maxtime,float(Data.maxsteps)]

     input_utils.ParamCheck(Data.startaltitude,Data.endparams)
 
     Data.IOPT = tsy_params_utils.IOPTprocess(Data.kp)

     Data.maxsteppercent = Data.gyropercent/100
 
     if ServerData == 1:
          if int(EventDate.year) >= 1981:
               server.DownloadServerFile(int(EventDate.year),Data.g,Data.h)
          elif int(EventDate.year) < 1981 and int(EventDate.year) > 1963:
               server.DownloadServerFileLowRes(int(EventDate.year))
          else:
               print("Server data only valid for 1963 to present, please enter a valid date.")
          BxS, ByS, BzS, VS, DensityS, PdynS, KpS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, ByAvgS, BzAvgS, NIndexS, BIndexS, SymHCorrectedS, External = server.GetServerData(EventDate,External,Data.AdaptiveExternalModel)
          Data.IOPT = tsy_params_utils.IOPTprocess(KpS)
          if External == 100:
               Data.IOPT = tsy_params_utils.IOPTprocess_refit(KpS)
          Data.Kp = KpS
          WindCreate = solar_wind.SolarWind(VS, Data.vy, Data.vz, BxS, ByS, BzS, DensityS, PdynS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, KpS, ByAvgS, BzAvgS, NIndexS, BIndexS, SymHCorrectedS)
          Data.windarray = WindCreate.GetWind()
          
     if LiveData == 1:
          if External == 7 or External == 11:
               print("LIVE DATA NOT SUPPORTED FOR TSY04 OR TA16 MAGNETOSPHERIC MODELS. PLEASE SELECT ANOTHER EXTERNAL MAGNETIC FIELD MODEL.")
               exit()
          input_utils.DateCheck(EventDate)
          DstLive, VxLive, DensityLive, ByLive, BzLive, IOPTLive, G1Live, G2Live, G3Live, KpLive, By_avgLive, Bx_avgLive, Bz_avgLive, NIndexLive, BIndexLive = pull_live_data.Get_Data(EventDate)
          PdynLive = tsy_params_utils.Pdyn_comp(DensityLive,VxLive)
          Data.IOPT = IOPTLive
          if External == 100:
               Data.IOPT = tsy_params_utils.IOPTprocess_refit(KpLive)
          Data.Kp = KpLive
          W1 = W2 = W3 = W4 = W5 = W6 = sym_h_corrected = 0
          WindCreate = solar_wind.SolarWind(VxLive, Data.vy, Data.vz, Bx_avgLive, ByLive, BzLive, DensityLive, PdynLive, DstLive, G1Live, G2Live, G3Live, W1, W2, W3, W4, W5, W6, KpLive, By_avgLive, Bz_avgLive, NIndexLive, BIndexLive, sym_h_corrected)
          Data.windarray = WindCreate.GetWind()
 
     if ServerData == 0 and LiveData == 0:
          if Data.vx > 0:
                Data.vx = -1*Data.vx
          Data.Kp = Data.kp
          if External == 100:
               Data.IOPT = tsy_params_utils.IOPTprocess_refit(Data.kp)
          WindCreate = solar_wind.SolarWind(Data.vx, Data.vy, Data.vz, Data.bx, Data.by, Data.bz, Data.density, Data.pdyn, Data.Dst, Data.G1, Data.G2, Data.G3, Data.W1, Data.W2, Data.W3, Data.W4, Data.W5, Data.W6, Data.kp, Data.by_avg, Data.bz_avg, Data.n_index, Data.b_index, Data.sym_h_corrected)
          Data.windarray = WindCreate.GetWind()
 
     Data.rigidityarray = [Data.startrigidity,Data.endrigidity,Data.rigiditystep]
 
     Data.model = np.array([Internal,External,Bobon,bobtype])
 
     CreateStations = stations.Stations(Data.Stations, Data.startaltitude, Data.zenith, Data.azimuth)
     InputtedStations = CreateStations
     InputtedStations.find_non_matching_stations()
     if hasattr(Data, 'customlocations') and Data.customlocations:
           CreateStations.AddLocation(Data.customlocations)
     Used_Stations_Temp = CreateStations.GetStations()
     temp = list(Used_Stations_Temp)
     Data.station_array = temp
 
     Data.particlearray = [Data.Anum,AntiCheck]
 
     return