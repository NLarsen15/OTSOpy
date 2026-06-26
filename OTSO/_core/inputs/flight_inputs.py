import numpy as np

from ..utils import input_utils, tsy_params_utils
from ..custom_classes import date as d
from ..custom_classes import solar_wind
from ..livedata import pull_live_data
from ..serverdata import server
from ..data_classes.flight_data import FlightData

def FlightInputs(Data: FlightData) -> None:
    
     DateArrayList = []
     for x in Data.dates:
          DateCreate = d.Date(x)
          DateArray = DateCreate.GetDate()
          DateArrayList.append(DateArray)

     Data.datearraylist = DateArrayList

     variablelist = [Data.vx, Data.vy, Data.vz, Data.by, Data.bz, Data.density, Data.pdyn, Data.Dst,
           Data.G1, Data.G2, Data.G3, Data.W1, Data.W2, Data.W3, Data.W4, Data.W5, Data.W6, Data.kp]
     variablelist2 = [Data.latitudes, Data.longitudes, Data.dates, Data.vx, Data.vy, Data.vz, Data.by, Data.bz, Data.density, Data.pdyn, Data.Dst,
           Data.G1, Data.G2, Data.G3, Data.W1, Data.W2, Data.W3, Data.W4, Data.W5, Data.W6, Data.kp]
     
     input_utils.asymptotic_check(Data.asymptotic, Data.unit)

     AntiCheck = input_utils.anti_check(Data.anti)

     Data.magnetopause = input_utils.magnetopause_check(Data.mpause)

     Data.integrationmodel = input_utils.intmodel_check(Data.intmodel)

     Data.Rcomp, Data.zenith, Data.azimuth = input_utils.cutoff_comp_check(Data.cutoff_comp, Data.zenith, Data.azimuth)

     ServerData = input_utils.serverdata_check(Data.serverdata)

     LiveData = input_utils.livedata_check(Data.livedata)

     tempGlist = []
     tempHList = []

     for date in Data.datearraylist:
          Internal, g, h = input_utils.internalmag_check(Data.internalmag, date, Data.g, Data.h)
          tempGlist.append(g)
          tempHList.append(h)
     Data.glist = tempGlist
     Data.hlist = tempHList

     External = input_utils.externalmag_check(Data.externalmag, Data.MHDfile)

     Bobon, bobtype = input_utils.BobergCheck(Data.boberg, Data.bobergtype)
 
     input_utils.coordsystem_check(Data.coordsystem, Data.inputcoord)
     
     Data.Rscan = input_utils.rigidityscan_check(Data.rigidityscan)
    
     input_utils.flight_server_live_check(variablelist, variablelist2, Data.serverdata, Data.livedata)
 
     Data.maxsteppercent = Data.gyropercent/100

     KpList = []
 
     WindArrayList = []
     IOPTList = []
     i = 0

     for x, g_coeffs, h_coeffs in zip(Data.dates, Data.glist, Data.hlist):
          
        if ServerData == 1:
           if int(x.year) >= 1981:
                server.DownloadServerFile(int(x.year), g_coeffs, h_coeffs)
           elif int(x.year) < 1981 and int(x.year) > 1963:
                server.DownloadServerFileLowRes(int(x.year))
           else:
                print("Server data only valid for 1963 to present, please enter a valid date.")
           BxS, ByS, BzS, VS, DensityS, PdynS, KpS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS, External = server.GetServerData(x,External,Data.AdaptiveExternalModel)
           KpList.append(KpS)
           IOPTinput = tsy_params_utils.IOPTprocess(KpS)
           if External == 100:
               IOPTinput = tsy_params_utils.IOPTprocess_refit(KpS)
           IOPTList.append(IOPTinput)
           vytemp = 0
           vztemp = 0
           WindCreate = solar_wind.SolarWind(VS, vytemp, vztemp, BxS, ByS, BzS, DensityS, PdynS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, KpS, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS)
           WindArray = WindCreate.GetWind()
           WindArrayList.append(WindArray)
           
           if LiveData == 1:
                if External == 7 or External == 11:
                    print("LIVE DATA NOT SUPPORTED FOR TSY04 OR TA16 MAGNETOSPHERIC MODELS. PLEASE SELECT ANOTHER EXTERNAL MAGNETIC FIELD MODEL.")
                    exit()
                input_utils.DateCheck(x)
                DstLive, VxLive, DensityLive, ByLive, BzLive, IOPTLive, G1Live, G2Live, G3Live, KpLive, Bx_avgLive, By_avgLive, Bz_avgLive, NIndexLive, BIndexLive = pull_live_data.Get_Data(x)
                PdynLive = tsy_params_utils.Pdyn_comp(DensityLive,VxLive)
                KpList.append(KpLive)
                IOPTinput = IOPTLive
                if External == 100:
                    IOPTinput = tsy_params_utils.IOPTprocess_refit(KpLive)
                IOPTList.append(IOPTinput)
                vytemp = 0
                vztemp = 0
                W1 = W2 = W3 = W4 = W5 = W6 = sym_h_corrected = 0
                WindCreate = solar_wind.SolarWind(VxLive, vytemp, vztemp, Bx_avgLive, ByLive, BzLive, DensityLive, PdynLive, DstLive, G1Live, G2Live, G3Live, W1, W2, W3, W4, W5, W6, KpLive, By_avgLive, Bz_avgLive, NIndexLive, BIndexLive, sym_h_corrected)
                WindArray = WindCreate.GetWind()
                WindArrayList.append(WindArray)
 
        if ServerData == 0 and LiveData == 0:
           if Data.vx[i] > 0:
                Data.vx[i] = -1*Data.vx[i]
           WindCreate = solar_wind.SolarWind(Data.vx[i], Data.vy[i], Data.vz[i], Data.bx[i], Data.by[i], Data.bz[i], Data.density[i], Data.pdyn[i], Data.Dst[i], Data.G1[i], Data.G2[i],
                                              Data.G3[i], Data.W1[i], Data.W2[i], Data.W3[i], Data.W4[i], Data.W5[i], Data.W6[i], Data.kp[i], Data.by_avg[i], Data.bz_avg[i], Data.n_index[i], Data.b_index[i], Data.sym_h_corrected[i])
           WindArray = WindCreate.GetWind()
           KpList.append(Data.kp[i])
           IOPTinput = tsy_params_utils.IOPTprocess(Data.kp[i])
           if External == 100:
               IOPTinput = tsy_params_utils.IOPTprocess_refit(Data.kp[i])
           IOPTList.append(IOPTinput)
           WindArrayList.append(WindArray)
           i += 1
          
     Data.windarraylist = WindArrayList
     Data.IOPTlist = IOPTList
     Data.Kplist = KpList
 
     Data.rigidityarray = [Data.startrigidity,Data.endrigidity,Data.rigiditystep]
 
     Data.model = np.array([Internal,External,Bobon,bobtype])
 
     Data.endparams = [Data.minaltitude,Data.maxdistance,Data.maxtime,float(Data.maxsteps)]
     
     stationslist = []
     for lat,long,alt in zip(Data.latitudes,Data.longitudes,Data.altitudes):
          station = ["temp",lat,long,alt,Data.zenith,Data.azimuth]
          stationslist.append(station)
          input_utils.ParamCheck(alt,Data.endparams)
     
     Data.station_array = stationslist
          
     Data.particlearray = [Data.Anum,AntiCheck]
 
     return