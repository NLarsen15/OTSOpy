import numpy as np
from datetime import datetime,timedelta
import os
from . import date, solar_wind, stations
from . import misc, Request, Server
from .igrf_process import compute_gauss_coefficients

def TrajectoryInputs(Stations,rigidity,customlocations,startaltitude,minaltitude,zenith,azimuth,
           maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile, MHDcoordsys,inputcoord):
    
    EventDate = datetime(year,month,day,hour,minute,second)
    DateCreate = date.Date(EventDate)
    DateArray = DateCreate.GetDate()

    if anti == "YES":
         AntiCheck = 1
    elif anti == "NO":
         AntiCheck = 0
    else:
         print("Please enter a valid anti value: ""YES"" or ""NO"" ")
         exit()

    if magnetopause == "Sphere":
         Magnetopause = 0
    elif magnetopause == "aFormisano":
         Magnetopause = 1
    elif magnetopause == "Sibeck":
         Magnetopause = 2
    elif magnetopause == "Kobel":
         Magnetopause = 3
    elif magnetopause == "Lin":
         Magnetopause = 4
    elif magnetopause == "NONE":
         Magnetopause = 99
    else:
         print("Please enter a valid magnetopause model: ""Sphere"", ""aFormisano"", ""Sibeck"", ""Kobel"", ""Lin"", ""NONE"" ")


    if intmodel == "4RK":
         IntModel = 1
    elif intmodel == "Boris":
         IntModel = 2
    elif intmodel == "Vay":
         IntModel = 3
    elif intmodel == "HC":
         IntModel = 4
    else:
         print("Please enter a valid intmodel model: ""4RK"", ""Boris"", ""Vay"", ""HC"" ")
         exit()

    if serverdata == "ON":
         ServerData = 1
    elif serverdata == "OFF":
         ServerData = 0
    else:
         print("Please enter a valid serverdata value: ""ON"" or ""OFF"" ")
         exit()

    if livedata == "ON":
         LiveData = 1
    elif livedata == "OFF":
         LiveData = 0
    else:
         print("Please enter a valid livedata value: ""ON"" or ""OFF"" ")
         exit()
    
    if internalmag == "NONE":
         Internal = 0
         if g is None or h is None or len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
    elif internalmag == "IGRF":
         Internal = 1
         coeffs = compute_gauss_coefficients(DateArray)
         g = coeffs['G_coefficients']
         h = coeffs['H_coefficients']
         if len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
    elif internalmag == "Dipole":
         Internal = 2
         coeffs = compute_gauss_coefficients(DateArray)
         g = coeffs['G_coefficients']
         h = coeffs['H_coefficients']
         if len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
    elif internalmag == "Custom Gauss":
         Internal = 4
         if g is None or h is None:
              print("Please enter values for the g and h Gaussian coefficients to use the Custom Gauss option")
              exit()
         elif len(g) != 136:
              print(f"There should be 136 g coefficents in the inputted list, you have entered {len(g)}")
              exit()
         elif len(h) != 136:
              print(f"There should be 136 h coefficents in the inputted list, you have enetered {len(h)}")
    else:
         print("Please enter a valid internalmag model: ""NONE"",""IGRF"",""Dipole"", or ""Custom Gauss""")
         exit()
      
    if externalmag == "NONE":
         External = 0
    elif externalmag == "TSY87short":
         External = 1
    elif externalmag == "TSY87long":
         External = 2
    elif externalmag == "TSY89":
         External = 3
    elif externalmag == "TSY96":
         External = 4
    elif externalmag == "TSY01":
         External = 5
    elif externalmag == "TSY01S":
         External = 6
    elif externalmag == "TSY04":
         External = 7
    elif externalmag == "TSY89_BOBERG":
         External = 8
    elif externalmag == "TSY15N":
         External = 9
    elif externalmag == "TSY15B":
         External = 10
    elif externalmag == "TA16_RBF":
         External = 11
    elif externalmag == "MHD":
         External = 99
         if not os.path.exists(MHDfile):
            print(f"The file '{MHDfile}' does not exist.")
            exit()
    else:
         print("Please enter a valid externalmag model: ""NONE"", ""TSY87short"",""TSy87long"",""TSY89"",""TSY89_BOBERG"",""TSY96"",""TSY01"",""TSY01S"",""TSY04"",""TSY15N"",""TSY15B"",""TA16_RBF""")
         exit()
         
    if coordsystem not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid coordsystem: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    if inputcoord not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid inputcoord: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

    misc.DataCheck(ServerData,LiveData,EventDate)

    Azimuth = azimuth
    Zenith = zenith

    IOPTinput = misc.IOPTprocess(kp)
    KpS = 0

    if ServerData == 1:
         if int(EventDate.year) >= 1981:
              Server.DownloadServerFile(int(EventDate.year),g,h)
         elif int(EventDate.year) < 1981 and int(EventDate.year) > 1963:
              Server.DownloadServerFileLowRes(int(EventDate.year))
         else:
              print("Server data only valid for 1963 to present, please enter a valid date.")
         BxS, ByS, BzS, VS, DensityS, PdynS, KpS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS = Server.GetServerData(EventDate,External)
         IOPTinput = misc.IOPTprocess(KpS)
         WindCreate = solar_wind.Solar_Wind(VS, vy, vz, BxS, ByS, BzS, DensityS, PdynS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, KpS, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS)
         WindArray = WindCreate.GetWind()
         
    if LiveData == 1:
         if External == 7 or External == 11:
              print("LIVE DATA NOT SUPPORTED FOR TSY04 OR TA16 MAGNETOSPHERIC MODELS. PLEASE SELECT ANOTHER EXTERNAL MAGNETIC FIELD MODEL.")
              exit()
         misc.DateCheck(EventDate)
         DstLive, VxLive, DensityLive, ByLive, BzLive, IOPTLive, G1Live, G2Live, G3Live, KpLive, By_avgLive, Bx_avgLive, Bz_avgLive, NIndexLive, BIndexLive = Request.Get_Data(EventDate)
         PdynLive = misc.Pdyn_comp(DensityLive,VxLive)
         IOPTinput = IOPTLive
         kp = KpLive
         WindCreate = solar_wind.Solar_Wind(VxLive, vy, vz, Bx_avgLive, ByLive, BzLive, DensityLive, PdynLive, DstLive, G1Live, G2Live, G3Live, W1, W2, W3, W4, W5, W6, KpLive, By_avgLive, Bz_avgLive, NIndexLive, BIndexLive, sym_h_corrected)
         WindArray = WindCreate.GetWind()

    if ServerData == 0 and LiveData == 0:
          if vx > 0:
               vx = -1*vx
          WindCreate = solar_wind.Solar_Wind(vx, vy, vz, bx, by, bz, density, pdyn, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, kp, by_avg, bz_avg, n_index, b_index, sym_h_corrected, External)
          WindArray = WindCreate.GetWind()

    Rigidity = rigidity

    MagFieldModel = np.array([Internal,External])

    EndParams = [minaltitude,maxdistance,maxtime]

    CreateStations = stations.Stations(Stations, startaltitude, Zenith, Azimuth)
    InputtedStations = CreateStations
    if 'customlocations' in locals() or 'customlocations' in globals():
          CreateStations.AddLocation(customlocations)
    Used_Stations_Temp = CreateStations.GetStations()
    temp = list(Used_Stations_Temp)
    Station_Array = temp

    
    ParticleArray = [Anum,AntiCheck]

    misc.ParamCheck(startaltitude,year,Internal,EndParams)

    TrajectoryInputArray = [Rigidity,DateArray,MagFieldModel,IntModel,ParticleArray,IOPTinput,WindArray,Magnetopause,
    coordsystem,gyropercent,EndParams, Station_Array, InputtedStations, KpS, corenum, LiveData, ServerData, g, h]

    return TrajectoryInputArray