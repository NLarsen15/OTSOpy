import numpy as np
from datetime import datetime,timedelta
import os
from . import date, solar_wind, stations
from . import misc, Request, Server

def TraceInputs(startaltitude,Coordsys,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,year,
           month,day,hour,minute,second,internalmag,externalmag,
           gyropercent,magnetopause,corenum,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile, MHDcoordsys):
    
    EventDate = datetime(year,month,day,hour,minute,second)
    DateCreate = date.Date(EventDate)
    DateArray = DateCreate.GetDate()

    if Coordsys not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
      print("Please select a valid coordsystem: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
      exit()

    if magnetopause == "Sphere":
         Magnetopause = 0
    elif magnetopause == "aFormisano":
         Magnetopause = 1
    elif magnetopause == "Sibeck":
         Magnetopause = 2
    elif magnetopause == "Kobel":
         Magnetopause = 3
    else:
         print("Please enter a valid magnetopause model: ""Sphere"", ""aFormisano"", ""Sibeck"", ""Kobel"" ")
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
         if not g or not h: 
            g = [0] * 105
            h = [0] * 105
    elif internalmag == "IGRF":
         Internal = 1
         if not g or not h: 
            g = [0] * 105
            h = [0] * 105
    elif internalmag == "Dipole":
         Internal = 2
         if not g or not h: 
            g = [0] * 105
            h = [0] * 105
    elif internalmag == "Custom Gauss":
         Internal = 4
         if not g or not h:
              print("Please enter values for the g and h Gaussian coefficients to use the Custom Gauss option")
              exit()
         elif len(g) != 105:
              print(f"There should be 105 g coefficents in the inputted list, you have entered {len(g)}")
              exit()
         elif len(h) != 105:
              print(f"There should be 105 h coefficents in the inputted list, you have enetered {len(h)}")
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
    elif externalmag == "MHD":
         External = 99
         if not os.path.exists(MHDfile):
            print(f"The file '{MHDfile}' does not exist.")
            exit()
    else:
         print("Please enter a valid externalmag model: ""NONE"", ""TSY87short"",""TSy87long"",""TSY89"",""TSY89_BOBERG"",""TSY96"",""TSY01"",""TSY01S"",""TSY04""")
         exit()
    
    MinAlt = 0
    MaxDist = 1000
    MaxTime = 0
    AtomicNum = 1
    AntiCheck = 1
    IntModel = 2
    gyropercent = 20
    Rigidity = 1
    Zenith = 0
    Azimuth = 0

    misc.DataCheck(ServerData,LiveData,EventDate)

    IOPTinput = misc.IOPTprocess(kp)
    KpS = 0

    if ServerData == 1:
         if int(EventDate.year) > 1981:
              Server.DownloadServerFile(int(EventDate.year))
         elif int(EventDate.year) < 1981 and int(EventDate.year) > 1963:
              Server.DownloadServerFileLowRes(int(EventDate.year))
         else:
              print("Server data only valid for 1963 to present, please enter a valid date.")
         ByS, BzS, VS, DensityS, PdynS, KpS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S = Server.GetServerData(EventDate,External)
         IOPTinput = misc.IOPTprocess(KpS)
         WindCreate = solar_wind.Solar_Wind(VS, vy, vz, ByS, BzS, DensityS, PdynS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, KpS)
         WindArray = WindCreate.GetWind()
         
    if LiveData == 1:
         misc.DateCheck(EventDate)
         DstLive, VxLive, DensityLive, ByLive, BzLive, IOPTLive, G1Live, G2Live, G3Live, KpLive = Request.Get_Data(EventDate)
         PdynLive = misc.Pdyn_comp(DensityLive,VxLive)
         IOPTinput = IOPTLive
         WindCreate = solar_wind.Solar_Wind(VxLive, vy, vz, ByLive, BzLive, DensityLive, PdynLive, DstLive, G1Live, G2Live, G3Live, W1, W2, W3, W4, W5, W6, KpLive)
         WindArray = WindCreate.GetWind()

    if ServerData == 0 and LiveData == 0:
          if vx > 0:
               vx = -1*vx
          WindCreate = solar_wind.Solar_Wind(vx, vy, vz, by, bz, density, pdyn, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, kp)
          WindArray = WindCreate.GetWind()

    MagFieldModel = np.array([Internal,External])
    ParticleArray = [AtomicNum,AntiCheck]

    EndParams = [MinAlt,MaxDist,MaxTime]

    if latstep > 0:
         latstep = -1*latstep

    if not ((maxlat-minlat)/abs(latstep)).is_integer():
          print(f"Latitude range selected can not be split into integer number of steps using inputted latstep. \nlat range:{maxlat-minlat} \nlatstep:{abs(latstep)} \nlatrange/latstep:{(maxlat-minlat)/abs(latstep)}\nPlease edit the inputs appropriately.")
          exit()

    if not ((maxlong-minlong)/abs(longstep)).is_integer():
          print(f"Longitude range selected can not be split into integer number of steps using inputted longstep. \nlong range:{maxlong-minlong} \nlongstep:{abs(longstep)} \nlongrange/longstep:{(maxlong-minlong)/abs(longstep)}\nPlease edit the inputs appropriately.")
          exit()

    LatitudeList = np.arange(maxlat,minlat + latstep,latstep)
    LongitudeList = np.arange(minlong,maxlong + longstep,longstep)

    misc.ParamCheck(startaltitude,year,Internal,EndParams)

    TraceInputArray = [LongitudeList,LatitudeList,Rigidity,ParticleArray,DateArray,MagFieldModel,IntModel,IOPTinput,WindArray,Magnetopause,gyropercent,EndParams, Zenith, Azimuth, corenum, startaltitude, LiveData, AntiCheck, g, h]

    return TraceInputArray