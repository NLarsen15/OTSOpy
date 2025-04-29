import numpy as np
from datetime import datetime,timedelta
import os
from . import date, solar_wind, stations
from . import misc, Request, Server

def MagFieldInputs(Locations,
           serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,year,
           month,day,hour,minute,second,internalmag,externalmag,
           coordsystemIN,g,h,corenum,MHDfile, MHDcoordsys):
    
    EventDate = datetime(year,month,day,hour,minute,second)
    DateCreate = date.Date(EventDate)
    DateArray = DateCreate.GetDate()

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

    if coordsystemIN not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid coordsystem: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

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

    startaltitude = 100
    EndParams = [20,100,0]
    misc.ParamCheck(startaltitude,year,Internal,EndParams)

    MagfieldInputArray = [Locations,DateArray,MagFieldModel,IOPTinput,WindArray,coordsystemIN,KpS,corenum,LiveData,serverdata,g,h]

    return MagfieldInputArray