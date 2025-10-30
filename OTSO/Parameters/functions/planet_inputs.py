import numpy as np
from datetime import datetime,timedelta
import os
from . import date, solar_wind, stations
from . import misc, Request, Server
import warnings # Import warnings
from .igrf_process import compute_gauss_coefficients, schmidt_normalize

def PlanetInputs(startaltitude,cutoff_comp,minaltitude,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,rigidityscan,
           gyropercent,magnetopause,corenum,azimuth,zenith, asymptotic,asymlevels,unit,
           latstep,longstep,maxlat,minlat,maxlong,minlong,g,h,MHDfile, MHDcoordsys,inputcoord,
           array_of_lats_and_longs=None,
           grid_params_user_set=False): # Add flag
    
    EventDate = datetime(year,month,day,hour,minute,second)
    DateCreate = date.Date(EventDate)
    DateArray = DateCreate.GetDate()

    if unit not in ["GeV", "GV"]:
         print("Please enter a valid asymlevel unit: ""GeV"" or ""GV"" ")
         exit()

    if anti == "YES":
         AntiCheck = 1
    elif anti == "NO":
         AntiCheck = 0
    else:
         print("Please enter a valid anti value: ""YES"" or ""NO"" ")
         exit()

    if asymptotic not in ["YES","NO"]:
         print("Please enter a valid asymptotic value: ""YES"" or ""NO"" ")
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

    if cutoff_comp == "Vertical":
         CutoffComputation = 0
         Zenith = 0
         Azimuth = 0
    elif cutoff_comp == "Apparent":
         CutoffComputation = 1
         Zenith = 0
         Azimuth = 0
    elif cutoff_comp == "Custom":
         CutoffComputation = 2
         Zenith = zenith
         Azimuth = azimuth
    else:
         print("Please enter a valid cutoff_comp: ""Vertical"", ""Apparent"", or ""Custom"" ")
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
         g, h = schmidt_normalize(g, h, 15)
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

    if inputcoord not in ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]:
         print("Please select a valid inputcoord: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()

    
    if rigidityscan == "ON":
         Rscan = 1
    elif rigidityscan == "OFF":
         Rscan = 0
    else:
         print("Please enter a valid rigidityscan value: ""ON"" or ""OFF"" ")
         exit()


    misc.DataCheck(ServerData,LiveData,EventDate)

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
          WindCreate = solar_wind.Solar_Wind(vx, vy, vz, bx, by, bz, density, pdyn, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, kp, by_avg, bz_avg, n_index, b_index, sym_h_corrected)
          WindArray = WindCreate.GetWind()

    RigidityArray = [startrigidity,endrigidity,rigiditystep]

    MagFieldModel = np.array([Internal,External])

    EndParams = [minaltitude,maxdistance,maxtime]

    ParticleArray = [Anum,AntiCheck]

    misc.ParamCheck(startaltitude,year,Internal,EndParams)

    # --- Start: Coordinate Generation Logic --- 
    coordinate_pairs = []
    LatitudeList_meta = [] # For metadata/README
    LongitudeList_meta = [] # For metadata/README

    # Check and handle coordinate input
    if array_of_lats_and_longs is not None:
        try:
            coord_array = np.array(array_of_lats_and_longs)
            if coord_array.ndim != 2 or coord_array.shape[1] != 2:
                raise ValueError("array_of_lats_and_longs must be a list of pairs or a Nx2 array.")
            if not np.issubdtype(coord_array.dtype, np.number):
                 raise ValueError("Coordinates in array_of_lats_and_longs must be numeric.")
            coordinate_pairs = coord_array.tolist()
        except Exception as e:
            raise ValueError(f"Error processing array_of_lats_and_longs: {e}")

        if not coordinate_pairs:
             raise ValueError("Provided array_of_lats_and_longs is empty.")

        lats = coord_array[:, 0]
        lons = coord_array[:, 1]
        LatitudeList_meta = sorted(np.unique(lats), reverse=True)
        LongitudeList_meta = sorted(np.unique(lons))

        # --- Conditional Warning --- 
        # Only warn if grid params were explicitly set by the user
        if grid_params_user_set:
        # if latstep is not None or longstep is not None or maxlat is not None or minlat is not None or maxlong is not None or minlong is not None: # Old check
            warnings.warn("Both array_of_lats_and_longs and step/range parameters were provided. The array_of_lats_and_longs will be used.", UserWarning)
        # --- End Conditional Warning ---
    else:
        # Generate grid using step/range parameters (which now have defaults)
        # Remove the check for missing parameters, as they have defaults now.
        # required_params = [latstep, longstep, maxlat, minlat, maxlong, minlong]
        # if any(p is None for p in required_params):
        #      raise ValueError("Missing latitude/longitude step or range parameters when array_of_lats_and_longs is not provided.")
        
        # Keep the check for zero step
        if latstep == 0 or longstep == 0:
            raise ValueError("latstep and longstep must be non-zero for grid generation.")

        if latstep > 0:
             warnings.warn("latstep is positive. Latitude generation usually expects a negative step to go from maxlat down to minlat. Proceeding, but check your parameters.", UserWarning)

        lat_start, lat_end = (maxlat, minlat) if latstep < 0 else (minlat, maxlat)
        lon_start, lon_end = (minlong, maxlong)

        epsilon_lat = abs(latstep / 1000.0)
        epsilon_lon = abs(longstep / 1000.0)
        
        _LatitudeList_np = np.arange(lat_start, lat_end + (np.sign(latstep) * epsilon_lat), latstep)
        _LongitudeList_np = np.arange(lon_start, lon_end + epsilon_lon, longstep)

        if _LatitudeList_np.size == 0:
            raise ValueError(f"Generated LatitudeList is empty. Check maxlat ({maxlat}), minlat ({minlat}), and latstep ({latstep}).")
        if _LongitudeList_np.size == 0:
            raise ValueError(f"Generated LongitudeList is empty. Check maxlong ({maxlong}), minlong ({minlong}), and longstep ({longstep}).")

        LatitudeList_meta = _LatitudeList_np.tolist()
        LongitudeList_meta = _LongitudeList_np.tolist()
        
        # Generate coordinate pairs, handling poles specially to avoid duplicates
        coordinate_pairs = []
        for lat in LatitudeList_meta:
            if abs(lat) == 90.0:  # At the poles (90° or -90°)
                # Only add one point per pole (use first longitude)
                if not any(abs(coord[0]) == 90.0 and coord[0] == lat for coord in coordinate_pairs):
                    coordinate_pairs.append([lat, LongitudeList_meta[0]])
            else:
                # For non-polar latitudes, add all longitude combinations
                for lon in LongitudeList_meta:
                    coordinate_pairs.append([lat, lon])
    # --- End: Coordinate Generation Logic --- 

    PlanetInputArray = [coordinate_pairs, RigidityArray, DateArray, MagFieldModel, IntModel, ParticleArray, IOPTinput, WindArray, Magnetopause, gyropercent, EndParams, CutoffComputation, Rscan, Zenith, Azimuth, corenum, asymptotic, asymlevels, startaltitude, LiveData, AntiCheck, g, h, LatitudeList_meta, LongitudeList_meta]

    return PlanetInputArray