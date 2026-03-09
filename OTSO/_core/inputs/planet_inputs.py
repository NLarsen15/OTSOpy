import numpy as np
from datetime import datetime
import warnings # Import warnings
import psutil

from ..utils import input_utils, tsy_params_utils
from ..custom_classes import date, solar_wind
from ..livedata import pull_live_data
from ..serverdata import server
from ..data_classes.planet_data import PlanetData

def PlanetInputs(Data: 'PlanetData') -> None: # Add flag
    
    EventDate = datetime(Data.year,Data.month,Data.day,Data.hour,Data.minute,Data.second)
    DateCreate = date.Date(EventDate)
    Data.datearray = DateCreate.GetDate()

    input_utils.asymptotic_check(Data.asymptotic, Data.unit)

    AntiCheck = input_utils.anti_check(Data.anti)

    Data.magnetopause = input_utils.magnetopause_check(Data.mpause)

    Data.integrationmodel = input_utils.intmodel_check(Data.intmodel)

    Data.Rcomp, Data.zenith, Data.azimuth = input_utils.cutoff_comp_check(Data.cutoff_comp, Data.zenith, Data.azimuth)

    ServerData = input_utils.serverdata_check(Data.serverdata)

    LiveData = input_utils.livedata_check(Data.livedata)

    Internal, Data.g, Data.h = input_utils.internalmag_check(Data.internalmag, Data.datearray, Data.g, Data.h)
    
    External = input_utils.externalmag_check(Data.externalmag, Data.MHDfile)

    Bobon, bobtype = input_utils.BobergCheck(Data.boberg, Data.bobergtype)

    # !!!!!!!! NEED TO CHECK THIS FUNCTION NO COORDSYSTEM VARIABLE IN SCRIPT !!!!!!!!!
    input_utils.coordsystem_check("GDZ", Data.inputcoord) # Dummy value for coordsystem to use existing function
    #######################################
    
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
         BxS, ByS, BzS, VS, DensityS, PdynS, KpS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS, External = server.GetServerData(EventDate,External,Data.AdaptiveExternalModel)
         Data.IOPT = tsy_params_utils.IOPTprocess(KpS)
         if External == 100:
               Data.IOPT = tsy_params_utils.IOPTprocess_refit(KpS)
         Data.Kp = KpS
         WindCreate = solar_wind.SolarWind(VS, Data.vy, Data.vz, BxS, ByS, BzS, DensityS, PdynS, DstS, G1S, G2S, G3S, W1S, W2S, W3S, W4S, W5S, W6S, KpS, By_avgS, Bz_avgS, N_indexS, B_indexS, SYM_H_correctedS)
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

    Anum = 1
    Data.particlearray = [Anum,AntiCheck]

    # --- Start: Coordinate Generation Logic --- 
    coordinate_pairs = []
    LatitudeList_meta = [] # For metadata/README
    LongitudeList_meta = [] # For metadata/README

    # Check and handle coordinate input
    if Data.array_of_lats_and_longs is not None:
        try:
            coord_array = np.array(Data.array_of_lats_and_longs)
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
        Data.LatitudeList_meta = LatitudeList_meta
        Data.LongitudeList_meta = LongitudeList_meta

        # --- Conditional Warning --- 
        # Only warn if grid params were explicitly set by the user
        if Data.grid_params_user_set:
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
        if Data.latstep == 0 or Data.longstep == 0:
            raise ValueError("latstep and longstep must be non-zero for grid generation.")

        if Data.latstep > 0:
             warnings.warn("latstep is positive. Latitude generation usually expects a negative step to go from maxlat down to minlat. Proceeding, but check your parameters.", UserWarning)

        lat_start, lat_end = (Data.maxlat, Data.minlat) if Data.latstep < 0 else (Data.minlat, Data.maxlat)
        lon_start, lon_end = (Data.minlong, Data.maxlong)

        epsilon_lat = abs(Data.latstep / 1000.0)
        epsilon_lon = abs(Data.longstep / 1000.0)
        
        _LatitudeList_np = np.arange(lat_start, lat_end + (np.sign(Data.latstep) * epsilon_lat), Data.latstep)
        _LongitudeList_np = np.arange(lon_start, lon_end + epsilon_lon, Data.longstep)

        if _LatitudeList_np.size == 0:
            raise ValueError(f"Generated LatitudeList is empty. Check maxlat ({Data.maxlat}), minlat ({Data.minlat}), and latstep ({Data.latstep}).")
        if _LongitudeList_np.size == 0:
            raise ValueError(f"Generated LongitudeList is empty. Check maxlong ({Data.maxlong}), minlong ({Data.minlong}), and longstep ({Data.longstep}).")

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
    Data.coordinate_pairs = coordinate_pairs

    #PlanetInputArray = [coordinate_pairs, RigidityArray, DateArray, MagFieldModel, IntModel, ParticleArray, IOPTinput, WindArray, Magnetopause, gyropercent, EndParams, CutoffComputation, Rscan, Zenith, Azimuth, corenum, asymptotic, asymlevels, startaltitude, LiveData, AntiCheck, g, h, LatitudeList_meta, LongitudeList_meta]

    return

def CheckCoreNumPlanet(x: int) -> int:
  NewCore = x
  if(psutil.cpu_count(logical=True) < x):
    print("ERROR: You have entered an invalid number of cores")
    print("You have " + str(psutil.cpu_count(logical=True)) + " and have tried to use " + str(x) + " cores")
    print("To ensure operational integrity of your computer OTSO will run using 2 less than the max cores available, with a minumum value of 1.")
    NewCore = psutil.cpu_count(logical=True) - 2
    if NewCore <= 0:
      NewCore = 1
  return NewCore