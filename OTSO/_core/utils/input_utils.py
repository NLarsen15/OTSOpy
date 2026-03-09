import os
from datetime import datetime, timedelta

from .igrf_utils import compute_gauss_coefficients, schmidt_normalize, create_custom_coefficient_arrays


def anti_check(anti: str) -> int:
    if anti == "YES":
         return 1
    elif anti == "NO":
         return 0
    else:
         print("Please enter a valid anti value: ""YES"" or ""NO"" ")
         exit()

def magnetopause_check(magnetopause: str) -> int:
    if magnetopause == "Sphere":
         return 0
    elif magnetopause == "aFormisano":
         return 1
    elif magnetopause == "Sibeck":
         return 2
    elif magnetopause == "Kobel":
         return 3
    elif magnetopause == "Lin":
         return 4
    elif magnetopause == "NONE":
         return 99
    else:
         print("Please enter a valid magnetopause model: ""Sphere"", ""aFormisano"", ""Sibeck"", ""Kobel"", ""Lin"", ""NONE"" ")
         exit()

def intmodel_check(intmodel: str) -> int:
    if intmodel == "4RK":
         return 1
    elif intmodel == "Boris":
         return 2
    elif intmodel == "Vay":
         return 3
    elif intmodel == "HC":
         return 4
    else:
         print("Please enter a valid intmodel model: ""4RK"", ""Boris"", ""Vay"", ""HC"" ")
         exit()

def serverdata_check(serverdata: str) -> int:
    if serverdata == "ON":
         return 1
    elif serverdata == "OFF":
         return 0
    else:
         print("Please enter a valid serverdata value: ""ON"" or ""OFF"" ")
         exit()

def livedata_check(livedata: str) -> int:
    if livedata == "ON":
         return 1
    elif livedata == "OFF":
         return 0
    else:
         print("Please enter a valid livedata value: ""ON"" or ""OFF"" ")
         exit()

def internalmag_check(internalmag: str, DateArray, g, h) -> tuple[int, list[float], list[float]]:
    if internalmag == "NONE":
        coeffs = compute_gauss_coefficients(DateArray)
        g = coeffs['G_coefficients']
        h = coeffs['H_coefficients']
        if g is None or h is None or len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
        return 0, g, h
    elif internalmag == "IGRF":
        coeffs = compute_gauss_coefficients(DateArray)
        g = coeffs['G_coefficients']
        h = coeffs['H_coefficients']
        if len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
        return 1, g, h
    elif internalmag == "Dipole":
        coeffs = compute_gauss_coefficients(DateArray)
        g = coeffs['G_coefficients']
        h = coeffs['H_coefficients']
        if len(g) == 0 or len(h) == 0: 
            g = [0] * 136
            h = [0] * 136
        return 2, g, h
    elif internalmag == "Custom Gauss":
        G, H = create_custom_coefficient_arrays(g, h, max_degree=15)
        G_norm, H_norm = schmidt_normalize(G, H, 15)
        g = G_norm[1:137]  # Extract indices 2-137 to match IGRF format (skip indices 0,1)
        h = H_norm[1:137]  # Extract indices 2-137 to match IGRF format (skip indices 0,1)
        if g is None or h is None:
                print("Please enter values for the g and h Gaussian coefficients to use the Custom Gauss option")
                exit()
        elif len(g) != 136:
                print(f"There should be 136 g coefficents in the inputted list, you have entered {len(g)}")
                exit()
        elif len(h) != 136:
                print(f"There should be 136 h coefficents in the inputted list, you have enetered {len(h)}")
        return 4, g, h
    else:
        print("Please enter a valid internalmag model: ""NONE"",""IGRF"",""Dipole"", or ""Custom Gauss""")
        exit()

def externalmag_check(externalmag: str, MHDfile: str) -> int:
    if externalmag == "NONE":
         return 0
    elif externalmag == "TSY87short":
         return 1
    elif externalmag == "TSY87long":
         return 2
    elif externalmag == "TSY89a":
         return 3
    elif externalmag == "TSY96":
         return 4
    elif externalmag == "TSY01":
         return 5
    elif externalmag == "TSY01S":
         return 6
    elif externalmag == "TSY04":
         return 7
    elif externalmag == "TSY89c":
         return 8
    elif externalmag == "TSY15N":
         return 9
    elif externalmag == "TSY15B":
         return 10
    elif externalmag == "TA16_RBF":
         return 11
    elif externalmag == "TSY89_refit":
         return 100
    elif externalmag == "MHD":
         if not os.path.exists(MHDfile):
            print(f"The file '{MHDfile}' does not exist.")
            exit()
         return 99
    else:
          print("Please enter a valid externalmag model: ""NONE"", ""TSY87short"",""TSy87long"",""TSY89a"",""TSY89c"",""TSY89_refit"",""TSY96"",""TSY01"",""TSY01S"",""TSY04"",""TSY15N"",""TSY15B"",""TA16_RBF""")
          exit()

def coordsystem_check(coordsystem: str, inputcoord: str) -> None:
    valid_coordsystems = ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]
    if coordsystem not in valid_coordsystems:
         print("Please select a valid coordsystem: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    if inputcoord not in valid_coordsystems:
         print("Please select a valid inputcoord: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    return

def coordsystem_in_check(coordsystemIN: str) -> None:
    valid_coordsystems = ["GDZ","GEO","GSM","GSE","SM","GEI","MAG","SPH","RLL"]
    if coordsystemIN not in valid_coordsystems:
         print("Please select a valid coordsystemIN: ""GDZ"", ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""SPH"", ""RLL""")
         exit()
    return

def coordsystem_out_check(coordsystemOUT: str) -> None:
    valid_coordsystems = ["GEO","GSM","GSE","SM","GEI","MAG","RLL"]
    if coordsystemOUT not in valid_coordsystems:
         print("Magfield currently only supports cartesian outputs.")
         print("Please select a valid coordout: ""GEO"", ""GSM"", ""GSE"", ""SM"", ""GEI"", ""MAG"", ""RLL""")
         exit()
    return

def rigidityscan_check(rigidityscan: str) -> int:
    if rigidityscan == "ON":
         return 1
    elif rigidityscan == "OFF":
         return 0
    else:
         print("Please enter a valid rigidityscan value: ""ON"" or ""OFF"" ")
         exit()

def cutoff_comp_check(cutoff_comp: str, zenith: float = None, azimuth: float = None) -> tuple[int, float, float]:
    if cutoff_comp == "Vertical":
         return 0, 0, 0
    elif cutoff_comp == "Apparent":
         return 1, 0, 0
    elif cutoff_comp == "Custom":
         return 2, zenith, azimuth
    else:
         print("Please enter a valid cutoff_comp: ""Vertical"", ""Apparent"", or ""Custom"" ")
         exit()

def flight_server_live_check(variablelist: list, variablelist2: list, serverdata: str, livedata: str) -> None:
     if serverdata == "OFF" and livedata == "OFF":
          if any(not lst for lst in variablelist):
              print("If not using livedata or server data you must provide full lists for the input variables\n"
              "vx,vy,vz,by,bz,density,pdyn,Dst,G1,G2,G3,W1,W2,W3,W4,W5,W6,kp")
              exit()
          if any(len(lst) != len(variablelist2[0]) for lst in variablelist2):
              print("If not using livedata or server data all provided variable lists must be the same length")
              exit()

def asymptotic_check(asymptotic: str, unit: str) -> None:
    if unit not in ["GeV", "GV"]:
        print("Please enter a valid asymlevel unit: ""GeV"" or ""GV"" ")
        exit()
    if asymptotic not in ["YES","NO"]:
        print("Please enter a valid asymptotic value: ""YES"" or ""NO"" ")
        exit()

def DataCheck(ServerData: int, LiveData: int, EventDate: datetime):
       current_date = datetime.utcnow()
       if (ServerData == 1 or LiveData == 1) and EventDate > current_date:
        print("ERROR: Future date entered. No valid data available. \nOTSO program will now terminate.")
        exit()
       if ServerData == 1 and LiveData == 1:
        print("ERROR: You have requested live and server data. Only one can be selected. \nOTSO program will now terminate.")
        exit()

def ParamCheck(Alt: float, EndParams: list[float]) -> None:
       if EndParams[0] > Alt:
        print("ERROR: Inputted minimum altitude is larger than starting altitude. Value must be less than or equal to the starting altitude. Please check inputs. \nOTSO program will now terminate.")
        exit()

def DateCheck(Date):
     current_date = datetime.utcnow()
     seven_days_ago = current_date - timedelta(days=7)

     if Date < seven_days_ago:
         print("ERROR: Inputed date is over 7 days ago from current time (" + str(current_date) + ").\n live data only available for the last week. \nOTSO will now terminate.")
         exit()

     return

def CustomLocationsCheck(customlocations: list) -> None:
    if isinstance(customlocations, list) and all(
        isinstance(item, list) and len(item) == 3 and isinstance(item[0], str) and isinstance(item[1], (int, float)) and isinstance(item[2], (int, float))
        for item in customlocations
    ):
        return
    else:
        print("customlocations must be a list of a list with the structure [[\"NAME1\", lat1, lon1], [\"NAME2\", lat2, lon2], ...]")
        exit()

def BobergCheck(boberg: bool, bobergtype: str) -> tuple | int:
     if boberg:
          Bobon = 1
     else:
          Bobon = 0

     if bobergtype == "EXTENSION":
          bobtype = 1
     elif bobergtype == "CONTINUOUS":
          bobtype = 2
     elif bobergtype == "DST_DEPENDENT":
          bobtype = 3
     else:
          print("Please enter a valid bobergtype: ""EXTENSION"", ""CONTINUOUS"", or ""DST_DEPENDENT"" ")
          exit()

     return Bobon, bobtype