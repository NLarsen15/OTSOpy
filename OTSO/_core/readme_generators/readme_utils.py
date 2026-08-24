from importlib.metadata import version, PackageNotFoundError


def OTSOVersion() -> str:
    try:
        return version("OTSO")
    except PackageNotFoundError:
        return "Unknown"

def OnlineDataStatus(LiveData, serverdata) -> str:
    if LiveData == 1:
       OnlineData = "Online Space Weather Data Used"
    elif serverdata == "ON":
       OnlineData = "Server Data Generated Using OMNI Database Used"
    else:
       OnlineData = "User Inputted Data Used"
    return OnlineData

def ParticleCheck(AntiCheck: int) -> str:
    return "anti-particle" if AntiCheck == 1 else "Normal Particle"

def IntegrationMethodCheck(IntModel: int) -> str:
    IntegrationMethods = [
        "4th Order Runge-Kutta Method",
        "Vay Method",
        "Higuera-Cary Method",
        "6th Order Runge-Kutta Method",
        "5th Order Runge-Kutta Method",
        "Boris-Buneman Method",
    ]
    return IntegrationMethods[IntModel-1] if 1 <= IntModel <= len(IntegrationMethods) else "Unknown Integration Method"

def InternalModelCheck(model: str) -> str:
    InternalModels = [
        "NONE",
        "IGRF",
        "Dipole",
        "Custom Gauss",
        "CHAOS"
    ]

    return model if model in InternalModels else "Unknown Internal Model"
    
def ExternalModelCheck(model: list) -> str:
        ExternalModels = [
        "No External Field", "Tsyganenko 87 Short", "Tsyganenko 87 Long", "Tsyganenko 89a",
        "Tsyganenko 96", "Tsyganenko 01", "Tsyganenko 01 Storm", "Tsyganenko 04", "Tsyganenko 89c",
        "Tsyganenko 15 N Model", "Tsyganenko 15 B Model", "TA16 RBF Model"
    ]
        External = ExternalModels[model[1]] if 0 <= model[1] <= 11 else "Unknown External Model"
        if model[1] == 99:
             External = "MHD file"
        elif model[1] == 100:
             External = "Tsyganenko 89 Refit"
        return External

def MagnetopauseModelCheck(Magnetopause: int, spheresize: float) -> str:
    PauseModels = {
        0: f"{spheresize}Re Sphere",
        1: "Aberrated Formisano Model",
        2: "Sibeck Model",
        3: "Kobel Model",
        4: "Lin 2010 Model",
        99: "No Magnetopause Model",
    }
    return PauseModels.get(Magnetopause, "Unknown Magnetopause Model")

def CutoffCompCheck(Rcomp: str) -> str:
    if Rcomp == "Vertical" or Rcomp == 0:
      return "Vertical Cutoff Rigidity"
    elif Rcomp == "Apparent" or Rcomp == 1:
      return "Apparent Cutoff Rigidity"
    else:
      return "Custom Cutoff"
    
def RigidityScanCheck(Rscan: int) -> str:
    return "Rigidity Scan Used" if Rscan != 0 else "No Rigidity Scan"

def solar_wind_readme_section(result: str, WindArray: list) -> str:
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(f"Vx = {abs(round(WindArray[0], 3))}\n")
    result.append(f"Vy = {round(WindArray[1], 3)}\n")
    result.append(f"Vz = {round(WindArray[2], 3)}\n\n")
    result.append(f"IMF [nT]:\n")
    result.append(f"Bx = {round(WindArray[3], 3)}\n")
    result.append(f"By = {round(WindArray[4], 3)}\n")
    result.append(f"Bz = {round(WindArray[5], 3)}\n\n")
    result.append(f"By Average (30mins) = {round(WindArray[20], 3)}\n")
    result.append(f"Bz Average (30mins) = {round(WindArray[21], 3)}\n\n")
    result.append(f"Density = {round(WindArray[6], 3)} cm^-3\n")
    result.append(f"Pdyn = {round(WindArray[7], 3)} nPa\n\n")
    result.append(f"Dst = {round(WindArray[8], 3)} nT\n\n")
    result.append(f"G1 = {round(WindArray[9], 3)}\n")
    result.append(f"G2 = {round(WindArray[10], 3)}\n")
    result.append(f"G3 = {round(WindArray[11], 3)}\n\n")
    result.append(f"W1 = {round(WindArray[12], 3)}\n")
    result.append(f"W2 = {round(WindArray[13], 3)}\n")
    result.append(f"W3 = {round(WindArray[14], 3)}\n")
    result.append(f"W4 = {round(WindArray[15], 3)}\n")
    result.append(f"W5 = {round(WindArray[16], 3)}\n")
    result.append(f"W6 = {round(WindArray[17], 3)}\n\n")
    result.append(f"N = {round(WindArray[22], 3)}\n")
    result.append(f"B = {round(WindArray[23], 3)}\n\n")
    result.append(f"SYM-H = {round(WindArray[24], 3)}\n\n")
    return result

def field_models_readme_section(result: str, Internal: str, External: str, PauseModel: str, max_degree: int) -> str:
    result.append(f"Magnetic Field Models:\n")
    result.append(f"Internal Model = {Internal}\n")
    result.append(f"Max Spherical Harmonics Order/Degree = {max_degree}\n")
    result.append(f"External Model = {External}\n\n")
    result.append(f"Magnetopause Model = {PauseModel}\n\n")
    return result

def boberg_readme_section(result: str, Boberg: bool, BobergType: str) -> str:
    if Boberg:
        result.append(f"Boberg Extension Applied: {BobergType}\n\n")
    else:
        result.append(f"Boberg Extension Applied: No Boberg Extension Applied\n\n")
    return result


def beta_readme_section(result: str, TotalBetaCheck: bool, betaerror: float,
                        adaptivestep: bool, fixedstep: float = 0) -> str:
    if TotalBetaCheck:
        result.append(f"Total Beta Error Check Applied: Yes\n")
    else:
        result.append(f"Total Beta Error Check Applied: No\n")

    if adaptivestep:
        result.append(f"Adaptive Step Size Used\n")
    elif fixedstep and fixedstep > 0:
        result.append(f"Fixed Step Size Used: {fixedstep}s\n")
    else:
        result.append(f"Fixed Step Size Used\n")

    result.append(f"Beta Error Tolerance: {betaerror}%\n\n")
    return result

def end_conditions_readme_section(result: str, EndParams: list) -> str:
    result.append(f"Minimum Altitude: {EndParams[0]}km\n")
    result.append(f"Max Distance: {EndParams[1]}Re\n")
    result.append(f"Max Time: {EndParams[2]}s\n")
    result.append(f"Max Steps: {EndParams[3]}\n\n")
    return result