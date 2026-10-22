from importlib.metadata import version, PackageNotFoundError

from ..custom_classes.parameter_types import DEFAULTS_SOURCE


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
    # Codes from input_utils.intmodel_check (2 was the removed Boris method)
    IntegrationMethods = {
        1: "4th Order Runge-Kutta Method",
        3: "Vay Method",
        4: "Higuera-Cary Method",
        5: "6th Order Runge-Kutta Method",
        6: "5th Order Runge-Kutta Method",
        7: "Boris-Buneman Method",
    }
    return IntegrationMethods.get(IntModel, "Unknown Integration Method")

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

def solar_wind_readme_section(result: str, WindArray: list, defaulted_inputs: list = None) -> str:
    defaulted = set(defaulted_inputs or [])

    def line(label, key, index, unit=""):
        value = abs(round(WindArray[index], 3)) if key == "vx" else round(WindArray[index], 3)
        mark = " (default)" if key in defaulted else ""
        return f"{label} = {value}{unit}{mark}\n"

    if defaulted:
        result.append(f"Inputs not provided, so OTSO's defaults were used ({DEFAULTS_SOURCE}):\n")
        result.append(f"{', '.join(defaulted_inputs)}\n\n")
    result.append(f"Solar Wind Speed [km/s]:\n")
    result.append(line("Vx", "vx", 0))
    result.append(line("Vy", "vy", 1))
    result.append(line("Vz", "vz", 2) + "\n")
    result.append(f"IMF [nT]:\n")
    result.append(line("Bx", "bx", 3))
    result.append(line("By", "by", 4))
    result.append(line("Bz", "bz", 5) + "\n")
    result.append(line("By Average (30mins)", "by_avg", 20))
    result.append(line("Bz Average (30mins)", "bz_avg", 21) + "\n")
    result.append(line("Density", "density", 6, " cm^-3"))
    if WindArray[7] == 0:
        result.append(f"Pdyn = computed from density and solar wind speed\n\n")
    else:
        result.append(line("Pdyn", "pdyn", 7, " nPa") + "\n")
    result.append(line("Dst", "Dst", 8, " nT") + "\n")
    for i, key in enumerate(("G1", "G2", "G3")):
        result.append(line(key, key, 9 + i))
    result.append("\n")
    for i, key in enumerate(("W1", "W2", "W3", "W4", "W5", "W6")):
        result.append(line(key, key, 12 + i))
    result.append("\n")
    result.append(line("N", "n_index", 22))
    result.append(line("B", "b_index", 23) + "\n")
    result.append(line("SYM-H", "sym_h_corrected", 24, " nT") + "\n")
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