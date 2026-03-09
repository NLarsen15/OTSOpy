"""
Parameter type definitions for OTSO functions.

This module contains TypedDict classes that define structured parameter groups
used across OTSO functions to improve code organization and maintainability.
"""

from __future__ import annotations
from typing import TypedDict, Sequence, Optional


class SolarWindParams(TypedDict, total=False):
    """Solar wind velocity and magnetic field parameters.
    
    Available keys (all optional):
    • vx: float = -500      # Solar wind velocity x-component (km/s)
    • vy: float = 0         # Solar wind velocity y-component (km/s) 
    • vz: float = 0         # Solar wind velocity z-component (km/s)
    • bx: float = 0         # IMF x-component (nT)
    • by: float = 5         # IMF y-component (nT)
    • bz: float = 5         # IMF z-component (nT)
    • by_avg: float = 0     # Averaged IMF By (nT)
    • bz_avg: float = 0     # Averaged IMF Bz (nT)
    • density: float = 1    # Solar wind density (particles/cm³)
    • pdyn: float = 0       # Solar wind dynamic pressure (nPa)
    
    Example: solar_wind = {"vx": -400, "by": 3.0, "bz": -2.0}
    """
    vx: float  # Solar wind velocity x-component (km/s)
    vy: float  # Solar wind velocity y-component (km/s) 
    vz: float  # Solar wind velocity z-component (km/s)
    bx: float  # IMF x-component (nT)
    by: float  # IMF y-component (nT)
    bz: float  # IMF z-component (nT)
    density: float  # Solar wind density (particles/cm³)
    pdyn: float  # Solar wind dynamic pressure (nPa)
    by_avg: float  # Averaged IMF By (nT)
    bz_avg: float  # Averaged IMF Bz (nT)
    
    # Default values
    DEFAULTS = {
        "vx": -500, "vy": 0, "vz": 0, "bx": 0, "by": 5, "bz": 5, 
        "density": 1, "pdyn": 0, "by_avg": 0, "bz_avg": 0
    }


class GeomagneticParams(TypedDict, total=False):
    """Geomagnetic activity indices and coupling functions.
    
    Available keys (all optional):
    • Dst: float = 0              # Dst index (nT)
    • kp: float = 0               # Kp index (0-9)
    • n_index: float = 0          # Newell coupling function
    • b_index: float = 0          # Boynton coupling function  
    • sym_h_corrected: float = 0  # Corrected SYM-H index (nT)
    
    Example: geomagnetic = {"Dst": -50, "kp": 4.0}
    """
    Dst: float  # Dst index (nT)
    kp: float   # Kp index
    n_index: float  # Newell coupling function
    b_index: float  # Boynton coupling function
    sym_h_corrected: float  # Corrected SYM-H index
    
    # Default values
    DEFAULTS = {
        "Dst": 0, "kp": 0, "n_index": 0, "b_index": 0, "sym_h_corrected": 0
    }


class TsyganenkoParams(TypedDict, total=False):
    """Tsyganenko model coefficients for external field models.
    
    Available keys (all optional):
    • G1: float = 0    # Tsyganenko G1 coefficient
    • G2: float = 0    # Tsyganenko G2 coefficient  
    • G3: float = 0    # Tsyganenko G3 coefficient
    • W1: float = 0    # Tsyganenko W1 coefficient
    • W2: float = 0    # Tsyganenko W2 coefficient
    • W3: float = 0    # Tsyganenko W3 coefficient
    • W4: float = 0    # Tsyganenko W4 coefficient
    • W5: float = 0    # Tsyganenko W5 coefficient
    • W6: float = 0    # Tsyganenko W6 coefficient
    
    Example: tsyganenko = {"G1": 2.5, "G2": -1.8, "W1": 0.5}
    Used with TSY89, TSY96, TSY01, TSY04, and newer models.
    """
    G1: float  # Tsyganenko G1 coefficient
    G2: float  # Tsyganenko G2 coefficient
    G3: float  # Tsyganenko G3 coefficient
    W1: float  # Tsyganenko W1 coefficient
    W2: float  # Tsyganenko W2 coefficient
    W3: float  # Tsyganenko W3 coefficient
    W4: float  # Tsyganenko W4 coefficient
    W5: float  # Tsyganenko W5 coefficient
    W6: float  # Tsyganenko W6 coefficient
    
    # Default values
    DEFAULTS = {
        "G1": 0, "G2": 0, "G3": 0, "W1": 0, "W2": 0, "W3": 0, "W4": 0, "W5": 0, "W6": 0
    }


class DateTimeParams(TypedDict, total=False):
    """Date and time parameters for magnetic field configuration.
    
    Available keys (all optional):
    • year: int = 2024    # Year (e.g., 2023)
    • month: int = 1      # Month (1-12)
    • day: int = 1        # Day (1-31) 
    • hour: int = 12      # Hour (0-23)
    • minute: int = 0     # Minute (0-59)
    • second: int = 0     # Second (0-59)
    
    Example: datetime_params = {"year": 2023, "month": 6, "day": 15}
    """
    year: int    # Year
    month: int   # Month (1-12)
    day: int     # Day (1-31)
    hour: int    # Hour (0-23)
    minute: int  # Minute (0-59)
    second: int  # Second (0-59)
    
    # Default values
    DEFAULTS = {
        "year": 2024, "month": 1, "day": 1, "hour": 12, "minute": 0, "second": 0
    }


class MagFieldParams(TypedDict, total=False):
    """Magnetic field model configuration parameters.
    
    Available keys (all optional):
    • internalmag: str = "IGRF"           # "IGRF", "Dipole", "NONE", "Custom Gauss"
    • externalmag: str = "TSY89c"         # "TSY89c", "TSY01", "TSY15B", "NONE", etc.
    • boberg: bool = False                # Enable Boberg modification  
    • bobergtype: str = "EXTENSION"       # Boberg modification type
    • magnetopause: str = "Kobel"         # "Kobel", "Shue", "Lin", "None", "Sphere"
    • spheresize: float = 25              # Spherical boundary radius (Re)
    • AdaptiveExternalModel: bool = False # Auto-select external model
    
    Example: magfield_params = {"externalmag": "TSY01", "magnetopause": "Shue"}
    """
    internalmag: str  # Internal field model
    externalmag: str  # External field model
    boberg: bool     # Enable Boberg modification
    bobergtype: str  # Boberg modification type
    magnetopause: str # Magnetopause model
    spheresize: float # Spherical boundary size
    AdaptiveExternalModel: bool # Adaptive model selection
    
    # Default values
    DEFAULTS = {
        "internalmag": "IGRF", "externalmag": "TSY89c", "boberg": False, "bobergtype": "EXTENSION",
        "magnetopause": "Kobel", "spheresize": 25, "AdaptiveExternalModel": False
    }


class IntegrationParams(TypedDict, total=False):
    """Particle integration and tracing parameters.
    
    Example usage:
        integration = {"intmodel": "Boris", "gyropercent": 10, "maxdistance": 50}
        cutoff("DOMC", integration_params=integration)
        
    Controls particle tracing algorithm and boundary conditions.
    """
    intmodel: str  # Integration method (Boris, 4RK, Vay, HC)
    gyropercent: float  # Gyration period percentage for time step
    minaltitude: float  # Minimum altitude (km or Re)
    maxdistance: float  # Maximum tracing distance (Re)
    maxtime: float  # Maximum integration time
    mintrapdist: float  # Minimum trapping distance
    startaltitude: float  # Starting altitude (km or Re)
    betaerror: float  # Maximum allowed beta error for integration steps %
    totalbetacheck: bool  # Enable total beta check for whole trace
    adaptivestep: bool  # Enable adaptive time steps
    maxsteps: int  # Maximum number of integration steps
    
    # Default values
    DEFAULTS = {
        "intmodel": "Boris", "gyropercent": 15, "minaltitude": 20, 
        "maxdistance": 100, "maxtime": 0, "mintrapdist": 0, "startaltitude": 20,
        "betaerror": 0.001, "totalbetacheck": True, "adaptivestep": True, "maxsteps": 0}
    


class CoordinateParams(TypedDict, total=False):
    """Coordinate system specifications.
    
    Example usage:
        coordinates = {"coordsystem": "GSM", "inputcoord": "GEO"}
        cutoff("DOMC", coordinate_params=coordinates)
        
    Available systems: GEO, GSM, GSE, SM, GEI, MAG, SPH, RLL.
    """
    coordsystem: str  # Output coordinate system
    inputcoord: str   # Input coordinate system
    coordout: str     # Output coordinate system
    
    # Default values
    DEFAULTS = {
        "coordsystem": "GEO", "inputcoord": "GDZ", "coordout": "GSM"
    }


class ParticleParams(TypedDict, total=False):
    """Particle type and charge specifications.
    
    Example usage:
        particle = {"Anum": 1, "anti": "NO"}  # Protons
        cutoff("DOMC", particle_params=particle)
        
    Anum: 0=electron, 1=proton, 2=alpha particle.
    """
    Anum: int   # Atomic number
    anti: str   # Particle vs anti-particle
    azimuth: float  # Particle azimuth angle (degrees)
    zenith: float   # Particle zenith angle (degrees)
    rigidity: float   # Particle rigidity
    
    # Default values
    DEFAULTS = {
        "Anum": 1, "anti": "YES","azimuth": 0, "zenith": 0, "rigidity": 1
    }


class RigidityParams(TypedDict, total=False):
    """Rigidity scanning parameters for cutoff calculations.
    
    Example usage:
        rigidity = {"startrigidity": 15.0, "endrigidity": 5.0, "rigiditystep": 0.1}
        cutoff("DOMC", rigidity_params=rigidity)
        
    Controls the energy range and resolution of cutoff scanning.
    """
    startrigidity: float  # Initial rigidity (GV)
    endrigidity: float    # Final rigidity (GV)
    rigiditystep: float   # Rigidity step size (GV)
    rigidityscan: str     # Enable/disable scanning
    
    # Default values
    DEFAULTS = {
        "startrigidity": 20, "endrigidity": 0, "rigiditystep": 0.01, "rigidityscan": "ON"
    }


class ComputationParams(TypedDict, total=False):
    """General computation settings.
    
    Example usage:
        computation = {"corenum": 4, "Verbose": False}
        cutoff("DOMC", computation_params=computation)
        
    Controls parallel processing and output verbosity.
    """
    corenum: Optional[int]  # Number of CPU cores
    Verbose: bool           # Enable verbose output
    delim: str              # Delimiter for output formatting
    
    # Default values
    DEFAULTS = {
        "corenum": None, "Verbose": True, "delim": ";"
    }


class DataRetrievalParams(TypedDict, total=False):
    """Server and live data retrieval settings.
    
    Example usage:
        data_retrieval = {"serverdata": "ON", "livedata": "OFF"}
        cutoff("DOMC", data_retrieval_params=data_retrieval)
        
    Enables automatic downloading of space weather data.
    """
    serverdata: str  # Server data retrieval (ON/OFF)
    livedata: str    # Live data retrieval (ON/OFF)
    
    # Default values
    DEFAULTS = {
        "serverdata": "OFF", "livedata": "OFF"
    }


class CustomFieldParams(TypedDict, total=False):
    """Parameters for custom magnetic field models.
    
    Example usage:
        custom = {"g": [1.0, 2.0, 3.0], "MHDfile": "/path/to/mhd_data.txt"}
        cutoff("DOMC", custom_field_params=custom)
        
    For advanced users requiring custom field configurations.
    """
    g: Optional[Sequence[float]]  # Gauss coefficients g
    h: Optional[Sequence[float]]  # Gauss coefficients h
    MHDfile: Optional[str]        # MHD simulation file path
    MHDcoordsys: Optional[str]    # MHD coordinate system
    
    # Default values
    DEFAULTS = {
        "g": None, "h": None, "MHDfile": None, "MHDcoordsys": None
    }


class GridParams(TypedDict, total=False):
    """Grid configuration parameters for field mapping.
    
    Example usage:
        grid = {"latstep": -2.5, "longstep": 2.5, "maxlat": 90, "minlat": -90,
                "maxlong": 360, "minlong": 0}
        planet("DOMC", grid_params=grid)
        
    Defines the spatial grid for magnetic field calculations.
    """
    latstep: float
    longstep: float
    maxlat: float
    minlat: float
    maxlong: float
    minlong: float
    array_of_lats_and_longs: Optional[list]
    
    # Default values
    DEFAULTS = {
        "latstep": -5, "longstep": 5, "maxlat": 90, "minlat": -90,
        "maxlong": 360, "minlong": 0, "array_of_lats_and_longs": None
    }

class AsymptoticParams(TypedDict, total=False):
    """Asymptotic direction calculation parameters.
    
    Example usage:
        asymptotic = {"altitude": 10.0}
        trace("GEO", asymptotic_params=asymptotic)
        
    Controls settings for asymptotic direction computations.
    """
    
    unit: str
    asymptotic: str
    asymlevels: list
    
    # Default values
    DEFAULTS = {
        "unit": "GeV","asymptotic": "NO", 
        "asymlevels": [0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]
    }