"""
OTSO (Open-source geomagneToSphere prOpagation tool)
============================================

A comprehensive Python package for charged particle trajectory and 
geomagnetic cut-off computations in the Earth’s magnetosphere.
 Main features include:

- Cosmic ray geomagnetic cut-off rigidity calculations
- Asymptotic cone of acceptance computations
- Transmission functions
- geomagnetic cut-off sky map generation  
- Magnetic field line tracing
- Coordinate system transformations
- Particle trajectory tracing in geomagnetic fields

For detailed documentation and example uses, visit: https://github.com/NLarsen15/OTSOpy
"""

from __future__ import annotations
import difflib
import functools
import inspect
import warnings
from typing import TYPE_CHECKING, Sequence, Optional, Union

# Import parameter type definitions from custom classes
from ._core.custom_classes.parameter_types import (
    SolarWindParams,
    GeomagneticParams, 
    TsyganenkoParams,
    DateTimeParams,
    MagFieldParams,
    IntegrationParams,
    CoordinateParams,
    ParticleParams,
    RigidityParams,
    ComputationParams,
    DataRetrievalParams,
    CustomFieldParams,
    GridParams,
    AsymptoticParams,
    TransmissionParams,
    SkymapParams,
    DEFAULTS_SOURCE
)

_PARAM_GROUPS = {cls.__name__: cls for cls in (
    SolarWindParams, GeomagneticParams, TsyganenkoParams, DateTimeParams, MagFieldParams,
    IntegrationParams, CoordinateParams, ParticleParams, RigidityParams, ComputationParams,
    DataRetrievalParams, CustomFieldParams, GridParams, AsymptoticParams, TransmissionParams, SkymapParams)}

# Inputs each solar-wind-driven Tsyganenko model reads (see SolarWindModule.f95). The Kp-driven
# models (TSY87, TSY89, TSY89c, TSY89_refit) are left out: TSY89c is the default external model,
# and its quiet-time default Kp is reported in the metadata instead.
_TSYGANENKO_INPUTS = {
    "TSY96": ("pdyn", "Dst", "by", "bz"),
    "TSY01": ("pdyn", "Dst", "by", "bz", "G1", "G2"),
    "TSY01S": ("pdyn", "Dst", "by", "bz", "G2", "G3"),
    "TSY04": ("pdyn", "Dst", "by", "bz", "W1", "W2", "W3", "W4", "W5", "W6"),
    "TSY15N": ("pdyn", "by_avg", "bz_avg", "n_index"),
    "TSY15B": ("pdyn", "by_avg", "bz_avg", "b_index"),
    "TA16_RBF": ("pdyn", "sym_h_corrected", "n_index", "by_avg"),
}

# Giving any of these without pdyn means pdyn is computed from density and velocity
_PDYN_SOURCES = {"density", "vx", "vy", "vz"}

def _group_keys(cls):
    return list(cls.__annotations__) + [key for key in cls.DEFAULTS if key not in cls.__annotations__]

def _closest(name, candidates):
    """Return the candidate closest to name, ignoring case, or None."""
    lowered = {candidate.lower(): candidate for candidate in candidates}
    match = difflib.get_close_matches(name.lower(), list(lowered), n=1)
    return lowered[match[0]] if match else None

def _check_group_keys(func_name, groups, bound):
    """Raise a TypeError for any key in a *_params dictionary that does not belong to that group."""
    for group, cls in groups.items():
        params = bound.get(group)
        if not isinstance(params, dict):
            continue
        keys = _group_keys(cls)
        for key in params:
            if key in keys:
                continue
            owner = next((other for other, other_cls in groups.items()
                          if other != group and key in _group_keys(other_cls)), None)
            if owner:
                hint = f"it belongs in {owner}"
            else:
                match = _closest(key, keys)
                hint = f"did you mean '{match}'?" if match else f"valid keys are: {', '.join(keys)}"
            raise TypeError(f"{func_name}(): '{key}' is not a valid key for {group} ({hint})")

def _merge_wind_inputs(solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params):
    """Merge the solar wind, geomagnetic and Tsyganenko inputs with OTSO's quiet-time defaults.

    Returns the three merged dictionaries and the names of the inputs that fell back to defaults
    (empty when server or live data supply the inputs), for reporting in the metadata.
    """
    solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
    given = {**solar_wind_params, **geomagnetic_params, **tsyganenko_params}
    if "pdyn" not in given and _PDYN_SOURCES & set(given):
        solar_wind["pdyn"] = 0  # computed from the given density and velocity
        given["pdyn"] = 0
    retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    if str(retrieval["serverdata"]).upper() == "ON" or str(retrieval["livedata"]).upper() == "ON":
        return solar_wind, geomagnetic, tsyganenko, []
    defaulted = [key for cls in (SolarWindParams, GeomagneticParams, TsyganenkoParams)
                 for key in cls.DEFAULTS if key not in given]
    return solar_wind, geomagnetic, tsyganenko, defaulted

def _warn_default_tsyganenko_inputs(func_name, groups, bound):
    """Warn when a solar-wind-driven Tsyganenko model will silently run on default inputs."""
    magfield = {**MagFieldParams.DEFAULTS, **(bound.get("magfield_params") or {})}
    retrieval = {**DataRetrievalParams.DEFAULTS, **(bound.get("data_retrieval_params") or {})}
    required = _TSYGANENKO_INPUTS.get(magfield["externalmag"])
    if (required is None or magfield["AdaptiveExternalModel"]
            or str(retrieval["serverdata"]).upper() == "ON" or str(retrieval["livedata"]).upper() == "ON"):
        return
    input_groups = ("solar_wind_params", "geomagnetic_params", "tsyganenko_params")
    given = {key for group in input_groups for key in (bound.get(group) or {})}
    if _PDYN_SOURCES & given:
        given.add("pdyn")  # computed from the given density and velocity
    missing = [key for key in required if key not in given]
    if not missing:
        return
    defaults = []
    for key in missing:
        group = next(group for group in input_groups if key in groups[group].DEFAULTS)
        defaults.append(f"{key}={groups[group].DEFAULTS[key]} ({group})")
    warnings.warn(
        f"{func_name}(): externalmag='{magfield['externalmag']}' is running with server and live data off, "
        f"but no value was given for {', '.join(missing)}. Using OTSO's defaults ({DEFAULTS_SOURCE}): "
        f"{', '.join(defaults)}. "
        f"Set these explicitly, or use data_retrieval_params={{'serverdata': 'ON'}}.",
        UserWarning, stacklevel=3)

def _check_arguments(func):
    """Check a public function's arguments before calling it.

    Raises a TypeError for unknown argument names and for keys that do not belong in a *_params
    dictionary, suggesting the closest valid name, and warns when a Tsyganenko model would run
    on default solar wind and geomagnetic inputs.
    """
    signature = inspect.signature(func)
    valid = list(signature.parameters)
    groups = {name: _PARAM_GROUPS[param.annotation] for name, param in signature.parameters.items()
              if param.annotation in _PARAM_GROUPS}
    warn_defaults = {"magfield_params", "solar_wind_params"} <= set(groups)

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        unknown = [key for key in kwargs if key not in valid]
        if unknown:
            messages = []
            for key in unknown:
                match = _closest(key, valid)
                messages.append(f"'{key}'" + (f" (did you mean '{match}'?)" if match else ""))
            accepted = f"Valid arguments are: {', '.join(valid)}" if valid else "It takes no arguments"
            raise TypeError(f"{func.__name__}() got unexpected keyword argument(s): {', '.join(messages)}. {accepted}")
        if groups:
            bound = signature.bind(*args, **kwargs).arguments
            _check_group_keys(func.__name__, groups, bound)
            if warn_defaults:
                _warn_default_tsyganenko_inputs(func.__name__, groups, bound)
        return func(*args, **kwargs)

    return wrapper


from ._cutoff import cutoff as cutoff_func
from ._cone import cone as cone_func
from ._planet import planet as planet_func
from ._trajectory import trajectory as trajectory_func
from ._magfield import magfield as magfield_func
from ._coordtrans import coordtrans as coordtrans_func
from ._flight import flight as flight_func
from ._trace import trace as trace_func
from ._transmission import transmission as transmission_func
from ._skymap import skymap as skymap_func

###############################################################################################################################
@_check_arguments
def cutoff(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # asymptotic parameters grouped
    asymptotic_params: AsymptoticParams = {},  # unit,asymptotic,asymlevels
    # transmission parameters grouped
    transmission_params: TransmissionParams = {} # transmission,transmissionsamples,transmissionRstep
) -> list:
    """
    Compute geomagnetic cutoff rigidities for given neutron monitor stations,
    or user-defined locations.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the cutoff
    rigidities and related metadata. Asymptotic viewing directions and transmission 
    functions can also be computed on request. 

    Args:
        Stations (str | list): Station name(s) or identifiers used for cutoff calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].
        cutoff_comp (str): Cutoff computation method ("Vertical", "Apparent", "Custom").

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ADAPTIVE"): Penumbra scanning ("ADAPTIVE"/"ON"/"OFF").
              "ADAPTIVE" extends the fine scan until the penumbra is bounded by clean guard bands;
              "ON" uses a fixed ±0.75 GV window around a 0.25 GV coarse scan; "OFF" scans the full range

        asymptotic_params (AsymptoticParams): Asymptotic computation parameters.

            Available keys:

            - `unit` (`str`, default="GeV"): level unit ("GeV", "GV")
            - `asymptotic` (`str`, default="NO"): Enable asymptotic cone computation ("YES"/"NO")
            - `asymlevels` (`list`, default=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]): level values for asymptotic computation

        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmission` (`bool`, default=False): Bool to enable or disable transmission computation
            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Output coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `threadnum` (`int`, default=None): Number of threads per CPU core for Fortran computations
            - `Verbose` (`bool`, default=True): Enable verbose output
            - `delim` (`str`, default=";"): Delimiter for asymptotic output formatting

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [cutoff_dataframe, asymptotic_dataframe, transmission_dataframe, readme_text]

            - cutoff_dataframe: Cutoff data values per input station in a pandas dataframe. (Ru, Rc, Rl, PTF)
            - asymptotic_dataframe: Global asymptotic viewing directions for the energy/rigidity levels input.
            - transmission_dataframe: Transmission pandas dataframe with transmission function values per input station as a function of rigidity 
            - readme_text: OTSO computation summary

    Examples:
    ```python
       from OTSO import cutoff
 
       if __name__ == '__main__':
 
            stations_list = ["OULU"] # list of neutron monitor stations (using their abbreviations)
 
            # Example using grouped parameters
            cutoff_results = cutoff(
               Stations=stations_list,
               cutoff_comp="Vertical",
               datetime_params={"year": 2005, "month": 5, "day": 1, "hour": 0},
               magfield_params={"internalmag": "IGRF", "externalmag": "NONE"},
               rigidity_params={"rigiditystep": 0.01}
            )
        
            # Access the results
            cutoff_data, asymptotic_data, transmission_data, metadata = result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_asymptotic = {**AsymptoticParams.DEFAULTS, **asymptotic_params}
    final_transmission = {**TransmissionParams.DEFAULTS, **transmission_params}

    #remove unused variables from cutoff
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    
    # Call the original function with expanded parameters
    return cutoff_func(
        Stations,
        customlocations=customlocations,
        cutoff_comp=cutoff_comp,
        # Solar wind parameters
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        # Geomagnetic parameters
        **final_geomagnetic,
        # Tsyganenko coefficients  
        **final_tsyganenko,
        # Date/time parameters
        **final_datetime,
        # Magnetic field parameters
        **final_magfield,
        # Integration parameters
        **final_integration,
        # Particle parameters
        **final_particle,
        # Rigidity parameters
        **final_rigidity,
        # Coordinate parameters
        **final_coordinate,
        # Computation parameters
        **final_computation,
        # Data retrieval parameters
        **final_data_retrieval,
        # Custom field parameters
        **final_custom_field,
        # Asymptotic parameters
        **final_asymptotic,
        # Transmission parameters
        **final_transmission
    )

#######################################################################################################################################
@_check_arguments
def cone(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # transmission parameters grouped
    transmission_params: TransmissionParams = {} # transmission,transmissionsamples,transmissionRstep
) -> list:
    """
    Compute asymptotic viewing cones for given neutron monitor stations,
    or user-defined locations. Also computes cutoff rigidities. Transmission
    functions can be computed on request.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the asymptotic 
    viewing cones, rigidities, and related metadata. Transmission functions will are
    computed only at the the request of the user.

    Args:
        Stations (str | list): Station name(s) or identifiers used for cone calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)

        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmission` (`bool`, default=False): Bool to enable or disable transmission computation
            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Output coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `threadnum` (`int`, default=None): Number of threads per CPU core for Fortran computations
            - `Verbose` (`bool`, default=True): Enable verbose output
            - `delim` (`str`, default=";"): Delimiter for asymptotic output formatting

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [cone_dataframe, cutoff_dataframe, readme_text]

            - cone_dataframe: Asymptotic directions per rigidity (filter;lat;lon) if transmission is enabled the format is (filter;lat;lon;TF)
            - cutoff_dataframe: Cutoff rigidity table (Ru, Rc, Rl, PTF)
            - readme_text: OTSO computation summary

    Examples:
    ```python
        from OTSO import cone

        if __name__ == '__main__':

            stations_list = ["OULU", "ROME"] # list of neutron monitor stations (using their abbreviations)

            cone_results = cone(Stations=stations_list,
                                computation_params={"corenum": 1, "threadnum": 8},
                                datetime_params={"year": 2005, "month": 5, "day": 1, "hour": 0},
                                integration_params={"gyropercent": 1},
                                magfield_params={"internalmag": "IGRF", "externalmag": "TSY89c"},
                                rigidity_params={"rigiditystep":0.001},
                                transmission_params={"transmission": True, "transmissionRstep": 0.0001, "transmissionsamples": 20})
        
            # Access the results
            cone_df, cutoff_df, metadata = cone_result
    ```
    """

    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_transmission = {**TransmissionParams.DEFAULTS, **transmission_params}

    #remove unused variables from cone
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    final_rigidity.pop("rigidityscan", None)

    
    # Call the original function with expanded parameters
    return cone_func(
        Stations,
        customlocations=customlocations,
        # Solar wind parameters
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        # Geomagnetic parameters
        **final_geomagnetic,
        # Tsyganenko coefficients  
        **final_tsyganenko,
        # Date/time parameters
        **final_datetime,
        # Magnetic field parameters
        **final_magfield,
        # Integration parameters (filtered)
        **final_integration,
        # Particle parameters
        **final_particle,
        # Rigidity parameters (filtered)
        **final_rigidity,
        # Coordinate parameters
        **final_coordinate,
        # Computation parameters
        **final_computation,
        # Data retrieval parameters
        **final_data_retrieval,
        # Custom field parameters
        **final_custom_field,
        # Transmission parameters
        **final_transmission
    )
#########################################################################################################################
@_check_arguments
def planet(
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # grid parameters grouped
    grid_params: GridParams = {},  # latstep,longstep,maxlat,minlat,maxlong,minlong,array_of_lats_and_longs
    # Asymptotic parameters grouped
    asymptotic_params: AsymptoticParams = {},  # unit,asymptotic,asymlevels
    # transmission parameters grouped
    transmission_params: TransmissionParams = {} # transmission,transmissionsamples,transmissionRstep
) -> list:
    """
    Compute planetary cutoff grid using the OTSO framework.
    
    Generates a global grid of geomagnetic cutoff rigidities across the Earth's
    surface. Users can also request the computation of asymptotic viewing directions
    and transmission functions by request. Useful for studying global cosmic ray
    accessibility patterns.

    Args:
        cutoff_comp (str): Cutoff computation method ("Vertical", "Apparent", "Custom").

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ADAPTIVE"): Penumbra scanning ("ADAPTIVE"/"ON"/"OFF").
              "ADAPTIVE" extends the fine scan until the penumbra is bounded by clean guard bands;
              "ON" uses a fixed ±0.75 GV window around a 0.25 GV coarse scan; "OFF" scans the full range

        asymptotic_params (AsymptoticParams): Asymptotic computation parameters.

            Available keys:

            - `unit` (`str`, default="GeV"): level unit ("GeV", "GV")
            - `asymptotic` (`str`, default="NO"): Enable asymptotic cone computation ("YES"/"NO")
            - `asymlevels` (`list`, default=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]): level values for asymptotic computation

        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmission` (`bool`, default=False): Bool to enable or disable transmission computation
            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range

        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        grid_params (GridParams): Grid configuration parameters. 
        
            Available keys:

            - `latstep` (`float`, default=-5):  Latitude step size for grid
            - `longstep` (`float`, default=5):  Longitude step size for grid
            - `maxlat` (`float`, default=90):  Maximum latitude for grid
            - `minlat` (`float`, default=-90):  Minimum latitude for grid
            - `maxlong` (`float`, default=360):  Maximum longitude for grid
            - `minlong` (`float`, default=0):  Minimum longitude for grid
            - `array_of_lats_and_longs` (`list`, default=None):  Custom grid points
            - `grid_type` (`str`, default="regular"):  "regular" (latstep/longstep grid), "healpix" (equal-area HEALPix points), "healpix_geomag" (HEALPix rings laid out in centred-dipole coordinates for the date; most accurate per point) or "fibonacci" (golden-angle spiral)
            - `nside` (`int`, default=12):  Resolution for the equal-area grid types: 12 * nside**2 points (12 -> 1728 points, ~550 km spacing, like a 5° grid at the equator)

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Output coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `threadnum` (`int`, default=None): Number of threads per CPU core for Fortran computations
            - `Verbose` (`bool`, default=True): Enable verbose output
            - `delim` (`str`, default=";"): Delimiter for asymptotic output formatting

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [planet_dataframe, asymptotic_dataframe, transmission_dataframe, readme_text]

            - planet_dataframe: Global cutoff rigidity grid
            - asymptotic_dataframe: Global asymptotic viewing directions for the energy/rigidity levels input.
            - transmission_dataframe: Global tranmssion values for set rigidity values.
            - readme_text: OTSO computation summary

    Examples:
    ```python
        from OTSO import planet

        if __name__ == '__main__':
        
            # Basic global grid
            planet_results = planet(
                cutoff_comp="Vertical",
                grid_params={"latstep": -10, "longstep": 15},
                computation_params={"corenum": 8, "threadnum": 1},
                datetime_params={"year": 2000},
                rigidity_params={"rigiditystep": 0.01},
                integration_params={"gyropercent":10}
            )

            
            # Access the results
            planet_df, asymptotic_df, transmission_df, metadata = planet_result
    ```
    """
    
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_grid = {**GridParams.DEFAULTS, **grid_params}
    final_asymptotic = {**AsymptoticParams.DEFAULTS, **asymptotic_params}
    final_transmission = {**TransmissionParams.DEFAULTS, **transmission_params}


    #remove unused variables from planet
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    final_particle.pop("Anum", None)

    planet_args = {
        'cutoff_comp': cutoff_comp,
        'defaulted_inputs': defaulted_inputs,
        **final_solar_wind,
        **final_geomagnetic,
        **final_tsyganenko,
        **final_datetime,
        **final_magfield,
        **final_integration,
        **final_rigidity,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field,
        **final_particle,
        **final_coordinate,
        **final_grid,
        **final_asymptotic,
        **final_transmission
    }
    
    return planet_func(**planet_args)
################################################################################################################
@_check_arguments
def trajectory(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
) -> list:
    """
    Compute cosmic-ray particle trajectories using the OTSO framework.
    
    Traces individual particles through the magnetosphere to determine their
    trajectories from given starting locations and rigidity. Returns a dictionary
    of trajectories and readme metadata for the computations.

    Args:
        Stations (str | list): Station name(s) or identifiers for trajectory calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Launch zenith angle
            - `azimuth` (`float`, default=0): Launch azimuth angle
            - `rigidity` (`float`, default=1): Particle rigidity (GV)

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Output coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `Verbose` (`bool`, default=True): Enable verbose output

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [trajectory_data, readme_text]

            - trajectory_data: Dictionary with positional, velocity, filter and asymptotic viewing directions information for trajectories
              (the "Filter" value gives why the trace stopped: 1 = escaped (magnetopause / outer boundary),
              -1 = reached minaltitude (Earth/atmosphere), -2 = maxdistance, -3 = maxsteps, -4 = maxtime,
              -5 = numerical failure (non-finite positions; the trace is not reliable))
            - readme_text: OTSO computation summary

    Examples:
    ```python
        from OTSO import trajectory

    if __name__ == '__main__':

        stations_list = ["OULU"] # list of neutron monitor stations (using their abbreviations)

        # Example using grouped parameters
        trajectory_results = trajectory(
            Stations=stations_list,
            particle_params={"rigidity": 0.6},
            computation_params={"corenum": 1},
        )

        trajectory_dict, metadata = trajectory_results

        # Accessing the Oulu trajectory value.
        Oulu_trajectory = next(
            (traj for traj in trajectory_dict if traj["station"] == "OULU"),
            None
            )
    ```
    """

    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}

    #remove unused variables from trajectory
    final_coordinate.pop("coordout", None)
    final_computation.pop("delim", None)
    final_rigidity.pop("endrigidity", None)
    final_rigidity.pop("rigidityscan", None)
    final_rigidity.pop("rigiditystep", None)
    final_rigidity.pop("startrigidity", None)

    return trajectory_func(
        Stations,
        customlocations=customlocations,
        # Expanded parameters
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        **final_geomagnetic,
        **final_tsyganenko,
        **final_datetime,
        **final_magfield,
        **final_integration,
        **final_particle,
        **final_rigidity,
        **final_coordinate,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field
    )
#################################################################################################################
@_check_arguments
def flight(
    latitudes: Sequence[float],
    longitudes: Sequence[float],
    dates: Sequence,
    altitudes: Sequence[float],
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped (optional - can be None for automatic retrieval)
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped (optional - can be None for automatic retrieval)
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped (optional - can be None for automatic retrieval)
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # asymptotic parameters grouped
    asymptotic_params: AsymptoticParams = {},
    # transmission parameters grouped
    transmission_params: TransmissionParams = {}, # transmission,transmissionsamples,transmissionRstep
) -> list:
    """
    Compute cosmic-ray cutoff rigidities along a flight path using OTSO.
    
    Calculates cutoff rigidities at specified time-varying locations, typically
    used for aircraft or satellite trajectory analysis. Supports automatic
    space weather data retrieval based on flight times. On request the asymptotic
    viewing directions and transmission functions can be computed as well.

    Args:
        latitudes (list): Latitude coordinates along flight path.
        longitudes (list): Longitude coordinates along flight path.
        dates (list): Date/time stamps for each location (datetime objects).
        altitudes (list): Altitude coordinates in km along flight path.
        cutoff_comp (str): Cutoff computation method ("Vertical", "Apparent", "Custom").

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ADAPTIVE"): Penumbra scanning ("ADAPTIVE"/"ON"/"OFF").
              "ADAPTIVE" extends the fine scan until the penumbra is bounded by clean guard bands;
              "ON" uses a fixed ±0.75 GV window around a 0.25 GV coarse scan; "OFF" scans the full range

        asymptotic_params (AsymptoticParams): Asymptotic computation parameters.

            Available keys:

            - `unit` (`str`, default="GeV"): level unit ("GeV", "GV")
            - `asymptotic` (`str`, default="NO"): Enable asymptotic cone computation ("YES"/"NO")
            - `asymlevels` (`list`, default=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]): level values for asymptotic computation
        
        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmission` (`bool`, default=False): Bool to enable or disable transmission computation
            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range
        
        solar_wind_params (SolarWindParams): Solar wind parameters (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            All values should be provided as lists of floats, one per flight point.

            Available keys:

            - `vx` (list of float, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (list of float, default=0): Solar wind velocity y-component (km/s)
            - `vz` (list of float, default=0): Solar wind velocity z-component (km/s)
            - `bx` (list of float, default=0.11): IMF x-component (nT)
            - `by` (list of float, default=-0.16): IMF y-component (nT)
            - `bz` (list of float, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (list of float, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (list of float, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (list of float,, default=-4): Dst index (nT)
            - `kp` (list of float,, default=0.7): Kp index (0-9)
            - `n_index` (list of float,, default=0.145): Newell coupling function
            - `b_index` (list of float,, default=0.032): Boynton coupling function
            - `sym_h_corrected` (list of float,, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:
            
            - `G1` (list of float,, default=0.38): Tsyganenko G1 coefficient
            - `G2` (list of float,, default=0.12): Tsyganenko G2 coefficient
            - `G3` (list of float,, default=0.47): Tsyganenko G3 coefficient
            - `W1` (list of float,, default=0.07): Tsyganenko W1 coefficient
            - `W2` (list of float,, default=0.06): Tsyganenko W2 coefficient
            - `W3` (list of float,, default=0.13): Tsyganenko W3 coefficient
            - `W4` (list of float,, default=0.03): Tsyganenko W4 coefficient
            - `W5` (list of float,, default=0.05): Tsyganenko W5 coefficient
            - `W6` (list of float,, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation


        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Output coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `Verbose` (`bool`, default=True): Enable verbose output
            - `delim` (`str`, default=";"): Delimiter for asymptotic output formatting

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [flight_df, asymptotic_df, transmission_df, readme_text, input_dataframe]

            - flight_df: Cutoff rigidities along flight path
            - asymptotic_df: Asymptotic viewing directions for input energy/rigidiy levels along flight path
            - transmission_df: Transmission functions as a function of rigidity along the flight path
            - readme_text: OTSO computation summary
            - input_dataframe: Flight path input data

    Examples:
    ```python
        from OTSO import flight
        import datetime

        if __name__ == "__main__":

            latitude_list = [10, 15, 20, 25, 30]
            longitude_list = [10, 15, 20, 25, 30]
            altitude_list = [30, 40, 50, 60, 80]
            date_list = [
                datetime.datetime(2000, 10, 12, 8),
                datetime.datetime(2000, 10, 12, 9),
                datetime.datetime(2000, 10, 12, 10),
                datetime.datetime(2000, 10, 12, 11),
                datetime.datetime(2000, 10, 12, 12),
            ]

            # Example using grouped parameters
            flight_results = flight(
                latitudes=latitude_list,
                longitudes=longitude_list,
                dates=date_list,
                altitudes=altitude_list,
                cutoff_comp="Vertical",
                computation_params={"corenum": 1, "threadnum": 8},
                asymptotic_params={"asymptotic": "YES"},
                transmission_params={"transmission": True},
                integration_params={"gyropercent": 1}
            )
        
            # Access the results
            flight_df, asymptotic_df, transmission_df, metadata, input_df = flight_result
    ```
    """
    
    # For flight function, we need special handling since some parameters might be None
    # for automatic data retrieval. Only merge with defaults if user provided values.
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params} if solar_wind_params else {}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params} if geomagnetic_params else {}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params} if tsyganenko_params else {}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_asymptotic = {**AsymptoticParams.DEFAULTS, **asymptotic_params}
    final_transmission = {**TransmissionParams.DEFAULTS, **transmission_params}

    #remove unused variables from trajectory
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    final_integration.pop("startaltitude", None)

    
    # Build arguments with conditional inclusion for optional solar wind parameters
    flight_args = {
        'latitudes': latitudes,
        'longitudes': longitudes,
        'dates': dates,
        'altitudes': altitudes,
        'cutoff_comp': cutoff_comp,
        **final_magfield,
        **final_integration,
        **final_particle,
        **final_rigidity,
        **final_coordinate,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field,
        **final_asymptotic,
        **final_transmission
    }
    
    # Only include solar wind/geomagnetic parameters if user provided them
    if final_solar_wind:
        flight_args.update(final_solar_wind)
    if final_geomagnetic:
        flight_args.update(final_geomagnetic)
    if final_tsyganenko:
        flight_args.update(final_tsyganenko)
    
    return flight_func(**flight_args)

#################################################################################################################
@_check_arguments
def trace(
    coordsys: str = "GEO",
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # integration parameters grouped (only minaltitude,maxdistance,maxtime,maxsteps,startaltitude,fixedstep apply to trace)
    integration_params: IntegrationParams = {},  # minaltitude,maxdistance,maxtime,maxsteps,startaltitude,fixedstep
    # coordinate parameters grouped (only inputcoord applies to trace)
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # grid parameters grouped (not used in trace but included for consistency)
    grid_params: GridParams = {},  # latstep,longstep,maxlat,minlat,maxlong,minlong
) -> list:
    """
    Trace magnetic field lines using the OTSO framework.
    
    Traces magnetic field lines across a global grid to visualize
    the magnetosphere structure. Useful for understanding magnetic connectivity
    and field line topology. Computation also returns the McIlwain L value for the
    field line and invariant latitude for the origin location of the field line. Users
    should be cautious when interpreting low- or high-latitude origin invariant latitudes 
    as they are ill-defined in such locations.

    Args:
        Coordsys (str): Coordinate system for field line positions.

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        grid_params (GridParams): Grid configuration parameters. 
        
            Available keys:

            - `latstep` (`float`, default=-5):  Latitude step size for grid
            - `longstep` (`float`, default=5):  Longitude step size for grid
            - `maxlat` (`float`, default=90):  Maximum latitude for grid
            - `minlat` (`float`, default=-90):  Minimum latitude for grid
            - `maxlong` (`float`, default=360):  Maximum longitude for grid
            - `minlong` (`float`, default=0):  Minimum longitude for grid

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `Verbose` (`bool`, default=True): Enable verbose output

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.
        
            Available keys:
            
            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [trace_data, readme_text]

            - trace_data: Dictionary with magnetic field line positions
            - readme_text: OTSO computation summary

    Examples:
    ```python
       from OTSO import trace

       if __name__ == '__main__':

           # Example using grouped parameters
           trace_results = trace(
               computation_params={"corenum": 4},
               magfield_params={"externalmag":"TSY89c"},
               grid_params={"latstep": -30, "longstep": 180, "maxlat": 60, "minlat": 0},
           )
        
            # Access the results
            trace_dict, metadata = trace_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_grid = {**GridParams.DEFAULTS, **grid_params}
    
    # Remove parameters not supported by trace function
    trace_magfield = {k: v for k, v in final_magfield.items() if k != 'AdaptiveExternalModel'}

    #remove unused variables from trajectory
    final_coordinate.pop("coordout", None)
    final_coordinate.pop("coordsystem", None)
    final_integration.pop("adaptivestep", None)
    final_integration.pop("betaerror", None)
    final_integration.pop("gyropercent", None)
    final_integration.pop("mintrapdist", None)
    final_integration.pop("totalbetacheck", None)
    final_integration.pop("intmodel", None)
    final_computation.pop("threadnum", None)
    final_computation.pop("delim", None)
    final_grid.pop("array_of_lats_and_longs", None)
    final_grid.pop("grid_type", None)
    final_grid.pop("nside", None)


    return trace_func(
        coordsys=coordsys,
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        **final_geomagnetic,
        **final_tsyganenko,
        **final_datetime,
        **trace_magfield,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field,
        **final_integration,
        **final_coordinate,
        **final_grid
    )
#################################################################################################################
@_check_arguments
def magfield(
    Locations: Sequence[float],
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
) -> list:
    """
    Compute magnetic field vectors at specified locations using OTSO models.
    
    Evaluates the total magnetic field (internal + external) at given locations
    using various geomagnetic field models.

    Args:
        Locations (list): Locations as [[x, y, z]] in specified coordinate system.

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        solar_wind_params (SolarWindParams): Solar wind parameters (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval). Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"
            - `coordout` (`str`, default="GSM"): Output coordinate system, cartesian only; "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `Verbose` (`bool`, default=True): Enable verbose output

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [magfield_dataframe, readme_text]

            - magfield_dataframe: Magnetic field vectors (components in `coordout`) at input locations. The index is each point's position in `Locations` (rows are sorted by position; use `.sort_index()` to restore input order)
            - readme_text: OTSO computation summary

    Examples:
    ```python
       from OTSO import magfield
 
       if __name__ == '__main__':
 
            location_list = [[10,10,10]] # [[X,Y,Z]] Earth radii Geocentric coordinates in this instance
 
            # Example using grouped parameters
            magfield_results = magfield(
               Locations=location_list,
               coordinate_params={"inputcoord": "GEO"},
               computation_params={"corenum": 1}
            )

            # Access the results
            field_df, metadata = field_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}

    #remove unused variables from planet
    final_coordinate.pop("coordsystem", None)
    final_computation.pop("threadnum", None)
    final_computation.pop("delim", None)
    
    # Remove parameters not supported by magfield function
    magfield_params_filtered = {k: v for k, v in final_magfield.items() 
                               if k not in ['magnetopause', 'spheresize', 'AdaptiveExternalModel']}
    
    return magfield_func(
        Locations,
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        **final_geomagnetic,
        **final_tsyganenko,
        **final_datetime,
        **magfield_params_filtered,
        **final_coordinate,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field
    )
############################################################################################################
#######################################################################################################################################
@_check_arguments
def transmission(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # transmission parameters grouped
    transmission_params: TransmissionParams = {} # transmission,transmissionsamples,transmissionRstep
) -> list:
    """
    Compute transmission function values as a function of particle rigidity 
    for given neutron monitor stations, or user-defined locations.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the transmission 
    function values and related metadata.

    Args:
        Stations (str | list): Station name(s) or identifiers used for cone calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)

        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `threadnum` (`int`, default=None): Number of threads per CPU core for Fortran computations
            - `Verbose` (`bool`, default=True): Enable verbose output

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [transmission_df, readme_text]

            - transmission_df: Transmission functions as a function of rigidity in a dataframe
            - readme_text: OTSO computation summary

    Examples:
    ```python
       from OTSO import transmission
 
       if __name__ == '__main__':
 
            stations_list = ["ROME"] # list of neutron monitor stations (using their abbreviations)
 
            transmission_results = transmission(Stations=stations_list,
                               computation_params={"corenum": 1, "threadnum": 8},
                               datetime_params={"year": 2005, "month": 5, "day": 1, "hour": 0},
                               integration_params={"gyropercent": 10},
                               magfield_params={"internalmag": "IGRF", "externalmag": "NONE"},
                               rigidity_params={"rigiditystep":0.01},
                               transmission_params={"transmission": True, "transmissionRstep": 0.0001, "transmissionsamples": 25})
        
            # Access the results
            transmission_df, metadata = transmission_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_transmission = {**TransmissionParams.DEFAULTS, **transmission_params}

    #remove unused variables from cone
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    final_rigidity.pop("rigidityscan", None)
    final_coordinate.pop("coordsystem", None)

    
    # Call the original function with expanded parameters
    return transmission_func(
        Stations,
        customlocations=customlocations,
        # Solar wind parameters
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        # Geomagnetic parameters
        **final_geomagnetic,
        # Tsyganenko coefficients  
        **final_tsyganenko,
        # Date/time parameters
        **final_datetime,
        # Magnetic field parameters
        **final_magfield,
        # Integration parameters (filtered)
        **final_integration,
        # Particle parameters
        **final_particle,
        # Rigidity parameters (filtered)
        **final_rigidity,
        # Coordinate parameters
        **final_coordinate,
        # Computation parameters
        **final_computation,
        # Data retrieval parameters
        **final_data_retrieval,
        # Custom field parameters
        **final_custom_field,
        # Transmission parameters
        **final_transmission
    )
#########################################################################################################################
###############################################################################################################################
@_check_arguments
def skymap(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind_params: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic_params: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko_params: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
    # Date/time parameters grouped
    datetime_params: DateTimeParams = {},  # year,month,day,hour,minute,second
    # Magnetic field model parameters grouped
    magfield_params: MagFieldParams = {},  # internalmag,externalmag,boberg,magnetopause,etc
    # integration parameters grouped
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude,maxdistance,etc
    # particle parameters grouped
    particle_params: ParticleParams = {},  # Anum,anti
    # rigidity parameters grouped
    rigidity_params: RigidityParams = {},  # startrigidity,endrigidity,rigiditystep,rigidityscan
    # coordinate parameters grouped
    coordinate_params: CoordinateParams = {},  # coordsystem,inputcoord
    # computation parameters grouped
    computation_params: ComputationParams = {},  # corenum,Verbose
    # data retrieval parameters grouped
    data_retrieval_params: DataRetrievalParams = {},  # serverdata,livedata
    # custom field parameters grouped
    custom_field_params: CustomFieldParams = {},  # g,h,MHDfile,MHDcoordsys
    # skymap parameters grouped
    skymap_params: SkymapParams = {},  # 

) -> list:
    """
    Compute geomagnetic cutoff rigidities over a range of zenith and azimuth
    combinations for given neutron monitor stations, or user-defined locations.
    Useful for showing the East-West asymmetry of cosmic-ray arrival.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the cutoff
    rigidities as functions of zenith and azimuth, and related metadata.

    Args:
        Stations (str | list): Station name(s) or identifiers used for cutoff calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].
        cutoff_comp (str): Cutoff computation method ("Vertical", "Apparent", "Custom").

        datetime_params (DateTimeParams): Date/time parameters.

            Available keys:

            - `year` (`int`, default=2024): Year (e.g., 2023)
            - `month` (`int`, default=1): Month (1–12)
            - `day` (`int`, default=1): Day (1–31)
            - `hour` (`int`, default=12): Hour (0–23)
            - `minute` (`int`, default=0): Minute (0–59)
            - `second` (`int`, default=0): Second (0–59)

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss", "CHAOS"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY87short", "TSY87long", "TSY89a", "TSY96", "TSY01", "TSY01S", "TSY04", "TSY89c", "TSY15N", "TSY15B", "TA16_RBF", "TSY89_refit", "MHD"
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): "EXTENSION", "CONTINUOUS", "DST_DEPENDENT", "DST_MIDPOINT"
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere", "aFormisano"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
            - `optimise_tsy` (`bool`, default=False): Use the faster implementation of the selected Tsyganenko model (affects externalmag="TSY01", "TSY01S", "TSY04", or "TSY96"; other models are unaffected by this flag). Both settings give numerically identical results (all known bugs are fixed in both); True is just faster

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ADAPTIVE"): Penumbra scanning ("ADAPTIVE"/"ON"/"OFF").
              "ADAPTIVE" extends the fine scan until the penumbra is bounded by clean guard bands;
              "ON" uses a fixed ±0.75 GV window around a 0.25 GV coarse scan; "OFF" scans the full range

        skymap_params (SkymapParams): Sky map zenith and azimuth resolution parameters

            Available keys:

            - `zenithstep` (`float`, default=15): Zenith resolution
            - `azimuthstep` (`float`, default=45): Azimuth resolution
            - `maxzenith` (`float`, default=75):  Maximum deviation from the zenith
            - `minzenith` (`float`, default=0):  Minimum deviation from the zenith
            - `maxazimuth` (`float`, default=360): Maximum azimuth value
            - `minazimuth` (`float`, default=0): Minimum azimuth value
        
        solar_wind_params (SolarWindParams): Solar wind parameters. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `vx` (`float`, default=-359.4): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0.11): IMF x-component (nT)
            - `by` (`float`, default=-0.16): IMF y-component (nT)
            - `bz` (`float`, default=0.58): IMF z-component (nT)
            - `by_avg` (`float`, default=-0.17): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0.59): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=4.94): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=1.31): Solar wind dynamic pressure (nPa); computed from density and velocity when 0, or when density or velocity is given without pdyn
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `Dst` (`float`, default=-4): Dst index (nT)
            - `kp` (`float`, default=0.7): Kp index (0-9)
            - `n_index` (`float`, default=0.145): Newell coupling function
            - `b_index` (`float`, default=0.032): Boynton coupling function
            - `sym_h_corrected` (`float`, default=-18.07): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients. Defaults are quiet-time medians (Kp <= 1) of OMNI data, 1996-2025.

            Available keys:

            - `G1` (`float`, default=0.38): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0.12): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0.47): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0.07): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0.06): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0.13): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0.03): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0.05): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0.08): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
            - `gyropercent` (`float`, default=1): Step size as a percentage of the gyration period. With adaptivestep=True (default) it caps each step using the local field; with adaptivestep=False the step is fixed for the whole trace at this percentage of the gyration period at the launch point
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum path length travelled (Re), not distance from Earth; 0 disables this check (use 0 for long trapped-particle runs: a 0.6c proton travels 100 Re in ~4 s)
            - `maxtime` (`float`, default=0): Maximum flight time in seconds (lab frame); 0 = no limit
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=0): Maximum number of integration steps; 0 = no limit

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (-1=muon, 0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `inputcoord` (`str`, default="GDZ"): Input coordinate system; "GDZ", "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL"

        computation_params (ComputationParams): Computation settings.

            Available keys:

            - `corenum` (`int`, default=None): Number of CPU cores for multicore processing
            - `threadnum` (`int`, default=None): Number of threads per CPU core for Fortran computations
            - `Verbose` (`bool`, default=True): Enable verbose output

        data_retrieval_params (DataRetrievalParams): Data retrieval.

            Available keys:

            - `serverdata` (`str`, default="OFF"): Server data retrieval from OMNI
            - `livedata` (`str`, default="OFF"): real-time data retrieval from NOAA

        custom_field_params (CustomFieldParams): Custom fields.

            Available keys:

            - `g` (`list`, default=None): Gauss g coefficients
            - `h` (`list`, default=None): Gauss h coefficients
            - `max_degree` (`int`, default=13): Max degree of spherical harmonic expansion
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system
            - `MHDgridtype` (`str`, default="auto"): How the MHD grid is looked up - "auto" detects per-axis spacing automatically, "uniform" forces the fast fixed-spacing path (only correct if every axis really is evenly spaced), "stretched" forces the general path needed for a grid whose resolution varies with position (e.g. finer near Earth)
            - `MHDinterpolation` (`str`, default="trilinear"): How the field is reconstructed between grid points - "trilinear" is the classic fast corner-blend (kinked gradients at cell boundaries), "tricubic" is a smooth C^1 fit that's far more accurate in smooth field regions but can overshoot near sharp gradients/boundaries, "monotonic" is the same tricubic fit with overshoot limiting (recommended for general use), "divfree" returns the curl of an interpolated vector potential so the field is exactly divergence-free (for long quasi-trapped trajectories; less point-wise accurate than "tricubic"; the gridded field must itself be divergence-free, e.g. an external-only field - a warning is issued otherwise)

    Returns:
        list: [skymap_dataframe, readme_text]

            - skymap_dataframe: Cutoff data values per input station and zenith and azimuth combination in a pandas dataframe. (Ru, Rc, Rl, PTF)
            - readme_text: OTSO computation summary

    Examples:
    ```python
        import OTSO

        if __name__ == "__main__":
        
            stations = ["ROME"] # list of neutron monitor stations (using their abbreviations)
        
            # Using parameter groups
            skymap_results = skymap(
                Stations=stations,
                computation_params={"corenum": 1,"threadnum": 8},
                datetime_params={"year": 2005,"month": 5,"day": 1,"hour": 0},
                integration_params={"gyropercent": 10},
                magfield_params={"internalmag": "IGRF","externalmag": "NONE"},
                rigidity_params={"rigiditystep": 0.01},
                skymap_params={"zenithstep": 5,"azimuthstep": 15}
            )
            
            # Access the results
            skymap_data, metadata = result
    ```
    """

    if rigidity_params.get("startrigidity") is None:
        rigidity_params["startrigidity"] = 50

    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind, final_geomagnetic, final_tsyganenko, defaulted_inputs = _merge_wind_inputs(
        solar_wind_params, geomagnetic_params, tsyganenko_params, data_retrieval_params)
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_skymap = {**SkymapParams.DEFAULTS, **skymap_params}

    #remove unused variables from cutoff
    final_coordinate.pop("coordout", None)
    final_particle.pop("rigidity", None)
    final_computation.pop("delim", None)
    
    # Call the original function with expanded parameters
    return skymap_func(
        Stations,
        customlocations=customlocations,
        # Solar wind parameters
        defaulted_inputs=defaulted_inputs,
        **final_solar_wind,
        # Geomagnetic parameters
        **final_geomagnetic,
        # Tsyganenko coefficients  
        **final_tsyganenko,
        # Date/time parameters
        **final_datetime,
        # Magnetic field parameters
        **final_magfield,
        # Integration parameters
        **final_integration,
        # Particle parameters
        **final_particle,
        # Rigidity parameters
        **final_rigidity,
        # Coordinate parameters
        **final_coordinate,
        # Computation parameters
        **final_computation,
        # Data retrieval parameters
        **final_data_retrieval,
        # Custom field parameters
        **final_custom_field,
        # skymap parameters
        **final_skymap
    )

#######################################################################################################################################

@_check_arguments
def coordtrans(
    Locations: Sequence[float],
    dates: Sequence,
    CoordIN: str = "GEO",
    CoordOUT: str = "GDZ",
    corenum: Optional[int] = None,
    Verbose: bool = True,
) -> list:
    """
    Transform coordinates between different coordinate systems using OTSO.
    
    Converts spatial coordinates between various reference frames used in
    space physics and geomagnetics. Supports time-dependent transformations
    for date-specific coordinate system orientations.

    Coordinate transformations are done via the IRBEM library.

    Args:
        Locations (list): Input coordinates as [[coord1, coord2, coord3]].
        dates (list): Date/time stamps for coordinate transformations.
        CoordIN (str): Input coordinate system.
            Options: "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL", "GDZ"
        CoordOUT (str): Output coordinate system.
            Options: "GEO", "GSM", "GSE", "SM", "GEI", "MAG", "SPH", "RLL", "GDZ"
        corenum (int, optional): Number of CPU cores for parallel processing.
        Verbose (bool): Enable verbose output.

    Returns:
        list: [coords_dataframe, summary_text]

            - coords_dataframe: Transformed coordinates. The index is each point's position in `Locations` (rows are sorted by date and position; use `.sort_index()` to restore input order)
            - summary_text: Transformation summary

    Examples:
    ```python
       from OTSO import coordtrans
        import datetime

        if __name__ == '__main__':

            lat_lon_alt_list = [[10,10,10]] # [[Latitude,Longitude,Altitude]]
            date_list = [datetime.datetime(2000,10,12,8)] # [dates]
            
            # Example using grouped parameters (coordtrans has fewer parameter groups)
            Coords = coordtrans(
                Locations=lat_lon_alt_list,
                dates=date_list,
                CoordIN="GEO",
                CoordOUT="GSM",
                corenum=1  # coordtrans uses individual parameters, not grouped ones
            )
        
            # Access the results
            coord_df, metadata = coord_result
    ```
    """
    return coordtrans_func(
        Locations = Locations, 
        dates = dates, 
        CoordIN = CoordIN, 
        CoordOUT = CoordOUT, 
        corenum = corenum, 
        Verbose = Verbose
    )
#######################################################################################################################
@_check_arguments
def clean():
    """
    Clean OTSO-generated files and temporary data.
    
    Removes temporary files, cached data, and output files created by OTSO
    calculations to free up disk space and reset the working environment.
    
    Usage:
        OTSO.clean()
        
    Can also be used as CLI command: OTSO.clean
    """
    from .otso_cli import clean as clean_func
    return clean_func()

@_check_arguments
def addstation(Name: str, Latitude: float, Longitude: float):
    """
    Add a new neutron monitor station to the OTSO database.
    
    Registers a new station location for use in OTSO calculations. Useful when
    new neutron monitor stations are established or for adding custom locations
    with specific names for repeated use.
    
    Args:
        Name (str): Name identifier for the new station.
        Latitude (float): Geographic latitude in degrees.
        Longitude (float): Geographic longitude in degrees.
    
    Usage:
        OTSO.addstation("NEWSTATION", 65.0, 25.0)
        
    CLI usage:
        OTSO.addstation NEWSTATION 65.0 25.0
        
    Note: If the station already exists, there is an option to overwrite it.
    """
    from .otso_cli import AddStation as addstation_func
    addstation_func(Name, float(Latitude), float(Longitude))
    return

@_check_arguments
def removestation(Name: str):
    """
    Remove a station from the OTSO database.
    
    Deletes a station entry from the OTSO database. Useful for correcting
    incorrectly entered station data or removing obsolete stations.
    
    Args:
        Name (str): Name of station to remove.
    
    Usage:
        OTSO.removestation("OLDSTATION")
        
    CLI usage:
        OTSO.removestation OLDSTATION
    """
    from .otso_cli import RemoveStation as removestation_func
    removestation_func(Name)
    return

@_check_arguments
def liststations():
    """
    List all available neutron monitor stations in the OTSO database.
    
    Displays all currently registered stations with their coordinates.
    Useful for finding available station names and verifying station data.
    
    Usage:
        OTSO.liststations()
        
    CLI usage:
        OTSO.liststations
        
    Returns:
        List of all available stations with coordinates.
    """
    from .otso_cli import ListStations as liststations_func
    return liststations_func()

@_check_arguments
def IGRFupdate():
    """
    CLI function to update the IGRF model data. Running this CLI function will download the latest
    IGRF coefficients from the official IGRF website and update the local OTSO database.
    You can specify older IGRF models by providing the desired model version as an argument 
    (e.g., 13 for IGRF-13).
    
    example usage:
        OTSO.IGRFupdate \n
        OTSO.IGRFupdate 13
    """
    from .otso_cli import IGRFupdate as IGRFupdate_func
    return IGRFupdate_func()

@_check_arguments
def serverdownload():
    """
    Download space physics data from online servers for offline use.
    
    Fetches essential solar wind and geomagnetic parameters from authoritative
    databases including NOAA Space Weather Prediction Center, NASA's Goddard
    Space Flight Center, and the World Data Center for Geomagnetism in Kyoto.
    
    Downloaded data includes:
        - Solar wind parameters (velocity, density, temperature, magnetic field)
        - Geomagnetic indices (Kp, Ap, Dst, AE, F10.7)
        - Interplanetary magnetic field conditions
        - Real-time and historical space weather data
    
    Usage:
        OTSO.serverdownload()
        
    CLI usage:
        OTSO.serverdownload
        
    Note: Requires internet connection. Downloaded data is cached locally
    to enable OTSO calculations in offline environments. This is particularly
    useful for field work or when internet connectivity is limited.
    
    Data Sources:
        - NOAA Space Weather Prediction Center
        - NASA/Goddard OMNI database 
        - World Data Center for Geomagnetism, Kyoto
        - Real-time space weather services
    """
    from .otso_cli import serverdownload as serverdownload_func
    return serverdownload_func()

@_check_arguments
def chaosdownload():
    """
    Download CHAOS files for offline use.
    
    Fetches essential mat files for chaosmagpy to use when extracting CHAOS gaussian
    coefficients. Will pull the most recent version from the zenodo repo and overwrite
    exisiting file if the same version.
    
    Usage:
        OTSO.chaosdownload()
        
    CLI usage:
        OTSO.chaosdownload
        
    Note: Requires internet connection. Downloaded data is cached locally
    to enable OTSO calculations in offline environments.

    Data Sources:
        - DTU, Denmark
    """
    from .otso_cli import CHAOSdownload as CHAOSdownload_func
    return CHAOSdownload_func()


__all__ = [
    "cutoff",
    "cone",
    "planet",
    "trajectory",
    "flight",
    "magfield",
    "coordtrans",
    "trace",
]