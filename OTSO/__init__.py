"""
OTSO (Open-source geomagneToSphere prOpagation tool)
============================================

A comprehensive Python package for charged particle trajectory and 
geomagnetic cut-off computations in the Earth’s magnetosphere.
 Main features include:

- Cosmic ray geomagnetic cut-off rigidity calculations
- Asymptotic cone of acceptance computations  
- Magnetic field line tracing
- Coordinate system transformations
- Particle trajectory tracing in geomagnetic fields

For detailed documentation and example uses visit: https://github.com/NLarsen15/OTSOpy
"""

from __future__ import annotations
import sys
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
    AsymptoticParams
)

from ._cutoff import cutoff as cutoff_func
from ._cone import cone as cone_func
from ._planet import planet as planet_func
from ._trajectory import trajectory as trajectory_func
from ._flight import flight as flight_func
from ._magfield import magfield as magfield_func
from ._coordtrans import coordtrans as coordtrans_func
from ._trace import trace as trace_func


def cutoff(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
) -> list:
    """
    Compute geomagnetic cutoff rigidities for given neutron monitor stations,
    or user-defined locations.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the cutoff
    rigidities and related metadata.

    Args:
        Stations (str | list): Station name(s) or identifiers used for cutoff calculations.
        customlocations (list, optional): Custom locations as [["NAME", lat, lon]].
        cutoff_comp (str): Cutoff computation method ("Vertical", "Custom", etc.).

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

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")
        
        solar_wind (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients.

            Available keys:

            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris"): "Boris", "4RK", "Vay", "HC"
            - `gyropercent` (`float`, default=15): Gyration period percentage
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum distance (Re)
            - `maxtime` (`float`, default=0): Maximum time
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=None): Maximum number of integration steps

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Coordinate system used for calculations
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Examples:
    ```python
        import OTSO
        
        # Using parameter groups
        result = OTSO.cutoff(["OULU"], 
                            datetime_params={"year": 2023, "month": 6},
                            solar_wind={"vx": -400, "by": 3.0},
                            rigidity_params={"startrigidity": 15})
        
        # Access the results
        cutoff_data, metadata = result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
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
    
    # Call the original function with expanded parameters
    return cutoff_func(
        Stations,
        customlocations=customlocations,
        cutoff_comp=cutoff_comp,
        # Solar wind parameters
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
        **final_asymptotic
    )

def cone(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
    # Individual parameters for backward compatibility
) -> list:
    """
    Compute asymptotic viewing cones for given neutron monitor stations,
    or user-defined locations. Also computes cutoff rigidities.

    Upon calling this function, OTSO will perform particle tracing
    simulations based on the specified parameters, returning the asymptotic 
    viewing cones, rigidities, and related metadata.

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

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")
        
        solar_wind (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients.

            Available keys:

            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris"): "Boris", "4RK", "Vay", "HC"
            - `gyropercent` (`float`, default=15): Gyration period percentage
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum distance (Re)
            - `maxtime` (`float`, default=0): Maximum time
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=None): Maximum number of integration steps

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Coordinate system used for calculations
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [cone_dataframe, cutoff_dataframe, summary_text]
            - cone_dataframe: Asymptotic directions per rigidity (filter;lat;lon format)
            - cutoff_dataframe: Cutoff rigidity table (Ru, Rc, Rl)
            - summary_text: Input parameter summary

    Examples:
    ```python
        import OTSO
        
        # Basic usage with parameter groups
        cone_result = OTSO.cone(["OULU", "ROME"], 
                               computation_params={"corenum": 4},
                               datetime_params={"year": 2023, "month": 6})
        
        # With solar wind conditions
        cone_data = OTSO.cone(["DOMC"], 
                             solar_wind={"vx": -400, "by": 3.0, "bz": -2.0},
                             rigidity_params={"startrigidity": 15, "rigiditystep": 0.1})
        
        # Access the results
        cone_df, cutoff_df, metadata = cone_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}

    
    # Call the original function with expanded parameters
    return cone_func(
        Stations,
        customlocations=customlocations,
        # Solar wind parameters
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
        **final_custom_field
    )

def planet(
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
):
    """
    Compute planetary cutoff grid using the OTSO framework.
    
    Generates a global grid of geomagnetic cutoff rigidities across the Earth's
    surface. Useful for studying global cosmic ray accessibility patterns.

    Args:
        cutoff_comp (str): Cutoff computation method ("Vertical", "Custom", etc.).

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

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")
        
        solar_wind (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients.

            Available keys:

            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris"): "Boris", "4RK", "Vay", "HC"
            - `gyropercent` (`float`, default=15): Gyration period percentage
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum distance (Re)
            - `maxtime` (`float`, default=0): Maximum time
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=None): Maximum number of integration steps

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation

        asymptotic_params (AsymptoticParams): Asymptotic computation parameters.

            Available keys:

            - `unit` (`str`, default="GeV"): Rigidity unit ("GeV", "GV")
            - `asymptotic` (`str`, default="NO"): Enable asymptotic cone computation ("YES"/"NO")
            - `asymlevels` (`list`, default=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]): Rigidity levels for asymptotic computation

        grid_params (GridParams): Grid configuration parameters. 
        
            Available keys:

            - `latstep` (`float`, default=-5):  Latitude step size for grid
            - `longstep` (`float`, default=5):  Longitude step size for grid
            - `maxlat` (`float`, default=90):  Maximum latitude for grid
            - `minlat` (`float`, default=-90):  Minimum latitude for grid
            - `maxlong` (`float`, default=360):  Maximum longitude for grid
            - `minlong` (`float`, default=0):  Minimum longitude for grid
            - `array_of_lats_and_longs` (`list`, default=None):  Custom grid points

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Coordinate system used for calculations
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [planet_dataframe, summary_text]
            - planet_dataframe: Global cutoff rigidity grid
            - summary_text: Input parameter summary

    Examples:
    ```python
        import OTSO
        
        # Basic global grid
        planet_result = OTSO.planet(
            computation_params={"corenum": 4},
            grid_params={"latstep": 10, "longstep": 20}
        )
        
        # High-resolution polar region
        polar_data = OTSO.planet(
            grid_params={"latstep": 2, "longstep": 5, "maxlat": 90, "minlat": 60},
            datetime_params={"year": 2023, "month": 3},
            rigidity_params={"rigiditystep": 0.05}
        )
        
        # Access the results
        planet_df, metadata = planet_result
    ```
    """
    
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
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
    
    # Map coordinate system parameter correctly for planet function
    planet_args = {
        'cutoff_comp': cutoff_comp,
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
        **final_asymptotic
    }
    
    return planet_func(**planet_args)

def trajectory(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
    *args, **kwargs
):
    """
    Compute cosmic-ray particle trajectories using the OTSO framework.
    
    Traces individual particles through the magnetosphere to determine their
    trajectories from given starting locations and rigidity.

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

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")
        
        solar_wind (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients.

            Available keys:

            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris"): "Boris", "4RK", "Vay", "HC"
            - `gyropercent` (`float`, default=15): Gyration period percentage
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum distance (Re)
            - `maxtime` (`float`, default=0): Maximum time
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=None): Maximum number of integration steps

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Launch zenith angle
            - `azimuth` (`float`, default=0): Launch azimuth angle

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Coordinate system used for calculations
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [trajectory_data, summary_text]
            - trajectory_data: Dictionary with positional information for trajectories
            - summary_text: Input parameter summary

    Examples:
    ```python
        import OTSO
        
        # Single particle trajectory
        traj_result = OTSO.trajectory(["OULU"], 
                                     rigidity_params={"startrigidity": 5.0},
                                     computation_params={"corenum": 2})
        
        # With specific magnetic field conditions
        traj_data = OTSO.trajectory(["DOMC"], 
                                   rigidity_params={"startrigidity": 10.0},
                                   datetime_params={"year": 2023, "month": 3},
                                   magfield_params={"externalmag": "TSY01"})
        
        # Access the results
        trajectory_dict, metadata = traj_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}

    
    return trajectory_func(
        Stations,
        customlocations=customlocations,
        # Expanded parameters
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

def flight(
    latitudes: Sequence[float],
    longitudes: Sequence[float],
    dates: Sequence,
    altitudes: Sequence[float],
    cutoff_comp: str = "Vertical",
    # Solar wind parameters grouped (optional - can be None for automatic retrieval)
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped (optional - can be None for automatic retrieval)
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped (optional - can be None for automatic retrieval)
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
    *args, **kwargs
):
    """
    Compute cosmic-ray cutoff rigidities along a flight path using OTSO.
    
    Calculates cutoff rigidities at specified time-varying locations, typically
    used for aircraft or satellite trajectory analysis. Supports automatic
    space weather data retrieval based on flight times.

    Args:
        latitudes (list): Latitude coordinates along flight path.
        longitudes (list): Longitude coordinates along flight path.
        dates (list): Date/time stamps for each location (datetime objects).
        altitudes (list): Altitude coordinates in km along flight path.
        cutoff_comp (str): Cutoff computation method ("Vertical", "Custom").

        magfield_params (MagFieldParams): Magnetic field models.

            Available keys:

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")
        
        solar_wind (SolarWindParams): Solar wind parameters (optional for auto-retrieval).

            All values should be provided as lists of floats, one per flight point.

            Available keys:

            - `vx` (list of float, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (list of float, default=0): Solar wind velocity y-component (km/s)
            - `vz` (list of float, default=0): Solar wind velocity z-component (km/s)
            - `bx` (list of float, default=0): IMF x-component (nT)
            - `by` (list of float, default=5): IMF y-component (nT)
            - `bz` (list of float, default=5): IMF z-component (nT)
            - `by_avg` (list of float, default=0): Averaged IMF By (nT)
            - `bz_avg` (list of float, default=0): Averaged IMF Bz (nT)
            - `density` (list of float, default=1): Solar wind density (particles/cm³)
            - `pdyn` (list of float, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval).

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval).

            Available keys:
            
            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris"): "Boris", "4RK", "Vay", "HC"
            - `gyropercent` (`float`, default=15): Gyration period percentage
            - `minaltitude` (`float`, default=20): Minimum altitude (GDZ = km or other = Re)
            - `maxdistance` (`float`, default=100): Maximum distance (Re)
            - `maxtime` (`float`, default=0): Maximum time
            - `mintrapdist` (`float`, default=0): Minimum trapping distance
            - `startaltitude` (`float`, default=20): Starting altitude (GDZ = km or other = Re)
            - `betaerror` (`float`, default=0.001): Maximum allowed beta error for integration steps %
            - `totalbetacheck` (`bool`, default=False): Enable cumulative beta check during integration
            - `adaptivestep` (`bool`, default=True): Enable adaptive time steps
            - `maxsteps` (`int`, default=None): Maximum number of integration steps

        particle_params (ParticleParams): Particle settings.

            Available keys:

            - `Anum` (`int`, default=1): Atomic number (0=electron, 1=proton, 2=alpha)
            - `anti` (`str`, default="YES"): YES = anti-particle, NO = particle
            - `zenith` (`float`, default=0): Zenith angle for Custom cutoff computation
            - `azimuth` (`float`, default=0): Azimuth angle for Custom cutoff computation


        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `coordsystem` (`str`, default="GEO"): Coordinate system used for calculations
            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

        asymptotic_params (AsymptoticParams): Asymptotic computation parameters.

            Available keys:

            - `unit` (`str`, default="GeV"): Rigidity unit ("GeV", "GV")
            - `asymptotic` (`str`, default="NO"): Enable asymptotic cone computation ("YES"/"NO")
            - `asymlevels` (`list`, default=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000]): Rigidity levels for asymptotic computation

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [flight_dataframe, summary_text, input_dataframe]
            - flight_dataframe: Cutoff rigidities along flight path
            - summary_text: Input parameter summary
            - input_dataframe: Flight path input data

    Examples:
    ```python
        import datetime
        import OTSO
        
        # Basic flight path analysis
        lats = [60.0, 65.0, 70.0]
        lons = [25.0, 20.0, 15.0] 
        alts = [35, 40, 45]  # km
        times = [datetime.datetime(2023, 6, 1, h) for h in [10, 11, 12]]
        
        flight_result = OTSO.flight(lats, lons, times, alts,
                                   computation_params={"corenum": 2})
        
        # With automatic space weather data retrieval
        flight_auto = OTSO.flight(lats, lons, times, alts,
                                 data_retrieval_params={"serverdata": "ON"})
        
        # Access the results
        flight_df, metadata, input_df = flight_result
    ```
    """
    
    
    # For flight function, we need special handling since some parameters might be None
    # for automatic data retrieval. Only merge with defaults if user provided values.
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind} if solar_wind else {}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic} if geomagnetic else {}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko} if tsyganenko else {}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_particle = {**ParticleParams.DEFAULTS, **particle_params}
    final_rigidity = {**RigidityParams.DEFAULTS, **rigidity_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_asymptotic = {**AsymptoticParams.DEFAULTS, **asymptotic_params}
    
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
        **final_asymptotic
    }
    
    # Only include solar wind/geomagnetic parameters if user provided them
    if final_solar_wind:
        flight_args.update(final_solar_wind)
    if final_geomagnetic:
        flight_args.update(final_geomagnetic)
    if final_tsyganenko:
        flight_args.update(final_tsyganenko)
    
    return flight_func(**flight_args)

def magfield(
    Locations: Sequence[float],
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
    *args, **kwargs
):
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

            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
        
        solar_wind (SolarWindParams): Solar wind parameters (optional for auto-retrieval).

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval).

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval).

            Available keys:

            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        coordinate_params (CoordinateParams): Coordinate systems.

            Available keys:

            - `inputcoord` (`str`, default="GDZ"): Input coordinate system

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [magfield_dataframe, summary_text]
            - magfield_dataframe: Magnetic field vectors at input locations
            - summary_text: Input parameter summary

    Examples:
    ```python
        import OTSO
        
        # Basic field evaluation
        field_result = OTSO.magfield([[10, 10, 10]], 
                                    coordinate_params={"coordsystem": "GEO"},
                                    computation_params={"corenum": 1})
        
        # With specific field models
        field_data = OTSO.magfield([[5, 0, 0]], 
                                  datetime_params={"year": 2023, "month": 6},
                                  magfield_params={"externalmag": "TSY01", "internalmag": "IGRF"})
        
        # Access the results
        field_df, metadata = field_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coordinate_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    
    # Remove parameters not supported by magfield function
    magfield_params_filtered = {k: v for k, v in final_magfield.items() 
                               if k not in ['magnetopause', 'spheresize', 'AdaptiveExternalModel']}
    
    return magfield_func(
        Locations,
        **final_solar_wind,
        **final_geomagnetic,
        **final_tsyganenko,
        **final_datetime,
        **magfield_params_filtered,
        **final_coordinate,
        **final_computation,
        **final_data_retrieval,
        **final_custom_field,
        **kwargs
    )

def coordtrans(
    Locations: Sequence[float],
    dates: Sequence,
    CoordIN: str = "GEO",
    CoordOUT: str = "GDZ",
    corenum: Optional[int] = None,
    Verbose: bool = True,
    *args, **kwargs
):
    """
    Transform coordinates between different coordinate systems using OTSO.
    
    Converts spatial coordinates between various reference frames used in
    space physics and geomagnetics. Supports time-dependent transformations
    for date-specific coordinate system orientations.

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
            - coords_dataframe: Transformed coordinates
            - summary_text: Transformation summary

    Examples:
    ```python
        import datetime
        import OTSO
        
        # GEO to GSM transformation
        coord_result = OTSO.coordtrans([[10, 10, 10]], 
                                      [datetime.datetime(2023, 6, 1)],
                                      CoordIN="GEO", CoordOUT="GSM")
        
        # Multiple locations and times
        locs = [[1, 0, 0], [0, 1, 0]]
        times = [datetime.datetime(2023, 6, 1, h) for h in [10, 14]]
        multi_coord = OTSO.coordtrans(locs, times, "GEO", "GSM", corenum=2)
        
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
        Verbose = Verbose, 
        *args, 
        **kwargs
    )

def trace(
    Coordsys: str = "GEO",
    # Solar wind parameters grouped
    solar_wind: SolarWindParams = {},  # vx,vy,vz,bx,by,bz,by_avg,bz_avg,density,pdyn
    # Geomagnetic parameters grouped  
    geomagnetic: GeomagneticParams = {},  # Dst,kp,n_index,b_index,sym_h_corrected
    # Tsyganenko coefficients grouped
    tsyganenko: TsyganenkoParams = {},  # G1,G2,G3,W1,W2,W3,W4,W5,W6
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
    # integration parameters grouped (not used in trace but included for consistency)
    integration_params: IntegrationParams = {},  # intmodel,gyropercent,minaltitude
    # coordinate parameters grouped (not used in trace but included for consistency)
    coord_params: CoordinateParams = {},  # coordsystem,inputcoord
    # grid parameters grouped (not used in trace but included for consistency)
    grid_params: GridParams = {},  # latstep,longstep,maxlat,minlat,maxlong,minlong
    *args, **kwargs
):
    """
    Trace magnetic field lines using the OTSO framework.
    
    Computes magnetic field line trajectories across a global grid to visualize
    the magnetosphere structure. Useful for understanding magnetic connectivity
    and field line topology.

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
            - `internalmag` (`str`, default="IGRF"): "NONE", "IGRF", "Dipole", "Custom Gauss"
            - `externalmag` (`str`, default="TSY89c"): "NONE", "TSY89c", "TSY01", "TSY15B", etc.
            - `boberg` (`bool`, default=False): Enable Boberg extension
            - `bobergtype` (`str`, default="EXTENSION"): Boberg extension type
            - `magnetopause` (`str`, default="Kobel"): "NONE", "Kobel", "Sibeck", "Lin", "Sphere"
            - `spheresize` (`float`, default=25): Spherical boundary radius (Re)
            - `AdaptiveExternalModel` (`bool`, default=False): Auto-select external model
        
        solar_wind (SolarWindParams): Solar wind parameters.
            Available keys:
            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic (GeomagneticParams): Geomagnetic indices.
            Available keys:
            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko (TsyganenkoParams): Tsyganenko model coefficients.
            Available keys:
            - `G1` (`float`, default=0): Tsyganenko G1 coefficient
            - `G2` (`float`, default=0): Tsyganenko G2 coefficient
            - `G3` (`float`, default=0): Tsyganenko G3 coefficient
            - `W1` (`float`, default=0): Tsyganenko W1 coefficient
            - `W2` (`float`, default=0): Tsyganenko W2 coefficient
            - `W3` (`float`, default=0): Tsyganenko W3 coefficient
            - `W4` (`float`, default=0): Tsyganenko W4 coefficient
            - `W5` (`float`, default=0): Tsyganenko W5 coefficient
            - `W6` (`float`, default=0): Tsyganenko W6 coefficient

        grid_params (GridParams): Grid configuration parameters. 
        
            Available keys:

            - `latstep` (`float`, default=-5):  Latitude step size for grid
            - `longstep` (`float`, default=5):  Longitude step size for grid
            - `maxlat` (`float`, default=90):  Maximum latitude for grid
            - `minlat` (`float`, default=-90):  Minimum latitude for grid
            - `maxlong` (`float`, default=360):  Maximum longitude for grid
            - `minlong` (`float`, default=0):  Minimum longitude for grid
            - `array_of_lats_and_longs` (`list`, default=None):  Custom grid points

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
            - `MHDfile` (`str`, default=None): MHD simulation file
            - `MHDcoordsys` (`str`, default=None): MHD coordinate system

    Returns:
        list: [trace_data, summary_text]
            - trace_data: Dictionary with magnetic field line positions
            - summary_text: Input parameter summary

    Examples:
    ```python
        import OTSO
        
        # Global field line tracing
        trace_result = OTSO.trace(
            grid_params={"latstep": 30, "longstep": 60},
            computation_params={"corenum": 4}
        )
        
        # High-resolution polar region tracing
        polar_trace = OTSO.trace(
            grid_params={"latstep": 5, "longstep": 10, "maxlat": 90, "minlat": 60},
            datetime_params={"year": 2023, "month": 3},
            magfield_params={"externalmag": "TSY01"}
        )
        
        # Access the results
        trace_dict, metadata = trace_result
    ```
    """
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko}
    final_datetime = {**DateTimeParams.DEFAULTS, **datetime_params}
    final_magfield = {**MagFieldParams.DEFAULTS, **magfield_params}
    final_computation = {**ComputationParams.DEFAULTS, **computation_params}
    final_data_retrieval = {**DataRetrievalParams.DEFAULTS, **data_retrieval_params}
    final_custom_field = {**CustomFieldParams.DEFAULTS, **custom_field_params}
    final_integration = {**IntegrationParams.DEFAULTS, **integration_params}
    final_coordinate = {**CoordinateParams.DEFAULTS, **coord_params}
    final_grid = {**GridParams.DEFAULTS, **grid_params}
    
    # Remove parameters not supported by trace function
    trace_magfield = {k: v for k, v in final_magfield.items() if k != 'AdaptiveExternalModel'}
    
    return trace_func(
        Coordsys=Coordsys,
        *args,
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
        **final_grid,
        **kwargs
    )

def clean(*args, **kwargs):
    """
    Clean OTSO-generated files and temporary data.
    
    Removes temporary files, cached data, and output files created by OTSO
    calculations to free up disk space and reset the working environment.
    
    Usage:
        OTSO.clean()
        
    Can also be used as CLI command: OTSO.clean
    """
    from .otso_cli import clean as clean_func
    return clean_func(*args, **kwargs)

def addstation(*args, **kwargs):
    """
    Add a new neutron monitor station to the OTSO database.
    
    Registers a new station location for use in OTSO calculations. Useful when
    new neutron monitor stations are established or for adding custom locations
    with specific names for repeated use.
    
    Args:
        StationName (str): Name identifier for the new station.
        Latitude (float): Geographic latitude in degrees.
        Longitude (float): Geographic longitude in degrees.
    
    Usage:
        OTSO.addstation("NEWSTATION", 65.0, 25.0)
        
    CLI usage:
        OTSO.addstation NEWSTATION 65.0 25.0
        
    Note: If station already exists, you will have option to overwrite.
    """
    from .otso_cli import AddStation as addstation_func
    Name, Latitude, Longitude = sys.argv[1], float(sys.argv[2]), float(sys.argv[3])
    addstation_func(Name, Latitude, Longitude)
    return

def removestation(*args, **kwargs):
    """
    Remove a station from the OTSO database.
    
    Deletes a station entry from the OTSO database. Useful for correcting
    incorrectly entered station data or removing obsolete stations.
    
    Args:
        StationName (str): Name of station to remove.
    
    Usage:
        OTSO.removestation("OLDSTATION")
        
    CLI usage:
        OTSO.removestation OLDSTATION
    """
    from .otso_cli import RemoveStation as removestation_func
    Name = sys.argv[1]
    removestation_func(Name)
    return

def liststations(*args, **kwargs):
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
    return liststations_func(*args, **kwargs)

def IGRFupdate(*args, **kwargs):
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
    return IGRFupdate_func(*args, **kwargs)

def serverdownload(*args, **kwargs):
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
    from .otso_cli import ServerDownload as serverdownload_func
    return serverdownload_func(*args, **kwargs)


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