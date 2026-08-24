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
    AsymptoticParams,
    TransmissionParams,
    SkymapParams
)

###############################################################################################################################

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

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")

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
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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

    if integration_params.get("gyropercent") is None:
        integration_params["gyropercent"] = 1

    if integration_params.get("adaptivestep") is None:
        integration_params["adaptivestep"] = False
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")

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

        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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
    *args, **kwargs
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
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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

    Returns:
        list: [trajectory_data, readme_text]

            - trajectory_data: Dictionary with positional, velocity, filter and asymptotic viewing directions information for trajectories
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

    if integration_params.get("gyropercent") is None:
        integration_params["gyropercent"] = 1

    if integration_params.get("adaptivestep") is None:
        integration_params["adaptivestep"] = False
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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
    *args, **kwargs
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

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")

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
        
        solar_wind_params (SolarWindParams): Solar wind parameters (optional for auto-retrieval).

            All values should be provided as lists of floats, one per flight point.

            Available keys:

            - `vx` (list of float, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (list of float, default=0): Solar wind velocity y-component (km/s)
            - `vz` (list of float, default=0): Solar wind velocity z-component (km/s)
            - `bx` (list of float, default=0): IMF x-component (nT)
            - `by` (list of float, default=5): IMF y-component (nT)
            - `bz` (list of float, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (list of float, default=1): Solar wind density (particles/cm³)
            - `pdyn` (list of float, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval).

            Available keys:

            - `Dst` (list of float,, default=0): Dst index (nT)
            - `kp` (list of float,, default=0): Kp index (0-9)
            - `n_index` (list of float,, default=0): Newell coupling function
            - `b_index` (list of float,, default=0): Boynton coupling function
            - `sym_h_corrected` (list of float,, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval).

            Available keys:
            
            - `G1` (list of float,, default=0): Tsyganenko G1 coefficient
            - `G2` (list of float,, default=0): Tsyganenko G2 coefficient
            - `G3` (list of float,, default=0): Tsyganenko G3 coefficient
            - `W1` (list of float,, default=0): Tsyganenko W1 coefficient
            - `W2` (list of float,, default=0): Tsyganenko W2 coefficient
            - `W3` (list of float,, default=0): Tsyganenko W3 coefficient
            - `W4` (list of float,, default=0): Tsyganenko W4 coefficient
            - `W5` (list of float,, default=0): Tsyganenko W5 coefficient
            - `W6` (list of float,, default=0): Tsyganenko W6 coefficient

        integration_params (IntegrationParams): Integration settings.

            Available keys:

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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
    *args, **kwargs
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
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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


    return trace_func(
        coordsys=coordsys,
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
#################################################################################################################
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
    *args, **kwargs
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
        
        solar_wind_params (SolarWindParams): Solar wind parameters (optional for auto-retrieval).

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices (optional for auto-retrieval).

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients (optional for auto-retrieval).

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

    Returns:
        list: [magfield_dataframe, readme_text]

            - magfield_dataframe: Magnetic field vectors at input locations
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
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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
############################################################################################################
#######################################################################################################################################
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

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)

        transmission_params (TransmissionParams): Transmission function computation parameters

            Available keys:

            - `transmissionRstep` (`float`, default=0.001): R ± transmissionRstep defines the transmission rigidity sampling range 
            - `transmissionsamples` (`int`, default=20): the number of sample rigidities to test within the sampling range
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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

        rigidity_params (RigidityParams): Rigidity scanning.

            Available keys:

            - `startrigidity` (`float`, default=20): Initial rigidity (GV)
            - `endrigidity` (`float`, default=0): Final rigidity (GV)
            - `rigiditystep` (`float`, default=0.01): Step size (GV)
            - `rigidityscan` (`str`, default="ON"): Enable scanning ("ON"/"OFF")

        skymap_params (SkymapParams): Sky map zenith and azimuth resolution parameters

            Available keys:

            - `zenithstep` (`float`, default=15): Zenith resolution
            - `azimuthstep` (`float`, default=45): Azimuth resolution
            - `maxzenith` (`float`, default=75):  Maximum deviation from the zenith
            - `minzenith` (`float`, default=0):  Minimum deviation from the zenith
            - `maxazimuth` (`float`, default=360): Maximum azimuth value
            - `minazimuth` (`float`, default=0): Minimum azimuth value
        
        solar_wind_params (SolarWindParams): Solar wind parameters.

            Available keys:

            - `vx` (`float`, default=-500): Solar wind velocity x-component (km/s)
            - `vy` (`float`, default=0): Solar wind velocity y-component (km/s)
            - `vz` (`float`, default=0): Solar wind velocity z-component (km/s)
            - `bx` (`float`, default=0): IMF x-component (nT)
            - `by` (`float`, default=5): IMF y-component (nT)
            - `bz` (`float`, default=5): IMF z-component (nT)
            - `by_avg` (`float`, default=0): Averaged IMF By over last 30mins (nT)
            - `bz_avg` (`float`, default=0): Averaged IMF Bz over last 30mins (nT)
            - `density` (`float`, default=1): Solar wind density (particles/cm³)
            - `pdyn` (`float`, default=0): Solar wind dynamic pressure (nPa)
            
        geomagnetic_params (GeomagneticParams): Geomagnetic indices.

            Available keys:

            - `Dst` (`float`, default=0): Dst index (nT)
            - `kp` (`float`, default=0): Kp index (0-9)
            - `n_index` (`float`, default=0): Newell coupling function
            - `b_index` (`float`, default=0): Boynton coupling function
            - `sym_h_corrected` (`float`, default=0): Corrected SYM-H index (nT)
            
        tsyganenko_params (TsyganenkoParams): Tsyganenko model coefficients.

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

            - `intmodel` (`str`, default="Boris-Buneman"): "4RK", "5RK", "6RK", "Vay", "HC", "Boris-Buneman"
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

    if integration_params.get("gyropercent") is None:
        integration_params["gyropercent"] = 1
    
    # Merge user parameters with defaults from parameter classes
    final_solar_wind = {**SolarWindParams.DEFAULTS, **solar_wind_params}
    final_geomagnetic = {**GeomagneticParams.DEFAULTS, **geomagnetic_params}
    final_tsyganenko = {**TsyganenkoParams.DEFAULTS, **tsyganenko_params}
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

def coordtrans(
    Locations: Sequence[float],
    dates: Sequence,
    CoordIN: str = "GEO",
    CoordOUT: str = "GDZ",
    corenum: Optional[int] = None,
    Verbose: bool = True,
    *args, **kwargs
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

            - coords_dataframe: Transformed coordinates
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
        Verbose = Verbose, 
        *args, 
        **kwargs
    )
#######################################################################################################################
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

def chaosdownload(*args, **kwargs):
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
    return CHAOSdownload_func(*args, **kwargs)


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