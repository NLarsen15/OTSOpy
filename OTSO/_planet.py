from __future__ import annotations
from typing import TYPE_CHECKING, Sequence, Optional, Union

if TYPE_CHECKING:
    from ._core.otso_functions.otso_planet import OTSO_planet

def planet(
    **kwargs
) -> "OTSO_planet":
    
    import psutil
    from dataclasses import fields
    from ._core.otso_functions import otso_planet
    from ._core.data_classes.planet_data import PlanetData

    allowed = {f.name for f in fields(PlanetData)}

    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(
            f"Unexpected keyword arguments: {', '.join(sorted(unknown))}"
    )

    # Set corenum if not provided in kwargs
    if 'corenum' not in kwargs or kwargs['corenum'] is None:
        kwargs['corenum'] = psutil.cpu_count(logical=True) - 2
        if kwargs['corenum'] <= 0:
            kwargs['corenum'] = 1

    # Check if grid parameters were explicitly set by the user  
    grid_param_keys = {'latstep', 'longstep', 'maxlat', 'minlat', 'maxlong', 'minlong'}
    defaults = {
        'latstep': -5, 'longstep': 5, 'maxlat': 90, 'minlat': -90, 'maxlong': 360, 'minlong': 0
    }
    grid_params_user_set = any(key in kwargs and kwargs[key] != defaults[key] for key in grid_param_keys)
    
    # Handle None values for list parameters
    for key in ['g', 'h', 'MHDfile', 'MHDcoordsys']:
        if key in kwargs and kwargs[key] is None:
            kwargs[key] = []

    # Handle array_of_lats_and_longs specially
    if 'array_of_lats_and_longs' not in kwargs:
        kwargs['array_of_lats_and_longs'] = None 

    # Create CutoffData instance - all parameters already provided by wrapper
    PlanetDataInstance = PlanetData(
        startaltitude = kwargs['startaltitude'], 
        cutoff_comp = kwargs['cutoff_comp'], 
        minaltitude = kwargs['minaltitude'], 
        maxdistance = kwargs['maxdistance'], 
        maxtime = kwargs['maxtime'],
        serverdata = kwargs['serverdata'], 
        livedata = kwargs['livedata'],
        vx = kwargs['vx'], 
        vy = kwargs['vy'], 
        vz = kwargs['vz'], 
        bx = kwargs['bx'], 
        by = kwargs['by'], 
        bz = kwargs['bz'],
        density = kwargs['density'], 
        pdyn = kwargs['pdyn'], 
        Dst = kwargs['Dst'], 
        G1 = kwargs['G1'], 
        G2 = kwargs['G2'], 
        G3 = kwargs['G3'],
        W1 = kwargs['W1'], 
        W2 = kwargs['W2'], 
        W3 = kwargs['W3'], 
        W4 = kwargs['W4'], 
        W5 = kwargs['W5'], 
        W6 = kwargs['W6'], 
        kp = kwargs['kp'], 
        by_avg = kwargs['by_avg'], 
        bz_avg = kwargs['bz_avg'], 
        n_index = kwargs['n_index'], 
        b_index = kwargs['b_index'], 
        sym_h_corrected = kwargs['sym_h_corrected'], 
        anti = kwargs['anti'], 
        year = kwargs['year'], 
        month = kwargs['month'],
        day = kwargs['day'], 
        hour = kwargs['hour'], 
        minute = kwargs['minute'], 
        second = kwargs['second'], 
        internalmag = kwargs['internalmag'], 
        externalmag = kwargs['externalmag'], 
        boberg = kwargs['boberg'], 
        bobergtype = kwargs['bobergtype'], 
        intmodel = kwargs['intmodel'],
        startrigidity = kwargs['startrigidity'], 
        endrigidity = kwargs['endrigidity'], 
        rigiditystep = kwargs['rigiditystep'], 
        rigidityscan = kwargs['rigidityscan'],
        coordsystem = kwargs['coordsystem'],
        gyropercent = kwargs['gyropercent'],
        fixedstep = kwargs['fixedstep'],
        magnetopause = kwargs['magnetopause'],
        corenum = kwargs['corenum'],
        azimuth = kwargs['azimuth'], 
        zenith = kwargs['zenith'], 
        g = kwargs['g'], 
        h = kwargs['h'], 
        MHDfile = kwargs['MHDfile'], 
        MHDcoordsys = kwargs['MHDcoordsys'],
        spheresize = kwargs['spheresize'], 
        inputcoord = kwargs['inputcoord'], 
        Verbose = kwargs['Verbose'], 
        AdaptiveExternalModel = kwargs['AdaptiveExternalModel'], 
        array_of_lats_and_longs = kwargs['array_of_lats_and_longs'],
        grid_params_user_set = grid_params_user_set,  # Pass the flag indicating user-set grid params
        mintrapdist = kwargs['mintrapdist'], 
        unit = kwargs['unit'], 
        asymptotic = kwargs['asymptotic'],
        asymlevels = kwargs['asymlevels'], 
        delim = kwargs['delim'],
        adaptivestep = kwargs['adaptivestep'], 
        betaerror = kwargs['betaerror'], 
        totalbetacheck = kwargs['totalbetacheck'],
        maxsteps = kwargs['maxsteps'],
        threadnum = kwargs['threadnum'],
        latstep = kwargs['latstep'],
        longstep = kwargs['longstep'],
        maxlat = kwargs['maxlat'],
        minlat = kwargs['minlat'],
        maxlong = kwargs['maxlong'],
        minlong = kwargs['minlong'],
        transmission = kwargs['transmission'],
        transmissionsamples = kwargs['transmissionsamples'],
        transmissionRstep = kwargs['transmissionRstep'],
        max_degree = kwargs["max_degree"]
    )

    planet = otso_planet.OTSO_planet(PlanetDataInstance)
    
    return planet