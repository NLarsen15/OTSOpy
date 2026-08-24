from __future__ import annotations
from typing import TYPE_CHECKING, Sequence, Optional, Union

if TYPE_CHECKING:
    from ._core.otso_functions.otso_skymap import OTSO_skymap

def skymap(
    Stations: Union[str, Sequence[str]],
    customlocations: Optional[list] = None,
    **kwargs
) -> OTSO_skymap:
    
    import psutil
    from dataclasses import fields
    from ._core.otso_functions import otso_skymap
    from ._core.data_classes.skymap_data import SkymapData

    allowed = {f.name for f in fields(SkymapData)}

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
    
    # Handle None values for list parameters
    for key in ['g', 'h', 'MHDfile', 'MHDcoordsys']:
        if key in kwargs and kwargs[key] is None:
            kwargs[key] = []
    if customlocations is None:
        customlocations = []

    # Create CutoffData instance - all parameters already provided by wrapper
    SkymapDataInstance = SkymapData(
        Stations = Stations, 
        customlocations = customlocations, 
        startaltitude = kwargs['startaltitude'],
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
        Anum = kwargs['Anum'], 
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
        mintrapdist = kwargs['mintrapdist'], 
        adaptivestep = kwargs['adaptivestep'], 
        betaerror = kwargs['betaerror'], 
        totalbetacheck = kwargs['totalbetacheck'],
        maxsteps = kwargs['maxsteps'],
        threadnum = kwargs['threadnum'],
        zenithstep = kwargs["zenithstep"],
        azimuthstep = kwargs["azimuthstep"],
        maxzenith = kwargs["maxzenith"],
        minzenith = kwargs["minzenith"],
        maxazimuth = kwargs["maxazimuth"],
        minazimuth = kwargs["minazimuth"],
        max_degree = kwargs["max_degree"]
    )

    skymap = otso_skymap.OTSO_skymap(SkymapDataInstance)
    
    return skymap