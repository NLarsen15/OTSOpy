def cone(Stations, customlocations=None, **kwargs):
    import psutil
    from ._core.otso_functions import otso_cone
    from ._core.data_classes.cone_data import ConeData

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

    ConeDataInstance = ConeData(
        Stations, customlocations, kwargs['startaltitude'], kwargs['minaltitude'],
        kwargs['maxdistance'], kwargs['maxtime'], kwargs['zenith'], kwargs['azimuth'],
        kwargs['serverdata'], kwargs['livedata'],
        kwargs['vx'], kwargs['vy'], kwargs['vz'],
        kwargs['bx'], kwargs['by'], kwargs['bz'],
        kwargs['density'], kwargs['pdyn'], kwargs['Dst'], 
        kwargs['G1'], kwargs['G2'], kwargs['G3'],
        kwargs['W1'], kwargs['W2'], kwargs['W3'], kwargs['W4'], 
        kwargs['W5'], kwargs['W6'], kwargs['kp'], 
        kwargs['by_avg'], kwargs['bz_avg'], kwargs['n_index'], 
        kwargs['b_index'], kwargs['sym_h_corrected'], kwargs['Anum'], 
        kwargs['anti'], kwargs['year'], kwargs['month'],
        kwargs['day'], kwargs['hour'], kwargs['minute'], 
        kwargs['second'], kwargs['internalmag'], 
        kwargs['externalmag'], kwargs['boberg'], 
        kwargs['bobergtype'], kwargs['intmodel'],
        kwargs['startrigidity'], kwargs['endrigidity'], 
        kwargs['rigiditystep'], kwargs['coordsystem'],
        kwargs['gyropercent'], kwargs['magnetopause'], kwargs['corenum'],
        kwargs['g'], kwargs['h'], kwargs['MHDfile'], kwargs['MHDcoordsys'],
        kwargs['spheresize'], kwargs['inputcoord'], kwargs['Verbose'],
        kwargs['AdaptiveExternalModel'], kwargs['mintrapdist'], kwargs['delim'],
        kwargs['adaptivestep'], kwargs['betaerror'], kwargs['totalbetacheck'],
        kwargs['maxsteps']
    )
    
    cone = otso_cone.OTSO_cone(ConeDataInstance)
    
    return cone