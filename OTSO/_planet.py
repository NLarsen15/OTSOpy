def planet(**kwargs):
    import psutil    
    from ._core.otso_functions import otso_planet
    from ._core.data_classes.planet_data import PlanetData

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

    PlanetDataInstance = PlanetData(
        kwargs['startaltitude'],
        kwargs['cutoff_comp'],
        kwargs['minaltitude'],
        kwargs['maxdistance'],
        kwargs['maxtime'],
        kwargs['serverdata'],
        kwargs['livedata'], 
        kwargs['vx'], 
        kwargs['vy'], 
        kwargs['vz'], 
        kwargs['bx'], 
        kwargs['by'], 
        kwargs['bz'], 
        kwargs['density'], 
        kwargs['pdyn'], 
        kwargs['Dst'],
        kwargs['G1'], 
        kwargs['G2'], 
        kwargs['G3'], 
        kwargs['W1'], 
        kwargs['W2'], 
        kwargs['W3'], 
        kwargs['W4'], 
        kwargs['W5'], 
        kwargs['W6'], 
        kwargs['kp'], 
        kwargs['by_avg'], 
        kwargs['bz_avg'], 
        kwargs['n_index'], 
        kwargs['b_index'], 
        kwargs['sym_h_corrected'], 
        kwargs['anti'], 
        kwargs['year'],
        kwargs['month'], 
        kwargs['day'], 
        kwargs['hour'], 
        kwargs['minute'], 
        kwargs['second'], 
        kwargs['internalmag'], 
        kwargs['externalmag'],
        kwargs['boberg'],
        kwargs['bobergtype'],
        kwargs['intmodel'], 
        kwargs['startrigidity'], 
        kwargs['endrigidity'], 
        kwargs['rigiditystep'], 
        kwargs['rigidityscan'],
        kwargs['gyropercent'], 
        kwargs['magnetopause'], 
        kwargs['corenum'], 
        kwargs['azimuth'], 
        kwargs['zenith'], 
        kwargs['asymptotic'], 
        kwargs['asymlevels'], 
        kwargs['unit'],
        kwargs['latstep'], 
        kwargs['longstep'], 
        kwargs['maxlat'], 
        kwargs['minlat'], 
        kwargs['maxlong'], 
        kwargs['minlong'], 
        kwargs['g'], 
        kwargs['h'], 
        kwargs['MHDfile'], 
        kwargs['MHDcoordsys'],
        kwargs['spheresize'],
        kwargs['inputcoord'], 
        kwargs['Verbose'],
        kwargs['AdaptiveExternalModel'],
        kwargs['array_of_lats_and_longs'],
        grid_params_user_set,  # Pass the flag indicating user-set grid params
        kwargs['mintrapdist'],
        kwargs['delim'],
        kwargs['adaptivestep'], 
        kwargs['betaerror'],
        kwargs['totalbetacheck'],
        kwargs['maxsteps']
    )

    # Pass the flag indicating user-set grid params
    planet_result = otso_planet.OTSO_planet(PlanetDataInstance)
    
    return planet_result