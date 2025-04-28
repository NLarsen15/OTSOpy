def planet(startaltitude=20,cutoff_comp="Vertical",minaltitude=20,maxdistance=100,maxtime=0,
           serverdata="OFF",livedata="OFF",vx=-500,vy=0,vz=0,by=5,bz=5,density=1,pdyn=0,Dst=0,
           G1=0,G2=0,G3=0,W1=0,W2=0,W3=0,W4=0,W5=0,W6=0,kp=0,anti="YES",year=2024,
           month=1,day=1,hour=12,minute=0,second=0,internalmag="IGRF",externalmag="TSY89",
           intmodel="Boris",startrigidity=20,endrigidity=0,rigiditystep=0.01,rigidityscan="ON",
           gyropercent=15,magnetopause="Kobel",corenum=1, azimuth=0, zenith=0, asymptotic="NO",
           asymlevels = [0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,10,15,20,30,50,70,100,300,500,700,1000],unit="GeV",
           latstep=-5,longstep=5, maxlat=90,minlat=-90,maxlong=360,minlong=0,g=None,h=None,MHDfile=None,MHDcoordsys=None,
           array_of_lats_and_longs=None):
    from .Parameters.functions import otso_planet
    
    # Capture all local variables, including defaults and user-provided kwargs
    arguments = locals()
    
    # --- Check if grid parameters were explicitly set by the user --- 
    grid_param_keys = {'latstep', 'longstep', 'maxlat', 'minlat', 'maxlong', 'minlong'}
    # We need a way to know what was *actually* passed in the function call.
    # inspect.signature might be too complex here. A common pattern is to use a 
    # sentinel default value (like object()) to distinguish, but that requires changing defaults.
    # Alternative: Check if the values in `arguments` differ from the defaults defined here.
    # This assumes defaults are static, which they are in this case.
    defaults = {
        'latstep': -5, 'longstep': 5, 'maxlat': 90, 'minlat': -90, 'maxlong': 360, 'minlong': 0
    }
    grid_params_user_set = False
    for key in grid_param_keys:
        if arguments[key] != defaults[key]:
            grid_params_user_set = True
            break
    # --- End Check --- 

    # Exclude specific keys where None is a valid/handled input
    excluded_keys = {'g', 'h', 'array_of_lats_and_longs'}
    for arg, value in arguments.items():
       if arg not in excluded_keys and value is None:
          # Defaulting None to [] might still be needed for some args?
          arguments[arg] = [] 

    # Pass the flag indicating user-set grid params
    planet_result = otso_planet.OTSO_planet(
        startaltitude, cutoff_comp, minaltitude, maxdistance, maxtime,
        serverdata, livedata, vx, vy, vz, by, bz, density, pdyn, Dst,
        G1, G2, G3, W1, W2, W3, W4, W5, W6, kp, anti, year,
        month, day, hour, minute, second, internalmag, externalmag,
        intmodel, startrigidity, endrigidity, rigiditystep, rigidityscan,
        gyropercent, magnetopause, corenum, azimuth, zenith, asymptotic, asymlevels, unit,
        latstep, longstep, maxlat, minlat, maxlong, minlong, g, h, MHDfile, MHDcoordsys,
        array_of_lats_and_longs=array_of_lats_and_longs,
        grid_params_user_set=grid_params_user_set # Pass the flag
    )
    
    return planet_result