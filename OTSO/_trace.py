def trace(**kwargs):
   import psutil
   from ._core.otso_functions import otso_trace
   from ._core.data_classes.trace_data import TraceData

   # Set corenum if not provided in kwargs
   if 'corenum' not in kwargs or kwargs['corenum'] is None:
       kwargs['corenum'] = psutil.cpu_count(logical=True) - 2
       if kwargs['corenum'] <= 0:
           kwargs['corenum'] = 1
   
   # Handle None values for list parameters
   for key in ['g', 'h', 'MHDfile', 'MHDcoordsys']:
       if key in kwargs and kwargs[key] is None:
           kwargs[key] = []

   TraceDataInstance = TraceData(
       kwargs['startaltitude'], kwargs.get('Coordsys', 'GEO'),
       kwargs['serverdata'], kwargs['livedata'],
       kwargs['vx'], kwargs['vy'], kwargs['vz'],
       kwargs['bx'], kwargs['by'], kwargs['bz'],
       kwargs['density'], kwargs['pdyn'], kwargs['Dst'],
       kwargs['G1'], kwargs['G2'], kwargs['G3'],
       kwargs['W1'], kwargs['W2'], kwargs['W3'], kwargs['W4'],
       kwargs['W5'], kwargs['W6'], kwargs['kp'],
       kwargs['by_avg'], kwargs['bz_avg'], kwargs['n_index'],
       kwargs['b_index'], kwargs['sym_h_corrected'],
       kwargs['year'], kwargs['month'], kwargs['day'],
       kwargs['hour'], kwargs['minute'], kwargs['second'],
       kwargs['internalmag'], kwargs['externalmag'],
       kwargs['boberg'], kwargs['bobergtype'],
       kwargs['gyropercent'], kwargs['magnetopause'], kwargs['corenum'],
       kwargs['latstep'], kwargs['longstep'], kwargs['maxlat'],
       kwargs['minlat'], kwargs['maxlong'], kwargs['minlong'],
       kwargs['g'], kwargs['h'], kwargs['MHDfile'], kwargs['MHDcoordsys'],
       kwargs['spheresize'], kwargs['inputcoord'], kwargs['Verbose'],
       kwargs['mintrapdist'], kwargs["maxsteps"]
   )
   
   trace = otso_trace.OTSO_trace(TraceDataInstance)
   
   return trace