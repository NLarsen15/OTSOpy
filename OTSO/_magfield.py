from __future__ import annotations
from typing import TYPE_CHECKING, Sequence, Optional, Union

if TYPE_CHECKING:
    from ._core.otso_functions.otso_magfield import OTSO_magfield

def magfield(Locations, **kwargs):
    import psutil
    from dataclasses import fields
    from ._core.otso_functions import otso_magfield
    from ._core.data_classes.magfield_data import MagfieldData

    # Set corenum if not provided in kwargs
    if 'corenum' not in kwargs or kwargs['corenum'] is None:
       kwargs['corenum'] = psutil.cpu_count(logical=True) - 2
       if kwargs['corenum'] <= 0:
           kwargs['corenum'] = 1

    allowed = {f.name for f in fields(MagfieldData)}

    unknown = set(kwargs) - allowed
    if unknown:
        raise TypeError(
            f"Unexpected keyword arguments: {', '.join(sorted(unknown))}"
    )
   
   #print(kwargs['g'])

   # Handle None values for list parameters
   #for key in ['g', 'h', 'MHDfile', 'MHDcoordsys']:
   #    if key in kwargs and kwargs[key] is None:
   #        kwargs[key] = []

    MagfieldDataInstance = MagfieldData(
       locations = Locations, 
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
       inputcoord = kwargs['inputcoord'],
       coordout = kwargs['coordout'],
       g = kwargs['g'],
       h = kwargs['h'],
       corenum = kwargs['corenum'],
       MHDfile = kwargs['MHDfile'],
       MHDcoordsys = kwargs['MHDcoordsys'],
       Verbose = kwargs['Verbose'],
       max_degree = kwargs["max_degree"]
   )

    magfield = otso_magfield.OTSO_magfield(MagfieldDataInstance)

   
   
    return magfield