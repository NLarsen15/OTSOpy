def coordtrans(Locations,dates,CoordIN,CoordOUT,corenum=None):
    from .Parameters.functions import otso_coordtrans
    import psutil

    if corenum is None:
       corenum = psutil.cpu_count(logical=False) - 2
       if corenum <= 0:
          corenum = 1
    
    arguments = locals()
    for arg in arguments:
       if arguments[arg] is None:
          arguments[arg] = []

    coordtrans = otso_coordtrans.OTSO_coordtrans(Locations,dates,CoordIN,CoordOUT,corenum)
    
    return coordtrans