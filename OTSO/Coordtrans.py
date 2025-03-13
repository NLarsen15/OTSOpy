def coordtrans(Locations,Dates,CoordIN,CoordOUT,corenum=1):
    from .Parameters.functions import otso_coordtrans
    
    arguments = locals()
    for arg in arguments:
       if arguments[arg] is None:
          arguments[arg] = []

    coordtrans = otso_coordtrans.OTSO_coordtrans(Locations,Dates,CoordIN,CoordOUT,corenum)
    
    return coordtrans