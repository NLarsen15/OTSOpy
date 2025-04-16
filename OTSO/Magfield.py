def magfield(Locations,serverdata="OFF",livedata="OFF",vx=-500,vy=0,vz=0,by=5,bz=5,density=1,pdyn=0,Dst=0,
           G1=0,G2=0,G3=0,W1=0,W2=0,W3=0,W4=0,W5=0,W6=0,kp=0,year=2024,
           month=1,day=1,hour=12,minute=0,second=0,internalmag="IGRF",externalmag="TSY89",
           coordsystem="GEO", g=None,h=None,corenum=1,MHDfile=None,MHDcoordsys=None):
    from .Parameters.functions import otso_magfield
    
    arguments = locals()
    for arg in arguments:
       if arguments[arg] is None:
          arguments[arg] = []

    magfield = otso_magfield.OTSO_magfield(Locations,serverdata,livedata,vx,vy,vz,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,year,
           month,day,hour,minute,second,internalmag,externalmag,
           coordsystem,arguments["g"],arguments["h"],corenum,MHDfile,MHDcoordsys)
    
    return magfield