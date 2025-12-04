def cone(Stations, customlocations=[], startaltitude=20,minaltitude=20,zenith=0,azimuth=0,maxdistance=100,maxtime=0,
           serverdata="OFF",livedata="OFF",vx=-500,vy=0,vz=0,bx=0,by=5,bz=5,density=1,pdyn=0,Dst=0,
           G1=0,G2=0,G3=0,W1=0,W2=0,W3=0,W4=0,W5=0,W6=0,kp=0,by_avg=0,bz_avg=0,n_index=0,b_index=0,sym_h_corrected=0,Anum=1,anti="YES",year=2024,
           month=1,day=1,hour=12,minute=0,second=0,internalmag="IGRF",externalmag="TSY89",
           intmodel="Boris",startrigidity=20,endrigidity=0,rigiditystep=0.01,
           coordsystem="GEO",gyropercent=15,magnetopause="Kobel",corenum=None,g=None,h=None,MHDfile=None,MHDcoordsys=None,
           spheresize=25, inputcoord="GDZ", Verbose=True, AdaptiveExternalModel=False):
    from .Parameters.functions import otso_cone
    import psutil

    if corenum is None:
       corenum = psutil.cpu_count(logical=False) - 2
       if corenum <= 0:
          corenum = 1
    
    arguments = locals()
    for arg in arguments:
       if arguments[arg] is None:
          arguments[arg] = []
    
    cone = otso_cone.OTSO_cone(Stations,customlocations,startaltitude,minaltitude,zenith,azimuth,maxdistance,maxtime,
           serverdata,livedata,vx,vy,vz,bx,by,bz,density,pdyn,Dst,
           G1,G2,G3,W1,W2,W3,W4,W5,W6,kp,by_avg,bz_avg,n_index,b_index,sym_h_corrected,Anum,anti,year,
           month,day,hour,minute,second,internalmag,externalmag,
           intmodel,startrigidity,endrigidity,rigiditystep,
           coordsystem,gyropercent,magnetopause,corenum,g,h,MHDfile,MHDcoordsys,spheresize, inputcoord, Verbose, AdaptiveExternalModel)
    
    return cone