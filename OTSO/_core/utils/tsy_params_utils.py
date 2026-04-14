import math
import pandas
import numpy as np

def IOPTprocess(Kp:float) -> int:
    if Kp >= 6:
        IOPT = 7
    else:
        IOPT = Kp + 1

    return IOPT

def IOPTprocess_refit(Kp:float) -> int:
    # Sub model uses fractional Kp values indexed as integers
    # Kp=0 -> IOPT=1 (0.0), Kp=1 -> IOPT=2 (0.333), Kp=2 -> IOPT=3 (0.666),
    # Kp=3 -> IOPT=4 (1.0), Kp=4 -> IOPT=5 (1.333), etc.
    # Goes up to Kp=7 -> IOPT=22 (7.0) and beyond
    if Kp >= 7:
        IOPT = 22  # Corresponds to Kp=7.0
    else:
        # Map Kp (0-9) to integer IOPT values (1-22)
        # Each Kp unit has 3 sub-divisions (0, 0.333, 0.666)
        IOPT = int(Kp * 3) + 1

    return IOPT


def TSY01_Constants(By:float,Bz:float,V:float,N:float) -> tuple[float,float,float]:
    G1 = 0
    G2 = 0
    G3 = 0

    if V == 9999.0:
        V = np.nan

    By = By
    Bz = Bz
    V = V
    W = 1

    B = (By*By + Bz*Bz)**(0.5)
    h = (((B/40)**(2))/(1 + B/40))

    if(By == 0 and Bz == 0):
        phi = 0
    else:
        phi = math.atan2(By,Bz)
        if(phi <= 0):
            phi = phi + 2*math.pi

    if(Bz < 0):
        Bs = abs(Bz)
    elif Bz >= 0:
        Bs = 0
    
    G1 = G1 + (W)*V*h*(math.sin(phi/2)**3)
    G2 = G2 + (0.005)*(W)*(V)*Bs
    G3 = G3 + (N*V*Bs)/2000

    return round(G1, 2), round(G2, 2), round(G3, 2)

def OMNI_TSY01_Constants(data: pandas.DataFrame) -> tuple[float,float,float]:
    IMFy = data["By"]
    IMFz = data["Bz"]
    Speed = data["V"]
    Density = data["Density"]

    G1 = 0
    G2 = 0
    G3 = 0

    for (By, Bz, V, N) in zip(IMFy, IMFz, Speed, Density):
        By = By
        Bz = Bz
        if V == 9999.0:
            V = np.nan
        if V < 0:
            V = -1*V
        W = 1/len(IMFy)
    
    
        B = (By*By + Bz*Bz)**(0.5)
        h = (((B/40)**(2))/(1 + B/40))
    
        if(By == 0 and Bz == 0):
            phi = 0
        else:
            phi = math.atan2(By,Bz)
            if(phi <= 0):
                phi = phi + 2*math.pi
    
        if(Bz < 0):
            Bs = abs(Bz)
        elif Bz >= 0:
            Bs = 0
        
        G1 = G1 + (W)*V*h*(math.sin(phi/2)**3)
        G2 = G2 + (0.005)*(W)*(V)*Bs
        G3 = G3 + (N*V*Bs)/2000

    return round(G1, 2), round(G2, 2), round(G3, 2)

def N_index(V:float, Bt:float, thetac:float) -> float:

    N = (10**(-4))*(V**(4/3))*(Bt**(2/3))*(math.sin(thetac/2)**(8/3))

    return N

def N_index_normalized(V:float, Bt:float, thetac:float) -> float:

    N = 0.86*((V/400)**(4/3))*((Bt/5)**(2/3))*(math.sin(thetac/2)**(8/3))

    return N

def B_index(Np:float, V:float, Bt:float, thetac:float) -> float:

    B = ((Np/5)**(1/2))*((V/400)**(5/2))*(Bt/5)*(math.sin(thetac/2)**6)

    return B

def Pdyn_comp(Density:float, Vx:float) -> float:
    Pressure = ((((Density))*10**6)*(1.672621898e-27)*((Vx*1000)**2))*10**(9)
    Pressure = round(Pressure,3)
    return Pressure

