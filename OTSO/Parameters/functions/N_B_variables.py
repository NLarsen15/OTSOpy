import math

# V is solar wind speed in km/s
# Bt is the tranverse magnetic field component in nT
# thetac is the IMF clock angle in radians
# Np is the solar wind proton density in cm^-3

def N_index(V, Bt, thetac):

    N = (10**(-4))*(V**(4/3))*(Bt**(2/3))*(math.sin(thetac/2)**(8/3))

    return N

def N_index_normalized(V, Bt, thetac):

    N = 0.86*((V/400)**(4/3))*((Bt/5)**(2/3))*(math.sin(thetac/2)**(8/3))

    return N

def B_index(Np, V, Bt, thetac):

    B = ((Np/5)**(1/2))*((V/400)**(5/2))*(Bt/5)*(math.sin(thetac/2)**6)

    return B
