import numpy as np

class SolarWind:
 def __init__(self, Vx, Vy, Vz, Bx, By, Bz, Density, Pdyn, Dst, G1, G2, G3, W1, W2, W3, W4, W5, W6, Kp, by_avg, bz_avg, n_index, b_index, sym_h_corrected):
    self.Vx = Vx
    self.Vy = Vy
    self.Vz = Vz
    self.Bx = Bx
    self.By = By
    self.Bz = Bz
    self.Density = Density
    # Pdyn = 0 means "compute from density and velocity". Done here, with OMNI's convention
    # P = 2e-6 Np V^2 (nPa; see tsy_params_utils.Pdyn_comp), so the Fortran never falls back
    # to its own proton-only formula (~17% lower than OMNI's pressure).
    if Pdyn == 0 and Density > 0:
      Pdyn = 2e-6 * Density * (Vx**2 + Vy**2 + Vz**2)
    self.Pdyn = Pdyn
    self.Dst = Dst
    self.G1 = G1
    self.G2 = G2
    self.G3 = G3
    self.W1 = W1
    self.W2 = W2
    self.W3 = W3
    self.W4 = W4
    self.W5 = W5
    self.W6 = W6
    self.KpIndex = Kp
    self.by_avg = by_avg
    self.bz_avg = bz_avg
    self.n_index = n_index
    self.b_index = b_index
    self.sym_h_corrected = sym_h_corrected


    self.DstBob = 0

    if Kp > 6:
      self.DstBob = (Kp-6)*(-50)

    if Dst != 0:
      self.DstBob = Dst

    self.WindArray = np.array([self.Vx, self.Vy, self.Vz, self.Bx, self.By, 
                               self.Bz, self.Density, self.Pdyn, self.Dst, 
                               self.G1, self.G2, self.G3, self.W1, self.W2, 
                               self.W3, self.W4, self.W5, 
                               self.W6, self.KpIndex, self.DstBob, 
                               self.by_avg, self.bz_avg, self.n_index, 
                               self.b_index, self.sym_h_corrected])

 def GetWind(self):
    return self.WindArray