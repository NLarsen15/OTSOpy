import pandas as pd
import numpy as np
import multiprocessing as mp

from ..libs import MiddleMan as OTSOLib
from ..utils import mhd_utils
from ..data_classes.magfield_data import MagfieldData

def FortranMagfield(Data: list, MagfieldDataInstance: MagfieldData, queue: mp.Queue) -> None:
    
  if MagfieldDataInstance.model[1] == 99:
    mhd_utils.MHDinitialise(MagfieldDataInstance.MHDfile)

  for x in Data:
    Position = x

    DateArray = MagfieldDataInstance.datearray
    model = MagfieldDataInstance.model
    IOPT = MagfieldDataInstance.IOPT
    WindArray = MagfieldDataInstance.windarray
    CoordinateSystem = MagfieldDataInstance.inputcoord
    CoordOUT = MagfieldDataInstance.coordout
    MHDCoordSys = MagfieldDataInstance.MHDcoordsys
    g = MagfieldDataInstance.g
    h = MagfieldDataInstance.h
    
    Bfield = OTSOLib.magstrength(Position, DateArray, model, IOPT, WindArray, CoordinateSystem, CoordOUT, MHDCoordSys, g, h)
    Bfield = Bfield*10**9
    combined_array = np.concatenate((Position, Bfield))
    coord_suffix = f"_{CoordinateSystem}"  # Add the coordinate system suffix
    columns = [f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]", f"{CoordOUT}_Bx [nT]", f"{CoordOUT}_By [nT]", f"{CoordOUT}_Bz [nT]"]
    df = pd.DataFrame(columns=columns)
    df.loc[len(df)] = combined_array
    queue.put(df)

  return