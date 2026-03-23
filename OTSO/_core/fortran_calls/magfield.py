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
    old_altitude = None
    if MagfieldDataInstance.inputcoord == "GDZ":
          old_altitude = x[2]  # Save original altitude
          Position = [1+(x[2]/6371.0),x[0],x[1]]
    else:
          Position = x  # Use original position for non-GDZ coordinates

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
    
    # For GDZ coordinates, replace the transformed altitude with original altitude
    if CoordinateSystem == "GDZ" and old_altitude is not None:
        display_position = [old_altitude, Position[1], Position[2]]  # Use original altitude
    else:
        display_position = Position
    
    combined_array = np.concatenate((display_position, Bfield))
    coord_suffix = f"_{CoordinateSystem}"  # Add the coordinate system suffix
    
    # Define position columns based on coordinate system
    if CoordinateSystem == "GDZ":
        position_columns = [f"altitude{coord_suffix} [km]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
    elif CoordinateSystem == "SPH":
        position_columns = [f"altitude{coord_suffix} [Re]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
    else:
        # Default to X, Y, Z for other coordinate systems
        position_columns = [f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]"]
    
    columns = position_columns + [f"{CoordOUT}_Bx [nT]", f"{CoordOUT}_By [nT]", f"{CoordOUT}_Bz [nT]"]
    df = pd.DataFrame(columns=columns)
    df.loc[len(df)] = combined_array
    queue.put(df)

  return