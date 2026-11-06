import pandas as pd
import numpy as np
import multiprocessing as mp

from ..libs.MiddleMan import Middleman as OTSOLib
from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.magfield_data import MagfieldData


def magfield_columns(CoordinateSystem: str, CoordOUT: str) -> list:
    coord_suffix = f"_{CoordinateSystem}"
    if CoordinateSystem == "GDZ":
        position_columns = [f"altitude{coord_suffix} [km]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
    elif CoordinateSystem == "SPH":
        position_columns = [f"altitude{coord_suffix} [Re]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
    else:
        position_columns = [f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]"]
    return position_columns + [f"{CoordOUT}_Bx [nT]", f"{CoordOUT}_By [nT]", f"{CoordOUT}_Bz [nT]"]


def FortranMagfield(Data: list, MagfieldDataInstance: MagfieldData, queue: mp.Queue, Index: list) -> None:
    """Field at every location in Data; puts ONE DataFrame (all of this worker's points) on the queue.

    The Fortran inputs are prepared once, and results go into a NumPy array, instead of a
    one-row DataFrame and a queue message per point (which cost ~100x the Fortran call)."""
    if MagfieldDataInstance.model[1] == 99:
        mhd_utils.MHDinitialise(MagfieldDataInstance.MHDfile, MagfieldDataInstance.MHDgridtype,
                                MagfieldDataInstance.MHDinterpolation)

    CoordinateSystem = MagfieldDataInstance.inputcoord
    CoordOUT = MagfieldDataInstance.coordout
    FortranData = fortran_data_utils.prepare_fortran_magfield(MagfieldDataInstance)
    g, h = MagfieldDataInstance.g, MagfieldDataInstance.h

    Data = np.asarray(Data, dtype=np.float64).reshape(-1, 3)
    out = np.empty((len(Data), 6), dtype=np.float64)
    Bfield = np.zeros(3, dtype=np.float64)
    for k, x in enumerate(Data):
        if CoordinateSystem == "GDZ":
            Position = [1 + (x[2] / 6371.0), x[0], x[1]]
            display_position = [x[2], Position[1], Position[2]]  # original altitude, as before
        else:
            Position = x
            display_position = x
        OTSOLib.magstrength(Position, FortranData, CoordinateSystem, CoordOUT, g, h, Bfield)
        out[k, :3] = display_position
        out[k, 3:] = Bfield * 10**9

    df = pd.DataFrame(out, columns=magfield_columns(CoordinateSystem, CoordOUT), index=list(Index))
    queue.put(df)
    return
