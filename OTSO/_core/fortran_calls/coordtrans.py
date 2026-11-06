import pandas as pd
import numpy as np

from ..custom_classes import date
from ..utils import fortran_data_utils
from ..libs.MiddleMan import Middleman as OTSOLib


def coordtrans_columns(CoordIN: str, CoordOUT: str) -> list:
    def cols(c):
        if c == "GDZ":
            return [f"altitude_{c} [km]", f"latitude_{c}", f"longitude_{c}"]
        if c == "SPH":
            return [f"altitude_{c} [Re]", f"latitude_{c}", f"longitude_{c}"]
        return [f"X_{c} [Re]", f"Y_{c} [Re]", f"Z_{c} [Re]"]
    return ["Date"] + cols(CoordIN) + cols(CoordOUT)


def FortranCoordtrans(Data, DateArray, CoordIN, CoordOUT, queue, maxdegree, g, h, Index):
    """Transform every location in Data; puts ONE DataFrame (all of this worker's points) on the queue.

    Results go into arrays instead of a one-row DataFrame and a queue message per point, and the
    Fortran inputs are only rebuilt when the date changes from one point to the next."""
    n = len(Data)
    positions = np.empty((n, 3), dtype=np.float64)
    coords_out = np.empty((n, 3), dtype=np.float64)
    dates = []
    Coords = np.zeros(3, dtype=np.float64)
    last_key, FortranData = None, None
    for k, (x, y, g_point, h_point) in enumerate(zip(Data, DateArray, g, h)):
        # GDZ input is given as [lat, lon, alt_km]; the Fortran side expects [alt_km, lat, lon]
        Position = [x[2], x[0], x[1]] if CoordIN == "GDZ" else x
        key = tuple(np.asarray(y).ravel().tolist())
        if key != last_key:
            FortranData = fortran_data_utils.prepare_fortran_Coordtrans(y, maxdegree, g_point, h_point)
            datetimeobj = date.convert_to_datetime(y)
            last_key = key
        OTSOLib.coordtrans(Position, FortranData, CoordIN, CoordOUT, g_point, h_point, Coords)
        positions[k] = Position
        coords_out[k] = Coords
        dates.append(datetimeobj)

    columns = coordtrans_columns(CoordIN, CoordOUT)
    df = pd.DataFrame(np.hstack([positions, coords_out]), columns=columns[1:], index=list(Index))
    df.insert(0, "Date", pd.Series(dates, index=df.index, dtype=object))
    queue.put(df)
    return
