import pandas as pd
import numpy as np

from ..custom_classes import date
from ..libs import MiddleMan as OTSOLib

def FortranCoordtrans(Data, DateArray, CoordIN, CoordOUT, queue, g, h):
    for x,y in zip(Data,DateArray):
      # Only apply position correction for GDZ input coordinates
      if CoordIN == "GDZ":
          Position = [1+(x[2]/6371.0),x[0],x[1]]
      else:
          Position = x  # Use original position for non-GDZ coordinates

      datetimeobj = date.convert_to_datetime(y)

      year = y[0]
      day = y[1]
      hour = y[2]
      minute = y[3]
      secint = y[4]
      sectot = y[5]
      
      Coords = OTSOLib.coordtrans(Position,year,day,hour,minute,secint,sectot,CoordIN,CoordOUT, g, h)
      combined_array = np.concatenate(([datetimeobj], Position, Coords))

      coord_suffix = f"_{CoordIN}"
      coord_suffix2 = f"_{CoordOUT}"
      
      # Define input coordinate columns based on CoordIN
      if CoordIN == "GDZ":
          input_columns = [f"altitude{coord_suffix} [km]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
      elif CoordIN == "SPH":
          input_columns = [f"altitude{coord_suffix} [Re]", f"latitude{coord_suffix}", f"longitude{coord_suffix}"]
      else:
          input_columns = [f"X{coord_suffix} [Re]", f"Y{coord_suffix} [Re]", f"Z{coord_suffix} [Re]"]
      
      # Define output coordinate columns based on CoordOUT
      if CoordOUT == "GDZ":
          output_columns = [f"altitude{coord_suffix2} [km]", f"latitude{coord_suffix2}", f"longitude{coord_suffix2}"]
      elif CoordOUT == "SPH":
          output_columns = [f"altitude{coord_suffix2} [Re]", f"latitude{coord_suffix2}", f"longitude{coord_suffix2}"]
      else:
          output_columns = [f"X{coord_suffix2} [Re]", f"Y{coord_suffix2} [Re]", f"Z{coord_suffix2} [Re]"]
      
      columns = ["Date"] + input_columns + output_columns
      df = pd.DataFrame(columns=columns)
      df.loc[len(df)] = combined_array
      queue.put(df)

    return