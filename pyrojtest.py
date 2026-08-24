from pyproj import Transformer
import numpy as np

lla_to_geo = Transformer.from_crs(
    {"proj": "longlat", "ellps": "WGS84"},
    {"proj": "geocent", "ellps": "WGS84"}
)

lat = 0
long = -50
alt = 382

alt2 = 382*1e3

RE = 6371.2e3
RE = 6378.137e3

r = lla_to_geo.transform(long, lat, alt * 1e3, radians=False)

print(r)

earth_radius = (np.linalg.norm(r) - alt2) / 1e3 # m -> km
print(earth_radius)

rnew = [r[0] / RE, r[1] / RE, r[2] / RE]

print(rnew)

rmag = (rnew[0]**2 + rnew[1]**2 +rnew[2]**2)**0.5

print(rmag)