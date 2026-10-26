# Latitude-longitude grid

Dictionary: `grid_params` · Class: `GridParams`

**Used by:** `planet`, `trace`

Defines the latitude/longitude grid over which a global field/cutoff map is computed.

| Key | Type | Default | Description |
|---|---|---|---|
| `maxlat` | float | 90 | Northern latitude bound (degrees). |
| `minlat` | float | -90 | Southern latitude bound (degrees). |
| `latstep` | float | -5 | Latitude step (degrees). **Must be negative** to iterate from `maxlat` down to `minlat` — a positive value triggers a warning and will likely produce an empty or unexpected grid. |
| `maxlong` | float | 360 | Eastern longitude bound (degrees). |
| `minlong` | float | 0 | Western longitude bound (degrees). |
| `longstep` | float | 5 | Longitude step (degrees), from `minlong` up to `maxlong`. |
| `array_of_lats_and_longs` | list or `None` | `None` | Optional explicit list of `(latitude, longitude)` pairs to use instead of generating a grid from the bounds/step above. |
| `grid_type` | str | `"regular"` | `planet` only. `"regular"` builds the latitude/longitude grid above. `"healpix"` uses equal-area HEALPix pixel centres. `"healpix_geomag"` uses the same points with the rings laid out in centred-dipole (geomagnetic) coordinates for the requested date. `"fibonacci"` uses a golden-angle (Fibonacci) spiral. The equal-area types ignore the steps and bounds. |
| `nside` | int | 12 | Resolution for the equal-area grid types: `12 * nside**2` points (the HEALPix count, also used for `"fibonacci"`). |

## Equal-area grids

A regular grid puts the same number of points in every latitude row, so the points crowd together towards the poles. A 5° grid has 29% of its points poleward of ±60°, which is only 13% of the Earth's surface, and the cutoff there is near zero and barely changes. HEALPix points each cover the same area.

| `nside` | Points | Spacing |
|---|---|---|
| 10 | 1200 | ~650 km |
| 12 | 1728 | ~550 km (like 5° at the equator) |
| 16 | 3072 | ~410 km |
| 23 | 6348 | ~290 km |

The points are built as equal-area directions from the Earth's centre (geocentric). With `inputcoord` `"GDZ"` (the default) their latitudes are converted to WGS84 geodetic latitude, which moves them up to ~0.19° poleward at mid-latitudes. With `"GEO"` or `"SPH"` they stay geocentric.

The points are not on a latitude/longitude grid. To map them, interpolate on the sphere. A smooth interpolant, e.g. `scipy.interpolate.RBFInterpolator` with a cubic kernel on the 3-D unit vectors, worked best for Ru and Rc, and linear interpolation worked best for Rl and the transmission.

```python
import OTSO

df, *_ = OTSO.planet(grid_params={"grid_type": "healpix_geomag", "nside": 12})
```
