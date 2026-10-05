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
