# Coordinate systems

Dictionary: `coordinate_params` · Class: `CoordinateParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `magfield`, `transmission`, `skymap`, `trace` (`trace` only uses `inputcoord`; `coordsystem`/`coordout` don't apply to field-line tracing)

Input/output coordinate systems for station locations and results.

| Key | Type | Default | Description |
|---|---|---|---|
| `inputcoord` | str | `"GDZ"` | Coordinate system station/location inputs are given in. One of `"GDZ"`, `"GEO"`, `"GSM"`, `"GSE"`, `"SM"`, `"GEI"`, `"MAG"`, `"SPH"`, `"RLL"`. `"GDZ"` (geodetic) is the natural choice when giving latitude/longitude/altitude. |
| `coordsystem` | str | `"GEO"`/`"GSM"` (varies by function — check its docstring) | Coordinate system results (e.g. trajectory positions) are reported in. Same 9 options as `inputcoord`. |
| `coordout` | str | `"GSM"` | `magfield`-specific output coordinate system. Cartesian only: `"GEO"`, `"GSM"`, `"GSE"`, `"SM"`, `"GEI"`, `"MAG"`, `"RLL"` (no `"GDZ"`/`"SPH"`). |
