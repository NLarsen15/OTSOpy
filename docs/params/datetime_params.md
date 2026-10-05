# Date and time

Dictionary: `datetime_params` · Class: `DateTimeParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `trace`, `magfield`, `transmission`, `skymap` (not `flight`, which takes a list of dates directly since a flight path moves through time)

The UTC date/time of the simulation. This drives which IGRF/CHAOS coefficients are interpolated for the internal field and which OMNI/NOAA record is pulled when `serverdata`/`livedata` are enabled.

| Key | Type | Default | Description |
|---|---|---|---|
| `year` | int | 2024 | Year |
| `month` | int | 1 | Month (1–12) |
| `day` | int | 1 | Day (1–31) |
| `hour` | int | 12 | Hour (0–23) |
| `minute` | int | 0 | Minute (0–59) |
| `second` | int | 0 | Second (0–59) |
