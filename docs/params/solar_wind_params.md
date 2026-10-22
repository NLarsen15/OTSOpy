# Solar wind

Dictionary: `solar_wind_params` · Class: `SolarWindParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Manual solar wind and IMF values, used when not pulling data from the OMNI server or live NOAA feed. Several external field models (see `magfield_params["externalmag"]`) need specific fields here to be non-default.

Defaults are quiet-time values: the median of each input over OMNI data from 1996 to 2025 when Kp ≤ 1. OTSO uses them for any key you leave out, unless `serverdata` or `livedata` is on. See [Default space weather conditions](../Parameters.md#default-space-weather-conditions).


| Key | Type | Default | Description |
|---|---|---|---|
| `vx` | float | -359.4 | Solar wind velocity x-component (km/s). OTSO expects this pointed sunward-to-Earth; a positive value is auto-flipped negative internally. |
| `vy` | float | 0 | Solar wind velocity y-component (km/s). |
| `vz` | float | 0 | Solar wind velocity z-component (km/s). |
| `bx` | float | 0.11 | IMF x-component (nT). |
| `by` | float | -0.16 | IMF y-component (nT). |
| `bz` | float | 0.58 | IMF z-component (nT). Southward (negative) IMF drives stronger geomagnetic activity. |
| `by_avg` | float | -0.17 | IMF By averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `bz_avg` | float | 0.59 | IMF Bz averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `density` | float | 4.94 | Solar wind proton density (particles/cm³). |
| `pdyn` | float | 1.31 | Solar wind dynamic pressure (nPa). If you set it to 0, or give `density` or a velocity without `pdyn`, OTSO computes it from the density and velocity. |
