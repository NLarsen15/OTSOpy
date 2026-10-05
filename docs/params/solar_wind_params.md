# Solar wind

Dictionary: `solar_wind_params` · Class: `SolarWindParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Manual solar wind and IMF values, used when not pulling data from the OMNI server or live NOAA feed. Several external field models (see `magfield_params["externalmag"]`) need specific fields here to be non-default.

| Key | Type | Default | Description |
|---|---|---|---|
| `vx` | float | -500 | Solar wind velocity x-component (km/s). OTSO expects this pointed sunward-to-Earth; a positive value is auto-flipped negative internally. |
| `vy` | float | 0 | Solar wind velocity y-component (km/s). |
| `vz` | float | 0 | Solar wind velocity z-component (km/s). |
| `bx` | float | 0 | IMF x-component (nT). |
| `by` | float | 5 | IMF y-component (nT). |
| `bz` | float | 5 | IMF z-component (nT). Southward (negative) IMF drives stronger geomagnetic activity. |
| `by_avg` | float | 0 | IMF By averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `bz_avg` | float | 0 | IMF Bz averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `density` | float | 1 | Solar wind proton density (particles/cm³). |
| `pdyn` | float | 0 | Solar wind dynamic pressure (nPa). |
