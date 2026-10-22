# Geomagnetic indices

Dictionary: `geomagnetic_params` · Class: `GeomagneticParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Geomagnetic activity indices, required by specific external field models.

Defaults are quiet-time values: the median of each input over OMNI data from 1996 to 2025 when Kp ≤ 1. OTSO uses them for any key you leave out, unless `serverdata` or `livedata` is on. See [Default space weather conditions](../Parameters.md#default-space-weather-conditions).


| Key | Type | Default | Description |
|---|---|---|---|
| `Dst` | float | -4 | Dst index (nT) — ring current strength. Used by the Boberg extension's Dst-dependent variants and some Tsyganenko models. |
| `kp` | float | 0.7 | Kp index (0–9). Drives `IOPT`, the storm-level input used by the older Tsyganenko models (TSY89 family). |
| `n_index` | float | 0.145 | Newell coupling function. Needed by TSY15N. |
| `b_index` | float | 0.032 | Boynton coupling function. Needed by TSY15B. |
| `sym_h_corrected` | float | -18.07 | Corrected SYM-H index (nT). Needed by TA16_RBF. |
