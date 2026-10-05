# Geomagnetic indices

Dictionary: `geomagnetic_params` · Class: `GeomagneticParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Geomagnetic activity indices, required by specific external field models.

| Key | Type | Default | Description |
|---|---|---|---|
| `Dst` | float | 0 | Dst index (nT) — ring current strength. Used by the Boberg extension's Dst-dependent variants and some Tsyganenko models. |
| `kp` | float | 0 | Kp index (0–9). Drives `IOPT`, the storm-level input used by the older Tsyganenko models (TSY89 family). |
| `n_index` | float | 0 | Newell coupling function. Needed by TSY15N. |
| `b_index` | float | 0 | Boynton coupling function. Needed by TSY15B. |
| `sym_h_corrected` | float | 0 | Corrected SYM-H index (nT). Needed by TA16_RBF. |
