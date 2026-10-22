# Tsyganenko coefficients

Dictionary: `tsyganenko_params` · Class: `TsyganenkoParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

The G1–G3 and W1–W6 driving coefficients used by the newer Tsyganenko external field models (TSY01 storm, TSY04, and later).

Defaults are quiet-time values: the median of each input over OMNI data from 1996 to 2025 when Kp ≤ 1. OTSO uses them for any key you leave out, unless `serverdata` or `livedata` is on. See [Default space weather conditions](../Parameters.md#default-space-weather-conditions).


| Key | Type | Default | Description |
|---|---|---|---|
| `G1` | float | 0.38 | Tsyganenko G1 coefficient (TSY01/TSY01S). |
| `G2` | float | 0.12 | Tsyganenko G2 coefficient (TSY01). |
| `G3` | float | 0.47 | Tsyganenko G3 coefficient (TSY01S). |
| `W1`–`W6` | float | 0.07, 0.06, 0.13, 0.03, 0.05, 0.08 | Tsyganenko W1–W6 coefficients (TSY04) — cumulative solar wind forcing terms. |
