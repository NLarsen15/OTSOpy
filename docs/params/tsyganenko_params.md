# Tsyganenko coefficients

Dictionary: `tsyganenko_params` · Class: `TsyganenkoParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

The G1–G3 and W1–W6 driving coefficients used by the newer Tsyganenko external field models (TSY01 storm, TSY04, and later).

| Key | Type | Default | Description |
|---|---|---|---|
| `G1` | float | 0 | Tsyganenko G1 coefficient (TSY01/TSY01S). |
| `G2` | float | 0 | Tsyganenko G2 coefficient (TSY01). |
| `G3` | float | 0 | Tsyganenko G3 coefficient (TSY01S). |
| `W1`–`W6` | float | 0 | Tsyganenko W1–W6 coefficients (TSY04) — cumulative solar wind forcing terms. |
