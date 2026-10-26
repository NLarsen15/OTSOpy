# Rigidity scan

Dictionary: `rigidity_params` · Class: `RigidityParams`

**Used by:** `cutoff`, `cone`, `planet`, `flight`, `transmission`, `skymap` (not `trajectory`, which traces a single particle at the fixed rigidity given in `particle_params["rigidity"]` rather than scanning a range)

Controls the descending rigidity scan used to locate the cutoff transitions (upper/effective/lower cutoff).

| Key | Type | Default | Description |
|---|---|---|---|
| `startrigidity` | float | 20 | Highest rigidity (GV) the scan starts from. |
| `endrigidity` | float | 0 | Lowest rigidity (GV) the scan stops at. |
| `rigiditystep` | float | 0.01 | Rigidity decrement (GV) between successive traced particles. Smaller = finer cutoff resolution but more traces. |
| `rigidityscan` | str | `"ADAPTIVE"` | `"ADAPTIVE"` performs a rough 0.25 GV scan to bracket the penumbra, scans the bracket ±0.25 GV at `rigiditystep`, then keeps extending that fine scan in 0.25 GV pieces until there is a guard of only forbidden rigidities below the lowest allowed one (1.25 GV, reduced to 0.25 × that rigidity below 5 GV, but at least 0.5 GV) and 0.75 GV of only allowed rigidities above the highest forbidden one. This catches narrow allowed bands far below the main penumbra. `"ON"` uses the same rough scan with a fixed ±0.75 GV fine-scan window (the pre-1.3.11 default; faster, but can miss such bands). `"OFF"` disables scanning. |
