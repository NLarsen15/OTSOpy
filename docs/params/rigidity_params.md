# Rigidity scan

Dictionary: `rigidity_params` · Class: `RigidityParams`

**Used by:** `cutoff`, `cone`, `planet`, `flight`, `transmission`, `skymap` (not `trajectory`, which traces a single particle at the fixed rigidity given in `particle_params["rigidity"]` rather than scanning a range)

Controls the descending rigidity scan used to locate the cutoff transitions (upper/effective/lower cutoff).

| Key | Type | Default | Description |
|---|---|---|---|
| `startrigidity` | float | 20 | Highest rigidity (GV) the scan starts from. |
| `endrigidity` | float | 0 | Lowest rigidity (GV) the scan stops at. |
| `rigiditystep` | float | 0.01 | Rigidity decrement (GV) between successive traced particles. Smaller = finer cutoff resolution but more traces. |
| `rigidityscan` | str | `"ON"` | `"ON"` performs a preliminary rough scan of rigidity range to find rough upper and lower bounds of rigidity range to test with user given resolution. `"OFF"` disables scanning. |
