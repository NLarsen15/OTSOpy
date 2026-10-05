# Computation

Dictionary: `computation_params` · Class: `ComputationParams`

**Used by:** all functions

Controls parallelism and console output. OTSO parallelizes over stations/grid points using multiple OS processes, and each of those processes can additionally use multiple CPU threads internally (OpenMP) for its own Fortran computation — so a rough rule of thumb is `corenum × threadnum` ≤ the machine's CPU count.

| Key | Type | Default | Description |
|---|---|---|---|
| `corenum` | int or `None` | `None` (auto: logical CPU count − 2, minimum 1) | Number of parallel OS worker processes. Each worker handles a chunk of the stations/grid points. |
| `threadnum` | int | 1 | Number of CPU threads each worker process is pinned to and allowed to use internally (OpenMP-parallelized Fortran routines) for its own computation. |
| `Verbose` | bool | `True` | Print progress (a progress bar where available) while computing. |
| `delim` | str | `";"` | Delimiter used to join multiple asymptotic-direction values into a single output field, when `asymptotic_params["asymptotic"]="YES"`. |
