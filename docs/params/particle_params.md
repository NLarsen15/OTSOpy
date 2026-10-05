# Particle

Dictionary: `particle_params` · Class: `ParticleParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `transmission`, `skymap` (not `trace` or `magfield`, which don't trace a physical particle)

The particle species and (for `trajectory`, or `cutoff_comp="Custom"`) its fixed direction/rigidity.

| Key | Type | Default | Description |
|---|---|---|---|
| `Anum` | int | 1 | Atomic number of the traced species: `-1` = muon, `0` = electron, `1` = hydrogen (proton), `2` = helium (alpha particle), `3` = lithium, `4` = beryllium. Values above 4 are not currently supported. |
| `anti` | str | `"YES"` | `"YES"` traces the anti-particle (standard for cosmic-ray cutoff work, since OTSO backtraces from Earth outward); `"NO"` traces the particle itself. |
| `zenith` | float | 0 | Zenith angle (degrees) of the initial direction. Only used when `cutoff_comp="Custom"`. |
| `azimuth` | float | 0 | Azimuth angle (degrees) of the initial direction. Only used when `cutoff_comp="Custom"`. |
| `rigidity` | float | 1 | Fixed rigidity (GV) to trace at. Used by `trajectory`, which traces one particle rather than scanning a rigidity range. |
