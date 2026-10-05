# Sky map grid

Dictionary: `skymap_params` · Class: `SkymapParams`

**Used by:** `skymap`

Defines the zenith/azimuth grid, centred on each station, that `skymap` computes a cutoff rigidity for.

| Key | Type | Default | Description |
|---|---|---|---|
| `minzenith` | float | 0 | Minimum zenith angle (degrees, 0 = straight up). |
| `maxzenith` | float | 75 | Maximum zenith angle (degrees). |
| `zenithstep` | float | 15 | Zenith step (degrees) between grid points. |
| `minazimuth` | float | 0 | Minimum azimuth angle (degrees, measured from geographic North). |
| `maxazimuth` | float | 360 | Maximum azimuth angle (degrees). |
| `azimuthstep` | float | 45 | Azimuth step (degrees) between grid points. |
