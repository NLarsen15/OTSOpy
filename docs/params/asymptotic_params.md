# Asymptotic directions

Dictionary: `asymptotic_params` · Class: `AsymptoticParams`

**Used by:** `cutoff`, `planet`, `flight`

Computes the asymptotic viewing direction (the direction a particle arrived from, once far from Earth) at a fixed set of energy/rigidity levels, independent of the cutoff rigidity scan.

| Key | Type | Default | Description |
|---|---|---|---|
| `asymptotic` | str | `"NO"` | `"YES"` enables asymptotic direction computation, `"NO"` disables it. |
| `unit` | str | `"GeV"` | Unit the `asymlevels` values are given in: `"GeV"` or `"GV"`. |
| `asymlevels` | list | `[0.1, 0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 70, 100, 300, 500, 700, 1000]` | The list of energy/rigidity levels (in `unit`) to compute an asymptotic direction for. |
