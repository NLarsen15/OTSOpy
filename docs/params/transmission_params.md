# Transmission function

Dictionary: `transmission_params` · Class: `TransmissionParams`

**Used by:** `cutoff`, `cone`, `planet`, `flight`, `transmission`, `skymap`

Computes the transmission function (probability that a given rigidity is geomagnetically allowed) by sampling a small band of rigidities around each scanned point rather than a single trace.

| Key | Type | Default | Description |
|---|---|---|---|
| `transmission` | bool | `False` | Enable transmission function computation. |
| `transmissionRstep` | float | 0.001 | Half-width (GV) of the sampling band around each tested rigidity: samples are drawn from `R ± transmissionRstep`. |
| `transmissionsamples` | int | 20 | Number of rigidities sampled within that band. |
