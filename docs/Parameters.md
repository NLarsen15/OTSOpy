# Input parameters

Each OTSO function takes its settings as a few grouped dictionaries, such as `magfield_params` or `integration_params`.
Every key is optional. Any key left out uses the default listed on its page.
Not every function accepts every group. Each page has a "Used by" line, and each function's page lists the groups it takes.

```python
from OTSO import cutoff

results = cutoff(
    Stations=["OULU"],
    magfield_params={"externalmag": "TSY89c"},
    integration_params={"gyropercent": 5},
)
```

## Field models

| Dictionary | Class | What it sets |
|---|---|---|
| [`magfield_params`](params/magfield_params.md) | `MagFieldParams` | Internal and external field models and the magnetopause. |
| [`solar_wind_params`](params/solar_wind_params.md) | `SolarWindParams` | Solar wind and IMF values for the external field models. |
| [`geomagnetic_params`](params/geomagnetic_params.md) | `GeomagneticParams` | Dst, Kp and other indices for the external field models. |
| [`tsyganenko_params`](params/tsyganenko_params.md) | `TsyganenkoParams` | G and W driving coefficients for the newer Tsyganenko models. |
| [`custom_field_params`](params/custom_field_params.md) | `CustomFieldParams` | Custom Gauss coefficients and MHD grid input. |

## Default space weather conditions

The external field models need solar wind, IMF and geomagnetic inputs. You can give them yourself, or have OTSO look them up for the date with `data_retrieval_params={"serverdata": "ON"}` (OMNI) or `{"livedata": "ON"}` (NOAA). Any input that is not supplied either way falls back to a default.

The defaults describe typical quiet-time conditions. Each one is the median of that input over every 5-minute interval from 1996 to 2025 with Kp ≤ 1 (38% of the record). The values come from OTSO's server data files, which hold OMNI high- and low-resolution data and the derived Tsyganenko inputs (G1–G3, W1–W6, N-index, B-index and corrected SYM-H). OMNI has no median for `vy` and `vz`, so they default to 0.

| Input | Default | Quiet-time 25th–75th percentile |
|---|---|---|
| `vx` (km/s) | -359.4 | -404.1 to -325.7 |
| `bx`, `by`, `bz` (nT) | 0.11, -0.16, 0.58 | Bz: -0.62 to 1.94 |
| `by_avg`, `bz_avg` (nT) | -0.17, 0.59 | Bz: -0.53 to 1.87 |
| `density` (cm⁻³) | 4.94 | 3.30 to 7.40 |
| `pdyn` (nPa) | 1.31 | 0.96 to 1.79 |
| `Dst` (nT) | -4 | -10 to 2 |
| `kp` | 0.7 | 0.3 to 1.0 |
| `n_index`, `b_index` | 0.145, 0.032 | N: 0.054 to 0.265 |
| `sym_h_corrected` (nT) | -18.07 | -22.35 to -14.40 |
| `G1`, `G2`, `G3` | 0.38, 0.12, 0.47 | G3: 0 to 5.61 |
| `W1`–`W6` | 0.07, 0.06, 0.13, 0.03, 0.05, 0.08 | W6: 0.02 to 0.24 |

Defaults are meant for quick, representative runs. For a specific event, set the inputs or turn on `serverdata`.

When a solar-wind-driven Tsyganenko model (TSY96, TSY01, TSY01S, TSY04, TSY15N, TSY15B, TA16_RBF) runs on any default input, OTSO prints a warning naming those inputs. The metadata output lists the inputs used and marks each default with `(default)` (for every function except `flight`, whose inputs can vary along the path).

## Particle and integration

| Dictionary | Class | What it sets |
|---|---|---|
| [`particle_params`](params/particle_params.md) | `ParticleParams` | Particle species, direction and rigidity. |
| [`rigidity_params`](params/rigidity_params.md) | `RigidityParams` | The rigidity range and step used to find cutoffs. |
| [`integration_params`](params/integration_params.md) | `IntegrationParams` | Integration method, step size and stopping conditions. |

## Outputs

| Dictionary | Class | What it sets |
|---|---|---|
| [`asymptotic_params`](params/asymptotic_params.md) | `AsymptoticParams` | Asymptotic viewing directions at fixed rigidity or energy levels. |
| [`transmission_params`](params/transmission_params.md) | `TransmissionParams` | Transmission function settings. |
| [`coordinate_params`](params/coordinate_params.md) | `CoordinateParams` | Coordinate systems for inputs and results. |

## Run settings

| Dictionary | Class | What it sets |
|---|---|---|
| [`datetime_params`](params/datetime_params.md) | `DateTimeParams` | The UTC date and time of the run. |
| [`computation_params`](params/computation_params.md) | `ComputationParams` | Number of processes and threads, and console output. |
| [`data_retrieval_params`](params/data_retrieval_params.md) | `DataRetrievalParams` | Automatic download of solar wind and index data. |

## Grids

| Dictionary | Class | What it sets |
|---|---|---|
| [`grid_params`](params/grid_params.md) | `GridParams` | The global grid used by planet and trace. |
| [`skymap_params`](params/skymap_params.md) | `SkymapParams` | The zenith and azimuth grid used by skymap. |
