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
