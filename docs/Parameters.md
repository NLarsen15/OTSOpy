# Parameters Reference

Every OTSO function (`cutoff`, `cone`, `planet`, `trajectory`, `flight`, `trace`, `magfield`, `transmission`, `skymap`) takes its settings as a handful of grouped dictionaries — `magfield_params`, `integration_params`, `rigidity_params`, and so on — instead of one long flat argument list. Each function's own docstring lists exactly which of these groups it accepts (not every function accepts every group — see the "Used by" line under each heading below), and this page explains what each key inside those groups actually does.

Every key is optional; anything you don't set falls back to the default shown.

---

## `datetime_params` — `DateTimeParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `trace`, `magfield`, `transmission`, `skymap` (not `flight`, which takes a list of dates directly since a flight path moves through time)

The UTC date/time of the simulation. This drives which IGRF/CHAOS coefficients are interpolated for the internal field and which OMNI/NOAA record is pulled when `serverdata`/`livedata` are enabled.

| Key | Type | Default | Description |
|---|---|---|---|
| `year` | int | 2024 | Year |
| `month` | int | 1 | Month (1–12) |
| `day` | int | 1 | Day (1–31) |
| `hour` | int | 12 | Hour (0–23) |
| `minute` | int | 0 | Minute (0–59) |
| `second` | int | 0 | Second (0–59) |

---

## `magfield_params` — `MagFieldParams`

**Used by:** all functions

Selects and configures the internal (Earth's core) and external (magnetospheric) magnetic field models, and the magnetopause boundary used to decide when a particle has escaped.

| Key | Type | Default | Description |
|---|---|---|---|
| `internalmag` | str | `"IGRF"` | Internal field model. One of `"NONE"`, `"IGRF"`, `"Dipole"`, `"Custom Gauss"`, `"CHAOS"`. `"NONE"` still computes Gauss coefficients but effectively disables the internal contribution; `"Custom Gauss"` requires `g`/`h` in `custom_field_params`. |
| `externalmag` | str | `"TSY89c"` | External (magnetospheric) field model. One of `"NONE"`, `"TSY87short"`, `"TSY87long"`, `"TSY89a"`, `"TSY96"`, `"TSY01"`, `"TSY01S"`, `"TSY04"`, `"TSY89c"`, `"TSY15N"`, `"TSY15B"`, `"TA16_RBF"`, `"TSY89_refit"`, `"MHD"` (requires `MHDfile`). Newer Tsyganenko models need more solar wind/index inputs — see `tsyganenko_params` and `geomagnetic_params`. |
| `boberg` | bool | `False` | Apply the Boberg extension on top of the TSY89 model to account for storm conditions, Boberg extension is TSY89 version specific. |
| `bobergtype` | str | `"EXTENSION"` | Boberg variant to use: `"EXTENSION"`, `"CONTINUOUS"`, `"DST_DEPENDENT"`, `"DST_MIDPOINT"`. Only relevant when `boberg=True`. |
| `magnetopause` | str | `"Kobel"` | Magnetopause shape model used as the outer escape boundary. One of `"NONE"` (disabled — particle never counted as escaped by this check), `"Sphere"` (radius set by `spheresize`), `"aFormisano"`, `"Sibeck"`, `"Kobel"`, `"Lin"`. |
| `spheresize` | float | 25 | Radius (Re) of the spherical magnetopause boundary. Only used when `magnetopause="Sphere"`. |
| `AdaptiveExternalModel` | bool | `False` | Only relevant with `serverdata="ON"`. If the requested external model's required OMNI inputs are missing for the given date, automatically fall back to the next-simplest Tsyganenko model that has valid inputs (e.g. TA16 → TSY15B → ... → TSY89) instead of raising an error. |

---

## `rigidity_params` — `RigidityParams`

**Used by:** `cutoff`, `cone`, `planet`, `flight`, `transmission`, `skymap` (not `trajectory`, which traces a single particle at the fixed rigidity given in `particle_params["rigidity"]` rather than scanning a range)

Controls the descending rigidity scan used to locate the cutoff transitions (upper/effective/lower cutoff).

| Key | Type | Default | Description |
|---|---|---|---|
| `startrigidity` | float | 20 | Highest rigidity (GV) the scan starts from. |
| `endrigidity` | float | 0 | Lowest rigidity (GV) the scan stops at. |
| `rigiditystep` | float | 0.01 | Rigidity decrement (GV) between successive traced particles. Smaller = finer cutoff resolution but more traces. |
| `rigidityscan` | str | `"ON"` | `"ON"` performs a preliminary rough scan of rigidity range to find rough upper and lower bounds of rigidity range to test with user given resolution. `"OFF"` disables scanning. |

---

## `asymptotic_params` — `AsymptoticParams`

**Used by:** `cutoff`, `planet`, `flight`

Computes the asymptotic viewing direction (the direction a particle arrived from, once far from Earth) at a fixed set of energy/rigidity levels, independent of the cutoff rigidity scan.

| Key | Type | Default | Description |
|---|---|---|---|
| `asymptotic` | str | `"NO"` | `"YES"` enables asymptotic direction computation, `"NO"` disables it. |
| `unit` | str | `"GeV"` | Unit the `asymlevels` values are given in: `"GeV"` or `"GV"`. |
| `asymlevels` | list | `[0.1, 0.3, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 50, 70, 100, 300, 500, 700, 1000]` | The list of energy/rigidity levels (in `unit`) to compute an asymptotic direction for. |

---

## `transmission_params` — `TransmissionParams`

**Used by:** `cutoff`, `cone`, `planet`, `flight`, `transmission`, `skymap`

Computes the transmission function (probability that a given rigidity is geomagnetically allowed) by sampling a small band of rigidities around each scanned point rather than a single trace.

| Key | Type | Default | Description |
|---|---|---|---|
| `transmission` | bool | `False` | Enable transmission function computation. |
| `transmissionRstep` | float | 0.001 | Half-width (GV) of the sampling band around each tested rigidity: samples are drawn from `R ± transmissionRstep`. |
| `transmissionsamples` | int | 20 | Number of rigidities sampled within that band. |

---

## `solar_wind_params` — `SolarWindParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Manual solar wind and IMF values, used when not pulling data from the OMNI server or live NOAA feed. Several external field models (see `magfield_params["externalmag"]`) need specific fields here to be non-default.

| Key | Type | Default | Description |
|---|---|---|---|
| `vx` | float | -500 | Solar wind velocity x-component (km/s). OTSO expects this pointed sunward-to-Earth; a positive value is auto-flipped negative internally. |
| `vy` | float | 0 | Solar wind velocity y-component (km/s). |
| `vz` | float | 0 | Solar wind velocity z-component (km/s). |
| `bx` | float | 0 | IMF x-component (nT). |
| `by` | float | 5 | IMF y-component (nT). |
| `bz` | float | 5 | IMF z-component (nT). Southward (negative) IMF drives stronger geomagnetic activity. |
| `by_avg` | float | 0 | IMF By averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `bz_avg` | float | 0 | IMF Bz averaged over the preceding 30 minutes (nT). Needed by TSY15N/TSY15B/TA16. |
| `density` | float | 1 | Solar wind proton density (particles/cm³). |
| `pdyn` | float | 0 | Solar wind dynamic pressure (nPa). |

---

## `geomagnetic_params` — `GeomagneticParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

Geomagnetic activity indices, required by specific external field models.

| Key | Type | Default | Description |
|---|---|---|---|
| `Dst` | float | 0 | Dst index (nT) — ring current strength. Used by the Boberg extension's Dst-dependent variants and some Tsyganenko models. |
| `kp` | float | 0 | Kp index (0–9). Drives `IOPT`, the storm-level input used by the older Tsyganenko models (TSY89 family). |
| `n_index` | float | 0 | Newell coupling function. Needed by TSY15N. |
| `b_index` | float | 0 | Boynton coupling function. Needed by TSY15B. |
| `sym_h_corrected` | float | 0 | Corrected SYM-H index (nT). Needed by TA16_RBF. |

---

## `tsyganenko_params` — `TsyganenkoParams`

**Used by:** all functions (ignored if `serverdata`/`livedata` supply these instead)

The G1–G3 and W1–W6 driving coefficients used by the newer Tsyganenko external field models (TSY01 storm, TSY04, and later).

| Key | Type | Default | Description |
|---|---|---|---|
| `G1` | float | 0 | Tsyganenko G1 coefficient (TSY01/TSY01S). |
| `G2` | float | 0 | Tsyganenko G2 coefficient (TSY01). |
| `G3` | float | 0 | Tsyganenko G3 coefficient (TSY01S). |
| `W1`–`W6` | float | 0 | Tsyganenko W1–W6 coefficients (TSY04) — cumulative solar wind forcing terms. |

---

## `integration_params` — `IntegrationParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `transmission`, `skymap` in full; `trace` accepts only `minaltitude`, `maxdistance`, `maxtime`, `maxsteps`, `startaltitude`, `fixedstep` (field-line tracing doesn't do particle-style adaptive/beta-checked integration, so `intmodel`/`gyropercent`/`betaerror`/`totalbetacheck`/`adaptivestep`/`mintrapdist` don't apply); not used by `magfield`

Controls how a particle's trajectory is numerically integrated through the field, and when a trace is stopped.

| Key | Type | Default | Description |
|---|---|---|---|
| `intmodel` | str | `"Boris-Buneman"` | Integration method: `"4RK"`, `"5RK"`, `"6RK"` (Runge-Kutta, 4th/5th/6th order), `"Boris-Buneman"`, `"Vay"`, `"HC"` (Higuera-Cary).|
| `gyropercent` | float | 5 | Adaptive step-size cap, as a percentage of the particle's local gyration (cyclotron) period. Smaller values give a finer, more accurate integration step but take longer; only applies when `adaptivestep=True`. |
| `minaltitude` | float | 20 | Termination boundary (km if `inputcoord="GDZ"`, otherwise Re): the trace stops once the particle's radial position drops below this altitude ("Earth encounter" — e.g. the particle is absorbed by the atmosphere). |
| `startaltitude` | float | 20 | Altitude (same units as `minaltitude`) the trace begins at, typically 20km is selected as the average altitude at which particle interact with the atmosphere. Distinct from `minaltitude`, which only governs when the trace *stops*. |
| `maxdistance` | float | 100 | Termination boundary (Re): the trace stops once the particle has travelled farther than this from Earth (it is assumed a forbidden trapped trajectory). |
| `maxtime` | float | 0 | Termination boundary (s) on total integration time. `0` disables this check. Acts as a safety net against traces that never hit another stopping condition. |
| `maxsteps` | int | 0 (`None` at the public API level) | Termination boundary on the number of integration steps taken. `0`/unset disables this check. Another safety net alongside `maxtime`. |
| `mintrapdist` | float | 0 | Reference radial distance (Re) used to detect trapped particles: if the particle's radial distance never grows past this value over the course of the trace, it's flagged as trapped rather than escaping or hitting the atmosphere. |
| `betaerror` | float | 0.001 | Maximum allowed fractional change in the particle's speed over a single integration step, as a percentage. Speed should be conserved in a static magnetic field (the Lorentz force does no work), so a larger-than-allowed change means the step was too coarse — the integrator rejects and retries it with a smaller step. |
| `totalbetacheck` | bool | `False` | Additionally checks the particle's *current* speed against its speed at the very start of the trace (not just step-to-step). If the cumulative drift exceeds `betaerror`, the integrator backs off and retries — catches slow numerical drift that per-step checks alone can miss. |
| `adaptivestep` | bool | `True` | Use an adaptive step size (capped by `gyropercent`, tightened by `betaerror`) instead of a fixed one. |
| `fixedstep` | float | 0.0 | Fixed step size in seconds, used only when `adaptivestep=False` and this is `> 0`. `0` (default) falls back to a `gyropercent`-derived step even with `adaptivestep=False`. |

---

## `particle_params` — `ParticleParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `transmission`, `skymap` (not `trace` or `magfield`, which don't trace a physical particle)

The particle species and (for `trajectory`, or `cutoff_comp="Custom"`) its fixed direction/rigidity.

| Key | Type | Default | Description |
|---|---|---|---|
| `Anum` | int | 1 | Atomic number of the traced species: `-1` = muon, `0` = electron, `1` = hydrogen (proton), `2` = helium (alpha particle), `3` = lithium, `4` = beryllium. Values above 4 are not currently supported. |
| `anti` | str | `"YES"` | `"YES"` traces the anti-particle (standard for cosmic-ray cutoff work, since OTSO backtraces from Earth outward); `"NO"` traces the particle itself. |
| `zenith` | float | 0 | Zenith angle (degrees) of the initial direction. Only used when `cutoff_comp="Custom"`. |
| `azimuth` | float | 0 | Azimuth angle (degrees) of the initial direction. Only used when `cutoff_comp="Custom"`. |
| `rigidity` | float | 1 | Fixed rigidity (GV) to trace at. Used by `trajectory`, which traces one particle rather than scanning a rigidity range. |

---

## `coordinate_params` — `CoordinateParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `magfield`, `transmission`, `skymap`, `trace` (`trace` only uses `inputcoord`; `coordsystem`/`coordout` don't apply to field-line tracing)

Input/output coordinate systems for station locations and results.

| Key | Type | Default | Description |
|---|---|---|---|
| `inputcoord` | str | `"GDZ"` | Coordinate system station/location inputs are given in. One of `"GDZ"`, `"GEO"`, `"GSM"`, `"GSE"`, `"SM"`, `"GEI"`, `"MAG"`, `"SPH"`, `"RLL"`. `"GDZ"` (geodetic) is the natural choice when giving latitude/longitude/altitude. |
| `coordsystem` | str | `"GEO"`/`"GSM"` (varies by function — check its docstring) | Coordinate system results (e.g. trajectory positions) are reported in. Same 9 options as `inputcoord`. |
| `coordout` | str | `"GSM"` | `magfield`-specific output coordinate system. Cartesian only: `"GEO"`, `"GSM"`, `"GSE"`, `"SM"`, `"GEI"`, `"MAG"`, `"RLL"` (no `"GDZ"`/`"SPH"`). |

---

## `computation_params` — `ComputationParams`

**Used by:** all functions

Controls parallelism and console output. OTSO parallelizes over stations/grid points using multiple OS processes, and each of those processes can additionally use multiple CPU threads internally (OpenMP) for its own Fortran computation — so a rough rule of thumb is `corenum × threadnum` ≤ your machine's CPU count.

| Key | Type | Default | Description |
|---|---|---|---|
| `corenum` | int or `None` | `None` (auto: logical CPU count − 2, minimum 1) | Number of parallel OS worker processes. Each worker handles a chunk of the stations/grid points. |
| `threadnum` | int | 1 | Number of CPU threads each worker process is pinned to and allowed to use internally (OpenMP-parallelized Fortran routines) for its own computation. |
| `Verbose` | bool | `True` | Print progress (a progress bar where available) while computing. |
| `delim` | str | `";"` | Delimiter used to join multiple asymptotic-direction values into a single output field, when `asymptotic_params["asymptotic"]="YES"`. |

---

## `data_retrieval_params` — `DataRetrievalParams`

**Used by:** all functions

Automatic space-weather data retrieval, as an alternative to manually filling in `solar_wind_params`/`geomagnetic_params`/`tsyganenko_params`.

| Key | Type | Default | Description |
|---|---|---|---|
| `serverdata` | str | `"OFF"` | `"ON"` downloads and uses the finalised OMNI database values for the given `datetime_params`, overriding manually-supplied solar wind/geomagnetic/Tsyganenko values. |
| `livedata` | str | `"OFF"` | `"ON"` pulls preliminary real-time data from NOAA instead. NOAA values are provisional and may later differ from the finalised OMNI database (`serverdata` takes precedence if both would apply). Not supported for `externalmag="TSY04"` or `"TA16_RBF"`. |

---

## `custom_field_params` — `CustomFieldParams`

**Used by:** all functions

Advanced/custom internal field configuration.

| Key | Type | Default | Description |
|---|---|---|---|
| `g` | list or `None` | `None` | Custom Gauss `g` coefficients. Required when `internalmag="Custom Gauss"`. |
| `h` | list or `None` | `None` | Custom Gauss `h` coefficients. Required when `internalmag="Custom Gauss"`. |
| `max_degree` | int | 13 | Maximum spherical harmonic degree/order used when expanding the internal field from Gauss coefficients. Higher = more accurate near Earth but more computation; capped at the underlying coefficient table's degree (IGRF: 13; `"Dipole"` is forced to degree 1 regardless of this value). |
| `MHDfile` | str or `None` | `None` | Path to an MHD simulation output file. Required when `externalmag="MHD"`; OTSO validates the file exists before running. |
| `MHDcoordsys` | str or `None` | `None` | Coordinate system the MHD file's field data is given in. |

---

## `grid_params` — `GridParams`

**Used by:** `planet`, `trace`

Defines the latitude/longitude grid over which a global field/cutoff map is computed.

| Key | Type | Default | Description |
|---|---|---|---|
| `maxlat` | float | 90 | Northern latitude bound (degrees). |
| `minlat` | float | -90 | Southern latitude bound (degrees). |
| `latstep` | float | -5 | Latitude step (degrees). **Must be negative** to iterate from `maxlat` down to `minlat` — a positive value triggers a warning and will likely produce an empty or unexpected grid. |
| `maxlong` | float | 360 | Eastern longitude bound (degrees). |
| `minlong` | float | 0 | Western longitude bound (degrees). |
| `longstep` | float | 5 | Longitude step (degrees), from `minlong` up to `maxlong`. |
| `array_of_lats_and_longs` | list or `None` | `None` | Optional explicit list of `(latitude, longitude)` pairs to use instead of generating a grid from the bounds/step above. |

---

## `skymap_params` — `SkymapParams`

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
