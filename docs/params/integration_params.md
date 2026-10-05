# Integration

Dictionary: `integration_params` · Class: `IntegrationParams`

**Used by:** `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `transmission`, `skymap` in full; `trace` accepts only `minaltitude`, `maxdistance`, `maxtime`, `maxsteps`, `startaltitude`, `fixedstep` (field-line tracing doesn't do particle-style adaptive/beta-checked integration, so `intmodel`/`gyropercent`/`betaerror`/`totalbetacheck`/`adaptivestep`/`mintrapdist` don't apply); not used by `magfield`

Controls how a particle's trajectory is numerically integrated through the field, and when a trace is stopped.

| Key | Type | Default | Description |
|---|---|---|---|
| `intmodel` | str | `"Boris-Buneman"` | Integration method: `"4RK"`, `"5RK"`, `"6RK"` (Runge-Kutta, 4th/5th/6th order; `"5RK"` is Nyström's classic 6-stage 5th-order scheme, with no embedded error estimate), `"Boris-Buneman"`, `"Vay"`, `"HC"` (Higuera-Cary). The last three are second-order leapfrog schemes. |
| `gyropercent` | float | 10 (1 for `cone`, `trajectory`, `skymap`) | Step size as a percentage of the particle's gyration (cyclotron) period. With `adaptivestep=False` the step is fixed at this percentage of the gyration period at the launch point; with `adaptivestep=True` it is the cap on each step, using the local field. Smaller values are more accurate but slower. |
| `minaltitude` | float | 20 | Termination boundary (km if `inputcoord="GDZ"`, otherwise Re): the trace stops once the particle's radial position drops below this altitude ("Earth encounter" — e.g. the particle is absorbed by the atmosphere). |
| `startaltitude` | float | 20 | Altitude (same units as `minaltitude`) the trace begins at, typically 20km is selected as the average altitude at which particle interact with the atmosphere. Distinct from `minaltitude`, which only governs when the trace *stops*. |
| `maxdistance` | float | 100 | Termination boundary (Re): the trace stops once the total path length travelled exceeds this (it is assumed a forbidden trapped trajectory). Note this is path length, not distance from Earth. |
| `maxtime` | float | 0 | Termination boundary on total flight time in seconds (lab frame). `0` disables this check. Acts as a safety net against traces that never hit another stopping condition. |
| `maxsteps` | int | 0 | Termination boundary on the number of integration steps taken. `0` disables this check. Another safety net alongside `maxtime`. |
| `mintrapdist` | float | 0 | Reference radial distance (Re) used to detect trapped particles: if the particle's radial distance never grows past this value over the course of the trace, it's flagged as trapped rather than escaping or hitting the atmosphere. |
| `betaerror` | float | 0.001 | Maximum allowed fractional change in the particle's speed over a single integration step, as a percentage. Speed should be conserved in a static magnetic field (the Lorentz force does no work), so a larger-than-allowed change means the step was too coarse — the integrator rejects and retries it with a smaller step. |
| `totalbetacheck` | bool | `False` | Additionally checks the particle's *current* speed against its speed at the very start of the trace (not just step-to-step). If the cumulative drift exceeds `betaerror`, the integrator backs off and retries — catches slow numerical drift that per-step checks alone can miss. |
| `adaptivestep` | bool | `False` | Use an adaptive step size (capped by `gyropercent`, tightened by `betaerror`) instead of a fixed one. |
| `fixedstep` | float | 0.0 | Fixed step size in seconds, used only when `adaptivestep=False` and this is `> 0`. `0` (default) falls back to a `gyropercent`-derived step even with `adaptivestep=False`. |
