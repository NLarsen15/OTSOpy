# Magnetic field models

Dictionary: `magfield_params` · Class: `MagFieldParams`

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
