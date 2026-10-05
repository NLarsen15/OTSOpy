# MHD input

OTSO normally uses analytic field models: IGRF or CHAOS for the internal field and the Tsyganenko models for the external field. The MHD option replaces the external model with a user-supplied field on a 3D grid. This is usually the output of a global MHD simulation, but any gridded vector field works. The field is interpolated between grid nodes at every step of a trace.

```python
magfield_params = {"internalmag": "IGRF", "externalmag": "MHD"}
custom_field_params = {
    "MHDfile": "path/to/grid.csv",
    "MHDcoordsys": "GSM",               # frame of the grid file
    "MHDgridtype": "auto",              # "auto" | "uniform" | "stretched"
    "MHDinterpolation": "monotonic",    # "trilinear" | "tricubic" | "monotonic" | "divfree"
}
```

All functions that take `custom_field_params` accept these keys: `cutoff`, `cone`, `planet`, `trajectory`, `flight`, `transmission`, `skymap`, `trace` and `magfield`.

## Interpolation methods

| Method | Description | More detail |
|---|---|---|
| `"trilinear"` (default) | Linear blend of the 8 corners of the grid cell. Fast, but the gradient kinks at every cell face. | [Trilinear interpolation](https://en.wikipedia.org/wiki/Trilinear_interpolation) |
| `"tricubic"` | Cubic fit through a 4×4×4 block of nodes, one axis at a time. Smooth and much more accurate. Can overshoot near sharp changes. | [Cubic Hermite spline](https://en.wikipedia.org/wiki/Cubic_Hermite_spline), [Tricubic interpolation](https://en.wikipedia.org/wiki/Tricubic_interpolation) |
| `"monotonic"` | Tricubic with limited slopes, so the fit never overshoots the data. Recommended for general use. | [Monotone cubic interpolation](https://en.wikipedia.org/wiki/Monotone_cubic_interpolation) |
| `"divfree"` | Interpolates the magnetic vector potential and returns its curl, so the field has zero divergence. Suited to long trapped orbits. | [Magnetic vector potential](https://en.wikipedia.org/wiki/Magnetic_vector_potential), [Mackay et al. (2006)](https://doi.org/10.1029/2005JA011382) |

## Grid file

- A CSV with columns `X, Y, Z` (Earth radii, in the `MHDcoordsys` frame) and `Bx, By, Bz` (nT).
- The grid must be rectilinear: every node is one value from each axis. The spacing does not have to be even.
- Points outside the grid get zero field.
- `MHDgridtype` sets how grid cells are found. `"uniform"` is fastest but only correct for evenly spaced axes. `"stretched"` works for any spacing. `"auto"` (default) checks the spacing and picks one.

## Notes

- In normal use the grid holds only the external (magnetospheric) field, and `internalmag` supplies the internal field. The two are added in a common frame, so the grid can be given in any of OTSO's Cartesian coordinate systems.
- `"divfree"` needs the gridded field to be divergence-free. This is true for an external-only field. A grid that samples a full internal field through the Earth is not, and OTSO gives a `RuntimeWarning` when the grid loads. `"monotonic"` is the better choice for such grids.
