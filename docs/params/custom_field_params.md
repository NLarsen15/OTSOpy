# Custom field and MHD grid

Dictionary: `custom_field_params` · Class: `CustomFieldParams`

**Used by:** all functions

Custom internal field (Gauss coefficients) and MHD grid configuration.

| Key | Type | Default | Description |
|---|---|---|---|
| `g` | list or `None` | `None` | Custom Gauss `g` coefficients. Required when `internalmag="Custom Gauss"`. |
| `h` | list or `None` | `None` | Custom Gauss `h` coefficients. Required when `internalmag="Custom Gauss"`. |
| `max_degree` | int | 13 | Maximum spherical harmonic degree/order used when expanding the internal field from Gauss coefficients. Higher = more accurate near Earth but more computation; capped at the underlying coefficient table's degree (IGRF: 13; `"Dipole"` is forced to degree 1 regardless of this value). |
| `MHDfile` | str or `None` | `None` | Path to an MHD simulation output file. Required when `externalmag="MHD"`; OTSO validates the file exists before running. |
| `MHDcoordsys` | str or `None` | `None` | Coordinate system the MHD file's grid positions and field components are given in (any of OTSO's Cartesian systems, e.g. `"GSM"`, `"GEO"`). OTSO rotates positions into this frame and the interpolated field back, so it can be combined with any `internalmag`. |
| `MHDgridtype` | str | `"auto"` | How grid cells are found: `"auto"` detects per-axis spacing, `"uniform"` forces the fast fixed-spacing lookup (only correct if every axis is evenly spaced), `"stretched"` forces the general lookup for grids whose spacing varies (e.g. finer near Earth). |
| `MHDinterpolation` | str | `"trilinear"` | How the field is reconstructed between grid points: `"trilinear"` (fast, kinked gradients at cell faces), `"tricubic"` (smooth, much more accurate in smooth regions, can overshoot near sharp features), `"monotonic"` (tricubic with overshoot limiting; recommended), `"divfree"` (curl of an interpolated vector potential; exactly divergence-free, for long quasi-trapped trajectories, less point-wise accurate). Works on uniform and stretched grids. See [MHD input](../MHD.md). |
