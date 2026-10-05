# MHD input examples

Scripts for building gridded fields for OTSO's MHD option and comparing the interpolation schemes. See the "MHD input" page of the
documentation (`docs/MHD.md`) for an overview of the methods.

Only `IGRF_grid_2013.csv` (IGRF 2013 on a 0.5 Re grid over ±10 Re, ~5 MB) is included, for `cutoff_grid_input.py`.
The larger test grids are not stored in the repository (45-570 MB each); generate them with:

| Script | Output |
|---|---|
| `gen_igrf_uniform.py` | `IGRF_grid_uniform_01.csv` (IGRF, uniform 0.1 Re) |
| `gen_igrf_stretched.py` | `IGRF_grid_stretched_r1_<ds0>_ext66.csv` (IGRF, stretched, finer near Earth) |
| `gen_lsmod2_uniform.py` | `LSMOD2_grid_uniform_01.csv` (LSMOD.2, uniform 0.1 Re) |
| `gen_lsmod2_stretched.py` | `LSMOD2_grid_stretched_r1_<ds0>_ext66.csv` (LSMOD.2, stretched) |

`planet_plot.py` / `planet_plot_lsmod2.py` run the cutoff maps for each scheme on those grids and write the
`planet_*.csv` files (the included ones are reference results); `plot_planet.py` plots them.

Note: these grids sample a full internal field through the Earth, which is not divergence-free near the centre, so
`MHDinterpolation="divfree"` warns on them. In normal use the MHD file holds the external (magnetospheric) field and
the internal field comes from `internalmag`.
