import time
import numpy as np
import pandas as pd
import OTSO

from pymagglobal import Model, coefficients
from transform_to_g_and_h import transform_to_g_and_h


def build_uniform_axis(extent=6.6, resolution=0.01):
    """Uniform axis covering [-extent, +extent] with fixed spacing, landing
    exactly on the endpoints."""
    n_half = round(extent / resolution)
    axis = np.linspace(-n_half * resolution, n_half * resolution, 2 * n_half + 1)
    return np.round(axis, 8)


def main():
    # 1. Obtain LSMOD.2 g/h coefficients, same as custom_gaussian_planet.py
    epoch = 40950  # BP
    lsmod = Model('LSMOD.2')
    _, _, coeffs = coefficients(1950 - epoch, lsmod)
    _, _, (glist, hlist) = transform_to_g_and_h(coeffs, lmax_out=13)

    axis = build_uniform_axis(extent=6.6, resolution=0.01)
    print(f"axis points per dimension: {len(axis)}  -> total grid points: {len(axis)**3}")
    print(f"spacing: {axis[1] - axis[0]:.5f} Re")

    X, Y, Z = np.meshgrid(axis, axis, axis, indexing="ij")
    X = X.ravel(); Y = Y.ravel(); Z = Z.ravel()
    n = len(X)
    locations = np.stack([X, Y, Z], axis=-1).tolist()

    chunk = 100000
    chunks_out = []

    t0 = time.perf_counter()
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        res = OTSO.magfield(
            Locations=locations[start:end],
            magfield_params={"internalmag": "Custom Gauss", "externalmag": "NONE"},
            custom_field_params={
                "g": glist,
                "h": hlist,
                "max_degree": 10,
                "MHDfile": None,
                "MHDcoordsys": None
            },
            datetime_params={"year": 2020, "month": 1, "day": 1, "hour": 1, "minute": 0, "second": 0},
            solar_wind_params={"vx": -500, "by": 1, "bz": 1, "density": 1},
            coordinate_params={"coordout": "GEO", "inputcoord": "GEO"},
            computation_params={"corenum": 7, "threadnum": 1, "Verbose": True},
        )
        # OTSO.magfield() sorts its output rows by position and, with
        # corenum > 1, concatenates worker-process results in a
        # nondeterministic arrival order - the row order (and thus a plain
        # .sort_index()) does NOT match the order `locations` was submitted
        # in. So read X/Y/Z back from the result's own position columns
        # instead of assuming alignment with `locations`, keeping each field
        # vector paired with the location it was actually computed at.
        df = res[0]
        chunks_out.append(pd.DataFrame({
            "X": df["X_GEO [Re]"].to_numpy(),
            "Y": df["Y_GEO [Re]"].to_numpy(),
            "Z": df["Z_GEO [Re]"].to_numpy(),
            "Bx": df["GEO_Bx [nT]"].to_numpy(),
            "By": df["GEO_By [nT]"].to_numpy(),
            "Bz": df["GEO_Bz [nT]"].to_numpy(),
        }))
        elapsed = time.perf_counter() - t0
        print(f"  {end}/{n} done ({elapsed:.1f}s elapsed)")

    out = pd.concat(chunks_out, ignore_index=True)

    # The Gauss expansion isn't valid inside Earth's surface - zero any node
    # with r < 1 Re, same convention as the IGRF grids.
    r2 = out["X"] ** 2 + out["Y"] ** 2 + out["Z"] ** 2
    inside = r2 < 1.0
    out.loc[inside, ["Bx", "By", "Bz"]] = 0.0
    print(f"zeroed {int(inside.sum())} points inside r < 1 Re")

    outpath = "LSMOD2_grid_uniform_01.csv"
    out.to_csv(outpath, index=False)
    t1 = time.perf_counter()
    print(f"Wrote {outpath} with {n} rows in {t1 - t0:.1f}s total")


if __name__ == "__main__":
    main()
