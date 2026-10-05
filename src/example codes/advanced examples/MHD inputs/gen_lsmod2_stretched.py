import time
import numpy as np
import pandas as pd
import OTSO

from pymagglobal import Model, coefficients
from transform_to_g_and_h import transform_to_g_and_h


def _solve_growth_ratio(dist, ds0, n, tol=1e-13, maxiter=400):
    """Bisect for the geometric ratio r such that n steps of spacing
    ds0, ds0*r, ds0*r**2, ... sum to exactly `dist`."""
    if n <= 0:
        return 1.0
    lo, hi = 1.0 + 1e-14, 10.0
    total = lambda r: ds0 * (r**n - 1.0) / (r - 1.0)
    r = hi
    for _ in range(maxiter):
        r = 0.5 * (lo + hi)
        if total(r) < dist:
            lo = r
        else:
            hi = r
        if hi - lo < tol:
            break
    return r


def build_stretched_axis(extent=6.6, ds0=0.001, n_in=12, n_out=80):
    """Geometric-stretched axis with its FINEST spacing anchored at +/-1 Re
    (Earth's surface) rather than at the coordinate origin - same
    construction as gen_igrf_stretched.py's build_stretched_axis (see there
    for the full rationale): two geometric segments per side, both starting
    at ds0 right at x=+-1 and growing coarser in *both* directions away from
    it (n_in steps inward to the mostly-unused 0-to-1 Re interior, n_out
    steps outward to `extent`), then mirrored onto the negative axis.
    """
    r_in = _solve_growth_ratio(1.0, ds0, n_in)
    k = np.arange(0, n_in + 1)
    cum_in = ds0 * (r_in**k - 1.0) / (r_in - 1.0)
    inner = 1.0 - cum_in[::-1]  # ascending 0 -> 1, coarse -> fine

    r_out = _solve_growth_ratio(extent - 1.0, ds0, n_out)
    k = np.arange(0, n_out + 1)
    cum_out = ds0 * (r_out**k - 1.0) / (r_out - 1.0)
    outer = 1.0 + cum_out  # ascending 1 -> extent, fine -> coarse

    positive_half = np.concatenate([inner[:-1], outer])  # 0 ... 1 ... extent
    axis = np.concatenate([-positive_half[::-1], positive_half[1:]])
    return np.round(axis, 8)


def main():
    # 1. Obtain LSMOD.2 g/h coefficients, same as gen_lsmod2_uniform.py
    epoch = 40950  # BP
    lsmod = Model('LSMOD.2')
    _, _, coeffs = coefficients(1950 - epoch, lsmod)
    _, _, (glist, hlist) = transform_to_g_and_h(coeffs, lmax_out=13)

    # Same surface-anchored geometry as gen_igrf_stretched.py: ds0=0.001 Re
    # right at Earth's surface (x=+-1 Re), n_out=80 to keep the outward
    # growth slow out to extent=6.6 Re (matching the magnetopause
    # spheresize/mintrapdist used in planet_plot_lsmod2.py), n_in=12 to
    # bridge the short, largely-unused 0-to-1 Re interior.
    # Cost: 185 points/axis, ~6.33M points total.
    ds0 = 0.001
    axis = build_stretched_axis(extent=6.6, ds0=ds0, n_in=12, n_out=80)
    diffs = np.diff(axis)
    idx1 = np.searchsorted(axis, 1.0)
    print(f"axis points per dimension: {len(axis)}  -> total grid points: {len(axis)**3}")
    print(f"spacing at Earth's surface (x=1 Re): {diffs[idx1-1]:.5f} Re, "
          f"at center (x=0): {diffs[len(diffs)//2]:.5f} Re, "
          f"at outer edge: {diffs[-1]:.5f} Re")

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

    # "_r1_" marks this as the surface-anchored design (finest spacing at
    # x=+-1 Re), matching the naming convention of the IGRF stretched grid.
    outpath = f"LSMOD2_grid_stretched_r1_{str(ds0).replace('.', '')}_ext66.csv"
    out.to_csv(outpath, index=False)
    t1 = time.perf_counter()
    print(f"Wrote {outpath} with {n} rows in {t1 - t0:.1f}s total")


if __name__ == "__main__":
    main()
