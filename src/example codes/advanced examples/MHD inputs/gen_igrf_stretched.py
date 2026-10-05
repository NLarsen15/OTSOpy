import time
import numpy as np
import pandas as pd
import OTSO


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


def build_stretched_axis(extent=10.0, ds0=0.001, n_in=12, n_out=80):
    """Geometric-stretched axis with its FINEST spacing anchored at +/-1 Re
    (Earth's surface) rather than at the coordinate origin.

    An earlier version of this axis grew spacing outward from x=0 - which
    seemed reasonable ("finest near Earth"), but a Cartesian axis's x=0 is
    the grid's coordinate origin, not the location that actually matters:
    real trajectories start just above r=1 Re (Earth's surface,
    e.g. startaltitude=20km -> r~1.003 Re) and the interior r<1 Re region is
    zeroed out and never queried by any trajectory. Refining spacing at x=0
    therefore refined the wrong place - by the time the geometric growth
    reached x=1 it had already grown back to roughly the old uniform grid's
    spacing, which is why that version didn't actually reduce interpolation
    noise near Earth.

    This version instead builds two geometric segments per side, both
    starting at ds0 right at x=+-1 and growing coarser in *both* directions
    away from it:
      - n_in steps from x=1 inward to x=0 (short, mostly-unused interior -
        coarsens quickly, doesn't need many points)
      - n_out steps from x=1 outward to x=extent (the actual magnetosphere
        region - most of the resolution budget goes here)
    then mirrors the same construction onto the negative axis.
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
    # extent=6.6 Re instead of 10 - matches the magnetopause spheresize/
    # mintrapdist (6.6 Re) already used in planet_plot.py's magfield/
    # integration params, so there's no point gridding further out than
    # trajectories are ever allowed to go.
    #
    # ds0=0.001 at x=+-1 Re (Earth's surface). n_out=80 slows the outward
    # growth ratio to ~1.07x/step over the shorter 1-to-6.6 Re span, so
    # spacing stays fine basically the whole way out: edge spacing is now
    # ~0.41 Re (vs. ~0.7 Re when this same n_out covered 1-to-10 Re).
    # n_in=12 just bridges the short, largely-unused 0-to-1 Re interior.
    # Cost: 185 points/axis, ~6.33M points total - raise/lower n_out to
    # trade far-field resolution against cost (roughly n_out**3).
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
            magfield_params={"internalmag": "IGRF", "externalmag": "NONE"},
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

    # IGRF isn't valid inside Earth's surface - zero any node with r < 1 Re,
    # same convention as IGRF_grid_2013.csv.
    r2 = out["X"] ** 2 + out["Y"] ** 2 + out["Z"] ** 2
    inside = r2 < 1.0
    out.loc[inside, ["Bx", "By", "Bz"]] = 0.0
    print(f"zeroed {int(inside.sum())} points inside r < 1 Re")

    # "_r1_" marks this as the surface-anchored design (finest spacing at
    # x=+-1 Re) - distinct from earlier origin-anchored files of the same
    # ds0 (e.g. IGRF_grid_stretched_001.csv), which had much coarser actual
    # resolution at Earth's surface despite the same nominal ds0.
    outpath = f"IGRF_grid_stretched_r1_{str(ds0).replace('.', '')}_ext66.csv"
    out.to_csv(outpath, index=False)
    t1 = time.perf_counter()
    print(f"Wrote {outpath} with {n} rows in {t1 - t0:.1f}s total")


if __name__ == "__main__":
    main()
