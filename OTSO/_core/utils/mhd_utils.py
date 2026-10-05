import warnings
import pandas as pd
import numpy as np

from ..libs.MiddleMan import Middleman as OTSOLib

def _is_uniformly_spaced(axis, rtol=1e-6):
      if len(axis) < 3:
            return True
      diffs = np.diff(axis)
      return bool(np.allclose(diffs, diffs[0], rtol=rtol, atol=rtol * max(abs(diffs[0]), 1e-12)))

_INTERP_METHOD_CODES = {"trilinear": 0, "tricubic": 1, "monotonic": 2, "divfree": 3}


def _cumulative_hermite_integral(values, coords, axis):
      """Cumulative integral of `values` along `axis`, from the first node, of the
      same piecewise cubic Hermite reconstruction the Fortran "tricubic"/"divfree"
      interpolators use (tangent at each node = mean of the two adjacent secants,
      one-sided at the ends). Exact for that interpolant, so it works for uniform
      and stretched axes alike. Per interval: h*(p1 + p2)/2 + h**2*(m1 - m2)/12."""
      v = np.moveaxis(values, axis, 0)
      h = np.diff(coords).reshape((-1,) + (1,) * (v.ndim - 1))
      d = np.diff(v, axis=0) / h                      # secant slopes, one per interval
      m = np.empty_like(v)
      m[0] = d[0]
      m[-1] = d[-1]
      m[1:-1] = 0.5 * (d[:-1] + d[1:])
      seg = h * (v[:-1] + v[1:]) / 2.0 + h**2 * (m[:-1] - m[1:]) / 12.0
      out = np.zeros_like(v)
      out[1:] = np.cumsum(seg, axis=0)
      return np.moveaxis(out, 0, axis)


def _build_vector_potential(XU, YU, ZU, MHDposition, MHDB):
      """Build a magnetic vector potential A with curl(A) = B on the same grid as
      MHDB, for the "divfree" interpolation (the Fortran side interpolates A with
      the tricubic Hermite scheme and returns B = curl(A) analytically, so the
      interpolated field is divergence-free by construction).

      Gauge Ax = 0, integrating from reference planes x = xr, y = yr:

        Ay(x,y,z) =  int_{xr}^{x} Bz(x',y,z) dx'
        Az(x,y,z) = -int_{xr}^{x} By(x',y,z) dx' + int_{yr}^{y} Bx(xr,y',z) dy'

      Then curl(A) = (Bx(xr,y,z) - int div_perp(B) dx', By, Bz), which equals B
      exactly when div(B) = 0. Any divergence in the input data shows up only in
      Bx, accumulated along x from the reference plane, so the reference planes
      are the grid faces with the weakest field (normally an outer face, well away
      from a strong or singular inner region).

      The input must represent a divergence-free field. A grid that contains a
      singular or zero-filled inner region (e.g. a full internal field sampled
      through r < 1 Re) is not, and the reconstruction is then checked and a
      warning issued; give only the external field in the MHD file (with the
      internal field from `internalmag`) or use "monotonic" in that case.
      """
      B = np.asarray(MHDB, dtype=np.float64)
      if not np.all(np.isfinite(B)):
            n_bad = int(np.sum(~np.isfinite(B).all(axis=-1)))
            warnings.warn(f"MHDinterpolation='divfree': {n_bad} grid node(s) have non-finite field values; "
                          "they are treated as zero field when building the vector potential.", RuntimeWarning)
            B = np.nan_to_num(B, nan=0.0, posinf=0.0, neginf=0.0)

      Bmag2 = np.sum(B**2, axis=-1)
      ir = int(np.argmin(np.sqrt(Bmag2.mean(axis=(1, 2)))))     # weakest x-plane
      jr = int(np.argmin(np.sqrt(Bmag2.mean(axis=(0, 2)))))     # weakest y-plane

      Ay = _cumulative_hermite_integral(B[..., 2], XU, 0)
      Ay -= Ay[ir:ir + 1, :, :]
      Az_x = -_cumulative_hermite_integral(B[..., 1], XU, 0)
      Az_x -= Az_x[ir:ir + 1, :, :]
      Az_y = _cumulative_hermite_integral(B[ir, :, :, 0], YU, 0)
      Az_y -= Az_y[jr:jr + 1, :]

      MHDA = np.zeros_like(B)
      MHDA[..., 1] = Ay
      MHDA[..., 2] = Az_x + Az_y[np.newaxis, :, :]

      _check_vector_potential(XU, YU, ZU, B, MHDA)
      return MHDA


def _check_vector_potential(XU, YU, ZU, B, A, tol=0.05):
      """Warn if curl(A) does not reproduce B at the grid nodes (away from the
      outermost two layers, where one-sided differences are less accurate)."""
      if min(len(XU), len(YU), len(ZU)) < 6:
            return
      g = lambda F, c, ax: np.gradient(F, c, axis=ax)
      cBx = g(A[..., 2], YU, 1) - g(A[..., 1], ZU, 2)
      cBy = -g(A[..., 2], XU, 0)
      cBz = g(A[..., 1], XU, 0)
      Bmag = np.linalg.norm(B, axis=-1)
      err = np.sqrt((cBx - B[..., 0])**2 + (cBy - B[..., 1])**2 + (cBz - B[..., 2])**2)
      inner = (slice(2, -2), slice(2, -2), slice(2, -2))
      scale = max(float(np.median(Bmag[inner])), 1e-30)
      rel = err[inner] / np.maximum(Bmag[inner], scale * 1e-3)
      bad = float(np.mean(rel > tol))
      if bad > 0.01:
            warnings.warn(
                  f"MHDinterpolation='divfree': the vector potential reproduces the grid field to within "
                  f"{100 * tol:.0f}% at only {100 * (1 - bad):.1f}% of grid nodes (median error "
                  f"{100 * float(np.median(rel)):.2g}%). The gridded field is probably not divergence-free "
                  "(e.g. a singular or zero-filled inner region). Consider giving only the external field in "
                  "the MHD file with the internal field from `internalmag`, or MHDinterpolation='monotonic'.",
                  RuntimeWarning)


def MHDinitialise(MHDfile, MHDgridtype="auto", MHDinterpolation="trilinear"):
      data = pd.read_csv(MHDfile)
      x1 = data["X"].values.astype(np.float64)
      y1 = data["Y"].values.astype(np.float64)
      z1 = data["Z"].values.astype(np.float64)
      bx = data["Bx"].values.astype(np.float64)
      by = data["By"].values.astype(np.float64)
      bz = data["Bz"].values.astype(np.float64)
  
      XU = np.unique(x1)
      YU = np.unique(y1)
      ZU = np.unique(z1)
      
      XUlen, YUlen, ZUlen = len(XU), len(YU), len(ZU)
  
      ix = np.searchsorted(XU, x1)
      iy = np.searchsorted(YU, y1)
      iz = np.searchsorted(ZU, z1)
  
      MHDposition = np.zeros((XUlen, YUlen, ZUlen, 3), dtype=np.float64)
      MHDB = np.zeros((XUlen, YUlen, ZUlen, 3), dtype=np.float64)

      #if x1**2 + y1**2 + z1**2 < 1.0:
      #  bx, by, bz = 0.0, 0.0, 0.0
      MHDB[ix, iy, iz, :] = np.stack([bx, by, bz], axis=-1)
      MHDposition[ix, iy, iz, :] = np.stack([x1, y1, z1], axis=-1)
  
      # Chunking
      min_chunk = 10
      n_x_split = max(1, XUlen // min_chunk)
      n_y_split = max(1, YUlen // min_chunk)
      n_z_split = max(1, ZUlen // min_chunk)
  
      chunk_x = XUlen // n_x_split
      chunk_y = YUlen // n_y_split
      chunk_z = ZUlen // n_z_split
  
      num_regions = n_x_split * n_y_split * n_z_split
      region_info = []
  
      for r in range(num_regions):
          rx = r % n_x_split
          ry = (r // n_x_split) % n_y_split
          rz = r // (n_x_split * n_y_split)
  
          sx = rx * chunk_x
          sy = ry * chunk_y
          sz = rz * chunk_z
  
          ex = (rx + 1) * chunk_x if rx < n_x_split - 1 else XUlen
          ey = (ry + 1) * chunk_y if ry < n_y_split - 1 else YUlen
          ez = (rz + 1) * chunk_z if rz < n_z_split - 1 else ZUlen
  
          cx = np.mean(XU[sx:ex])
          cy = np.mean(YU[sy:ey])
          cz = np.mean(ZU[sz:ez])
          dist = np.sqrt(cx**2 + cy**2 + cz**2)
  
          # ex/ey/ez are Python 0-indexed, slice-exclusive upper bounds (i.e. the
          # chunk covers indices sx..ex-1). The Fortran-side arrays are 1-indexed
          # and inclusive, so the last included index sx..ex-1 maps to Fortran
          # index ex directly - NOT ex+1, which pointed one past the end of the
          # grid's last chunk on each axis (an out-of-bounds Fortran array read).
          region_info.append((dist, r+1, sx+1, ex, sy+1, ey, sz+1, ez))
  
      region_info.sort()
      region_order = [int(r[1]) for r in region_info]
      start_x = [int(r[2]) for r in region_info]
      end_x   = [int(r[3]) for r in region_info]
      start_y = [int(r[4]) for r in region_info]
      end_y   = [int(r[5]) for r in region_info]
      start_z = [int(r[6]) for r in region_info]
      end_z   = [int(r[7]) for r in region_info]
  
      minX_global = np.min(MHDposition[:, :, :, 0])  # Min value in X axis for the whole grid
      maxX_global = np.max(MHDposition[:, :, :, 0])  # Max value in X axis for the whole grid
      minY_global = np.min(MHDposition[:, :, :, 1])  # Min value in Y axis for the whole grid
      maxY_global = np.max(MHDposition[:, :, :, 1])  # Max value in Y axis for the whole grid
      minZ_global = np.min(MHDposition[:, :, :, 2])  # Min value in Z axis for the whole grid
      maxZ_global = np.max(MHDposition[:, :, :, 2])  # Max value in Z axis for the whole grid

      del data

      # MHDgridtype selects how the Fortran interpolator looks up the
      # bracketing grid point for a query position:
      #   "uniform"   - force the fast, fixed-spacing index-arithmetic path.
      #                 Only correct if every axis really is evenly spaced -
      #                 forcing this on a non-uniform grid silently gives
      #                 wrong results, since it only ever looks at the
      #                 spacing between the first two grid points.
      #   "stretched" - force the general binary-search path, correct for
      #                 both uniform and non-uniform axis spacing (e.g. a
      #                 grid that gets coarser with distance from Earth).
      #   "auto" (default) - detect per-axis spacing from the data itself.
      if MHDgridtype == "uniform":
            uniform_grid = True
      elif MHDgridtype == "stretched":
            uniform_grid = False
      elif MHDgridtype == "auto":
            uniform_grid = (_is_uniformly_spaced(XU) and _is_uniformly_spaced(YU)
                             and _is_uniformly_spaced(ZU))
      else:
            raise ValueError(
                  f"MHDgridtype must be 'auto', 'uniform', or 'stretched', got {MHDgridtype!r}"
            )

      # MHDinterpolation selects the Fortran field reconstruction:
      #   "trilinear" (default) - classic 8-corner linear blend, C^0 only.
      #   "tricubic"  - smooth (C^1) cubic fit from a 4x4x4 stencil, far more
      #                 accurate in smooth regions but can overshoot near a
      #                 sharp gradient/discontinuity (bow shock, magnetopause,
      #                 hard inner boundary).
      #   "monotonic" - same tricubic fit with Fritsch-Carlson slope
      #                 limiting, so it never overshoots - safe everywhere.
      #   "divfree"   - interpolates the magnetic vector potential A (built
      #                 once by line integrals of B, see
      #                 _build_vector_potential) instead of B directly, and
      #                 returns B = curl(A) evaluated analytically from that
      #                 same interpolant. Divergence-free by construction.
      #                 Works on uniform and stretched grids.
      try:
            interp_method = _INTERP_METHOD_CODES[MHDinterpolation]
      except KeyError:
            raise ValueError(
                  "MHDinterpolation must be 'trilinear', 'tricubic', 'monotonic', "
                  f"or 'divfree', got {MHDinterpolation!r}"
            )

      if interp_method == 3:
            MHDA = _build_vector_potential(XU, YU, ZU, MHDposition, MHDB)
      else:
            MHDA = np.zeros_like(MHDB)

      # Call Fortran (numpy arrays are passed directly - f2py converts them to
      # Fortran-contiguous arrays itself, so there is no need to go through
      # Python lists first)
      OTSOLib.mhdstartupsorted(
          XU, YU, ZU,
          MHDposition, MHDB, MHDA, n_x_split, n_y_split, n_z_split,
          minX_global,maxX_global,minY_global,maxY_global,minZ_global,maxZ_global,
          region_order, start_x, end_x,
          start_y, end_y, start_z, end_z,
          num_regions,XUlen, YUlen, ZUlen, uniform_grid, interp_method)