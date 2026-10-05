module Interpolation
    implicit none

    !SHARED VARIABLES
    
    !MHDPosition = grid positions
    !MHDB = grid field values
    !MHDA = grid vector potential values (only populated/used when
    !interp_method == 3, "divfree" - see _build_vector_potential in
    !mhd_utils.py and InterpolateDivFree below)
    real(8), allocatable :: MHDposition(:,:,:,:), MHDB(:,:,:,:), MHDA(:,:,:,:)

    !n_xyz = length of arrays in xyz axis
    integer(4) :: n_x, n_y, n_z

    !regions = number of regions
    integer(4) :: regions

    !xyz_res = grid resolution in xyz axis
    real(8) :: x_res,y_res,z_res

    !start and end xyz region = grid region boundaries in xyz axis
    integer(4), allocatable :: start_idx_x_region(:), end_idx_x_region(:)
    integer(4), allocatable :: start_idx_y_region(:), end_idx_y_region(:)
    integer(4), allocatable :: start_idx_z_region(:), end_idx_z_region(:)

    !n_xyz_split = length of chunked regions array
    integer :: n_x_split, n_y_split, n_z_split

    !region_order = order in which to process/search regions
    integer, allocatable :: region_order(:)

    !Min/Max(XYZ) = min and max position values for the MHD grid
    real :: MinX, MaxX, MinY, MaxY, MinZ, MaxZ

    !is_uniform_grid = .true. for a fixed-spacing grid (fast index-arithmetic
    !lookup); .false. for a grid with varying spacing per axis (binary-search
    !lookup against the stored axis coordinate arrays below)
    logical :: is_uniform_grid = .true.
    real(8), allocatable :: XU_axis(:), YU_axis(:), ZU_axis(:)

    !interp_method = which reconstruction Interpolate uses between grid
    !points: 0 = trilinear (default), 1 = tricubic, 2 = monotonic tricubic,
    !3 = divergence-free (curl of an interpolated vector potential A).
    !Set once at grid load time (MHDstartupSorted), consumed every call.
    integer :: interp_method = 0

    !THREAD SPECIFIC VARIABLES
    logical :: first_region
    logical :: first_region_check = .true.
    integer :: last_region
    logical :: has_last_region = .false.
    integer :: first_region_val
    SAVE

    !$omp threadprivate(has_last_region, first_region, &
    !$omp               first_region_val, first_region_check, &
    !$omp               last_region)

contains

subroutine count_lines(filename, n_lines)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_lines
    integer :: unit, stat
    character(len=256) :: line
    logical :: is_header

    n_lines = 0
    is_header = .true.

    open(unit=99, file=filename, status='old', action='read', iostat=stat)
    if (stat /= 0) then
        write(*,*) "Error: Unable to open file", filename
        stop
    end if

    do
        read(99,'(A)', iostat=stat) line
        if (stat /= 0) exit
        if (is_header) then
            is_header = .false.
        else
            n_lines = n_lines + 1
        end if
    end do

    close(99)
end subroutine count_lines


function bracket_index(arr, n, val) result(idx)
    ! Returns idx (1 <= idx <= n-1) such that arr(idx) <= val <= arr(idx+1),
    ! for a sorted ascending array arr(1:n). Values outside the array's range
    ! are clamped to the nearest edge bracket. Used for the non-uniform
    ! ("stretched") grid lookup path, where spacing isn't constant so the
    ! O(1) index-arithmetic formula used for uniform grids doesn't apply -
    ! this is an O(log n) bisection instead.
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: arr(n)
    real(8), intent(in) :: val
    integer :: idx
    integer :: lo, hi, mid

    if (n <= 1) then
        idx = 1
        return
    end if

    if (val <= arr(1)) then
        idx = 1
        return
    end if
    if (val >= arr(n)) then
        idx = n - 1
        return
    end if

    lo = 1
    hi = n
    do while (hi - lo > 1)
        mid = (lo + hi) / 2
        if (arr(mid) <= val) then
            lo = mid
        else
            hi = mid
        end if
    end do
    idx = lo
end function bracket_index


pure function clamp_index(i, n) result(ci)
    implicit none
    integer, intent(in) :: i, n
    integer :: ci
    ci = min(max(i, 1), n)
end function clamp_index


pure subroutine hermite_tangents(x0, x1, x2, x3, p0, p1, p2, p3, monotonic, m1, m2)
    ! Shared tangent estimation used by both hermite_cubic1d (value) and
    ! hermite_cubic1d_deriv (derivative w.r.t. xt) below, so the two stay
    ! consistent by construction instead of duplicating this logic.
    !
    ! monotonic=.false. ("tricubic"): tangent at each interior knot is the
    ! average of its two adjacent secant slopes - the standard "natural"
    ! cubic reconstruction. Smoother (C^1) than trilinear, but like any
    ! unlimited cubic it can overshoot near a sharp gradient or
    ! discontinuity.
    !
    ! monotonic=.true. ("monotonic tricubic"): tangent is the Fritsch-Carlson
    ! weighted harmonic mean of the two secants, clamped to zero at a local
    ! extremum. This never overshoots the range of the input samples.
    !
    ! A duplicated outer knot (x0==x1 or x2==x3, from index clamping at a
    ! grid boundary) makes the corresponding secant singular; that case
    ! falls back to the one interior secant that's available instead of
    ! dividing by zero (a natural/not-a-knot edge condition).
    implicit none
    real(8), intent(in) :: x0, x1, x2, x3, p0, p1, p2, p3
    logical, intent(in) :: monotonic
    real(8), intent(out) :: m1, m2
    real(8) :: h0, h1, h2, d0, d1, d2, w1, w2

    h0 = x1 - x0
    h1 = x2 - x1
    h2 = x3 - x2
    d1 = (p2 - p1) / h1

    if (h0 > 0.0d0) then
        d0 = (p1 - p0) / h0
    else
        d0 = d1
    end if
    if (h2 > 0.0d0) then
        d2 = (p3 - p2) / h2
    else
        d2 = d1
    end if

    if (monotonic) then
        if (d0 * d1 <= 0.0d0) then
            m1 = 0.0d0
        else
            w1 = 2.0d0 * h1 + h0
            w2 = h1 + 2.0d0 * h0
            m1 = (w1 + w2) / (w1 / d0 + w2 / d1)
        end if
        if (d1 * d2 <= 0.0d0) then
            m2 = 0.0d0
        else
            w1 = 2.0d0 * h2 + h1
            w2 = h2 + 2.0d0 * h1
            m2 = (w1 + w2) / (w1 / d1 + w2 / d2)
        end if
    else
        m1 = 0.5d0 * (d0 + d1)
        m2 = 0.5d0 * (d1 + d2)
    end if
end subroutine hermite_tangents


pure function hermite_cubic1d(x0, x1, x2, x3, p0, p1, p2, p3, xt, monotonic) result(val)
    ! Cubic Hermite interpolation of p(x) between knots x1 and x2, using the
    ! outer knots x0/x3 only (via hermite_tangents) to estimate the tangents
    ! at x1 and x2. Knots need not be evenly spaced (the formulas reduce to
    ! the familiar Catmull-Rom / uniform Fritsch-Carlson case when they are).
    implicit none
    real(8), intent(in) :: x0, x1, x2, x3, p0, p1, p2, p3, xt
    logical, intent(in) :: monotonic
    real(8) :: val
    real(8) :: h1, m1, m2, t
    real(8) :: h00, h10, h01, h11

    h1 = x2 - x1
    call hermite_tangents(x0, x1, x2, x3, p0, p1, p2, p3, monotonic, m1, m2)

    t = (xt - x1) / h1
    h00 = (1.0d0 + 2.0d0 * t) * (1.0d0 - t)**2
    h10 = t * (1.0d0 - t)**2
    h01 = t**2 * (3.0d0 - 2.0d0 * t)
    h11 = t**2 * (t - 1.0d0)

    val = h00 * p1 + h10 * h1 * m1 + h01 * p2 + h11 * h1 * m2
end function hermite_cubic1d


pure function hermite_cubic1d_deriv(x0, x1, x2, x3, p0, p1, p2, p3, xt, monotonic) result(dval)
    ! d/d(xt) of hermite_cubic1d above, using the exact same tangents (so
    ! that a curl assembled from this and hermite_cubic1d is the analytic
    ! derivative of the identical interpolant - see InterpolateDivFree).
    implicit none
    real(8), intent(in) :: x0, x1, x2, x3, p0, p1, p2, p3, xt
    logical, intent(in) :: monotonic
    real(8) :: dval
    real(8) :: h1, m1, m2, t
    real(8) :: dh00, dh10, dh01, dh11

    h1 = x2 - x1
    call hermite_tangents(x0, x1, x2, x3, p0, p1, p2, p3, monotonic, m1, m2)

    t = (xt - x1) / h1
    dh00 = 6.0d0 * t**2 - 6.0d0 * t
    dh10 = 3.0d0 * t**2 - 4.0d0 * t + 1.0d0
    dh01 = -6.0d0 * t**2 + 6.0d0 * t
    dh11 = 3.0d0 * t**2 - 2.0d0 * t

    dval = (dh00 * p1 + dh01 * p2) / h1 + dh10 * m1 + dh11 * m2
end function hermite_cubic1d_deriv


subroutine InterpolateCubic(x_target, y_target, z_target, i0, j0, k0, n_x, n_y, n_z, monotonic, &
                             Bx_out, By_out, Bz_out)
    ! Tricubic (monotonic=.false.) or monotonic-tricubic (monotonic=.true.)
    ! interpolation, built from three nested passes of the 1D cubic Hermite
    ! reconstruction above: cubic along x for each of a 4x4 grid of (y,z)
    ! lines, then cubic along y, then cubic along z - the standard
    ! separable/tensor-product construction of a tricubic interpolant, using
    ! a 4x4x4 neighbor stencil per field component instead of trilinear's
    ! 2x2x2 corners. Knot coordinates come from XU_axis/YU_axis/ZU_axis,
    ! which are populated for both uniform and stretched grids (see
    ! MHDstartupSorted), so this same code serves both grid types without
    ! branching on is_uniform_grid. i0/j0/k0 (the lower corner of the cell
    ! containing the target, 1 <= i0 <= n_x-ish) must already be computed by
    ! the caller, same as for the trilinear path.
    implicit none
    real(8), intent(in) :: x_target, y_target, z_target
    integer, intent(in) :: i0, j0, k0, n_x, n_y, n_z
    logical, intent(in) :: monotonic
    real(8), intent(out) :: Bx_out, By_out, Bz_out
    integer :: ii(4), jj(4), kk(4)
    real(8) :: xc(4), yc(4), zc(4)
    real(8) :: vx(4,4), vy(4)
    integer :: a, b, c, comp
    real(8) :: Bvals(3)

    do a = 1, 4
        ii(a) = clamp_index(i0 - 2 + a, n_x)
        jj(a) = clamp_index(j0 - 2 + a, n_y)
        kk(a) = clamp_index(k0 - 2 + a, n_z)
    end do
    xc = XU_axis(ii)
    yc = YU_axis(jj)
    zc = ZU_axis(kk)

    do comp = 1, 3
        do b = 1, 4
            do c = 1, 4
                vx(b, c) = hermite_cubic1d(xc(1), xc(2), xc(3), xc(4), &
                    MHDB(ii(1), jj(b), kk(c), comp), MHDB(ii(2), jj(b), kk(c), comp), &
                    MHDB(ii(3), jj(b), kk(c), comp), MHDB(ii(4), jj(b), kk(c), comp), &
                    x_target, monotonic)
            end do
        end do
        do c = 1, 4
            vy(c) = hermite_cubic1d(yc(1), yc(2), yc(3), yc(4), &
                vx(1, c), vx(2, c), vx(3, c), vx(4, c), y_target, monotonic)
        end do
        Bvals(comp) = hermite_cubic1d(zc(1), zc(2), zc(3), zc(4), &
            vy(1), vy(2), vy(3), vy(4), z_target, monotonic)
    end do

    Bx_out = Bvals(1)
    By_out = Bvals(2)
    Bz_out = Bvals(3)
end subroutine InterpolateCubic


subroutine tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, comp, axis, dval)
    ! Evaluate d(A_comp)/d(axis) at the target point from a 4x4x4 MHDA
    ! stencil (comp: 1=Ax, 2=Ay, 3=Az; axis: 1=x, 2=y, 3=z), using the same
    ! three-pass (x, then y, then z) separable Hermite construction as
    ! InterpolateCubic, with exactly one pass - the requested axis - run in
    ! derivative mode (hermite_cubic1d_deriv) and the other two in ordinary
    ! value mode (hermite_cubic1d). Because each pass's tangent depends only
    ! on its own axis's knots and the 4 numbers handed to it (not on the
    ! other axes' target coordinates), differentiating exactly one pass and
    ! leaving the other two as plain value-interpolation gives the analytic
    ! partial derivative of the full tensor-product interpolant - see
    ! InterpolateDivFree.
    !
    ! Always uses the unlimited (non-monotonic) tangent estimate - monotone
    ! slope-limiting has no clear physical meaning applied to a vector
    ! potential component (unlike limiting B itself, which prevents
    ! overshoot of a physical quantity).
    implicit none
    integer, intent(in) :: ii(4), jj(4), kk(4), comp, axis
    real(8), intent(in) :: xc(4), yc(4), zc(4), x_target, y_target, z_target
    real(8), intent(out) :: dval
    real(8) :: vx(4,4), vy(4)
    integer :: b, c
    logical, parameter :: monotonic = .false.

    do b = 1, 4
        do c = 1, 4
            if (axis == 1) then
                vx(b, c) = hermite_cubic1d_deriv(xc(1), xc(2), xc(3), xc(4), &
                    MHDA(ii(1), jj(b), kk(c), comp), MHDA(ii(2), jj(b), kk(c), comp), &
                    MHDA(ii(3), jj(b), kk(c), comp), MHDA(ii(4), jj(b), kk(c), comp), &
                    x_target, monotonic)
            else
                vx(b, c) = hermite_cubic1d(xc(1), xc(2), xc(3), xc(4), &
                    MHDA(ii(1), jj(b), kk(c), comp), MHDA(ii(2), jj(b), kk(c), comp), &
                    MHDA(ii(3), jj(b), kk(c), comp), MHDA(ii(4), jj(b), kk(c), comp), &
                    x_target, monotonic)
            end if
        end do
    end do

    do c = 1, 4
        if (axis == 2) then
            vy(c) = hermite_cubic1d_deriv(yc(1), yc(2), yc(3), yc(4), &
                vx(1, c), vx(2, c), vx(3, c), vx(4, c), y_target, monotonic)
        else
            vy(c) = hermite_cubic1d(yc(1), yc(2), yc(3), yc(4), &
                vx(1, c), vx(2, c), vx(3, c), vx(4, c), y_target, monotonic)
        end if
    end do

    if (axis == 3) then
        dval = hermite_cubic1d_deriv(zc(1), zc(2), zc(3), zc(4), &
            vy(1), vy(2), vy(3), vy(4), z_target, monotonic)
    else
        dval = hermite_cubic1d(zc(1), zc(2), zc(3), zc(4), &
            vy(1), vy(2), vy(3), vy(4), z_target, monotonic)
    end if
end subroutine tensor_partial


subroutine InterpolateDivFree(x_target, y_target, z_target, i0, j0, k0, n_x, n_y, n_z, &
                               Bx_out, By_out, Bz_out)
    ! Divergence-free reconstruction (Mackay, Marchand & Kabin, 2006, JGR
    ! 111, A06208): interpolate the vector potential A - precomputed once in
    ! mhd_utils.py via an FFT solve of curl(A)=B, see
    ! _build_vector_potential/MHDA - with the same separable tensor-product
    ! Hermite construction used for B in InterpolateCubic, and return
    ! B = curl(A) evaluated analytically from that same interpolant (one
    ! pass swapped to derivative mode per partial, via tensor_partial).
    !
    ! Because curl is linear and the interpolant is linear in the stencil
    ! samples, this is exact by construction: div(curl(A)) = 0 within every
    ! cell, regardless of how smooth the reconstruction is across cell
    ! boundaries elsewhere.
    implicit none
    real(8), intent(in) :: x_target, y_target, z_target
    integer, intent(in) :: i0, j0, k0, n_x, n_y, n_z
    real(8), intent(out) :: Bx_out, By_out, Bz_out
    integer :: ii(4), jj(4), kk(4)
    real(8) :: xc(4), yc(4), zc(4)
    real(8) :: dAy_dz, dAz_dy, dAx_dz, dAz_dx, dAy_dx, dAx_dy
    integer :: a

    do a = 1, 4
        ii(a) = clamp_index(i0 - 2 + a, n_x)
        jj(a) = clamp_index(j0 - 2 + a, n_y)
        kk(a) = clamp_index(k0 - 2 + a, n_z)
    end do
    xc = XU_axis(ii)
    yc = YU_axis(jj)
    zc = ZU_axis(kk)

    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 3, 2, dAz_dy)
    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 2, 3, dAy_dz)
    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 1, 3, dAx_dz)
    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 3, 1, dAz_dx)
    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 2, 1, dAy_dx)
    call tensor_partial(ii, jj, kk, xc, yc, zc, x_target, y_target, z_target, 1, 2, dAx_dy)

    Bx_out = dAz_dy - dAy_dz
    By_out = dAx_dz - dAz_dx
    Bz_out = dAy_dx - dAx_dy
end subroutine InterpolateDivFree


subroutine Interpolate(x_target, y_target, z_target, n_x, n_y, n_z, Bx_out, By_out, Bz_out)
    implicit none
    real(8), intent(in) :: x_target, y_target, z_target  ! Target coordinates in Earth radii
    integer, intent(in) :: n_x, n_y, n_z  ! Grid dimensions
    integer :: i, j, k
    real(8) :: dist, min_dist, x_round, y_round, z_round
    integer :: i0, i1, j0, j1, k0, k1, dx, dy, dz
    real(8) :: diff_x, diff_y, diff_z
    real(8) :: xd, yd, zd
    real(8) :: c000, c100, c010, c110, c001, c101, c011, c111
    real(8) :: c00, c01, c10, c11, c0, c1
    real(8) :: Bx, By, Bz, Bx_out, By_out, Bz_out
    integer :: region
    logical :: found_region, found
    integer :: region_x,region_y,region_z
    integer :: neighbor_x, neighbor_y, neighbor_z
    integer :: chunk_x, chunk_y, chunk_z

    min_dist = 1.0E30

   !print*, x_res

    if (x_target > MaxX .or. y_target > MaxY .or. z_target > MaxZ) GOTO 100
    if (x_target < MinX .or. y_target < MinY .or. z_target < MinZ) GOTO 100

    if (is_uniform_grid) then

        ! The grid is a uniform rectilinear grid, and the region chunking is
        ! also a fixed, regular partition of it (see mhd_utils.py's chunking
        ! logic: each axis is split into n_*_split chunks of size n_*/n_*_split,
        ! with the remainder absorbed into the last chunk). So both the nearest
        ! grid point and the region it falls in can be computed directly by
        ! index arithmetic, with no per-call region search and no need for the
        ! region-cache bookkeeping (first_region/has_last_region/etc) that used
        ! to be required to make that search cheap.
        x_round = floor(x_target / x_res) * x_res
        y_round = floor(y_target / y_res) * y_res
        z_round = floor(z_target / z_res) * z_res

        i0 = nint((x_round - MHDposition(1,1,1,1)) / x_res) + 1
        j0 = nint((y_round - MHDposition(1,1,1,2)) / y_res) + 1
        k0 = nint((z_round - MHDposition(1,1,1,3)) / z_res) + 1

        i0 = min(max(i0, 1), n_x)
        j0 = min(max(j0, 1), n_y)
        k0 = min(max(k0, 1), n_z)

        chunk_x = n_x / n_x_split
        chunk_y = n_y / n_y_split
        chunk_z = n_z / n_z_split

        region_x = min((i0 - 1) / chunk_x, n_x_split - 1)
        region_y = min((j0 - 1) / chunk_y, n_y_split - 1)
        region_z = min((k0 - 1) / chunk_z, n_z_split - 1)

        region = region_z * (n_x_split * n_y_split) + region_y * n_x_split + region_x + 1

    else

        ! Non-uniform ("stretched") grid: spacing isn't constant, so the
        ! nearest grid point can't be found by a fixed-step formula. Instead,
        ! binary-search each axis's stored coordinate array directly for the
        ! bracketing index - O(log n) instead of the O(1) uniform-grid
        ! formula, but still far cheaper than the brute-force search this
        ! replaced originally. No region concept is needed here: the
        ! bisection already operates on the full axis in one step.
        i0 = bracket_index(XU_axis, n_x, x_target)
        j0 = bracket_index(YU_axis, n_y, y_target)
        k0 = bracket_index(ZU_axis, n_z, z_target)

    end if

    found_region = .true.

    i1 = min(i0 + 1, n_x)
    j1 = min(j0 + 1, n_y)
    k1 = min(k0 + 1, n_z)

    !print *, region

    !print *, "Target Position:", x_round, y_round, z_round
    !print *, "MinX:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 1)
    !print *, "MaxX:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 1)
    !print *, "MinY:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 2)
    !print *, "MaxY:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 2)
    !print *, "MinZ:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 3)
    !print *, "MaxZ:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 3)

    !print *, "(", MHDposition(i0,j0,k0,1), ",", MHDposition(i0,j0,k0,2), ",", MHDposition(i0,j0,k0,3), ")"
    !print *, "(", MHDposition(i1,j0,k0,1), ",", MHDposition(i1,j0,k0,2), ",", MHDposition(i1,j0,k0,3), ")"
    !print *, "(", MHDposition(i0,j1,k0,1), ",", MHDposition(i0,j1,k0,2), ",", MHDposition(i0,j1,k0,3), ")"
    !print *, "(", MHDposition(i1,j1,k0,1), ",", MHDposition(i1,j1,k0,2), ",", MHDposition(i1,j1,k0,3), ")"
    !print *, "(", MHDposition(i0,j0,k1,1), ",", MHDposition(i0,j0,k1,2), ",", MHDposition(i0,j0,k1,3), ")"
    !print *, "(", MHDposition(i1,j0,k1,1), ",", MHDposition(i1,j0,k1,2), ",", MHDposition(i1,j0,k1,3), ")"
    !print *, "(", MHDposition(i0,j1,k1,1), ",", MHDposition(i0,j1,k1,2), ",", MHDposition(i0,j1,k1,3), ")"
    !print *, "(", MHDposition(i1,j1,k1,1), ",", MHDposition(i1,j1,k1,2), ",", MHDposition(i1,j1,k1,3), ")"

    !print *, "(", MHDB(i0,j0,k0,1), ",", MHDB(i0,j0,k0,2), ",", MHDB(i0,j0,k0,3), ")"
    !print *, "(", MHDB(i1,j0,k0,1), ",", MHDB(i1,j0,k0,2), ",", MHDB(i1,j0,k0,3), ")"
    !print *, "(", MHDB(i0,j1,k0,1), ",", MHDB(i0,j1,k0,2), ",", MHDB(i0,j1,k0,3), ")"
    !print *, "(", MHDB(i1,j1,k0,1), ",", MHDB(i1,j1,k0,2), ",", MHDB(i1,j1,k0,3), ")"
    !print *, "(", MHDB(i0,j0,k1,1), ",", MHDB(i0,j0,k1,2), ",", MHDB(i0,j0,k1,3), ")"
    !print *, "(", MHDB(i1,j0,k1,1), ",", MHDB(i1,j0,k1,2), ",", MHDB(i1,j0,k1,3), ")"
    !print *, "(", MHDB(i0,j1,k1,1), ",", MHDB(i0,j1,k1,2), ",", MHDB(i0,j1,k1,3), ")"
    !print *, "(", MHDB(i1,j1,k1,1), ",", MHDB(i1,j1,k1,2), ",", MHDB(i1,j1,k1,3), ")"

    if (interp_method == 0) then

    if (MHDposition(i1,j0,k0,1) /= MHDposition(i0,j0,k0,1)) then
    xd = (x_target - MHDposition(i0,j0,k0,1)) / (MHDposition(i1,j0,k0,1) - MHDposition(i0,j0,k0,1))
    if (xd < 1.0E-10) then
    xd = 0.0
    end if
    else
        xd = 0.0
    end if

    if (MHDposition(i0,j1,k0,2) /= MHDposition(i0,j0,k0,2)) then
    yd = (y_target - MHDposition(i0,j0,k0,2)) / (MHDposition(i0,j1,k0,2) - MHDposition(i0,j0,k0,2))
    if (yd < 1.0E-10) then
    yd = 0.0
    end if
    else
        yd = 0.0
    end if

    if (MHDposition(i0,j0,k1,3) /= MHDposition(i0,j0,k0,3)) then
    zd = (z_target - MHDposition(i0,j0,k0,3)) / (MHDposition(i0,j0,k1,3) - MHDposition(i0,j0,k0,3))
    if (zd < 1.0E-10) then
    zd = 0.0
    end if
    else
        zd = 0.0
    end if

    c000 = MHDB(i0,j0,k0,1)
    c100 = MHDB(i1,j0,k0,1)
    c010 = MHDB(i0,j1,k0,1)
    c110 = MHDB(i1,j1,k0,1)
    c001 = MHDB(i0,j0,k1,1)
    c101 = MHDB(i1,j0,k1,1)
    c011 = MHDB(i0,j1,k1,1)
    c111 = MHDB(i1,j1,k1,1)

    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    Bx_out = c0 * (1 - zd) + c1 * zd

    c000 = MHDB(i0,j0,k0,2)
    c100 = MHDB(i1,j0,k0,2)
    c010 = MHDB(i0,j1,k0,2)
    c110 = MHDB(i1,j1,k0,2)
    c001 = MHDB(i0,j0,k1,2)
    c101 = MHDB(i1,j0,k1,2)
    c011 = MHDB(i0,j1,k1,2)
    c111 = MHDB(i1,j1,k1,2)

    c00 = c000 * (1 - xd) + c100 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    By_out = c0 * (1 - zd) + c1 * zd

    c000 = MHDB(i0,j0,k0,3)
    c100 = MHDB(i1,j0,k0,3)
    c010 = MHDB(i0,j1,k0,3)
    c110 = MHDB(i1,j1,k0,3)
    c001 = MHDB(i0,j0,k1,3)
    c101 = MHDB(i1,j0,k1,3)
    c011 = MHDB(i0,j1,k1,3)
    c111 = MHDB(i1,j1,k1,3)

    c00 = c000 * (1 - xd) + c100 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    Bz_out = c0 * (1 - zd) + c1 * zd

    else if (interp_method == 3) then

        ! interp_method 3 (divfree) - see InterpolateDivFree above.
        call InterpolateDivFree(x_target, y_target, z_target, i0, j0, k0, n_x, n_y, n_z, &
            Bx_out, By_out, Bz_out)

    else

        ! interp_method 1 (tricubic) or 2 (monotonic tricubic) - see
        ! InterpolateCubic above.
        call InterpolateCubic(x_target, y_target, z_target, i0, j0, k0, n_x, n_y, n_z, &
            (interp_method == 2), Bx_out, By_out, Bz_out)

    end if

    100 if (.not. found_region) then
        !print *, "no region"
        Bx_out = 0
        By_out = 0
        Bz_out = 0
    end if

    IF (ISNAN(Bx_out)) THEN
      Bx_out = 0.0
    END IF
    IF (ISNAN(By_out)) THEN
      By_out = 0.0
    END IF
    IF (ISNAN(Bz_out)) THEN
      Bz_out = 0.0
    END IF

    
    !print *, Bx_out, By_out, Bz_out
    !print *, first_region

end subroutine Interpolate

end module Interpolation