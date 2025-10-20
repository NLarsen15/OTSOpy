module Interpolation
    implicit none
    real(8), allocatable :: MHDposition(:,:,:,:), MHDB(:,:,:,:)
    real(8) :: x_res,y_res,z_res
    integer(4) :: n_x, n_y, n_z, regions, region
    integer(4), allocatable :: start_idx_x_region(:), end_idx_x_region(:)
    integer(4), allocatable :: start_idx_y_region(:), end_idx_y_region(:)
    integer(4), allocatable :: start_idx_z_region(:), end_idx_z_region(:)
    character(len=3) :: CoordINMHD, CoordOUTMHD
    integer :: last_region
    logical :: has_last_region
    integer :: n_x_split, n_y_split, n_z_split
    real(8), allocatable :: region_distance(:)
    integer, allocatable :: region_order(:)
    logical :: first_region, first_region_check
    integer :: first_region_val
    real :: MinX, MaxX, MinY, MaxY, MinZ, MaxZ
    integer :: search_range, num_combinations, idx, temp
    integer, allocatable :: dx_list(:), dy_list(:), dz_list(:), dist2(:), order(:)
    SAVE
contains

subroutine InitializeMHD(trimmed_filename, XU, YU, ZU, XUlen, YUlen, ZUlen)
    use grid_utils
    implicit none

    character(len=*), intent(in) :: trimmed_filename
    character(len=:), allocatable :: filename
    integer(4), intent(in) :: XUlen, YUlen, ZUlen
    real(8), intent(in) :: XU(XUlen), YU(YUlen), ZU(ZUlen)
    real(8), allocatable :: x(:), y(:), z(:), bx(:), by(:), bz(:)
    real(8) :: x1, y1, z1, x2, y2, z2
    integer :: n_points
    integer :: i, idx_x, idx_y, idx_z
    integer :: region_x, region_y, region_z
    integer :: chunk_x, chunk_y, chunk_z
    integer :: start_idx_x, end_idx_x
    integer :: start_idx_y, end_idx_y
    integer :: start_idx_z, end_idx_z
    integer :: min_chunk_size
    integer :: key, j
    real(8) :: key_dist

    filename = trimmed_filename
    has_last_region = .false.

    ! Count number of points
    call count_lines(filename, n_points)
    allocate(x(n_points), y(n_points), z(n_points))
    allocate(bx(n_points), by(n_points), bz(n_points))

    ! Read data and capture first two points for resolution
    open(unit=10, file=filename, status='old', action='read')
    read(10,*)  ! Skip header

    do i = 1, n_points
        read(10,*) x(i), y(i), z(i), bx(i), by(i), bz(i)
        if (i == 1) then
            x1 = x(i); y1 = y(i); z1 = z(i)
        else if (i == 2) then
            x2 = x(i); y2 = y(i); z2 = z(i)
        end if
    end do
    close(10)

    x_res = abs(x2 - x1)
    y_res = x_res
    z_res = x_res

    ! Grid dimensions passed from Python
    n_x = XUlen
    n_y = YUlen
    n_z = ZUlen

    ! Region splitting setup
    min_chunk_size = 4
    n_x_split = max(1, n_x / min_chunk_size)
    n_y_split = max(1, n_y / min_chunk_size)
    n_z_split = max(1, n_z / min_chunk_size)
    regions = n_x_split * n_y_split * n_z_split

    ! Allocate main data arrays
    allocate(MHDposition(n_x, n_y, n_z, 3))
    allocate(MHDB(n_x, n_y, n_z, 3))

    ! Fill in the MHDposition and MHDB grids
    do i = 1, n_points
        call find_index(x(i), XU, n_x, idx_x)
        call find_index(y(i), YU, n_y, idx_y)
        call find_index(z(i), ZU, n_z, idx_z)

        MHDposition(idx_x, idx_y, idx_z, :) = [x(i), y(i), z(i)]
        MHDB(idx_x, idx_y, idx_z, :) = [bx(i), by(i), bz(i)]
    end do

    ! Allocate region index bounds
    allocate(start_idx_x_region(regions), end_idx_x_region(regions))
    allocate(start_idx_y_region(regions), end_idx_y_region(regions))
    allocate(start_idx_z_region(regions), end_idx_z_region(regions))

    chunk_x = n_x / n_x_split
    chunk_y = n_y / n_y_split
    chunk_z = n_z / n_z_split

    do region = 1, regions
        region_x = mod(region - 1, n_x_split) + 1
        region_y = mod((region - 1) / n_x_split, n_y_split) + 1
        region_z = (region - 1) / (n_x_split * n_y_split) + 1
    
        start_idx_x = (region_x - 1) * chunk_x + 1
        if (region_x < n_x_split) then
            end_idx_x = region_x * chunk_x
        else
            end_idx_x = n_x
        end if
    
        start_idx_y = (region_y - 1) * chunk_y + 1
        if (region_y < n_y_split) then
            end_idx_y = region_y * chunk_y
        else
            end_idx_y = n_y
        end if
    
        start_idx_z = (region_z - 1) * chunk_z + 1
        if (region_z < n_z_split) then
            end_idx_z = region_z * chunk_z
        else
            end_idx_z = n_z
        end if
    
        start_idx_x_region(region) = start_idx_x
        end_idx_x_region(region) = end_idx_x
        start_idx_y_region(region) = start_idx_y
        end_idx_y_region(region) = end_idx_y
        start_idx_z_region(region) = start_idx_z
        end_idx_z_region(region) = end_idx_z

        region_distance(region) = sqrt( &
            (XU((start_idx_x + end_idx_x)/2))**2 + &
            (YU((start_idx_y + end_idx_y)/2))**2 + &
            (ZU((start_idx_z + end_idx_z)/2))**2 )
        region_order(region) = region

        do i = 2, regions
        key = region_order(i)
        key_dist = region_distance(key)
        j = i - 1

        do while (j > 0 .and. region_distance(region_order(j)) > key_dist)
            region_order(j+1) = region_order(j)
            j = j - 1
        end do
        region_order(j+1) = key
        end do

        !print *, "Region ", region
        !print *, "Start indices (i, j, k): ", start_idx_x, start_idx_y, start_idx_z
        !print *, "End indices (i, j, k): ", end_idx_x, end_idx_y, end_idx_z
        !print *, "Start (x, y, z): ", MHDposition(start_idx_x, start_idx_y, start_idx_z, 1), &
        !          MHDposition(start_idx_x, start_idx_y, start_idx_z, 2), &
        !          MHDposition(start_idx_x, start_idx_y, start_idx_z, 3)
        !print *, "End (x, y, z): ", MHDposition(end_idx_x, end_idx_y, end_idx_z, 1), &
        !          MHDposition(end_idx_x, end_idx_y, end_idx_z, 2), &
        !          MHDposition(end_idx_x, end_idx_y, end_idx_z, 3)
    
        !print *, "Start (bx, by, bz): ", MHDB(start_idx_x, start_idx_y, start_idx_z, 1), &
        !          MHDB(start_idx_x, start_idx_y, start_idx_z, 2), &
        !          MHDB(start_idx_x, start_idx_y, start_idx_z, 3)
        !print *, "End (bx, by, bz): ", MHDB(end_idx_x, end_idx_y, end_idx_z, 1), &
        !          MHDB(end_idx_x, end_idx_y, end_idx_z, 2), &
        !          MHDB(end_idx_x, end_idx_y, end_idx_z, 3)
    end do

    ! Clean-up: leave x,y,z,bx,by,bz allocated for optional later use

end subroutine InitializeMHD

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
    real(8), parameter :: step_size = 0.5
    integer :: region
    logical :: found_region, found
    integer :: region_x,region_y,region_z
    integer :: neighbor_x, neighbor_y, neighbor_z

    min_dist = 1.0E30

    !print *, x_res

    x_round = floor(x_target / x_res) * x_res
    y_round = floor(y_target / y_res) * y_res
    z_round = floor(z_target / z_res) * z_res

    !print *, "Target Position:", x_target, y_target, z_target
    !print *, "Target Position:", x_round, y_round, z_round

    if (x_target > MaxX .or. y_target > MaxY .or. z_target > MaxZ) GOTO 100
    if (x_target < MinX .or. y_target < MinY .or. z_target < MinZ) GOTO 100


    found_region = .false.
    
if (first_region_check .and. first_region) then
    region = first_region_val
    !print *, "searching first region: ", first_region_val
    !print *, "Target Position:", x_round, y_round, z_round
    !print *, "MinX:   ", MHDposition(start_idx_x_region(first_region_val), start_idx_y_region(first_region_val), &
    !start_idx_z_region(first_region_val), 1)
    !print *, "MaxX:   ", MHDposition(end_idx_x_region(first_region_val), end_idx_y_region(first_region_val), &
    !end_idx_z_region(first_region_val), 1)
    !print *, "MinY:   ", MHDposition(start_idx_x_region(first_region_val), start_idx_y_region(first_region_val), &
    !start_idx_z_region(first_region_val), 2)
    !print *, "MaxY:   ", MHDposition(end_idx_x_region(first_region_val), end_idx_y_region(first_region_val), &
    !end_idx_z_region(first_region_val), 2)
    !print *, "MinZ:   ", MHDposition(start_idx_x_region(first_region_val), start_idx_y_region(first_region_val), &
    !start_idx_z_region(first_region_val), 3)
    !print *, "MaxZ:   ", MHDposition(end_idx_x_region(first_region_val), end_idx_y_region(first_region_val), &
    !end_idx_z_region(first_region_val), 3)
    if (x_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 1) .and. &
        x_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 1) .and. &
        y_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 2) .and. &
        y_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 2) .and. &
        z_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 3) .and. &
        z_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 3)) then
        
        found_region = .true.
        last_region = region
        has_last_region = .true.
        first_region_check = .false.
    end if
end if


    if (has_last_region) then
    !print *, "searching last region: ", last_region

    region_x = mod(last_region - 1, n_x_split) + 1
    region_y = mod((last_region - 1) / n_x_split, n_y_split) + 1
    region_z = (last_region - 1) / (n_x_split * n_y_split) + 1

    found = .false.
        do dx = -1, 1
        do dy = -1, 1
            do dz = -1, 1
                neighbor_x = region_x + dx
                neighbor_y = region_y + dy
                neighbor_z = region_z + dz

                if (neighbor_x < 1 .or. neighbor_x > n_x_split) cycle
                if (neighbor_y < 1 .or. neighbor_y > n_y_split) cycle
                if (neighbor_z < 1 .or. neighbor_z > n_z_split) cycle

                region = (neighbor_z - 1) * (n_x_split * n_y_split) + (neighbor_y - 1) * n_x_split + neighbor_x

                !print *, "Target Position:", x_round, y_round, z_round
                !print *, "MinX:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                !start_idx_z_region(region), 1)
                !print *, "MaxX:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                !end_idx_z_region(region), 1)
                !print *, "MinY:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                !start_idx_z_region(region), 2)
                !print *, "MaxY:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                !end_idx_z_region(region), 2)
                !print *, "MinZ:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                !start_idx_z_region(region), 3)
                !print *, "MaxZ:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                !end_idx_z_region(region), 3)
                
                
                
                if (x_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                start_idx_z_region(region), 1) .and. &
                    x_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                    end_idx_z_region(region), 1) .and. &
                    y_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                    start_idx_z_region(region), 2) .and. &
                    y_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                    end_idx_z_region(region), 2) .and. &
                    z_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
                    start_idx_z_region(region), 3) .and. &
                    z_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
                    end_idx_z_region(region), 3)) then

                    last_region = region
                    has_last_region = .true.
                    found = .true.
                    found_region = .true.
                    exit
                end if
            end do
            if (found) then
            !print *, "In neigboring region"
            !print *, region
            found_region = .true.
            exit
            end if
        end do
        if (found) then
        !print *, "In neigboring region"
        !print *, region
        found_region = .true.
        exit
        end if
    end do
    if (.not. found) then
    has_last_region = .false.
    found_region = .false.
    !print *, "Not in neigboring regions"
    end if
end if


if (.not. found_region) then
    !print *, "Searching All regions"
    do i = 1, regions
        region = region_order(i)
        !print *, region
        !print *, "Target Position:", x_round, y_round, z_round
        !print *, "MinX:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 1)
        !print *, "MaxX:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 1)
        !print *, "MinY:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 2)
        !print *, "MaxY:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 2)
        !print *, "MinZ:   ", MHDposition(start_idx_x_region(region), start_idx_y_region(region), start_idx_z_region(region), 3)
        !print *, "MaxZ:   ", MHDposition(end_idx_x_region(region), end_idx_y_region(region), end_idx_z_region(region), 3)
        if (x_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
        start_idx_z_region(region), 1) .and. &
            x_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
            end_idx_z_region(region), 1) .and. &
            y_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region),&
            start_idx_z_region(region), 2) .and. &
            y_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region),&
            end_idx_z_region(region), 2) .and. &
            z_round >= MHDposition(start_idx_x_region(region), start_idx_y_region(region), &
            start_idx_z_region(region), 3) .and. &
            z_round <= MHDposition(end_idx_x_region(region), end_idx_y_region(region), &
            end_idx_z_region(region), 3)) then
            found_region = .true.
            last_region = region
            has_last_region = .true.
            exit
        end if
    end do
end if


if (.not. found_region) then
        !print *, "Error: Target position is out of bounds!"
        !print *, x_target, y_target, z_target
        !print *, x_round, y_round, z_round
    has_last_region = .false.
    GOTO 100
end if

if (.not. first_region) then
   first_region_val = region
   first_region = .true.
   first_region_check = .true.
   !print *, "assigning first region: ", first_region_val
end if


    !print *, "Target Position is in Region ", region

    i0 = start_idx_x_region(region)
    j0 = start_idx_y_region(region)
    k0 = start_idx_z_region(region)

    i1 = min(i0 + 1, end_idx_x_region(region))
    j1 = min(j0 + 1, end_idx_y_region(region))
    k1 = min(k0 + 1, end_idx_z_region(region))
    
    do i = start_idx_x_region(region), end_idx_x_region(region)
        do j = start_idx_y_region(region), end_idx_y_region(region)
            do k = start_idx_z_region(region), end_idx_z_region(region)
                diff_x = MHDposition(i,j,k,1) - x_round
                diff_y = MHDposition(i,j,k,2) - y_round
                diff_z = MHDposition(i,j,k,3) - z_round
    
                dist = sqrt(diff_x**2 + diff_y**2 + diff_z**2)
    
                if (dist < min_dist) then
                    min_dist = dist
                    i0 = i
                    j0 = j
                    k0 = k
                end if
            end do
        end do
    end do

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

subroutine swap(a, b)
    real(8), intent(inout) :: a, b
    real(8) :: temp
    temp = a
    a = b
    b = temp
end subroutine swap

end module Interpolation