module Interpolation
    implicit none
    real(8), allocatable :: MHDposition(:,:,:,:), MHDB(:,:,:,:)
    integer(4) :: n_x, n_y, n_z, regions
    integer(4), allocatable :: start_idx_x_region(:), end_idx_x_region(:)
    integer(4), allocatable :: start_idx_y_region(:), end_idx_y_region(:)
    integer(4), allocatable :: start_idx_z_region(:), end_idx_z_region(:)
    character(len=3) :: CoordINMHD, CoordOUTMHD
contains
subroutine InitializeMHD(trimmed_filename)
    implicit none
    real(8) :: x_target, y_target, z_target
    real(8) :: Bx_target, By_target, Bz_target
    integer :: i, j, k, n_points, n_unique,loop,region
    integer, allocatable :: x_values(:), y_values(:), z_values(:)
    character(len=:), allocatable :: filename
    character(len=:), allocatable :: trimmed_filename
    integer :: n_x_split, n_y_split, n_z_split,region_x,region_y,region_z
    integer :: start_idx_x, end_idx_x, start_idx_y, end_idx_y, start_idx_z, end_idx_z

    ! Define number of splits along each axis
    n_x_split = 10  ! Example: split x-axis into 4 regions
    n_y_split = 10  ! Example: split y-axis into 4 regions
    n_z_split = 10  ! Example: split z-axis into 4 regions

    regions = n_x_split * n_y_split * n_z_split

    allocate(start_idx_x_region(regions), end_idx_x_region(regions))
    allocate(start_idx_y_region(regions), end_idx_y_region(regions))
    allocate(start_idx_z_region(regions), end_idx_z_region(regions))

    ! Input file name
    filename = trimmed_filename

    ! Dynamically determine the number of points (lines in the file)
    call count_lines(filename, n_points)

    !print *, n_points

    allocate(x_values(n_points), y_values(n_points), z_values(n_points))

    call count_unique_x_values(filename, n_points, n_unique)

    !print *, n_unique

    ! Example grid size (adjust based on your data)
    n_x = n_unique  ! example size, adjust based on your data
    n_y = n_unique  ! example size, adjust based on your data
    n_z = n_unique  ! example size, adjust based on your data

    ! Allocate 4D arrays for positions (x, y, z for each point) and B fields (Bx, By, Bz for each point)
    allocate(MHDposition(n_x, n_y, n_z, 3), MHDB(n_x, n_y, n_z, 3))

    ! Open and read the CSV file
    open(unit=10, file=filename, status='old', action='read')
    read(10,*)  ! Skip header line

    ! Read in the data into the position and B arrays
    do i = 1, n_x
        do j = 1, n_y
            do k = 1, n_z
                ! Read position (x, y, z) and B field data (Bx, By, Bz) for each point
                read(10,'(3F20.6, 3F20.6, 3F20.6)') MHDposition(i,j,k,1), MHDposition(i,j,k,2), MHDposition(i,j,k,3), &
                                                 MHDB(i,j,k,1), MHDB(i,j,k,2), MHDB(i,j,k,3)
            end do
        end do
    end do
    close(10)

    do region = 1, regions
    ! Compute region indices along x, y, z
        region_x = mod(region - 1, n_x_split) + 1
        region_y = mod((region - 1) / n_x_split, n_y_split) + 1
        region_z = (region - 1) / (n_x_split * n_y_split) + 1
    
        ! Start and end indices for x-axis
        start_idx_x = (region_x - 1) * (n_x / n_x_split) + 1
        if (region_x < n_x_split) then
            end_idx_x = region_x * (n_x / n_x_split)
        else
            end_idx_x = n_x  ! Ensure the last region ends at the grid boundary
        end if
    
        ! Start and end indices for y-axis
        start_idx_y = (region_y - 1) * (n_y / n_y_split) + 1
        if (region_y < n_y_split) then
            end_idx_y = region_y * (n_y / n_y_split)
        else
            end_idx_y = n_y  ! Ensure the last region ends at the grid boundary
        end if
    
        ! Start and end indices for z-axis
        start_idx_z = (region_z - 1) * (n_z / n_z_split) + 1
        if (region_z < n_z_split) then
            end_idx_z = region_z * (n_z / n_z_split)
        else
            end_idx_z = n_z  ! Ensure the last region ends at the grid boundary
        end if
    
        ! Store the indices in the respective arrays
        start_idx_x_region(region) = start_idx_x
        end_idx_x_region(region) = end_idx_x
        start_idx_y_region(region) = start_idx_y
        end_idx_y_region(region) = end_idx_y
        start_idx_z_region(region) = start_idx_z
        end_idx_z_region(region) = end_idx_z

    end do

end subroutine InitializeMHD


subroutine count_lines(filename, n_lines)
    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(out) :: n_lines
    integer :: unit, stat
    character(len=256) :: line
    logical :: is_header

    n_lines = 0
    is_header = .true.  ! Assume the first line is a header

    open(unit=99, file=filename, status='old', action='read', iostat=stat)
    if (stat /= 0) then
        write(*,*) "Error: Unable to open file", filename
        stop
    end if

    ! Read lines, skipping the header
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

subroutine unique_values(input_array, n_points, n_unique)
    implicit none
    integer, dimension(:), intent(in) :: input_array
    integer, intent(in) :: n_points
    integer, intent(out) :: n_unique
    integer :: i, count

    ! Count unique values
    count = 1
    n_unique = 1
    do i = 2, n_points
        if (input_array(i) /= input_array(i-1)) then
            n_unique = n_unique + 1
        end if
    end do
end subroutine unique_values

subroutine sort(array)
    implicit none
    integer, dimension(:), intent(inout) :: array
    integer :: i, j, temp
    integer :: n

    n = size(array)
    do i = 1, n-1
        do j = i+1, n
            if (array(i) > array(j)) then
                ! Swap elements
                temp = array(i)
                array(i) = array(j)
                array(j) = temp
            end if
        end do
    end do
end subroutine sort

subroutine count_unique_x_values(filename, npoints, n_unique)
    implicit none
    character(len=*) :: filename
    integer :: npoints
    character(len=256) :: line
    real(8) :: x, y, z, bx, by, bz
    integer :: iunit, n_unique
    real(8), dimension(:), allocatable :: x_values
    integer :: i, unique_count, j
    logical :: is_unique

    ! Open the file for reading
    iunit = 10  ! Assign a valid unit number
    open(unit=iunit, file=filename, status='old', action='read')

    ! Skip the first line (header row)
    read(iunit, '(A)')

    ! Allocate x_values with an initial size
    allocate(x_values(npoints))

    unique_count = 0

    ! Read npoints lines and process each
    do i = 1, npoints
        ! Read a line into the string 'line'
        read(iunit, '(A)') line

        ! Parse the line into X, Y, Z, Bx, By, Bz
        read(line, *) x, y, z, bx, by, bz

        ! Check if X is already in the x_values array
        is_unique = .true.
        do j = 1, unique_count
            if (abs(x - x_values(j)) < 1.0e-6) then
                is_unique = .false.
                exit
            end if
        end do

        ! If X is unique, add it to the x_values array
        if (is_unique) then
            unique_count = unique_count + 1
            x_values(unique_count) = x
        end if
    end do

    ! Close the file
    close(iunit)

    ! Return the number of unique X values
    n_unique = unique_count

    ! Deallocate the array
    deallocate(x_values)

end subroutine count_unique_x_values


  subroutine Interpolate(x_target, y_target, z_target, position, B, n_x, n_y, n_z, Bx_out, By_out, Bz_out)
    implicit none
    real(8), intent(in) :: x_target, y_target, z_target  ! Target coordinates in Earth radii
    real(8), intent(in) :: position(:,:,:,:), B(:,:,:,:)  ! Position array (coordinates) in Earth radii
    integer, intent(in) :: n_x, n_y, n_z  ! Grid dimensions
    integer :: i, j, k
    real(8) :: dist, min_dist, x_round, y_round, z_round
    integer :: i0, i1, j0, j1, k0, k1
    real(8) :: diff_x, diff_y, diff_z
    real(8) :: xd, yd, zd
    real(8) :: c000, c100, c010, c110, c001, c101, c011, c111
    real(8) :: c00, c01, c10, c11, c0, c1
    real(8) :: Bx, By, Bz, Bx_out, By_out, Bz_out
    real(8), parameter :: step_size = 0.5
    integer :: region
    logical :: found_region

    ! Initialize minimum distance to a large value
    min_dist = 1.0E30

    x_round = floor(x_target * 2) / 2.0
    y_round = floor(y_target * 2) / 2.0
    z_round = floor(z_target * 2) / 2.0

    !print *, "Target Position:", x_round, y_round, z_round

    ! Initialize the flag to identify the region
    found_region = .false.

    !print *, start_idx_x_region(:)
    !print *, end_idx_x_region(:)

    ! Check which region the target position belongs to
    do region = 1, regions
        !print *, region
        !print *, "Target Position:", x_round, y_round, z_round
        !print *, "MinX:   ", MHDposition(1, 1, start_idx_x_region(region), 1)
        !print *, "MaxX:   ", MHDposition(1, 1, end_idx_x_region(region), 1)
        !print *, "MinY:   ", MHDposition(1, start_idx_y_region(region), 1, 2)
        !print *, "MaxY:   ", MHDposition(1, end_idx_y_region(region), 1, 2)
        !print *, "MinZ:   ", MHDposition(start_idx_z_region(region), 1, 1, 3)
        !print *, "MaxZ:   ", MHDposition(end_idx_z_region(region), 1, 1, 3)
        if (x_round >= MHDposition(1, 1, start_idx_x_region(region), 1) .and. &
            x_round <= MHDposition(1, 1, end_idx_x_region(region), 1) .and. &
            y_round >= MHDposition(1, start_idx_y_region(region), 1, 2) .and. &
            y_round <= MHDposition(1, end_idx_y_region(region), 1, 2) .and. &
            z_round >= MHDposition(start_idx_z_region(region), 1, 1, 3) .and. &
            z_round <= MHDposition(end_idx_z_region(region), 1, 1, 3)) then
            ! We found the region, now limit the search to this region
            found_region = .true.
            exit
        end if
    end do

    if (.not. found_region) then
        !print *, "Error: Target position is out of bounds!"
        GOTO 100
    end if

    !print *, "Target Position is in Region ", region

    ! Now that we know the region, set the bounds for the search within this region
    i0 = start_idx_x_region(region)
    j0 = start_idx_y_region(region)
    k0 = start_idx_z_region(region)

    ! Use the search bounds for the region
    i1 = min(i0 + 1, end_idx_x_region(region))
    j1 = min(j0 + 1, end_idx_y_region(region))
    k1 = min(k0 + 1, end_idx_z_region(region))
    
    do i = start_idx_z_region(region), end_idx_z_region(region)
        do j = start_idx_y_region(region), end_idx_y_region(region)
            do k = start_idx_x_region(region), end_idx_x_region(region)
                ! Calculate the difference between the target and the current position
                diff_x = position(i,j,k,1) - x_round
                diff_y = position(i,j,k,2) - y_round
                diff_z = position(i,j,k,3) - z_round
    
                ! Calculate the Euclidean distance to the target position
                dist = sqrt(diff_x**2 + diff_y**2 + diff_z**2)
    
                ! If this is the closest point so far, update min_dist and indices
                if (dist < min_dist) then
                    min_dist = dist
                    i0 = i
                    j0 = j
                    k0 = k
                end if
            end do
        end do
    end do

    ! Now we know the closest point (i0, j0, k0), let's print the 8 surrounding positions
    i1 = min(i0 + 1, n_x)
    j1 = min(j0 + 1, n_y)
    k1 = min(k0 + 1, n_z)

    !print *, "(", position(i0,j0,k0,1), ",", position(i0,j0,k0,2), ",", position(i0,j0,k0,3), ")"
    !print *, "(", position(i1,j0,k0,1), ",", position(i1,j0,k0,2), ",", position(i1,j0,k0,3), ")"
    !print *, "(", position(i0,j1,k0,1), ",", position(i0,j1,k0,2), ",", position(i0,j1,k0,3), ")"
    !print *, "(", position(i1,j1,k0,1), ",", position(i1,j1,k0,2), ",", position(i1,j1,k0,3), ")"
    !print *, "(", position(i0,j0,k1,1), ",", position(i0,j0,k1,2), ",", position(i0,j0,k1,3), ")"
    !print *, "(", position(i1,j0,k1,1), ",", position(i1,j0,k1,2), ",", position(i1,j0,k1,3), ")"
    !print *, "(", position(i0,j1,k1,1), ",", position(i0,j1,k1,2), ",", position(i0,j1,k1,3), ")"
    !print *, "(", position(i1,j1,k1,1), ",", position(i1,j1,k1,2), ",", position(i1,j1,k1,3), ")"

    !print *, "(", B(i0,j0,k0,1), ",", B(i0,j0,k0,2), ",", B(i0,j0,k0,3), ")"
    !print *, "(", B(i1,j0,k0,1), ",", B(i1,j0,k0,2), ",", B(i1,j0,k0,3), ")"
    !print *, "(", B(i0,j1,k0,1), ",", B(i0,j1,k0,2), ",", B(i0,j1,k0,3), ")"
    !print *, "(", B(i1,j1,k0,1), ",", B(i1,j1,k0,2), ",", B(i1,j1,k0,3), ")"
    !print *, "(", B(i0,j0,k1,1), ",", B(i0,j0,k1,2), ",", B(i0,j0,k1,3), ")"
    !print *, "(", B(i1,j0,k1,1), ",", B(i1,j0,k1,2), ",", B(i1,j0,k1,3), ")"
    !print *, "(", B(i0,j1,k1,1), ",", B(i0,j1,k1,2), ",", B(i0,j1,k1,3), ")"
    !print *, "(", B(i1,j1,k1,1), ",", B(i1,j1,k1,2), ",", B(i1,j1,k1,3), ")"

    ! Calculate the fractional distances xd, yd, zd along each axis
    xd = (x_target - position(i0,j0,k0,1)) / (position(i0,j0,k1,1) - position(i0,j0,k0,1))
    yd = (y_target - position(i0,j0,k0,2)) / (position(i0,j1,k0,2) - position(i0,j0,k0,2))
    zd = (z_target - position(i0,j0,k0,3)) / (position(i1,j0,k0,3) - position(i0,j0,k0,3))

    !print *, xd, yd, zd

    ! Interpolate Bx
    c000 = B(i0,j0,k0,1)
    c100 = B(i0,j0,k1,1)
    c010 = B(i0,j1,k0,1)
    c110 = B(i0,j1,k1,1)
    c001 = B(i1,j0,k0,1)
    c101 = B(i1,j0,k1,1)
    c011 = B(i1,j1,k0,1)
    c111 = B(i1,j1,k1,1)

    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    Bx_out = c0 * (1 - zd) + c1 * zd

    ! Interpolate By
    c000 = B(i0,j0,k0,2)
    c100 = B(i0,j0,k1,2)
    c010 = B(i0,j1,k0,2)
    c110 = B(i0,j1,k1,2)
    c001 = B(i1,j0,k0,2)
    c101 = B(i1,j0,k1,2)
    c011 = B(i1,j1,k0,2)
    c111 = B(i1,j1,k1,2)

    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    By_out = c0 * (1 - zd) + c1 * zd

    ! Interpolate Bz
    c000 = B(i0,j0,k0,3)
    c100 = B(i0,j0,k1,3)
    c010 = B(i0,j1,k0,3)
    c110 = B(i0,j1,k1,3)
    c001 = B(i1,j0,k0,3)
    c101 = B(i1,j0,k1,3)
    c011 = B(i1,j1,k0,3)
    c111 = B(i1,j1,k1,3)

    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd

    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd

    Bz_out = c0 * (1 - zd) + c1 * zd

    100 if (.not. found_region) then
        Bx_out = 0
        By_out = 0
        Bz_out = 0
    end if


end subroutine Interpolate

end module Interpolation