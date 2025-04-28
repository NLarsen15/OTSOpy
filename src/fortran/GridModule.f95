module grid_utils
    implicit none
    real(8), allocatable :: x_unique(:), y_unique(:), z_unique(:)
    save
contains

    subroutine extract_unique_sorted_x(input, n_unique)
        implicit none
        real(8), intent(in) :: input(:)
        integer, intent(out) :: n_unique

        real(8), allocatable :: temp(:)
        integer :: i, j, count, n

        n = size(input)
        allocate(temp(n))
        temp = input

        ! Simple bubble sort
        do i = 1, n - 1
            do j = i + 1, n
                if (temp(j) < temp(i)) then
                    call swap(temp(i), temp(j))
                end if
            end do
        end do

        ! Count and collect unique values
        count = 1
        do i = 2, n
            if (abs(temp(i) - temp(count)) > 1.0d-8) then
                count = count + 1
                temp(count) = temp(i)
            end if
        end do

        allocate(x_unique(count))
        x_unique = temp(1:count)
        n_unique = count

        deallocate(temp)
    end subroutine extract_unique_sorted_x

    subroutine extract_unique_sorted_y(input, n_unique)
        implicit none
        real(8), intent(in) :: input(:)
        integer, intent(out) :: n_unique

        real(8), allocatable :: temp(:)
        integer :: i, j, count, n

        n = size(input)
        allocate(temp(n))
        temp = input

        ! Simple bubble sort
        do i = 1, n - 1
            do j = i + 1, n
                if (temp(j) < temp(i)) then
                    call swap(temp(i), temp(j))
                end if
            end do
        end do

        ! Count and collect unique values
        count = 1
        do i = 2, n
            if (abs(temp(i) - temp(count)) > 1.0d-8) then
                count = count + 1
                temp(count) = temp(i)
            end if
        end do

        allocate(y_unique(count))
        y_unique = temp(1:count)
        n_unique = count

        deallocate(temp)
    end subroutine extract_unique_sorted_y

    subroutine extract_unique_sorted_z(input, n_unique)
        implicit none
        real(8), intent(in) :: input(:)
        integer, intent(out) :: n_unique

        real(8), allocatable :: temp(:)
        integer :: i, j, count, n

        n = size(input)
        allocate(temp(n))
        temp = input

        ! Simple bubble sort
        do i = 1, n - 1
            do j = i + 1, n
                if (temp(j) < temp(i)) then
                    call swap(temp(i), temp(j))
                end if
            end do
        end do

        ! Count and collect unique values
        count = 1
        do i = 2, n
            if (abs(temp(i) - temp(count)) > 1.0d-8) then
                count = count + 1
                temp(count) = temp(i)
            end if
        end do

        allocate(z_unique(count))
        z_unique = temp(1:count)
        n_unique = count

        deallocate(temp)
    end subroutine extract_unique_sorted_z

subroutine find_index(value, sorted_array, n_vals, index)
    implicit none
    real(8), intent(in) :: value
    real(8), intent(in) :: sorted_array(:)
    integer, intent(in) :: n_vals
    integer, intent(out) :: index
    integer :: i

    index = -1
    do i = 1, n_vals
        if (abs(value - sorted_array(i)) < 1.0d-8) then
            index = i
            exit
        end if
    end do

    if (index == -1) then
        print *, 'Value not found in sorted array: ', value
        stop
    end if
end subroutine find_index

    subroutine swap(a, b)
        implicit none
        real(8), intent(inout) :: a, b
        real(8) :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap

end module grid_utils
