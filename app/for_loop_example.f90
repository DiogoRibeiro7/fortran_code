program for_loop_example
    implicit none (type, external)

    integer :: i

    do i = 1, 5
        print '(a, 1x, i0)', "Iteration:", i
    end do
end program for_loop_example
