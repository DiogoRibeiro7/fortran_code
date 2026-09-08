module factorial_m
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none (type, external)
    private

    public :: factorial

contains

    function factorial(n) result(value)
        integer, intent(in) :: n
        real(real64) :: value
        integer :: i

        if (n < 0) error stop "factorial requires n >= 0"

        value = 1.0_real64
        do i = 2, n
            value = value * real(i, kind=real64)
        end do
    end function factorial

end module factorial_m
