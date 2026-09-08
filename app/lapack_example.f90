program lapack_example
    use, intrinsic :: iso_fortran_env, only : real64
    use lapack_solvers_m, only : solve_linear_system, solve_success
    implicit none (type, external)

    real(real64) :: a(3, 3)
    real(real64) :: b(3)
    real(real64) :: x(3)
    integer :: info

    a = reshape([ &
        2.0_real64, -3.0_real64, -2.0_real64, &
        1.0_real64, -1.0_real64, 1.0_real64, &
        -1.0_real64, 2.0_real64, 2.0_real64 &
    ], shape(a))
    b = [8.0_real64, -11.0_real64, -3.0_real64]

    call solve_linear_system(a, b, x, info)
    if (info /= solve_success) error stop "LAPACK failed to solve the system."

    print '(a, 3(f8.3, 1x))', "Solution: ", x
end program lapack_example
