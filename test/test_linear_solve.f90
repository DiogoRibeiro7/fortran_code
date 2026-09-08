program test_linear_solve
    use, intrinsic :: iso_fortran_env, only : real64
    use lapack_solvers_m, only : solve_linear_system, solve_success, solve_invalid_shape
    implicit none (type, external)

    real(real64), parameter :: tolerance = 1.0e-12_real64
    real(real64) :: a(3, 3)
    real(real64) :: b(3)
    real(real64) :: x(3)
    real(real64) :: bad_a(2, 3)
    real(real64) :: bad_b(2)
    real(real64) :: bad_x(2)
    integer :: info

    a = reshape([ &
        2.0_real64, -3.0_real64, -2.0_real64, &
        1.0_real64, -1.0_real64, 1.0_real64, &
        -1.0_real64, 2.0_real64, 2.0_real64 &
    ], shape(a))
    b = [8.0_real64, -11.0_real64, -3.0_real64]

    call solve_linear_system(a, b, x, info)
    if (info /= solve_success) error stop "DGESV returned a nonzero status"
    if (maxval(abs(x - [2.0_real64, 3.0_real64, -1.0_real64])) > tolerance) then
        error stop "linear-system solution mismatch"
    end if

    bad_a = 0.0_real64
    bad_b = 0.0_real64
    call solve_linear_system(bad_a, bad_b, bad_x, info)
    if (info /= solve_invalid_shape) error stop "shape validation failed"
end program test_linear_solve
