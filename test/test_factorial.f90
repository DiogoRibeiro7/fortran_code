program test_factorial
    use, intrinsic :: iso_fortran_env, only : real64
    use factorial_m, only : factorial
    implicit none (type, external)

    real(real64), parameter :: tolerance = 1.0e-12_real64

    if (abs(factorial(0) - 1.0_real64) > tolerance) error stop "0! failed"
    if (abs(factorial(1) - 1.0_real64) > tolerance) error stop "1! failed"
    if (abs(factorial(5) - 120.0_real64) > tolerance) error stop "5! failed"
end program test_factorial
