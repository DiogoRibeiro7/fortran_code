program factorial_example
    use, intrinsic :: iso_fortran_env, only : real64
    use factorial_m, only : factorial
    implicit none (type, external)

    integer, parameter :: n = 10
    real(real64) :: value

    value = factorial(n)
    print '(a, i0, a, es16.8)', "Factorial of ", n, " is ", value
end program factorial_example
