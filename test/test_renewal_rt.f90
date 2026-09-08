program test_renewal_rt
    use iso_fortran_env, only: real64
    use renewal_rt_m, only: rt_posterior_summary, infectiousness_from_incidence, &
        estimate_rt_window, estimate_rt_series
    implicit none

    real(real64), parameter :: tol = 1.0e-10_real64
    real(real64) :: incidence(8), weights(2), infectiousness(8)
    type(rt_posterior_summary) :: summary
    type(rt_posterior_summary) :: summaries(8)
    logical :: valid(8)

    incidence = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, &
        50.0_real64, 60.0_real64, 70.0_real64, 80.0_real64]
    weights = [1.0_real64, 0.0_real64]

    call infectiousness_from_incidence(incidence, weights, infectiousness)
    call assert_close(infectiousness(1), 0.0_real64, tol, "first infectiousness")
    call assert_close(infectiousness(2), 10.0_real64, tol, "lagged infectiousness")
    call assert_close(infectiousness(8), 70.0_real64, tol, "last infectiousness")

    call estimate_rt_window(incidence, infectiousness, 5, 8, 1.0_real64, 1.0_real64, &
        0.95_real64, summary)
    call assert_close(summary%shape, 261.0_real64, tol, "posterior shape")
    call assert_close(summary%rate, 221.0_real64, tol, "posterior rate")
    call assert_close(summary%mean, 261.0_real64 / 221.0_real64, tol, "posterior mean")
    if (summary%lower >= summary%mean .or. summary%upper <= summary%mean) then
        error stop "credible interval should contain posterior mean"
    end if
    if (summary%lower <= 0.0_real64) error stop "Rt lower bound should be positive"

    call estimate_rt_series(incidence, weights, 3, 1.0_real64, 1.0_real64, 0.95_real64, &
        summaries, valid)
    if (valid(1) .or. valid(2)) error stop "early windows should be unavailable"
    if (.not. all(valid(3:8))) error stop "later Rt windows should be available"
    if (summaries(8)%mean <= 1.0_real64) error stop "growing incidence should imply Rt above one"

    call test_constant_incidence()
    print *, "Renewal Rt tests passed"

contains

    subroutine test_constant_incidence()
        real(real64) :: counts(20), serial(20), gi(1)
        type(rt_posterior_summary) :: post

        counts = 100.0_real64
        gi = [1.0_real64]
        call infectiousness_from_incidence(counts, gi, serial)
        call estimate_rt_window(counts, serial, 10, 20, 1.0_real64, 1.0_real64, &
            0.95_real64, post)
        if (abs(post%mean - 1.0_real64) > 0.01_real64) then
            error stop "constant incidence should imply Rt near one under weak prior"
        end if
    end subroutine test_constant_incidence

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, "Assertion failed: ", trim(message)
            print *, "  actual:   ", actual
            print *, "  expected: ", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_renewal_rt
