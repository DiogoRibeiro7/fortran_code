program test_negative_binomial_rt
    use iso_fortran_env, only: real64
    use negative_binomial_rt_m, only: nb_rt_summary, negative_binomial_variance, &
        negative_binomial_log_probability, estimate_nb_rt_window, estimate_nb_rt_series
    implicit none

    call test_variance_and_poisson_limit()
    call test_constant_infectiousness_mle()
    call test_overdispersion_widens_interval()
    call test_zero_count_boundary()
    call test_series_behavior()
    call test_invalid_support()

    print *, "Negative-binomial Rt tests passed"

contains

    subroutine test_variance_and_poisson_limit()
        real(real64) :: nb_logp, poisson_logp

        call assert_close(negative_binomial_variance(10.0_real64, 5.0_real64), &
            30.0_real64, 1.0e-12_real64, "NB2 variance")
        nb_logp = negative_binomial_log_probability(7.0_real64, 6.2_real64, 1.0e8_real64)
        poisson_logp = -6.2_real64 + 7.0_real64 * log(6.2_real64) - log_gamma(8.0_real64)
        call assert_close(nb_logp, poisson_logp, 1.0e-6_real64, "Poisson limit")
    end subroutine test_variance_and_poisson_limit

    subroutine test_constant_infectiousness_mle()
        real(real64) :: incidence(4), infectiousness(4)
        type(nb_rt_summary) :: summary

        incidence = [8.0_real64, 12.0_real64, 10.0_real64, 10.0_real64]
        infectiousness = 5.0_real64
        call estimate_nb_rt_window(incidence, infectiousness, 1, 4, 3.0_real64, &
            0.95_real64, summary)
        if (.not. summary%valid) error stop "constant-infectiousness window should be valid"
        call assert_close(summary%mle, 2.0_real64, 1.0e-10_real64, "constant exposure MLE")
        if (summary%lower >= summary%mle .or. summary%upper <= summary%mle) then
            error stop "profile interval should contain the MLE"
        end if
    end subroutine test_constant_infectiousness_mle

    subroutine test_overdispersion_widens_interval()
        real(real64) :: incidence(8), infectiousness(8)
        type(nb_rt_summary) :: dispersed, near_poisson

        incidence = [8.0_real64, 12.0_real64, 10.0_real64, 10.0_real64, &
            9.0_real64, 11.0_real64, 8.0_real64, 12.0_real64]
        infectiousness = 5.0_real64
        call estimate_nb_rt_window(incidence, infectiousness, 1, 8, 2.0_real64, &
            0.95_real64, dispersed)
        call estimate_nb_rt_window(incidence, infectiousness, 1, 8, 1000.0_real64, &
            0.95_real64, near_poisson)
        if (dispersed%upper - dispersed%lower <= near_poisson%upper - near_poisson%lower) then
            error stop "greater overdispersion should widen the profile interval"
        end if
    end subroutine test_overdispersion_widens_interval

    subroutine test_zero_count_boundary()
        real(real64) :: incidence(1), infectiousness(1)
        type(nb_rt_summary) :: summary

        incidence = 0.0_real64
        infectiousness = 10.0_real64
        call estimate_nb_rt_window(incidence, infectiousness, 1, 1, 5.0_real64, &
            0.95_real64, summary)
        if (.not. summary%valid) error stop "zero-count window with exposure should be valid"
        call assert_close(summary%mle, 0.0_real64, 0.0_real64, "boundary MLE")
        call assert_close(summary%lower, 0.0_real64, 0.0_real64, "boundary lower limit")
        if (summary%upper <= 0.0_real64) error stop "boundary profile interval needs positive upper limit"
    end subroutine test_zero_count_boundary

    subroutine test_series_behavior()
        real(real64), parameter :: weights(1) = [1.0_real64]
        real(real64) :: constant_incidence(10), growing_incidence(10)
        type(nb_rt_summary) :: constant_summary(10), growing_summary(10)
        logical :: constant_valid(10), growing_valid(10)

        constant_incidence = 20.0_real64
        growing_incidence = [10.0_real64, 12.0_real64, 15.0_real64, 19.0_real64, &
            24.0_real64, 30.0_real64, 38.0_real64, 48.0_real64, 60.0_real64, 75.0_real64]

        call estimate_nb_rt_series(constant_incidence, weights, 4, 10.0_real64, 0.95_real64, &
            constant_summary, constant_valid)
        call estimate_nb_rt_series(growing_incidence, weights, 4, 10.0_real64, 0.95_real64, &
            growing_summary, growing_valid)

        if (.not. constant_valid(10) .or. .not. growing_valid(10)) then
            error stop "late Rt windows should be valid"
        end if
        call assert_close(constant_summary(10)%mle, 1.0_real64, 1.0e-10_real64, &
            "constant-incidence Rt")
        if (growing_summary(10)%mle <= 1.0_real64) then
            error stop "growing incidence should imply Rt above one"
        end if
    end subroutine test_series_behavior

    subroutine test_invalid_support()
        real(real64) :: incidence(2), infectiousness(2)
        type(nb_rt_summary) :: summary

        incidence = [5.0_real64, 0.0_real64]
        infectiousness = [0.0_real64, 5.0_real64]
        call estimate_nb_rt_window(incidence, infectiousness, 1, 2, 5.0_real64, &
            0.95_real64, summary)
        if (summary%valid) error stop "positive incidence with zero infectiousness should be invalid"
    end subroutine test_invalid_support

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

end program test_negative_binomial_rt
