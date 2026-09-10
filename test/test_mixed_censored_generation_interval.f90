program test_mixed_censored_generation_interval
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: regularized_gamma_q
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        fit_mixed_censored_gamma, discretize_interval_censored_gamma
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    real(real64) :: lower(20), upper(20), weights(40), tail, q_tail
    logical :: right_censored(20)
    type(interval_censored_gamma_fit) :: fit

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64
    right_censored = .false.
    right_censored(17:20) = .true.
    upper(17:20) = lower(17:20)

    call fit_mixed_censored_gamma(lower, upper, right_censored, fit)
    if (.not. fit%converged) error stop "mixed censored fit did not converge"
    call assert_close(fit%shape, 2.26534655_real64, 2.0e-6_real64, "mixed shape")
    call assert_close(fit%rate, 0.34202380_real64, 2.0e-6_real64, "mixed rate")
    call assert_close(fit%mean, 6.62335942_real64, 3.0e-6_real64, "mixed mean")
    if (fit%shape_se <= 0.0_real64 .or. fit%rate_se <= 0.0_real64 .or. fit%mean_se <= 0.0_real64) then
        error stop "mixed observed-information standard errors must be positive"
    end if

    q_tail = regularized_gamma_q(2.0_real64, 20.0_real64)
    call assert_close(q_tail, exp(-20.0_real64) * 21.0_real64, 1.0e-14_real64, &
        "stable Gamma upper tail")

    call discretize_interval_censored_gamma(fit, size(weights), weights, tail)
    call assert_close(sum(weights), 1.0_real64, 1.0e-12_real64, "renewal normalization")
    if (tail <= 0.0_real64) error stop "renewal tail should remain positive"

    print *, "Mixed-censored generation-interval tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_mixed_censored_generation_interval
