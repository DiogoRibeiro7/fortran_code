program test_interval_censored_generation_interval
    use iso_fortran_env, only: real64
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        fit_interval_censored_gamma, discretize_interval_censored_gamma
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    real(real64) :: lower(20), upper(20), weights(30), tail
    type(interval_censored_gamma_fit) :: fit

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64

    call fit_interval_censored_gamma(lower, upper, fit)
    if (.not. fit%converged) error stop "interval-censored Gamma fit did not converge"
    call assert_close(fit%shape, 2.996167998648517_real64, 2.0e-6_real64, "shape MLE")
    call assert_close(fit%rate, 0.503203531270069_real64, 2.0e-6_real64, "rate MLE")
    call assert_close(fit%mean, 5.954187147864975_real64, 2.0e-6_real64, "fitted mean")

    if (fit%shape_se <= 0.0_real64 .or. fit%rate_se <= 0.0_real64 .or. fit%mean_se <= 0.0_real64) then
        error stop "observed-information standard errors must be positive"
    end if

    call discretize_interval_censored_gamma(fit, size(weights), weights, tail)
    call assert_close(sum(weights), 1.0_real64, 1.0e-12_real64, "renewal weight normalization")
    if (any(weights < 0.0_real64)) error stop "renewal weights must be non-negative"
    if (tail >= 1.0e-3_real64) error stop "30-day renewal grid leaves too much fitted tail mass"

    print *, "Interval-censored generation-interval tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_interval_censored_generation_interval
