program test_generation_interval_fit
    use iso_fortran_env, only: real64
    use generation_interval_fit_m, only: gamma_generation_interval_fit, &
        fit_gamma_generation_interval, discretize_gamma_generation_interval
    implicit none

    real(real64), parameter :: intervals(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    type(gamma_generation_interval_fit) :: fit
    real(real64) :: weights(30), tail
    integer :: i

    call fit_gamma_generation_interval(intervals, fit)
    if (.not. fit%converged) error stop "Gamma generation-interval fit did not converge"
    call assert_close(fit%shape, 2.987568742129545_real64, 1.0e-7_real64, "shape MLE")
    call assert_close(fit%rate, 0.5016908047236851_real64, 1.0e-7_real64, "rate MLE")
    call assert_close(fit%mean, 5.955_real64, 1.0e-12_real64, "fitted mean")
    if (fit%shape_se <= 0.0_real64 .or. fit%rate_se <= 0.0_real64 .or. fit%mean_se <= 0.0_real64) then
        error stop "asymptotic standard errors must be positive"
    end if

    call discretize_gamma_generation_interval(fit, size(weights), weights, tail)
    call assert_close(sum(weights), 1.0_real64, 1.0e-12_real64, "discrete weights normalization")
    if (any(weights < 0.0_real64)) error stop "discrete generation weights must be non-negative"
    if (tail >= 1.0e-3_real64) error stop "30-day window should capture essentially all fitted mass"
    if (sum([(real(i, real64) * weights(i), i = 1, size(weights))]) < 5.0_real64) then
        error stop "discretized generation interval mean is implausibly small"
    end if

    print *, "Generation-interval fit tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_generation_interval_fit
