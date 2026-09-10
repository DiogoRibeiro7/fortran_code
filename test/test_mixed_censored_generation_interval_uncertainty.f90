program test_mixed_censored_generation_interval_uncertainty
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use interval_censored_generation_interval_uncertainty_m, only: &
        interval_censored_bootstrap_summary, interval_censored_rt_summary
    use mixed_censored_generation_interval_uncertainty_m, only: &
        propagate_mixed_censored_generation_interval_uncertainty_rt
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64,1.8_real64,2.2_real64,2.7_real64,3.1_real64,3.4_real64,3.8_real64,4.2_real64, &
        4.6_real64,5.0_real64,5.4_real64,5.9_real64,6.4_real64,7.0_real64,7.8_real64,8.6_real64, &
        9.5_real64,10.7_real64,12.0_real64,13.6_real64 ]
    real(real64) :: lower(20), upper(20), reported(40)
    logical :: right_censored(20)
    integer :: ages(40), i
    type(delay_distribution) :: delay
    type(interval_censored_rt_summary) :: rt_a(40), rt_b(40)
    type(interval_censored_bootstrap_summary) :: boot_a, boot_b

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64
    right_censored = .false.
    right_censored(17:20) = .true.
    upper(17:20) = lower(17:20)

    delay%probability = [0.50_real64,0.25_real64,0.12_real64,0.07_real64,0.04_real64,0.02_real64]
    do i = 1, 40
        reported(i) = 20.0_real64 * exp(0.025_real64 * real(i-1,real64))
        ages(i) = max(0, 40-i)
    end do

    call propagate_mixed_censored_generation_interval_uncertainty_rt(lower, upper, right_censored, &
        reported, ages, delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 97531_int64, 300, rt_a, boot_a)
    call propagate_mixed_censored_generation_interval_uncertainty_rt(lower, upper, right_censored, &
        reported, ages, delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 97531_int64, 300, rt_b, boot_b)

    if (.not. boot_a%valid) error stop "mixed-censoring bootstrap summary must be valid"
    if (boot_a%n_converged < 250) error stop "insufficient mixed-censoring bootstrap convergence"
    if (boot_a%generation_mean_sd <= 0.0_real64) error stop "GI bootstrap uncertainty must be positive"
    if (abs(boot_a%generation_mean - 6.6_real64) > 1.0_real64) then
        error stop "mixed-censoring bootstrap mean generation interval is implausible"
    end if
    if (.not. rt_a(40)%valid) error stop "final mixed-censoring Rt summary must be valid"
    if (rt_a(40)%lower >= rt_a(40)%upper) error stop "Rt interval must be ordered"

    call assert_close(boot_a%generation_mean, boot_b%generation_mean, 1.0e-12_real64, &
        "bootstrap reproducibility")
    call assert_close(boot_a%generation_mean_sd, boot_b%generation_mean_sd, 1.0e-12_real64, &
        "bootstrap uncertainty reproducibility")
    call assert_close(rt_a(40)%mean, rt_b(40)%mean, 1.0e-12_real64, "Rt reproducibility")

    print *, "Mixed-censored generation-interval uncertainty tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        if (abs(actual-expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_mixed_censored_generation_interval_uncertainty
