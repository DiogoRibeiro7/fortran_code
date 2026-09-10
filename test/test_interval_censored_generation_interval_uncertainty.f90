program test_interval_censored_generation_interval_uncertainty
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use interval_censored_generation_interval_uncertainty_m, only: &
        interval_censored_bootstrap_summary, interval_censored_rt_summary, &
        propagate_interval_censored_generation_interval_uncertainty_rt
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64,1.8_real64,2.2_real64,2.7_real64,3.1_real64,3.4_real64,3.8_real64,4.2_real64, &
        4.6_real64,5.0_real64,5.4_real64,5.9_real64,6.4_real64,7.0_real64,7.8_real64,8.6_real64, &
        9.5_real64,10.7_real64,12.0_real64,13.6_real64 ]
    real(real64) :: lower(20), upper(20), reported(40), lower_big(200), upper_big(200)
    integer :: ages(40), i, j, idx
    type(delay_distribution) :: delay
    type(interval_censored_rt_summary) :: rt_a(40), rt_b(40), rt_big(40)
    type(interval_censored_bootstrap_summary) :: boot_a, boot_b, boot_big

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64
    do i = 1, 10
        do j = 1, 20
            idx = (i-1)*20 + j
            lower_big(idx) = lower(j)
            upper_big(idx) = upper(j)
        end do
    end do

    delay%probability = [0.50_real64,0.25_real64,0.12_real64,0.07_real64,0.04_real64,0.02_real64]
    do i = 1, 40
        reported(i) = 20.0_real64 * exp(0.025_real64 * real(i-1,real64))
        ages(i) = max(0, 40-i)
    end do

    call propagate_interval_censored_generation_interval_uncertainty_rt(lower, upper, reported, ages, &
        delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
        24680_int64, 250, rt_a, boot_a)
    call propagate_interval_censored_generation_interval_uncertainty_rt(lower, upper, reported, ages, &
        delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
        24680_int64, 250, rt_b, boot_b)

    if (.not. boot_a%valid .or. boot_a%n_converged < 200) error stop "insufficient bootstrap convergence"
    call assert_close(boot_a%generation_mean, boot_b%generation_mean, 1.0e-12_real64, "bootstrap reproducibility")
    call assert_close(rt_a(40)%mean, rt_b(40)%mean, 1.0e-12_real64, "Rt reproducibility")
    if (.not. rt_a(40)%valid) error stop "final Rt summary must be valid"
    if (boot_a%generation_mean_sd <= 0.0_real64) error stop "bootstrap GI uncertainty must be positive"
    if (abs(boot_a%generation_mean - 5.95_real64) > 0.5_real64) error stop "bootstrap mean GI is implausible"

    call propagate_interval_censored_generation_interval_uncertainty_rt(lower_big, upper_big, reported, ages, &
        delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
        13579_int64, 250, rt_big, boot_big)
    if (.not. boot_big%valid) error stop "large-sample bootstrap must be valid"
    if (boot_big%generation_mean_sd >= boot_a%generation_mean_sd) then
        error stop "replicated pair dataset should reduce GI bootstrap uncertainty"
    end if

    print *, "Interval-censored generation-interval uncertainty tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        if (abs(actual-expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_interval_censored_generation_interval_uncertainty
