program test_fitted_generation_interval_uncertainty
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use fitted_generation_interval_uncertainty_m, only: generation_interval_bootstrap_summary, &
        fitted_rt_monte_carlo_summary, propagate_fitted_generation_interval_uncertainty_rt
    implicit none

    real(real64), parameter :: pair_intervals(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    real(real64), parameter :: reported(16) = [ &
        18.0_real64, 20.0_real64, 23.0_real64, 26.0_real64, &
        30.0_real64, 35.0_real64, 41.0_real64, 48.0_real64, &
        56.0_real64, 65.0_real64, 75.0_real64, 86.0_real64, &
        98.0_real64, 111.0_real64, 125.0_real64, 140.0_real64 ]
    integer, parameter :: ages(16) = [15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0]
    type(delay_distribution) :: delay
    type(fitted_rt_monte_carlo_summary) :: first(16), second(16)
    type(generation_interval_bootstrap_summary) :: boot_first, boot_second
    real(real64) :: replicated_pairs(200)
    type(fitted_rt_monte_carlo_summary) :: large_sample(16)
    type(generation_interval_bootstrap_summary) :: boot_large
    integer :: block

    allocate(delay%probability(5))
    delay%probability = [0.35_real64, 0.30_real64, 0.18_real64, 0.10_real64, 0.07_real64]

    call propagate_fitted_generation_interval_uncertainty_rt(pair_intervals, reported, ages, delay, &
        20, 0.001_real64, 0.001_real64, 5, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260909_int64, 400, first, boot_first)
    call propagate_fitted_generation_interval_uncertainty_rt(pair_intervals, reported, ages, delay, &
        20, 0.001_real64, 0.001_real64, 5, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260909_int64, 400, second, boot_second)

    if (.not. boot_first%valid .or. .not. boot_second%valid) then
        error stop "bootstrap generation-interval summaries must be valid"
    end if
    if (boot_first%n_converged /= boot_second%n_converged) then
        error stop "bootstrap convergence count is not reproducible"
    end if
    call assert_close(boot_first%shape_mean, boot_second%shape_mean, 0.0_real64, &
        "bootstrap shape mean reproducibility")
    call assert_close(boot_first%generation_mean, boot_second%generation_mean, 0.0_real64, &
        "bootstrap generation mean reproducibility")
    call assert_close(first(16)%mean, second(16)%mean, 0.0_real64, &
        "Rt Monte Carlo reproducibility")

    if (boot_first%n_converged < 350) error stop "too many bootstrap Gamma fits failed"
    if (abs(boot_first%generation_mean - 5.955_real64) > 0.35_real64) then
        error stop "bootstrap mean generation interval drifted too far from fitted reference"
    end if
    if (boot_first%generation_mean_sd <= 0.0_real64) then
        error stop "bootstrap generation-interval uncertainty must be positive"
    end if
    if (.not. first(16)%valid) error stop "final Rt summary must be valid"
    if (first(16)%lower >= first(16)%upper) error stop "final Rt interval must have positive width"

    do block = 1, 10
        replicated_pairs((block - 1) * 20 + 1:block * 20) = pair_intervals
    end do

    call propagate_fitted_generation_interval_uncertainty_rt(replicated_pairs, reported, ages, delay, &
        20, 0.001_real64, 0.001_real64, 5, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260910_int64, 400, large_sample, boot_large)

    if (.not. boot_large%valid) error stop "large-sample bootstrap summary must be valid"
    if (boot_large%generation_mean_sd >= boot_first%generation_mean_sd) then
        error stop "larger transmission-pair sample should reduce fitted generation-time uncertainty"
    end if

    print *, "Fitted generation-interval uncertainty tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_fitted_generation_interval_uncertainty
