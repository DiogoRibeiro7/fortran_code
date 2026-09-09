program fitted_generation_interval_uncertainty_example
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
    type(fitted_rt_monte_carlo_summary) :: rt(16)
    type(generation_interval_bootstrap_summary) :: boot
    integer :: day

    allocate(delay%probability(5))
    delay%probability = [0.35_real64, 0.30_real64, 0.18_real64, 0.10_real64, 0.07_real64]

    call propagate_fitted_generation_interval_uncertainty_rt(pair_intervals, reported, ages, delay, &
        20, 0.001_real64, 0.001_real64, 5, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260909_int64, 1000, rt, boot)

    print '(A,I0)', '# converged bootstrap fits = ', boot%n_converged
    print '(A,F10.4)', '# bootstrap mean generation interval = ', boot%generation_mean
    print '(A,F10.4)', '# bootstrap SD of mean generation interval = ', boot%generation_mean_sd
    print '(A,F10.6)', '# mean truncated tail probability = ', boot%mean_tail_probability
    print '(A)', 'day,reported,rt_mean,rt_lower,rt_upper'

    do day = 1, size(reported)
        if (rt(day)%valid) then
            write(*,'(I0,A,F10.3,A,F10.4,A,F10.4,A,F10.4)') day, ',', reported(day), ',', &
                rt(day)%mean, ',', rt(day)%lower, ',', rt(day)%upper
        end if
    end do
end program fitted_generation_interval_uncertainty_example
