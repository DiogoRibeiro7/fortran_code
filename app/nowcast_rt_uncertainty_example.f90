program nowcast_rt_uncertainty_example
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use nowcast_rt_uncertainty_m, only: rt_monte_carlo_summary, propagate_nowcast_uncertainty_rt
    implicit none

    integer, parameter :: n = 20
    real(real64) :: reported(n), generation_weights(1)
    integer :: ages(n), idx
    type(delay_distribution) :: delay
    type(rt_monte_carlo_summary) :: summaries(n)

    reported = 100.0_real64
    ages = 10
    ages(n - 2) = 2
    ages(n - 1) = 1
    ages(n) = 0
    reported(n - 2) = 80.0_real64
    reported(n - 1) = 50.0_real64
    reported(n) = 20.0_real64
    generation_weights = [1.0_real64]

    allocate(delay%probability(4))
    delay%probability = [0.20_real64, 0.30_real64, 0.30_real64, 0.20_real64]

    call propagate_nowcast_uncertainty_rt(reported, ages, delay, generation_weights, &
        0.01_real64, 0.0001_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260909_int64, 2000, summaries)

    print '(A)', 'day,reported,rt_mean,rt_lower,rt_upper'
    do idx = 1, n
        if (summaries(idx)%valid) then
            print '(I0,A,F8.2,A,F8.4,A,F8.4,A,F8.4)', idx, ',', reported(idx), ',', &
                summaries(idx)%mean, ',', summaries(idx)%lower, ',', summaries(idx)%upper
        end if
    end do
end program nowcast_rt_uncertainty_example
