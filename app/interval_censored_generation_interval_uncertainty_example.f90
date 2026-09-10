program interval_censored_generation_interval_uncertainty_example
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
    real(real64) :: lower(20), upper(20), reported(40)
    integer :: ages(40), day
    type(delay_distribution) :: delay
    type(interval_censored_rt_summary) :: rt(40)
    type(interval_censored_bootstrap_summary) :: boot

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64
    delay%probability = [0.50_real64,0.25_real64,0.12_real64,0.07_real64,0.04_real64,0.02_real64]

    do day = 1, 40
        reported(day) = 20.0_real64 * exp(0.025_real64 * real(day-1,real64))
        ages(day) = max(0, 40-day)
    end do

    call propagate_interval_censored_generation_interval_uncertainty_rt(lower, upper, reported, ages, &
        delay, 30, 0.01_real64, 0.01_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
        20260910_int64, 500, rt, boot)

    write(*,'(a,i0)') 'converged_bootstrap_fits=', boot%n_converged
    write(*,'(a,f10.5)') 'mean_generation_interval=', boot%generation_mean
    write(*,'(a,f10.5)') 'sd_generation_interval_mean=', boot%generation_mean_sd
    write(*,'(a,f12.8)') 'mean_truncated_tail=', boot%mean_tail_probability
    write(*,'(a)') 'day,reported,rt_mean,rt_lower,rt_upper'
    do day = 1, 40
        if (rt(day)%valid) then
            write(*,'(i0,a,f10.4,a,f10.4,a,f10.4,a,f10.4)') day, ',', reported(day), ',', &
                rt(day)%mean, ',', rt(day)%lower, ',', rt(day)%upper
        end if
    end do
end program interval_censored_generation_interval_uncertainty_example
