program mixed_censored_generation_interval_uncertainty_example
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
    type(interval_censored_rt_summary) :: rt(40)
    type(interval_censored_bootstrap_summary) :: boot

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
        0.95_real64, 86420_int64, 1000, rt, boot)

    write(*,'(a,i0)') 'bootstrap_requested=', boot%n_requested
    write(*,'(a,i0)') 'bootstrap_converged=', boot%n_converged
    write(*,'(a,f12.6)') 'generation_mean=', boot%generation_mean
    write(*,'(a,f12.6)') 'generation_mean_sd=', boot%generation_mean_sd
    write(*,'(a,f12.8)') 'mean_tail_probability=', boot%mean_tail_probability
    write(*,'(a)') 'day,reported,rt_mean,rt_lower,rt_upper'
    do i = 1, 40
        if (rt(i)%valid) then
            write(*,'(i0,a,f12.6,a,f12.6,a,f12.6,a,f12.6)') &
                i, ',', reported(i), ',', rt(i)%mean, ',', rt(i)%lower, ',', rt(i)%upper
        end if
    end do
end program mixed_censored_generation_interval_uncertainty_example
