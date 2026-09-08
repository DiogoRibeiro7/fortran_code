program reporting_delay_nowcast_example
    use iso_fortran_env, only: real64
    use reporting_delay_nowcast_m, only: delay_distribution, nowcast_summary, nowcast_series
    implicit none

    real(real64), parameter :: reported(7) = [120.0_real64, 118.0_real64, 115.0_real64, &
        100.0_real64, 82.0_real64, 55.0_real64, 28.0_real64]
    integer, parameter :: ages(7) = [6, 5, 4, 3, 2, 1, 0]
    type(delay_distribution) :: delay
    type(nowcast_summary) :: summaries(7)
    integer :: idx

    allocate(delay%probability(7))
    delay%probability = [0.25_real64, 0.22_real64, 0.18_real64, 0.14_real64, &
        0.10_real64, 0.07_real64, 0.04_real64]

    call nowcast_series(reported, ages, delay, 0.001_real64, 0.001_real64, &
        0.95_real64, summaries)

    print '(a)', 'age,reported,completeness,nowcast_mean,lower,upper,expected_future_reports'
    do idx = 1, size(reported)
        print '(i0,a,f8.2,a,f7.4,a,f10.2,a,f10.2,a,f10.2,a,f10.2)', &
            ages(idx), ',', reported(idx), ',', summaries(idx)%completeness, ',', &
            summaries(idx)%mean, ',', summaries(idx)%lower, ',', summaries(idx)%upper, ',', &
            summaries(idx)%expected_future_reports
    end do
end program reporting_delay_nowcast_example
