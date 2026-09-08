program renewal_rt_example
    use iso_fortran_env, only: real64
    use renewal_rt_m, only: rt_posterior_summary, estimate_rt_series
    implicit none

    real(real64), parameter :: incidence(21) = [ &
        5.0_real64, 6.0_real64, 8.0_real64, 10.0_real64, 13.0_real64, &
        17.0_real64, 22.0_real64, 28.0_real64, 35.0_real64, 43.0_real64, &
        52.0_real64, 61.0_real64, 69.0_real64, 75.0_real64, 78.0_real64, &
        76.0_real64, 71.0_real64, 64.0_real64, 56.0_real64, 48.0_real64, 40.0_real64]
    real(real64), parameter :: weights(7) = [ &
        0.05_real64, 0.12_real64, 0.20_real64, 0.24_real64, &
        0.19_real64, 0.12_real64, 0.08_real64]
    type(rt_posterior_summary) :: summaries(size(incidence))
    logical :: valid(size(incidence))
    integer :: day

    call estimate_rt_series(incidence, weights, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, summaries, valid)

    print '(a)', 'day,incidence,rt_mean,rt_lower,rt_upper'
    do day = 1, size(incidence)
        if (valid(day)) then
            print '(i0,a,f0.1,a,f0.4,a,f0.4,a,f0.4)', day, ',', incidence(day), ',', &
                summaries(day)%mean, ',', summaries(day)%lower, ',', summaries(day)%upper
        end if
    end do
end program renewal_rt_example
