program negative_binomial_rt_example
    use iso_fortran_env, only: real64
    use negative_binomial_rt_m, only: nb_rt_summary, estimate_nb_rt_series
    implicit none

    real(real64), parameter :: incidence(18) = [ &
        8.0_real64, 9.0_real64, 11.0_real64, 15.0_real64, 14.0_real64, 21.0_real64, &
        18.0_real64, 26.0_real64, 24.0_real64, 33.0_real64, 29.0_real64, 39.0_real64, &
        34.0_real64, 45.0_real64, 40.0_real64, 52.0_real64, 47.0_real64, 60.0_real64]
    real(real64), parameter :: weights(5) = [0.08_real64, 0.22_real64, 0.30_real64, &
        0.24_real64, 0.16_real64]
    integer, parameter :: window_size = 5
    type(nb_rt_summary) :: strong_overdispersion(size(incidence))
    type(nb_rt_summary) :: near_poisson(size(incidence))
    logical :: valid_overdispersed(size(incidence)), valid_near_poisson(size(incidence))
    integer :: day

    call estimate_nb_rt_series(incidence, weights, window_size, 2.0_real64, 0.95_real64, &
        strong_overdispersion, valid_overdispersed)
    call estimate_nb_rt_series(incidence, weights, window_size, 1000.0_real64, 0.95_real64, &
        near_poisson, valid_near_poisson)

    print '(a)', 'day,incidence,rt_nb,lower_nb,upper_nb,rt_near_poisson,lower_near_poisson,upper_near_poisson'
    do day = 1, size(incidence)
        if (.not. valid_overdispersed(day) .or. .not. valid_near_poisson(day)) cycle
        print '(i0,7(",",f10.4))', day, incidence(day), strong_overdispersion(day)%mle, &
            strong_overdispersion(day)%lower, strong_overdispersion(day)%upper, &
            near_poisson(day)%mle, near_poisson(day)%lower, near_poisson(day)%upper
    end do
end program negative_binomial_rt_example
