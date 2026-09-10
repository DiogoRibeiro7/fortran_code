program uncertain_infector_generation_interval_example
    use iso_fortran_env, only: real64
    use uncertain_infector_generation_interval_m, only: fit_uncertain_infector_gamma
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        discretize_interval_censored_gamma
    implicit none

    integer, parameter :: candidate_start(4) = [1,3,5,7]
    integer, parameter :: candidate_end(4) = [2,4,6,8]
    real(real64), parameter :: lower(8) = [ &
        1.0_real64,4.0_real64,2.0_real64,5.0_real64,3.0_real64,7.0_real64,4.0_real64,9.0_real64 ]
    real(real64), parameter :: upper(8) = [ &
        1.8_real64,4.8_real64,2.8_real64,5.8_real64,3.8_real64,7.8_real64,4.8_real64,9.8_real64 ]
    logical, parameter :: right_censored(8) = [ &
        .false.,.false.,.false.,.false.,.false.,.true.,.false.,.true. ]
    real(real64), parameter :: prior_weight(8) = [ &
        0.80_real64,0.20_real64,0.70_real64,0.30_real64,0.60_real64,0.40_real64,0.75_real64,0.25_real64 ]

    type(interval_censored_gamma_fit) :: fit
    real(real64) :: weights(30), tail
    integer :: lag

    call fit_uncertain_infector_gamma(candidate_start, candidate_end, lower, upper, right_censored, &
        prior_weight, fit)
    if (.not. fit%converged) error stop "uncertain-infector fit failed"
    call discretize_interval_censored_gamma(fit, size(weights), weights, tail)

    print '(A,F12.6)', 'shape=', fit%shape
    print '(A,F12.6)', 'rate=', fit%rate
    print '(A,F12.6)', 'mean_generation_interval=', fit%mean
    print '(A,F12.6)', 'sd_generation_interval=', fit%sd
    print '(A,ES14.6)', 'tail_probability_after_30_days=', tail
    print '(A)', 'lag,weight'
    do lag = 1, size(weights)
        print '(I0,A,ES14.6)', lag, ',', weights(lag)
    end do
end program uncertain_infector_generation_interval_example
