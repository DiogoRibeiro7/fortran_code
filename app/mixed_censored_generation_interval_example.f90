program mixed_censored_generation_interval_example
    use iso_fortran_env, only: real64
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        fit_mixed_censored_gamma, discretize_interval_censored_gamma
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    real(real64) :: lower(20), upper(20), weights(40), tail
    logical :: right_censored(20)
    type(interval_censored_gamma_fit) :: fit
    integer :: lag

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64
    right_censored = .false.
    right_censored(17:20) = .true.
    upper(17:20) = lower(17:20)

    call fit_mixed_censored_gamma(lower, upper, right_censored, fit)
    if (.not. fit%converged) error stop "mixed censored Gamma fit failed"
    call discretize_interval_censored_gamma(fit, size(weights), weights, tail)

    write(*,'(A,F12.6)') "shape=", fit%shape
    write(*,'(A,F12.6)') "rate=", fit%rate
    write(*,'(A,F12.6)') "mean_days=", fit%mean
    write(*,'(A,F12.6)') "sd_days=", fit%sd
    write(*,'(A,ES12.4)') "tail_probability=", tail
    write(*,'(A)') "lag,weight"
    do lag = 1, size(weights)
        write(*,'(I0,A,ES16.8)') lag, ",", weights(lag)
    end do

end program mixed_censored_generation_interval_example
