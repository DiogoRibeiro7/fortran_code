program interval_censored_generation_interval_example
    use iso_fortran_env, only: real64
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        fit_interval_censored_gamma, discretize_interval_censored_gamma
    implicit none

    real(real64), parameter :: centers(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    real(real64) :: lower(20), upper(20), weights(30), tail
    type(interval_censored_gamma_fit) :: fit
    integer :: lag

    lower = max(0.0_real64, centers - 0.4_real64)
    upper = centers + 0.4_real64

    call fit_interval_censored_gamma(lower, upper, fit)
    if (.not. fit%converged) error stop "interval-censored Gamma fit did not converge"

    call discretize_interval_censored_gamma(fit, size(weights), weights, tail)

    write(*,'(a,f12.6)') "shape = ", fit%shape
    write(*,'(a,f12.6)') "rate = ", fit%rate
    write(*,'(a,f12.6)') "mean_days = ", fit%mean
    write(*,'(a,f12.6)') "sd_days = ", fit%sd
    write(*,'(a,f12.6)') "mean_se_days = ", fit%mean_se
    write(*,'(a,es14.6)') "tail_probability_gt_30 = ", tail
    write(*,'(a)') "lag,weight"
    do lag = 1, size(weights)
        write(*,'(i0,a,es16.8)') lag, ",", weights(lag)
    end do

end program interval_censored_generation_interval_example
