program generation_interval_fit_example
    use iso_fortran_env, only: real64
    use generation_interval_fit_m, only: gamma_generation_interval_fit, &
        fit_gamma_generation_interval, discretize_gamma_generation_interval
    implicit none

    real(real64), parameter :: intervals(20) = [ &
        1.4_real64, 1.8_real64, 2.2_real64, 2.7_real64, 3.1_real64, &
        3.4_real64, 3.8_real64, 4.2_real64, 4.6_real64, 5.0_real64, &
        5.4_real64, 5.9_real64, 6.4_real64, 7.0_real64, 7.8_real64, &
        8.6_real64, 9.5_real64, 10.7_real64, 12.0_real64, 13.6_real64 ]
    type(gamma_generation_interval_fit) :: fit
    real(real64) :: weights(30), tail_probability
    integer :: lag

    call fit_gamma_generation_interval(intervals, fit)
    if (.not. fit%converged) error stop "generation-interval fit did not converge"
    call discretize_gamma_generation_interval(fit, size(weights), weights, tail_probability)

    write(*,'(A,F12.6)') '# shape = ', fit%shape
    write(*,'(A,F12.6)') '# rate = ', fit%rate
    write(*,'(A,F12.6)') '# mean_days = ', fit%mean
    write(*,'(A,F12.6)') '# sd_days = ', fit%sd
    write(*,'(A,F12.6)') '# mean_se = ', fit%mean_se
    write(*,'(A,ES12.4)') '# tail_probability_after_30_days = ', tail_probability
    write(*,'(A)') 'lag,weight'
    do lag = 1, size(weights)
        write(*,'(I0,A,F12.8)') lag, ',', weights(lag)
    end do
end program generation_interval_fit_example
