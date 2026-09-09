program generation_interval_uncertainty_example
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use generation_interval_uncertainty_m, only: rt_joint_uncertainty_summary, &
        propagate_joint_rt_uncertainty
    implicit none

    integer, parameter :: n_days = 24
    real(real64) :: incidence(n_days)
    integer :: ages(n_days)
    real(real64) :: generation_weights(7)
    type(delay_distribution) :: delay
    type(rt_joint_uncertainty_summary) :: uncertain(n_days), near_fixed(n_days)
    integer :: idx

    do idx = 1, n_days
        incidence(idx) = 20.0_real64 * exp(0.08_real64 * real(idx - 1, real64))
        ages(idx) = 30
    end do

    allocate(delay%probability(1))
    delay%probability = [1.0_real64]
    generation_weights = [0.05_real64, 0.10_real64, 0.20_real64, 0.30_real64, &
        0.20_real64, 0.10_real64, 0.05_real64]

    call propagate_joint_rt_uncertainty(incidence, ages, delay, generation_weights, &
        10.0_real64, 1.0e-6_real64, 1.0e-6_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 20260909_int64, 2000, uncertain)

    call propagate_joint_rt_uncertainty(incidence, ages, delay, generation_weights, &
        10000.0_real64, 1.0e-6_real64, 1.0e-6_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 20260909_int64, 2000, near_fixed)

    print '(A)', 'day,incidence,rt_mean_uncertain,rt_lower_uncertain,rt_upper_uncertain,' // &
        'rt_mean_near_fixed,rt_lower_near_fixed,rt_upper_near_fixed'
    do idx = 1, n_days
        if (uncertain(idx)%valid .and. near_fixed(idx)%valid) then
            write (*,'(I0,A,F10.3,A,6(F10.4,:,A))') idx, ',', incidence(idx), ',', &
                uncertain(idx)%mean, ',', uncertain(idx)%lower, ',', uncertain(idx)%upper, ',', &
                near_fixed(idx)%mean, ',', near_fixed(idx)%lower, ',', near_fixed(idx)%upper
        end if
    end do

    print '(A,F8.4)', '# sampled generation-time SD, kappa=10: ', uncertain(n_days)%generation_sd
    print '(A,F8.4)', '# sampled generation-time SD, kappa=10000: ', near_fixed(n_days)%generation_sd
end program generation_interval_uncertainty_example
