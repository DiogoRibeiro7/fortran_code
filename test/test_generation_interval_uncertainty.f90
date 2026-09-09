program test_generation_interval_uncertainty
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
    type(rt_joint_uncertainty_summary) :: low_uncertainty(n_days), high_uncertainty(n_days)
    type(rt_joint_uncertainty_summary) :: repeat_run(n_days)
    real(real64) :: low_width, high_width
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
        0.95_real64, 20260909_int64, 1000, low_uncertainty)

    call propagate_joint_rt_uncertainty(incidence, ages, delay, generation_weights, &
        10000.0_real64, 1.0e-6_real64, 1.0e-6_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 20260909_int64, 1000, high_uncertainty)

    call propagate_joint_rt_uncertainty(incidence, ages, delay, generation_weights, &
        10.0_real64, 1.0e-6_real64, 1.0e-6_real64, 7, 1.0_real64, 1.0_real64, &
        0.95_real64, 20260909_int64, 1000, repeat_run)

    if (.not. low_uncertainty(n_days)%valid .or. .not. high_uncertainty(n_days)%valid) then
        error stop "final Rt uncertainty summary should be valid"
    end if

    if (abs(low_uncertainty(n_days)%mean - repeat_run(n_days)%mean) > 0.0_real64 .or. &
            abs(low_uncertainty(n_days)%lower - repeat_run(n_days)%lower) > 0.0_real64 .or. &
            abs(low_uncertainty(n_days)%upper - repeat_run(n_days)%upper) > 0.0_real64) then
        error stop "fixed-seed joint uncertainty propagation is not reproducible"
    end if

    if (low_uncertainty(n_days)%generation_sd <= &
            3.0_real64 * high_uncertainty(n_days)%generation_sd) then
        error stop "lower Dirichlet concentration should materially broaden generation-time uncertainty"
    end if

    low_width = low_uncertainty(n_days)%upper - low_uncertainty(n_days)%lower
    high_width = high_uncertainty(n_days)%upper - high_uncertainty(n_days)%lower
    if (low_width <= high_width) then
        error stop "generation-interval uncertainty should widen Rt intervals for growing incidence"
    end if

    if (low_uncertainty(n_days)%generation_mean < 3.5_real64 .or. &
            low_uncertainty(n_days)%generation_mean > 4.5_real64) then
        error stop "sampled generation-interval mean is inconsistent with the baseline PMF"
    end if

    print *, "Generation-interval uncertainty tests passed"
end program test_generation_interval_uncertainty
