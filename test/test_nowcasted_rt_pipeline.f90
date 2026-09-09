program test_nowcasted_rt_pipeline
    use iso_fortran_env, only: real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use nowcasted_rt_pipeline_m, only: nowcasted_rt_pipeline_result, estimate_nowcasted_rt_pipeline
    implicit none

    integer, parameter :: n = 20
    type(delay_distribution) :: delay
    type(nowcasted_rt_pipeline_result) :: result
    real(real64) :: reported(n), weights(1), completeness
    integer :: ages(n), idx

    allocate(delay%probability(6))
    delay%probability = [0.25_real64, 0.25_real64, 0.20_real64, 0.15_real64, 0.10_real64, 0.05_real64]
    weights = [1.0_real64]

    do idx = 1, n
        ages(idx) = n - idx
        completeness = sum(delay%probability(1:min(ages(idx) + 1, size(delay%probability))))
        reported(idx) = 100.0_real64 * completeness
    end do

    call estimate_nowcasted_rt_pipeline(reported, ages, delay, weights, &
        1.0e-6_real64, 1.0e-6_real64, 3, 1.0e-6_real64, 1.0e-6_real64, &
        0.95_real64, result)

    if (.not. result%all_nowcasts_valid) error stop "all nowcasts should be valid"
    if (.not. result%raw_rt_valid(n)) error stop "raw Rt should be valid at final day"
    if (.not. result%corrected_rt_valid(n)) error stop "corrected Rt should be valid at final day"

    if (result%raw_rt(n)%mean >= 0.80_real64) then
        error stop "right-truncated reports should create a false recent Rt decline"
    end if

    if (abs(result%corrected_rt(n)%mean - 1.0_real64) > 0.05_real64) then
        error stop "nowcast-corrected Rt should recover approximately constant transmission"
    end if

    if (abs(result%nowcast_incidence(n) - 100.0_real64) > 1.0e-3_real64) then
        error stop "weak-prior nowcast should recover the latent recent incidence"
    end if

    if (result%corrected_rt(n)%mean <= result%raw_rt(n)%mean) then
        error stop "nowcasting should increase the biased recent Rt estimate"
    end if

    print *, "Nowcasted Rt pipeline tests passed"
end program test_nowcasted_rt_pipeline
