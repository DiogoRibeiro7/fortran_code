program nowcasted_rt_pipeline_example
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

    print '(a)', 'day,age,reported,nowcast_mean,raw_rt,corrected_rt'
    do idx = 1, n
        if (result%raw_rt_valid(idx) .and. result%corrected_rt_valid(idx)) then
            write(*,'(i0,a,i0,a,f10.4,a,f10.4,a,f10.4,a,f10.4)') &
                idx, ',', ages(idx), ',', reported(idx), ',', result%nowcast_incidence(idx), ',', &
                result%raw_rt(idx)%mean, ',', result%corrected_rt(idx)%mean
        end if
    end do
end program nowcasted_rt_pipeline_example
