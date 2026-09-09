module nowcasted_rt_pipeline_m
    use iso_fortran_env, only: real64
    use reporting_delay_nowcast_m, only: delay_distribution, nowcast_summary, nowcast_series
    use renewal_rt_m, only: rt_posterior_summary, estimate_rt_series
    implicit none
    private

    public :: nowcasted_rt_pipeline_result, estimate_nowcasted_rt_pipeline

    type :: nowcasted_rt_pipeline_result
        real(real64), allocatable :: nowcast_incidence(:)
        type(nowcast_summary), allocatable :: nowcasts(:)
        type(rt_posterior_summary), allocatable :: raw_rt(:)
        type(rt_posterior_summary), allocatable :: corrected_rt(:)
        logical, allocatable :: raw_rt_valid(:)
        logical, allocatable :: corrected_rt_valid(:)
        logical :: all_nowcasts_valid = .false.
    end type nowcasted_rt_pipeline_result

contains

    subroutine estimate_nowcasted_rt_pipeline(reported_counts, ages, delay, generation_weights, &
            nowcast_prior_shape, nowcast_prior_rate, rt_window_size, rt_prior_shape, rt_prior_rate, &
            credible_level, result)
        real(real64), intent(in) :: reported_counts(:)
        integer, intent(in) :: ages(:)
        type(delay_distribution), intent(in) :: delay
        real(real64), intent(in) :: generation_weights(:)
        real(real64), intent(in) :: nowcast_prior_shape, nowcast_prior_rate
        integer, intent(in) :: rt_window_size
        real(real64), intent(in) :: rt_prior_shape, rt_prior_rate, credible_level
        type(nowcasted_rt_pipeline_result), intent(out) :: result
        integer :: idx, n

        n = size(reported_counts)
        if (n < 1) error stop "reported counts must not be empty"
        if (size(ages) /= n) error stop "reported counts and ages must have the same length"
        if (any(reported_counts < 0.0_real64)) error stop "reported counts must be non-negative"
        if (any(ages < 0)) error stop "reporting ages must be non-negative"
        if (size(generation_weights) < 1) error stop "generation weights must not be empty"
        if (rt_window_size < 1) error stop "Rt window size must be positive"

        allocate(result%nowcast_incidence(n), result%nowcasts(n))
        allocate(result%raw_rt(n), result%corrected_rt(n))
        allocate(result%raw_rt_valid(n), result%corrected_rt_valid(n))

        call nowcast_series(reported_counts, ages, delay, nowcast_prior_shape, nowcast_prior_rate, &
            credible_level, result%nowcasts)

        result%all_nowcasts_valid = all([(result%nowcasts(idx)%valid, idx = 1, n)])
        result%nowcast_incidence = 0.0_real64
        do idx = 1, n
            if (result%nowcasts(idx)%valid) then
                result%nowcast_incidence(idx) = result%nowcasts(idx)%mean
            end if
        end do

        call estimate_rt_series(reported_counts, generation_weights, rt_window_size, rt_prior_shape, &
            rt_prior_rate, credible_level, result%raw_rt, result%raw_rt_valid)

        result%corrected_rt_valid = .false.
        if (result%all_nowcasts_valid) then
            call estimate_rt_series(result%nowcast_incidence, generation_weights, rt_window_size, &
                rt_prior_shape, rt_prior_rate, credible_level, result%corrected_rt, &
                result%corrected_rt_valid)
        end if
    end subroutine estimate_nowcasted_rt_pipeline

end module nowcasted_rt_pipeline_m
