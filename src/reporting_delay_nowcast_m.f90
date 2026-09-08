module reporting_delay_nowcast_m
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: gamma_quantile
    implicit none
    private

    public :: delay_distribution, nowcast_summary
    public :: validate_delay_distribution, reporting_completeness
    public :: nowcast_event_day, nowcast_series

    type :: delay_distribution
        real(real64), allocatable :: probability(:)
    end type delay_distribution

    type :: nowcast_summary
        real(real64) :: completeness = 0.0_real64
        real(real64) :: posterior_shape = 0.0_real64
        real(real64) :: posterior_rate = 0.0_real64
        real(real64) :: mean = 0.0_real64
        real(real64) :: sd = 0.0_real64
        real(real64) :: lower = 0.0_real64
        real(real64) :: upper = 0.0_real64
        real(real64) :: expected_future_reports = 0.0_real64
        logical :: valid = .false.
    end type nowcast_summary

contains

    logical function validate_delay_distribution(delay) result(valid)
        type(delay_distribution), intent(in) :: delay
        real(real64), parameter :: tolerance = 1.0e-12_real64

        valid = .false.
        if (.not. allocated(delay%probability)) return
        if (size(delay%probability) < 1) return
        if (any(delay%probability < 0.0_real64)) return
        if (abs(sum(delay%probability) - 1.0_real64) > tolerance) return
        valid = .true.
    end function validate_delay_distribution

    real(real64) function reporting_completeness(delay, age) result(completeness)
        type(delay_distribution), intent(in) :: delay
        integer, intent(in) :: age
        integer :: last_delay_index

        if (.not. validate_delay_distribution(delay)) error stop "invalid reporting-delay distribution"
        if (age < 0) error stop "reporting age must be non-negative"

        last_delay_index = min(age + 1, size(delay%probability))
        completeness = sum(delay%probability(1:last_delay_index))
    end function reporting_completeness

    subroutine nowcast_event_day(reported_count, age, delay, prior_shape, prior_rate, &
            credible_level, summary)
        real(real64), intent(in) :: reported_count
        integer, intent(in) :: age
        type(delay_distribution), intent(in) :: delay
        real(real64), intent(in) :: prior_shape, prior_rate, credible_level
        type(nowcast_summary), intent(out) :: summary
        real(real64) :: tail_probability

        summary = nowcast_summary()
        if (reported_count < 0.0_real64) error stop "reported count must be non-negative"
        if (prior_shape <= 0.0_real64 .or. prior_rate <= 0.0_real64) then
            error stop "Gamma prior shape and rate must be positive"
        end if
        if (credible_level <= 0.0_real64 .or. credible_level >= 1.0_real64) then
            error stop "credible level must lie in (0,1)"
        end if

        summary%completeness = reporting_completeness(delay, age)
        if (summary%completeness <= 0.0_real64) return

        summary%posterior_shape = prior_shape + reported_count
        summary%posterior_rate = prior_rate + summary%completeness
        summary%mean = summary%posterior_shape / summary%posterior_rate
        summary%sd = sqrt(summary%posterior_shape) / summary%posterior_rate
        summary%expected_future_reports = (1.0_real64 - summary%completeness) * summary%mean

        tail_probability = 0.5_real64 * (1.0_real64 - credible_level)
        summary%lower = gamma_quantile(tail_probability, summary%posterior_shape, summary%posterior_rate)
        summary%upper = gamma_quantile(1.0_real64 - tail_probability, &
            summary%posterior_shape, summary%posterior_rate)
        summary%valid = .true.
    end subroutine nowcast_event_day

    subroutine nowcast_series(reported_counts, ages, delay, prior_shape, prior_rate, &
            credible_level, summaries)
        real(real64), intent(in) :: reported_counts(:)
        integer, intent(in) :: ages(:)
        type(delay_distribution), intent(in) :: delay
        real(real64), intent(in) :: prior_shape, prior_rate, credible_level
        type(nowcast_summary), intent(out) :: summaries(:)
        integer :: idx

        if (size(reported_counts) /= size(ages) .or. size(summaries) /= size(reported_counts)) then
            error stop "nowcast input/output lengths must match"
        end if

        do idx = 1, size(reported_counts)
            call nowcast_event_day(reported_counts(idx), ages(idx), delay, prior_shape, prior_rate, &
                credible_level, summaries(idx))
        end do
    end subroutine nowcast_series

end module reporting_delay_nowcast_m
