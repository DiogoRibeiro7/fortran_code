module renewal_rt_m
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: gamma_quantile
    implicit none
    private

    public :: rt_posterior_summary
    public :: infectiousness_from_incidence, estimate_rt_window, estimate_rt_series

    type :: rt_posterior_summary
        real(real64) :: shape = 0.0_real64
        real(real64) :: rate = 0.0_real64
        real(real64) :: mean = 0.0_real64
        real(real64) :: sd = 0.0_real64
        real(real64) :: lower = 0.0_real64
        real(real64) :: upper = 0.0_real64
    end type rt_posterior_summary

contains

    subroutine infectiousness_from_incidence(incidence, weights, infectiousness)
        real(real64), intent(in) :: incidence(:)
        real(real64), intent(in) :: weights(:)
        real(real64), intent(out) :: infectiousness(:)
        integer :: t, lag, max_lag

        call validate_inputs(incidence, weights)
        if (size(infectiousness) /= size(incidence)) then
            error stop "infectiousness must have same length as incidence"
        end if

        infectiousness = 0.0_real64
        do t = 1, size(incidence)
            max_lag = min(size(weights), t - 1)
            do lag = 1, max_lag
                infectiousness(t) = infectiousness(t) + incidence(t - lag) * weights(lag)
            end do
        end do
    end subroutine infectiousness_from_incidence

    subroutine estimate_rt_window(incidence, infectiousness, start_idx, end_idx, &
            prior_shape, prior_rate, credible_level, summary)
        real(real64), intent(in) :: incidence(:), infectiousness(:)
        integer, intent(in) :: start_idx, end_idx
        real(real64), intent(in) :: prior_shape, prior_rate, credible_level
        type(rt_posterior_summary), intent(out) :: summary
        real(real64) :: sum_incidence, sum_infectiousness, tail_probability

        if (size(incidence) /= size(infectiousness)) then
            error stop "incidence and infectiousness lengths must match"
        end if
        if (start_idx < 1 .or. end_idx > size(incidence) .or. start_idx > end_idx) then
            error stop "invalid Rt estimation window"
        end if
        if (any(incidence < 0.0_real64) .or. any(infectiousness < 0.0_real64)) then
            error stop "incidence and infectiousness must be non-negative"
        end if
        if (prior_shape <= 0.0_real64 .or. prior_rate <= 0.0_real64) then
            error stop "Gamma prior shape and rate must be positive"
        end if
        if (credible_level <= 0.0_real64 .or. credible_level >= 1.0_real64) then
            error stop "credible level must lie in (0,1)"
        end if

        sum_incidence = sum(incidence(start_idx:end_idx))
        sum_infectiousness = sum(infectiousness(start_idx:end_idx))

        summary%shape = prior_shape + sum_incidence
        summary%rate = prior_rate + sum_infectiousness
        summary%mean = summary%shape / summary%rate
        summary%sd = sqrt(summary%shape) / summary%rate

        tail_probability = 0.5_real64 * (1.0_real64 - credible_level)
        summary%lower = gamma_quantile(tail_probability, summary%shape, summary%rate)
        summary%upper = gamma_quantile(1.0_real64 - tail_probability, summary%shape, summary%rate)
    end subroutine estimate_rt_window

    subroutine estimate_rt_series(incidence, weights, window_size, prior_shape, prior_rate, &
            credible_level, summaries, valid)
        real(real64), intent(in) :: incidence(:), weights(:)
        integer, intent(in) :: window_size
        real(real64), intent(in) :: prior_shape, prior_rate, credible_level
        type(rt_posterior_summary), intent(out) :: summaries(:)
        logical, intent(out) :: valid(:)
        real(real64), allocatable :: infectiousness(:)
        integer :: t, start_idx

        call validate_inputs(incidence, weights)
        if (size(summaries) /= size(incidence) .or. size(valid) /= size(incidence)) then
            error stop "Rt outputs must match incidence length"
        end if
        if (window_size < 1) error stop "window size must be positive"

        allocate(infectiousness(size(incidence)))
        call infectiousness_from_incidence(incidence, weights, infectiousness)
        valid = .false.

        do t = window_size, size(incidence)
            start_idx = t - window_size + 1
            if (sum(infectiousness(start_idx:t)) <= 0.0_real64) cycle
            call estimate_rt_window(incidence, infectiousness, start_idx, t, prior_shape, prior_rate, &
                credible_level, summaries(t))
            valid(t) = .true.
        end do
    end subroutine estimate_rt_series

    subroutine validate_inputs(incidence, weights)
        real(real64), intent(in) :: incidence(:), weights(:)
        real(real64), parameter :: tolerance = 1.0e-12_real64

        if (size(incidence) < 1) error stop "incidence must not be empty"
        if (size(weights) < 1) error stop "generation-interval weights must not be empty"
        if (any(incidence < 0.0_real64)) error stop "incidence must be non-negative"
        if (any(weights < 0.0_real64)) error stop "generation-interval weights must be non-negative"
        if (abs(sum(weights) - 1.0_real64) > tolerance) then
            error stop "generation-interval weights must sum to one"
        end if
    end subroutine validate_inputs

end module renewal_rt_m
