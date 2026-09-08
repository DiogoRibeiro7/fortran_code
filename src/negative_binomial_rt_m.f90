module negative_binomial_rt_m
    use iso_fortran_env, only: real64
    use renewal_rt_m, only: infectiousness_from_incidence
    implicit none
    private

    public :: nb_rt_summary
    public :: negative_binomial_variance, negative_binomial_log_probability
    public :: estimate_nb_rt_window, estimate_nb_rt_series

    type :: nb_rt_summary
        logical :: valid = .false.
        real(real64) :: dispersion = 0.0_real64
        real(real64) :: mle = 0.0_real64
        real(real64) :: lower = 0.0_real64
        real(real64) :: upper = 0.0_real64
        real(real64) :: log_likelihood = 0.0_real64
    end type nb_rt_summary

contains

    pure real(real64) function negative_binomial_variance(mean_count, dispersion) result(variance)
        real(real64), intent(in) :: mean_count, dispersion

        if (mean_count < 0.0_real64) error stop "mean count must be non-negative"
        if (dispersion <= 0.0_real64) error stop "dispersion must be positive"
        variance = mean_count + mean_count * mean_count / dispersion
    end function negative_binomial_variance

    pure real(real64) function negative_binomial_log_probability(count, mean_count, dispersion) result(logp)
        real(real64), intent(in) :: count, mean_count, dispersion

        if (count < 0.0_real64) error stop "count must be non-negative"
        if (mean_count < 0.0_real64) error stop "mean count must be non-negative"
        if (dispersion <= 0.0_real64) error stop "dispersion must be positive"
        if (abs(count - anint(count)) > 1.0e-8_real64) then
            error stop "negative-binomial observations must be integer-valued"
        end if

        if (mean_count <= 0.0_real64) then
            if (count <= 0.0_real64) then
                logp = 0.0_real64
            else
                logp = -huge(1.0_real64)
            end if
            return
        end if

        logp = log_gamma(count + dispersion) - log_gamma(dispersion) - log_gamma(count + 1.0_real64) &
            - dispersion * log(1.0_real64 + mean_count / dispersion) &
            + count * (log(mean_count) - log(dispersion + mean_count))
    end function negative_binomial_log_probability

    subroutine estimate_nb_rt_window(incidence, infectiousness, start_idx, end_idx, &
            dispersion, confidence_level, summary)
        real(real64), intent(in) :: incidence(:), infectiousness(:)
        integer, intent(in) :: start_idx, end_idx
        real(real64), intent(in) :: dispersion, confidence_level
        type(nb_rt_summary), intent(out) :: summary
        real(real64) :: mle, max_log_likelihood, cutoff

        call validate_window_inputs(incidence, infectiousness, start_idx, end_idx, &
            dispersion, confidence_level)
        summary = nb_rt_summary()
        summary%dispersion = dispersion

        if (.not. window_has_likelihood_support(incidence, infectiousness, start_idx, end_idx)) return

        mle = find_mle(incidence, infectiousness, start_idx, end_idx, dispersion)
        max_log_likelihood = window_log_likelihood(incidence, infectiousness, start_idx, end_idx, &
            mle, dispersion)
        cutoff = max_log_likelihood - likelihood_drop(confidence_level)

        summary%valid = .true.
        summary%mle = mle
        summary%log_likelihood = max_log_likelihood
        call profile_interval(incidence, infectiousness, start_idx, end_idx, dispersion, &
            mle, cutoff, summary%lower, summary%upper)
    end subroutine estimate_nb_rt_window

    subroutine estimate_nb_rt_series(incidence, weights, window_size, dispersion, confidence_level, &
            summaries, valid)
        real(real64), intent(in) :: incidence(:), weights(:)
        integer, intent(in) :: window_size
        real(real64), intent(in) :: dispersion, confidence_level
        type(nb_rt_summary), intent(out) :: summaries(:)
        logical, intent(out) :: valid(:)
        real(real64), allocatable :: infectiousness(:)
        integer :: t, start_idx

        if (size(summaries) /= size(incidence) .or. size(valid) /= size(incidence)) then
            error stop "negative-binomial Rt outputs must match incidence length"
        end if
        if (window_size < 1) error stop "window size must be positive"
        if (dispersion <= 0.0_real64) error stop "dispersion must be positive"
        if (confidence_level <= 0.0_real64 .or. confidence_level >= 1.0_real64) then
            error stop "confidence level must lie in (0,1)"
        end if
        call validate_count_vector(incidence)

        allocate(infectiousness(size(incidence)))
        call infectiousness_from_incidence(incidence, weights, infectiousness)
        summaries = nb_rt_summary()
        valid = .false.

        do t = window_size, size(incidence)
            start_idx = t - window_size + 1
            call estimate_nb_rt_window(incidence, infectiousness, start_idx, t, dispersion, &
                confidence_level, summaries(t))
            valid(t) = summaries(t)%valid
        end do
    end subroutine estimate_nb_rt_series

    pure real(real64) function window_log_likelihood(incidence, infectiousness, start_idx, end_idx, &
            reproduction_number, dispersion) result(log_likelihood)
        real(real64), intent(in) :: incidence(:), infectiousness(:)
        integer, intent(in) :: start_idx, end_idx
        real(real64), intent(in) :: reproduction_number, dispersion
        real(real64) :: mean_count
        integer :: idx

        if (reproduction_number < 0.0_real64) then
            log_likelihood = -huge(1.0_real64)
            return
        end if

        log_likelihood = 0.0_real64
        do idx = start_idx, end_idx
            mean_count = reproduction_number * infectiousness(idx)
            if (mean_count <= 0.0_real64 .and. incidence(idx) > 0.0_real64) then
                log_likelihood = -huge(1.0_real64)
                return
            end if
            log_likelihood = log_likelihood + &
                negative_binomial_log_probability(incidence(idx), mean_count, dispersion)
        end do
    end function window_log_likelihood

    pure real(real64) function score_rt(incidence, infectiousness, start_idx, end_idx, &
            reproduction_number, dispersion) result(score)
        real(real64), intent(in) :: incidence(:), infectiousness(:)
        integer, intent(in) :: start_idx, end_idx
        real(real64), intent(in) :: reproduction_number, dispersion
        integer :: idx
        real(real64) :: lambda

        if (reproduction_number <= 0.0_real64) then
            if (sum(incidence(start_idx:end_idx)) > 0.0_real64) then
                score = huge(1.0_real64)
            else
                score = -sum(infectiousness(start_idx:end_idx))
            end if
            return
        end if

        score = 0.0_real64
        do idx = start_idx, end_idx
            lambda = infectiousness(idx)
            score = score + dispersion * &
                (incidence(idx) - reproduction_number * lambda) / &
                (reproduction_number * (dispersion + reproduction_number * lambda))
        end do
    end function score_rt

    real(real64) function find_mle(incidence, infectiousness, start_idx, end_idx, dispersion) result(mle)
        real(real64), intent(in) :: incidence(:), infectiousness(:), dispersion
        integer, intent(in) :: start_idx, end_idx
        real(real64) :: lower, upper, midpoint, naive_ratio
        integer :: iteration

        if (sum(incidence(start_idx:end_idx)) <= 0.0_real64) then
            mle = 0.0_real64
            return
        end if

        naive_ratio = sum(incidence(start_idx:end_idx)) / &
            max(sum(infectiousness(start_idx:end_idx)), tiny(1.0_real64))
        lower = 0.0_real64
        upper = max(1.0_real64, 2.0_real64 * naive_ratio)
        do while (score_rt(incidence, infectiousness, start_idx, end_idx, upper, dispersion) > 0.0_real64)
            upper = 2.0_real64 * upper
        end do

        do iteration = 1, 100
            midpoint = 0.5_real64 * (lower + upper)
            if (score_rt(incidence, infectiousness, start_idx, end_idx, midpoint, dispersion) > 0.0_real64) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        mle = 0.5_real64 * (lower + upper)
    end function find_mle

    subroutine profile_interval(incidence, infectiousness, start_idx, end_idx, dispersion, &
            mle, cutoff, lower_limit, upper_limit)
        real(real64), intent(in) :: incidence(:), infectiousness(:), dispersion, mle, cutoff
        integer, intent(in) :: start_idx, end_idx
        real(real64), intent(out) :: lower_limit, upper_limit
        real(real64) :: lower, upper, midpoint
        integer :: iteration

        if (mle <= 0.0_real64) then
            lower_limit = 0.0_real64
        else
            lower = 0.0_real64
            upper = mle
            do iteration = 1, 100
                midpoint = 0.5_real64 * (lower + upper)
                if (window_log_likelihood(incidence, infectiousness, start_idx, end_idx, &
                        midpoint, dispersion) < cutoff) then
                    lower = midpoint
                else
                    upper = midpoint
                end if
            end do
            lower_limit = upper
        end if

        lower = mle
        upper = max(1.0_real64, 2.0_real64 * max(mle, 1.0e-6_real64))
        do while (window_log_likelihood(incidence, infectiousness, start_idx, end_idx, &
                upper, dispersion) > cutoff)
            upper = 2.0_real64 * upper
        end do
        do iteration = 1, 100
            midpoint = 0.5_real64 * (lower + upper)
            if (window_log_likelihood(incidence, infectiousness, start_idx, end_idx, &
                    midpoint, dispersion) > cutoff) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        upper_limit = upper
    end subroutine profile_interval

    pure real(real64) function likelihood_drop(confidence_level) result(drop)
        real(real64), intent(in) :: confidence_level
        real(real64) :: z

        z = normal_quantile(0.5_real64 * (1.0_real64 + confidence_level))
        drop = 0.5_real64 * z * z
    end function likelihood_drop

    pure real(real64) function normal_quantile(probability) result(quantile)
        real(real64), intent(in) :: probability
        real(real64) :: lower, upper, midpoint, cdf
        integer :: iteration

        if (probability <= 0.0_real64 .or. probability >= 1.0_real64) then
            error stop "normal quantile probability must lie in (0,1)"
        end if

        lower = -10.0_real64
        upper = 10.0_real64
        do iteration = 1, 100
            midpoint = 0.5_real64 * (lower + upper)
            cdf = 0.5_real64 * (1.0_real64 + erf(midpoint / sqrt(2.0_real64)))
            if (cdf < probability) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        quantile = 0.5_real64 * (lower + upper)
    end function normal_quantile

    pure logical function window_has_likelihood_support(incidence, infectiousness, start_idx, end_idx) result(valid)
        real(real64), intent(in) :: incidence(:), infectiousness(:)
        integer, intent(in) :: start_idx, end_idx
        integer :: idx

        valid = sum(infectiousness(start_idx:end_idx)) > 0.0_real64
        if (.not. valid) return
        do idx = start_idx, end_idx
            if (incidence(idx) > 0.0_real64 .and. infectiousness(idx) <= 0.0_real64) then
                valid = .false.
                return
            end if
        end do
    end function window_has_likelihood_support

    subroutine validate_window_inputs(incidence, infectiousness, start_idx, end_idx, &
            dispersion, confidence_level)
        real(real64), intent(in) :: incidence(:), infectiousness(:), dispersion, confidence_level
        integer, intent(in) :: start_idx, end_idx

        if (size(incidence) /= size(infectiousness)) then
            error stop "incidence and infectiousness lengths must match"
        end if
        if (start_idx < 1 .or. end_idx > size(incidence) .or. start_idx > end_idx) then
            error stop "invalid negative-binomial Rt estimation window"
        end if
        if (any(infectiousness < 0.0_real64)) error stop "infectiousness must be non-negative"
        if (dispersion <= 0.0_real64) error stop "dispersion must be positive"
        if (confidence_level <= 0.0_real64 .or. confidence_level >= 1.0_real64) then
            error stop "confidence level must lie in (0,1)"
        end if
        call validate_count_vector(incidence)
    end subroutine validate_window_inputs

    subroutine validate_count_vector(incidence)
        real(real64), intent(in) :: incidence(:)

        if (size(incidence) < 1) error stop "incidence must not be empty"
        if (any(incidence < 0.0_real64)) error stop "incidence must be non-negative"
        if (any(abs(incidence - anint(incidence)) > 1.0e-8_real64)) then
            error stop "negative-binomial incidence must be integer-valued"
        end if
    end subroutine validate_count_vector

end module negative_binomial_rt_m
