module mixed_censored_generation_interval_uncertainty_m
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng, uniform_open
    use random_distributions_m, only: gamma_random, poisson_random
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        fit_mixed_censored_gamma, discretize_interval_censored_gamma
    use reporting_delay_nowcast_m, only: delay_distribution, nowcast_summary, nowcast_series
    use renewal_rt_m, only: rt_posterior_summary, estimate_rt_series
    use interval_censored_generation_interval_uncertainty_m, only: &
        interval_censored_bootstrap_summary, interval_censored_rt_summary
    implicit none
    private

    public :: propagate_mixed_censored_generation_interval_uncertainty_rt

contains

    subroutine propagate_mixed_censored_generation_interval_uncertainty_rt(pair_lower, pair_upper, &
            right_censored, reported_counts, ages, delay, max_lag, nowcast_prior_shape, &
            nowcast_prior_rate, rt_window_size, rt_prior_shape, rt_prior_rate, credible_level, &
            seed, n_draws, rt_summaries, bootstrap_summary)
        real(real64), intent(in) :: pair_lower(:), pair_upper(:), reported_counts(:)
        logical, intent(in) :: right_censored(:)
        integer, intent(in) :: ages(:), max_lag, rt_window_size
        type(delay_distribution), intent(in) :: delay
        real(real64), intent(in) :: nowcast_prior_shape, nowcast_prior_rate
        real(real64), intent(in) :: rt_prior_shape, rt_prior_rate, credible_level
        integer(int64), intent(in) :: seed
        integer, intent(in) :: n_draws
        type(interval_censored_rt_summary), intent(out) :: rt_summaries(:)
        type(interval_censored_bootstrap_summary), intent(out) :: bootstrap_summary

        type(nowcast_summary), allocatable :: nowcasts(:)
        type(rt_posterior_summary), allocatable :: rt_conditional(:)
        type(interval_censored_gamma_fit) :: fit
        type(rng_state) :: rng
        logical, allocatable :: boot_right(:), rt_valid(:)
        real(real64), allocatable :: boot_lower(:), boot_upper(:), weights(:), complete_incidence(:)
        real(real64), allocatable :: rt_draws(:,:), shape_draws(:), rate_draws(:), mean_draws(:), tail_draws(:)
        integer, allocatable :: valid_counts(:)
        real(real64) :: latent_rate, future_mean, tail_probability
        integer(int64) :: future_reports
        integer :: n_days, n_pairs, draw_idx, day_idx, source_idx, fit_idx

        n_days = size(reported_counts)
        n_pairs = size(pair_lower)
        if (n_pairs < 2 .or. size(pair_upper) /= n_pairs .or. size(right_censored) /= n_pairs) then
            error stop "mixed-censoring arrays must have matching length >= 2"
        end if
        if (any(pair_lower < 0.0_real64)) error stop "censoring lower bounds must be non-negative"
        if (any((.not. right_censored) .and. pair_upper <= pair_lower)) then
            error stop "finite censoring upper bounds must exceed lower bounds"
        end if
        if (count(.not. right_censored) < 2) error stop "at least two finite records are required"
        if (n_days < 1 .or. size(ages) /= n_days .or. size(rt_summaries) /= n_days) then
            error stop "surveillance vectors must have matching non-empty lengths"
        end if
        if (any(reported_counts < 0.0_real64) .or. any(ages < 0)) error stop "invalid surveillance values"
        if (max_lag < 1 .or. rt_window_size < 1) error stop "lag and window sizes must be positive"
        if (n_draws < 2) error stop "at least two Monte Carlo draws are required"
        if (credible_level <= 0.0_real64 .or. credible_level >= 1.0_real64) then
            error stop "credible level must lie in (0,1)"
        end if

        allocate(nowcasts(n_days), rt_conditional(n_days), rt_valid(n_days), complete_incidence(n_days))
        allocate(boot_lower(n_pairs), boot_upper(n_pairs), boot_right(n_pairs), weights(max_lag))
        allocate(valid_counts(n_days), rt_draws(n_days,n_draws), shape_draws(n_draws), &
            rate_draws(n_draws), mean_draws(n_draws), tail_draws(n_draws))
        rt_summaries = interval_censored_rt_summary()
        bootstrap_summary = interval_censored_bootstrap_summary()
        bootstrap_summary%n_requested = n_draws
        rt_draws = 0.0_real64
        valid_counts = 0
        fit_idx = 0

        call nowcast_series(reported_counts, ages, delay, nowcast_prior_shape, nowcast_prior_rate, &
            credible_level, nowcasts)
        if (.not. all([(nowcasts(day_idx)%valid, day_idx=1,n_days)])) return

        call seed_rng(rng, seed)
        do draw_idx = 1, n_draws
            do day_idx = 1, n_pairs
                source_idx = 1 + int(floor(uniform_open(rng) * real(n_pairs, real64)))
                source_idx = min(n_pairs, max(1, source_idx))
                boot_lower(day_idx) = pair_lower(source_idx)
                boot_upper(day_idx) = pair_upper(source_idx)
                boot_right(day_idx) = right_censored(source_idx)
            end do

            if (count(.not. boot_right) < 2) cycle
            call fit_mixed_censored_gamma(boot_lower, boot_upper, boot_right, fit)
            if (.not. fit%converged) cycle
            call discretize_interval_censored_gamma(fit, max_lag, weights, tail_probability)

            fit_idx = fit_idx + 1
            shape_draws(fit_idx) = fit%shape
            rate_draws(fit_idx) = fit%rate
            mean_draws(fit_idx) = fit%mean
            tail_draws(fit_idx) = tail_probability

            do day_idx = 1, n_days
                if (nowcasts(day_idx)%completeness >= 1.0_real64 - 1.0e-15_real64) then
                    complete_incidence(day_idx) = reported_counts(day_idx)
                else
                    latent_rate = gamma_random(rng, nowcasts(day_idx)%posterior_shape, &
                        nowcasts(day_idx)%posterior_rate)
                    future_mean = latent_rate * (1.0_real64 - nowcasts(day_idx)%completeness)
                    future_reports = poisson_random(rng, future_mean)
                    complete_incidence(day_idx) = reported_counts(day_idx) + real(future_reports, real64)
                end if
            end do

            call estimate_rt_series(complete_incidence, weights, rt_window_size, rt_prior_shape, &
                rt_prior_rate, credible_level, rt_conditional, rt_valid)
            do day_idx = 1, n_days
                if (.not. rt_valid(day_idx)) cycle
                valid_counts(day_idx) = valid_counts(day_idx) + 1
                rt_draws(day_idx,valid_counts(day_idx)) = gamma_random(rng, rt_conditional(day_idx)%shape, &
                    rt_conditional(day_idx)%rate)
            end do
        end do

        bootstrap_summary%n_converged = fit_idx
        if (fit_idx >= 2) then
            call summarize_bootstrap(shape_draws(:fit_idx), rate_draws(:fit_idx), mean_draws(:fit_idx), &
                tail_draws(:fit_idx), bootstrap_summary)
        end if
        do day_idx = 1, n_days
            if (valid_counts(day_idx) >= 2) then
                call summarize_rt(rt_draws(day_idx,:valid_counts(day_idx)), credible_level, rt_summaries(day_idx))
            end if
        end do
    end subroutine propagate_mixed_censored_generation_interval_uncertainty_rt

    subroutine summarize_bootstrap(shape_draws, rate_draws, mean_draws, tail_draws, summary)
        real(real64), intent(in) :: shape_draws(:), rate_draws(:), mean_draws(:), tail_draws(:)
        type(interval_censored_bootstrap_summary), intent(inout) :: summary
        summary%shape_mean = sample_mean(shape_draws)
        summary%shape_sd = sample_sd(shape_draws)
        summary%rate_mean = sample_mean(rate_draws)
        summary%rate_sd = sample_sd(rate_draws)
        summary%generation_mean = sample_mean(mean_draws)
        summary%generation_mean_sd = sample_sd(mean_draws)
        summary%mean_tail_probability = sample_mean(tail_draws)
        summary%valid = .true.
    end subroutine summarize_bootstrap

    subroutine summarize_rt(draws, credible_level, summary)
        real(real64), intent(in) :: draws(:), credible_level
        type(interval_censored_rt_summary), intent(out) :: summary
        real(real64), allocatable :: ordered(:)
        real(real64) :: tail
        ordered = draws
        call sort_real(ordered)
        tail = 0.5_real64 * (1.0_real64 - credible_level)
        summary%mean = sample_mean(draws)
        summary%lower = empirical_quantile(ordered, tail)
        summary%upper = empirical_quantile(ordered, 1.0_real64-tail)
        summary%n_draws = size(draws)
        summary%valid = .true.
    end subroutine summarize_rt

    pure real(real64) function sample_mean(values) result(value)
        real(real64), intent(in) :: values(:)
        value = sum(values) / real(size(values), real64)
    end function sample_mean

    pure real(real64) function sample_sd(values) result(value)
        real(real64), intent(in) :: values(:)
        real(real64) :: m
        if (size(values) < 2) then
            value = 0.0_real64
            return
        end if
        m = sample_mean(values)
        value = sqrt(sum((values-m)**2) / real(size(values)-1, real64))
    end function sample_sd

    real(real64) function empirical_quantile(ordered, probability) result(value)
        real(real64), intent(in) :: ordered(:), probability
        real(real64) :: position, fraction
        integer :: lower_idx, upper_idx
        position = 1.0_real64 + min(1.0_real64,max(0.0_real64,probability)) * real(size(ordered)-1,real64)
        lower_idx = floor(position)
        upper_idx = min(lower_idx+1,size(ordered))
        fraction = position-real(lower_idx,real64)
        value = ordered(lower_idx) + fraction*(ordered(upper_idx)-ordered(lower_idx))
    end function empirical_quantile

    subroutine sort_real(values)
        real(real64), intent(inout) :: values(:)
        real(real64) :: key
        integer :: i,j
        do i=2,size(values)
            key=values(i)
            j=i-1
            do while (j>=1)
                if (values(j)<=key) exit
                values(j+1)=values(j)
                j=j-1
            end do
            values(j+1)=key
        end do
    end subroutine sort_real

end module mixed_censored_generation_interval_uncertainty_m
