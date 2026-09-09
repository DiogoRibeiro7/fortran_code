module nowcast_rt_uncertainty_m
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng
    use random_distributions_m, only: gamma_random, poisson_random
    use reporting_delay_nowcast_m, only: delay_distribution, nowcast_summary, nowcast_series
    use renewal_rt_m, only: rt_posterior_summary, estimate_rt_series
    implicit none
    private

    public :: rt_monte_carlo_summary, propagate_nowcast_uncertainty_rt

    type :: rt_monte_carlo_summary
        real(real64) :: mean = 0.0_real64
        real(real64) :: lower = 0.0_real64
        real(real64) :: upper = 0.0_real64
        integer :: n_draws = 0
        logical :: valid = .false.
    end type rt_monte_carlo_summary

contains

    subroutine propagate_nowcast_uncertainty_rt(reported_counts, ages, delay, generation_weights, &
            nowcast_prior_shape, nowcast_prior_rate, rt_window_size, rt_prior_shape, rt_prior_rate, &
            credible_level, seed, n_draws, summaries)
        real(real64), intent(in) :: reported_counts(:)
        integer, intent(in) :: ages(:)
        type(delay_distribution), intent(in) :: delay
        real(real64), intent(in) :: generation_weights(:)
        real(real64), intent(in) :: nowcast_prior_shape, nowcast_prior_rate
        integer, intent(in) :: rt_window_size
        real(real64), intent(in) :: rt_prior_shape, rt_prior_rate, credible_level
        integer(int64), intent(in) :: seed
        integer, intent(in) :: n_draws
        type(rt_monte_carlo_summary), intent(out) :: summaries(:)

        type(nowcast_summary), allocatable :: nowcasts(:)
        type(rt_posterior_summary), allocatable :: rt_conditional(:)
        logical, allocatable :: rt_valid(:)
        real(real64), allocatable :: complete_incidence(:), rt_draws(:,:)
        integer, allocatable :: valid_counts(:)
        type(rng_state) :: rng
        real(real64) :: latent_rate, future_mean
        integer(int64) :: future_reports
        integer :: n, draw_idx, day_idx

        n = size(reported_counts)
        if (n < 1) error stop "reported counts must not be empty"
        if (size(ages) /= n .or. size(summaries) /= n) then
            error stop "reported counts, ages, and summaries must have matching lengths"
        end if
        if (n_draws < 2) error stop "at least two Monte Carlo draws are required"
        if (credible_level <= 0.0_real64 .or. credible_level >= 1.0_real64) then
            error stop "credible level must lie in (0,1)"
        end if

        allocate(nowcasts(n), rt_conditional(n), rt_valid(n), complete_incidence(n))
        allocate(rt_draws(n, n_draws), valid_counts(n))
        rt_draws = 0.0_real64
        valid_counts = 0
        summaries = rt_monte_carlo_summary()

        call nowcast_series(reported_counts, ages, delay, nowcast_prior_shape, nowcast_prior_rate, &
            credible_level, nowcasts)
        if (.not. all([(nowcasts(day_idx)%valid, day_idx = 1, n)])) return

        call seed_rng(rng, seed)
        do draw_idx = 1, n_draws
            do day_idx = 1, n
                if (nowcasts(day_idx)%completeness >= 1.0_real64 - 1.0e-15_real64) then
                    complete_incidence(day_idx) = reported_counts(day_idx)
                else
                    latent_rate = gamma_random(rng, nowcasts(day_idx)%posterior_shape, &
                        nowcasts(day_idx)%posterior_rate)
                    future_mean = latent_rate * (1.0_real64 - nowcasts(day_idx)%completeness)
                    future_reports = poisson_random(rng, future_mean)
                    complete_incidence(day_idx) = reported_counts(day_idx) + &
                        real(future_reports, real64)
                end if
            end do

            call estimate_rt_series(complete_incidence, generation_weights, rt_window_size, &
                rt_prior_shape, rt_prior_rate, credible_level, rt_conditional, rt_valid)

            do day_idx = 1, n
                if (.not. rt_valid(day_idx)) cycle
                valid_counts(day_idx) = valid_counts(day_idx) + 1
                rt_draws(day_idx, valid_counts(day_idx)) = gamma_random(rng, &
                    rt_conditional(day_idx)%shape, rt_conditional(day_idx)%rate)
            end do
        end do

        do day_idx = 1, n
            if (valid_counts(day_idx) /= n_draws) cycle
            call summarize_draws(rt_draws(day_idx, :), credible_level, summaries(day_idx))
        end do
    end subroutine propagate_nowcast_uncertainty_rt

    subroutine summarize_draws(draws, credible_level, summary)
        real(real64), intent(in) :: draws(:)
        real(real64), intent(in) :: credible_level
        type(rt_monte_carlo_summary), intent(out) :: summary
        real(real64), allocatable :: ordered(:)
        real(real64) :: tail_probability

        ordered = draws
        call sort_real(ordered)
        tail_probability = 0.5_real64 * (1.0_real64 - credible_level)

        summary%mean = sum(draws) / real(size(draws), real64)
        summary%lower = empirical_quantile(ordered, tail_probability)
        summary%upper = empirical_quantile(ordered, 1.0_real64 - tail_probability)
        summary%n_draws = size(draws)
        summary%valid = .true.
    end subroutine summarize_draws

    real(real64) function empirical_quantile(ordered, probability) result(value)
        real(real64), intent(in) :: ordered(:)
        real(real64), intent(in) :: probability
        real(real64) :: position, fraction
        integer :: lower_idx, upper_idx

        if (size(ordered) < 1) error stop "quantile requires at least one draw"
        if (probability <= 0.0_real64) then
            value = ordered(1)
            return
        end if
        if (probability >= 1.0_real64) then
            value = ordered(size(ordered))
            return
        end if

        position = 1.0_real64 + probability * real(size(ordered) - 1, real64)
        lower_idx = floor(position)
        upper_idx = min(lower_idx + 1, size(ordered))
        fraction = position - real(lower_idx, real64)
        value = ordered(lower_idx) + fraction * (ordered(upper_idx) - ordered(lower_idx))
    end function empirical_quantile

    subroutine sort_real(values)
        real(real64), intent(inout) :: values(:)
        real(real64) :: key
        integer :: i, j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_real

end module nowcast_rt_uncertainty_m
