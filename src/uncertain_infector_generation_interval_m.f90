module uncertain_infector_generation_interval_m
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: regularized_gamma_p, regularized_gamma_q
    use interval_censored_generation_interval_m, only: interval_censored_gamma_fit, &
        discretize_interval_censored_gamma
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: fit_uncertain_infector_gamma
    public :: uncertain_infector_gamma_log_likelihood

contains

    subroutine fit_uncertain_infector_gamma(candidate_start, candidate_end, lower, upper, &
            right_censored, prior_weight, fit)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        real(real64), intent(in) :: lower(:), upper(:), prior_weight(:)
        logical, intent(in) :: right_censored(:)
        type(interval_censored_gamma_fit), intent(out) :: fit
        real(real64), parameter :: x_tolerance = 1.0e-10_real64
        real(real64), parameter :: f_tolerance = 1.0e-10_real64
        integer, parameter :: max_iterations = 1200
        real(real64) :: simplex(2,3), f(3), centroid(2), reflected(2), expanded(2), contracted(2)
        real(real64) :: f_reflected, f_expanded, f_contracted, span_x, span_f
        real(real64) :: mean0, var0, shape0, rate0
        integer :: iteration

        call validate_candidates(candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        fit = interval_censored_gamma_fit()
        fit%n_intervals = size(candidate_start)

        call weighted_initial_moments(candidate_start, candidate_end, lower, upper, right_censored, &
            prior_weight, mean0, var0)
        if (var0 <= tiny(1.0_real64)) return
        shape0 = max(1.0e-6_real64, mean0 * mean0 / var0)
        rate0 = max(1.0e-6_real64, mean0 / var0)

        simplex(:,1) = [log(shape0), log(rate0)]
        simplex(:,2) = simplex(:,1) + [0.15_real64, 0.0_real64]
        simplex(:,3) = simplex(:,1) + [0.0_real64, 0.15_real64]
        f(1) = objective(simplex(:,1), candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        f(2) = objective(simplex(:,2), candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        f(3) = objective(simplex(:,3), candidate_start, candidate_end, lower, upper, right_censored, prior_weight)

        do iteration = 1, max_iterations
            call sort_simplex(simplex, f)
            span_x = max(maxval(abs(simplex(:,2)-simplex(:,1))), maxval(abs(simplex(:,3)-simplex(:,1))))
            span_f = max(abs(f(2)-f(1)), abs(f(3)-f(1)))
            if (span_x <= x_tolerance .and. span_f <= f_tolerance) exit

            centroid = 0.5_real64 * (simplex(:,1) + simplex(:,2))
            reflected = centroid + (centroid - simplex(:,3))
            f_reflected = objective(reflected, candidate_start, candidate_end, lower, upper, &
                right_censored, prior_weight)

            if (f_reflected < f(1)) then
                expanded = centroid + 2.0_real64 * (reflected - centroid)
                f_expanded = objective(expanded, candidate_start, candidate_end, lower, upper, &
                    right_censored, prior_weight)
                if (f_expanded < f_reflected) then
                    simplex(:,3) = expanded
                    f(3) = f_expanded
                else
                    simplex(:,3) = reflected
                    f(3) = f_reflected
                end if
            else if (f_reflected < f(2)) then
                simplex(:,3) = reflected
                f(3) = f_reflected
            else
                if (f_reflected < f(3)) then
                    contracted = centroid + 0.5_real64 * (reflected - centroid)
                    f_contracted = objective(contracted, candidate_start, candidate_end, lower, upper, &
                        right_censored, prior_weight)
                    if (f_contracted <= f_reflected) then
                        simplex(:,3) = contracted
                        f(3) = f_contracted
                    else
                        call shrink(simplex, f, candidate_start, candidate_end, lower, upper, &
                            right_censored, prior_weight)
                    end if
                else
                    contracted = centroid + 0.5_real64 * (simplex(:,3) - centroid)
                    f_contracted = objective(contracted, candidate_start, candidate_end, lower, upper, &
                        right_censored, prior_weight)
                    if (f_contracted < f(3)) then
                        simplex(:,3) = contracted
                        f(3) = f_contracted
                    else
                        call shrink(simplex, f, candidate_start, candidate_end, lower, upper, &
                            right_censored, prior_weight)
                    end if
                end if
            end if
        end do

        call sort_simplex(simplex, f)
        fit%iterations = min(iteration, max_iterations)
        fit%shape = exp(simplex(1,1))
        fit%rate = exp(simplex(2,1))
        fit%mean = fit%shape / fit%rate
        fit%sd = sqrt(fit%shape) / fit%rate
        fit%log_likelihood = -f(1)
        fit%converged = iteration <= max_iterations .and. ieee_is_finite(fit%log_likelihood)
    end subroutine fit_uncertain_infector_gamma

    real(real64) function uncertain_infector_gamma_log_likelihood(candidate_start, candidate_end, &
            lower, upper, right_censored, prior_weight, shape, rate) result(value)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        real(real64), intent(in) :: lower(:), upper(:), prior_weight(:), shape, rate
        logical, intent(in) :: right_censored(:)
        real(real64) :: weight_sum, log_scale, term, mixture, probability
        integer :: infectee, candidate

        call validate_candidates(candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        if (shape <= 0.0_real64 .or. rate <= 0.0_real64) then
            value = -huge(1.0_real64)
            return
        end if

        value = 0.0_real64
        do infectee = 1, size(candidate_start)
            weight_sum = sum(prior_weight(candidate_start(infectee):candidate_end(infectee)))
            log_scale = -huge(1.0_real64)
            do candidate = candidate_start(infectee), candidate_end(infectee)
                if (prior_weight(candidate) <= 0.0_real64) cycle
                probability = candidate_probability(lower(candidate), upper(candidate), &
                    right_censored(candidate), shape, rate)
                if (probability <= 0.0_real64) cycle
                term = log(prior_weight(candidate) / weight_sum) + log(probability)
                log_scale = max(log_scale, term)
            end do
            if (.not. ieee_is_finite(log_scale)) then
                value = -huge(1.0_real64)
                return
            end if
            mixture = 0.0_real64
            do candidate = candidate_start(infectee), candidate_end(infectee)
                if (prior_weight(candidate) <= 0.0_real64) cycle
                probability = candidate_probability(lower(candidate), upper(candidate), &
                    right_censored(candidate), shape, rate)
                if (probability <= 0.0_real64) cycle
                term = log(prior_weight(candidate) / weight_sum) + log(probability)
                mixture = mixture + exp(term - log_scale)
            end do
            value = value + log_scale + log(mixture)
        end do
    end function uncertain_infector_gamma_log_likelihood

    real(real64) function candidate_probability(lower, upper, right_censored, shape, rate) result(probability)
        real(real64), intent(in) :: lower, upper, shape, rate
        logical, intent(in) :: right_censored
        real(real64) :: p_lower, p_upper, q_lower, q_upper

        if (right_censored) then
            probability = regularized_gamma_q(shape, rate * lower)
        else
            p_lower = regularized_gamma_p(shape, rate * lower)
            p_upper = regularized_gamma_p(shape, rate * upper)
            q_lower = regularized_gamma_q(shape, rate * lower)
            q_upper = regularized_gamma_q(shape, rate * upper)
            probability = max(p_upper-p_lower, q_lower-q_upper)
        end if
    end function candidate_probability

    subroutine weighted_initial_moments(candidate_start, candidate_end, lower, upper, right_censored, &
            prior_weight, mean_value, variance_value)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        real(real64), intent(in) :: lower(:), upper(:), prior_weight(:)
        logical, intent(in) :: right_censored(:)
        real(real64), intent(out) :: mean_value, variance_value
        real(real64) :: x, group_total, global_weight
        integer :: infectee, candidate

        mean_value = 0.0_real64
        variance_value = 0.0_real64
        global_weight = real(size(candidate_start), real64)
        do infectee = 1, size(candidate_start)
            group_total = sum(prior_weight(candidate_start(infectee):candidate_end(infectee)))
            do candidate = candidate_start(infectee), candidate_end(infectee)
                if (right_censored(candidate)) then
                    x = lower(candidate) + 1.0_real64
                else
                    x = 0.5_real64 * (lower(candidate) + upper(candidate))
                end if
                mean_value = mean_value + prior_weight(candidate) * x / group_total
            end do
        end do
        mean_value = mean_value / global_weight

        do infectee = 1, size(candidate_start)
            group_total = sum(prior_weight(candidate_start(infectee):candidate_end(infectee)))
            do candidate = candidate_start(infectee), candidate_end(infectee)
                if (right_censored(candidate)) then
                    x = lower(candidate) + 1.0_real64
                else
                    x = 0.5_real64 * (lower(candidate) + upper(candidate))
                end if
                variance_value = variance_value + prior_weight(candidate) * (x-mean_value)**2 / group_total
            end do
        end do
        variance_value = variance_value / global_weight
    end subroutine weighted_initial_moments

    real(real64) function objective(z, candidate_start, candidate_end, lower, upper, &
            right_censored, prior_weight) result(value)
        real(real64), intent(in) :: z(2), lower(:), upper(:), prior_weight(:)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        logical, intent(in) :: right_censored(:)
        value = -uncertain_infector_gamma_log_likelihood(candidate_start, candidate_end, lower, upper, &
            right_censored, prior_weight, exp(z(1)), exp(z(2)))
    end function objective

    subroutine shrink(simplex, f, candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        real(real64), intent(inout) :: simplex(2,3), f(3)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        real(real64), intent(in) :: lower(:), upper(:), prior_weight(:)
        logical, intent(in) :: right_censored(:)
        integer :: idx
        do idx = 2, 3
            simplex(:,idx) = simplex(:,1) + 0.5_real64 * (simplex(:,idx)-simplex(:,1))
            f(idx) = objective(simplex(:,idx), candidate_start, candidate_end, lower, upper, &
                right_censored, prior_weight)
        end do
    end subroutine shrink

    subroutine sort_simplex(simplex, f)
        real(real64), intent(inout) :: simplex(2,3), f(3)
        real(real64) :: ftmp, xtmp(2)
        integer :: i, j
        do i = 1, 2
            do j = i+1, 3
                if (f(j) < f(i)) then
                    ftmp=f(i); f(i)=f(j); f(j)=ftmp
                    xtmp=simplex(:,i); simplex(:,i)=simplex(:,j); simplex(:,j)=xtmp
                end if
            end do
        end do
    end subroutine sort_simplex

    subroutine validate_candidates(candidate_start, candidate_end, lower, upper, right_censored, prior_weight)
        integer, intent(in) :: candidate_start(:), candidate_end(:)
        real(real64), intent(in) :: lower(:), upper(:), prior_weight(:)
        logical, intent(in) :: right_censored(:)
        integer :: infectee, n_candidates

        n_candidates = size(lower)
        if (size(candidate_start) < 2 .or. size(candidate_end) /= size(candidate_start)) then
            error stop "at least two infectees with matching candidate ranges are required"
        end if
        if (size(upper) /= n_candidates .or. size(right_censored) /= n_candidates .or. &
                size(prior_weight) /= n_candidates) error stop "candidate arrays must have matching lengths"
        if (n_candidates < size(candidate_start)) error stop "candidate array is too short"
        if (any(lower < 0.0_real64) .or. any(prior_weight < 0.0_real64)) then
            error stop "candidate lower bounds and weights must be non-negative"
        end if
        if (any((.not. right_censored) .and. upper <= lower)) then
            error stop "finite candidate upper bounds must exceed lower bounds"
        end if
        do infectee = 1, size(candidate_start)
            if (candidate_start(infectee) < 1 .or. candidate_end(infectee) > n_candidates .or. &
                    candidate_end(infectee) < candidate_start(infectee)) error stop "invalid candidate range"
            if (sum(prior_weight(candidate_start(infectee):candidate_end(infectee))) <= 0.0_real64) then
                error stop "each infectee must have positive total candidate weight"
            end if
            if (infectee > 1) then
                if (candidate_start(infectee) /= candidate_end(infectee-1)+1) then
                    error stop "candidate ranges must be contiguous and non-overlapping"
                end if
            end if
        end do
        if (candidate_start(1) /= 1 .or. candidate_end(size(candidate_end)) /= n_candidates) then
            error stop "candidate ranges must cover the full candidate array"
        end if
    end subroutine validate_candidates

end module uncertain_infector_generation_interval_m
