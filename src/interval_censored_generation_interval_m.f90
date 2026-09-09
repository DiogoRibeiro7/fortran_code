module interval_censored_generation_interval_m
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: regularized_gamma_p
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    public :: interval_censored_gamma_fit
    public :: fit_interval_censored_gamma, interval_censored_gamma_log_likelihood
    public :: discretize_interval_censored_gamma

    type :: interval_censored_gamma_fit
        integer :: n_intervals = 0
        integer :: iterations = 0
        real(real64) :: shape = 0.0_real64
        real(real64) :: rate = 0.0_real64
        real(real64) :: mean = 0.0_real64
        real(real64) :: sd = 0.0_real64
        real(real64) :: shape_se = 0.0_real64
        real(real64) :: rate_se = 0.0_real64
        real(real64) :: mean_se = 0.0_real64
        real(real64) :: log_likelihood = -huge(1.0_real64)
        logical :: converged = .false.
    end type interval_censored_gamma_fit

contains

    subroutine fit_interval_censored_gamma(lower, upper, fit)
        real(real64), intent(in) :: lower(:), upper(:)
        type(interval_censored_gamma_fit), intent(out) :: fit
        real(real64), parameter :: x_tolerance = 1.0e-10_real64
        real(real64), parameter :: f_tolerance = 1.0e-10_real64
        integer, parameter :: max_iterations = 1000
        real(real64) :: simplex(2,3), f(3), centroid(2), reflected(2), expanded(2), contracted(2)
        real(real64) :: f_reflected, f_expanded, f_contracted, span_x, span_f
        real(real64) :: mean_mid, var_mid, shape0, rate0
        real(real64), allocatable :: midpoint(:)
        integer :: iteration

        call validate_intervals(lower, upper)
        fit = interval_censored_gamma_fit()
        fit%n_intervals = size(lower)

        allocate(midpoint(size(lower)))
        midpoint = 0.5_real64 * (lower + upper)
        mean_mid = sum(midpoint) / real(size(midpoint), real64)
        var_mid = sum((midpoint - mean_mid)**2) / real(max(1, size(midpoint)-1), real64)
        if (var_mid <= tiny(1.0_real64)) return
        shape0 = max(1.0e-6_real64, mean_mid * mean_mid / var_mid)
        rate0 = max(1.0e-6_real64, mean_mid / var_mid)

        simplex(:,1) = [log(shape0), log(rate0)]
        simplex(:,2) = simplex(:,1) + [0.15_real64, 0.0_real64]
        simplex(:,3) = simplex(:,1) + [0.0_real64, 0.15_real64]
        f(1) = objective(simplex(:,1), lower, upper)
        f(2) = objective(simplex(:,2), lower, upper)
        f(3) = objective(simplex(:,3), lower, upper)

        do iteration = 1, max_iterations
            call sort_simplex(simplex, f)
            span_x = max(maxval(abs(simplex(:,2) - simplex(:,1))), &
                maxval(abs(simplex(:,3) - simplex(:,1))))
            span_f = max(abs(f(2)-f(1)), abs(f(3)-f(1)))
            if (span_x <= x_tolerance .and. span_f <= f_tolerance) exit

            centroid = 0.5_real64 * (simplex(:,1) + simplex(:,2))
            reflected = centroid + (centroid - simplex(:,3))
            f_reflected = objective(reflected, lower, upper)

            if (f_reflected < f(1)) then
                expanded = centroid + 2.0_real64 * (reflected - centroid)
                f_expanded = objective(expanded, lower, upper)
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
                    f_contracted = objective(contracted, lower, upper)
                    if (f_contracted <= f_reflected) then
                        simplex(:,3) = contracted
                        f(3) = f_contracted
                    else
                        call shrink_simplex(simplex, f, lower, upper)
                    end if
                else
                    contracted = centroid + 0.5_real64 * (simplex(:,3) - centroid)
                    f_contracted = objective(contracted, lower, upper)
                    if (f_contracted < f(3)) then
                        simplex(:,3) = contracted
                        f(3) = f_contracted
                    else
                        call shrink_simplex(simplex, f, lower, upper)
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
        if (fit%converged) call observed_standard_errors(lower, upper, simplex(:,1), fit)
    end subroutine fit_interval_censored_gamma

    real(real64) function interval_censored_gamma_log_likelihood(lower, upper, shape, rate) result(value)
        real(real64), intent(in) :: lower(:), upper(:), shape, rate
        real(real64) :: cdf_lower, cdf_upper, probability
        integer :: idx

        call validate_intervals(lower, upper)
        if (shape <= 0.0_real64 .or. rate <= 0.0_real64) then
            value = -huge(1.0_real64)
            return
        end if

        value = 0.0_real64
        do idx = 1, size(lower)
            cdf_lower = regularized_gamma_p(shape, rate * lower(idx))
            cdf_upper = regularized_gamma_p(shape, rate * upper(idx))
            probability = cdf_upper - cdf_lower
            if (probability <= 0.0_real64) then
                value = -huge(1.0_real64)
                return
            end if
            value = value + log(probability)
        end do
    end function interval_censored_gamma_log_likelihood

    subroutine discretize_interval_censored_gamma(fit, max_lag, weights, tail_probability)
        type(interval_censored_gamma_fit), intent(in) :: fit
        integer, intent(in) :: max_lag
        real(real64), intent(out) :: weights(:), tail_probability
        real(real64) :: previous_cdf, current_cdf, captured_mass
        integer :: lag

        if (.not. fit%converged) error stop "interval-censored Gamma fit has not converged"
        if (max_lag < 1 .or. size(weights) /= max_lag) error stop "invalid renewal grid"
        previous_cdf = 0.0_real64
        do lag = 1, max_lag
            current_cdf = regularized_gamma_p(fit%shape, fit%rate * real(lag, real64))
            weights(lag) = max(0.0_real64, current_cdf - previous_cdf)
            previous_cdf = current_cdf
        end do
        tail_probability = max(0.0_real64, 1.0_real64 - previous_cdf)
        captured_mass = sum(weights)
        if (captured_mass <= 0.0_real64) error stop "no renewal probability captured"
        weights = weights / captured_mass
    end subroutine discretize_interval_censored_gamma

    real(real64) function objective(z, lower, upper) result(value)
        real(real64), intent(in) :: z(2), lower(:), upper(:)
        value = -interval_censored_gamma_log_likelihood(lower, upper, exp(z(1)), exp(z(2)))
    end function objective

    subroutine shrink_simplex(simplex, f, lower, upper)
        real(real64), intent(inout) :: simplex(2,3), f(3)
        real(real64), intent(in) :: lower(:), upper(:)
        integer :: idx

        do idx = 2, 3
            simplex(:,idx) = simplex(:,1) + 0.5_real64 * (simplex(:,idx) - simplex(:,1))
            f(idx) = objective(simplex(:,idx), lower, upper)
        end do
    end subroutine shrink_simplex

    subroutine sort_simplex(simplex, f)
        real(real64), intent(inout) :: simplex(2,3), f(3)
        real(real64) :: ftmp, xtmp(2)
        integer :: i, j

        do i = 1, 2
            do j = i + 1, 3
                if (f(j) < f(i)) then
                    ftmp = f(i)
                    f(i) = f(j)
                    f(j) = ftmp
                    xtmp = simplex(:,i)
                    simplex(:,i) = simplex(:,j)
                    simplex(:,j) = xtmp
                end if
            end do
        end do
    end subroutine sort_simplex

    subroutine observed_standard_errors(lower, upper, z, fit)
        real(real64), intent(in) :: lower(:), upper(:), z(2)
        type(interval_censored_gamma_fit), intent(inout) :: fit
        real(real64), parameter :: h = 1.0e-4_real64
        real(real64) :: f0, fp1, fm1, fp2, fm2, fpp, fpm, fmp, fmm
        real(real64) :: h11, h22, h12, determinant, cov11, cov22, cov12, mean_variance
        real(real64) :: dz1(2), dz2(2)

        dz1 = [h, 0.0_real64]
        dz2 = [0.0_real64, h]
        f0 = objective(z, lower, upper)
        fp1 = objective(z + dz1, lower, upper)
        fm1 = objective(z - dz1, lower, upper)
        fp2 = objective(z + dz2, lower, upper)
        fm2 = objective(z - dz2, lower, upper)
        fpp = objective(z + dz1 + dz2, lower, upper)
        fpm = objective(z + dz1 - dz2, lower, upper)
        fmp = objective(z - dz1 + dz2, lower, upper)
        fmm = objective(z - dz1 - dz2, lower, upper)
        h11 = (fp1 - 2.0_real64*f0 + fm1) / (h*h)
        h22 = (fp2 - 2.0_real64*f0 + fm2) / (h*h)
        h12 = (fpp - fpm - fmp + fmm) / (4.0_real64*h*h)
        determinant = h11*h22 - h12*h12
        if (h11 <= 0.0_real64 .or. h22 <= 0.0_real64 .or. determinant <= 0.0_real64) return
        cov11 = h22 / determinant
        cov22 = h11 / determinant
        cov12 = -h12 / determinant
        fit%shape_se = fit%shape * sqrt(cov11)
        fit%rate_se = fit%rate * sqrt(cov22)
        mean_variance = fit%mean**2 * max(0.0_real64, cov11 + cov22 - 2.0_real64*cov12)
        fit%mean_se = sqrt(mean_variance)
    end subroutine observed_standard_errors

    subroutine validate_intervals(lower, upper)
        real(real64), intent(in) :: lower(:), upper(:)

        if (size(lower) < 2 .or. size(lower) /= size(upper)) then
            error stop "interval bounds must have matching length >= 2"
        end if
        if (any(lower < 0.0_real64)) error stop "lower bounds must be non-negative"
        if (any(upper <= lower)) error stop "each upper bound must exceed its lower bound"
    end subroutine validate_intervals

end module interval_censored_generation_interval_m
