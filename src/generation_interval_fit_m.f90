module generation_interval_fit_m
    use iso_fortran_env, only: real64
    use gamma_distribution_m, only: regularized_gamma_p
    implicit none
    private

    public :: gamma_generation_interval_fit
    public :: fit_gamma_generation_interval, discretize_gamma_generation_interval
    public :: gamma_generation_interval_log_likelihood

    type :: gamma_generation_interval_fit
        integer :: n_observations = 0
        real(real64) :: shape = 0.0_real64
        real(real64) :: rate = 0.0_real64
        real(real64) :: mean = 0.0_real64
        real(real64) :: sd = 0.0_real64
        real(real64) :: shape_se = 0.0_real64
        real(real64) :: rate_se = 0.0_real64
        real(real64) :: mean_se = 0.0_real64
        real(real64) :: log_likelihood = -huge(1.0_real64)
        logical :: converged = .false.
    end type gamma_generation_interval_fit

contains

    subroutine fit_gamma_generation_interval(intervals, fit)
        real(real64), intent(in) :: intervals(:)
        type(gamma_generation_interval_fit), intent(out) :: fit
        real(real64), parameter :: root_tolerance = 1.0e-12_real64
        real(real64), parameter :: degeneracy_tolerance = 1.0e-14_real64
        real(real64) :: sample_mean, mean_log, target
        real(real64) :: lower, upper, midpoint, f_mid
        integer :: iteration

        fit = gamma_generation_interval_fit()
        if (size(intervals) < 2) error stop "at least two transmission-pair intervals are required"
        if (any(intervals <= 0.0_real64)) error stop "generation intervals must be strictly positive"

        fit%n_observations = size(intervals)
        sample_mean = sum(intervals) / real(size(intervals), real64)
        mean_log = sum(log(intervals)) / real(size(intervals), real64)
        target = log(sample_mean) - mean_log
        fit%mean = sample_mean

        if (target <= degeneracy_tolerance) return

        lower = 1.0e-10_real64
        upper = 1.0_real64
        do while (shape_equation(upper, target) > 0.0_real64)
            upper = 2.0_real64 * upper
            if (upper > 1.0e14_real64) return
        end do

        do iteration = 1, 200
            midpoint = 0.5_real64 * (lower + upper)
            f_mid = shape_equation(midpoint, target)
            if (f_mid > 0.0_real64) then
                lower = midpoint
            else
                upper = midpoint
            end if
            if (abs(upper - lower) <= root_tolerance * max(1.0_real64, midpoint)) exit
        end do

        fit%shape = 0.5_real64 * (lower + upper)
        fit%rate = fit%shape / sample_mean
        fit%sd = sqrt(fit%shape) / fit%rate
        fit%log_likelihood = gamma_generation_interval_log_likelihood(intervals, fit%shape, fit%rate)
        call asymptotic_standard_errors(fit)
        fit%converged = .true.
    end subroutine fit_gamma_generation_interval

    real(real64) function gamma_generation_interval_log_likelihood(intervals, shape, rate) result(value)
        real(real64), intent(in) :: intervals(:), shape, rate
        real(real64) :: n

        if (size(intervals) < 1) error stop "intervals must not be empty"
        if (any(intervals <= 0.0_real64)) error stop "generation intervals must be strictly positive"
        if (shape <= 0.0_real64 .or. rate <= 0.0_real64) then
            value = -huge(1.0_real64)
            return
        end if

        n = real(size(intervals), real64)
        value = n * (shape * log(rate) - log_gamma(shape)) + &
            (shape - 1.0_real64) * sum(log(intervals)) - rate * sum(intervals)
    end function gamma_generation_interval_log_likelihood

    subroutine discretize_gamma_generation_interval(fit, max_lag, weights, tail_probability)
        type(gamma_generation_interval_fit), intent(in) :: fit
        integer, intent(in) :: max_lag
        real(real64), intent(out) :: weights(:)
        real(real64), intent(out) :: tail_probability
        real(real64) :: previous_cdf, current_cdf, captured_mass
        integer :: lag

        if (.not. fit%converged) error stop "generation-interval fit has not converged"
        if (max_lag < 1) error stop "max_lag must be positive"
        if (size(weights) /= max_lag) error stop "weights size must equal max_lag"

        previous_cdf = 0.0_real64
        do lag = 1, max_lag
            current_cdf = regularized_gamma_p(fit%shape, fit%rate * real(lag, real64))
            weights(lag) = max(0.0_real64, current_cdf - previous_cdf)
            previous_cdf = current_cdf
        end do

        tail_probability = max(0.0_real64, 1.0_real64 - previous_cdf)
        captured_mass = sum(weights)
        if (captured_mass <= 0.0_real64) error stop "no generation-interval probability captured"
        weights = weights / captured_mass
    end subroutine discretize_gamma_generation_interval

    pure real(real64) function shape_equation(shape, target) result(value)
        real(real64), intent(in) :: shape, target
        value = log(shape) - digamma_positive(shape) - target
    end function shape_equation

    pure real(real64) function digamma_positive(x_input) result(value)
        real(real64), intent(in) :: x_input
        real(real64) :: x, inv, inv2

        x = x_input
        value = 0.0_real64
        do while (x < 8.0_real64)
            value = value - 1.0_real64 / x
            x = x + 1.0_real64
        end do

        inv = 1.0_real64 / x
        inv2 = inv * inv
        value = value + log(x) - 0.5_real64 * inv - inv2 * &
            (1.0_real64 / 12.0_real64 - inv2 * &
            (1.0_real64 / 120.0_real64 - inv2 * &
            (1.0_real64 / 252.0_real64 - inv2 * &
            (1.0_real64 / 240.0_real64 - inv2 * 5.0_real64 / 660.0_real64))))
    end function digamma_positive

    pure real(real64) function trigamma_positive(x_input) result(value)
        real(real64), intent(in) :: x_input
        real(real64) :: x, inv, inv2

        x = x_input
        value = 0.0_real64
        do while (x < 8.0_real64)
            value = value + 1.0_real64 / (x * x)
            x = x + 1.0_real64
        end do

        inv = 1.0_real64 / x
        inv2 = inv * inv
        value = value + inv + 0.5_real64 * inv2 + inv * inv2 / 6.0_real64 - &
            inv * inv2 * inv2 / 30.0_real64 + inv * inv2 * inv2 * inv2 / 42.0_real64 - &
            inv * inv2**4 / 30.0_real64 + 5.0_real64 * inv * inv2**5 / 66.0_real64
    end function trigamma_positive

    subroutine asymptotic_standard_errors(fit)
        type(gamma_generation_interval_fit), intent(inout) :: fit
        real(real64) :: n, info_shape, info_rate, info_cross, determinant
        real(real64) :: cov_shape, cov_rate, cov_cross
        real(real64) :: dmean_dshape, dmean_drate, mean_variance

        n = real(fit%n_observations, real64)
        info_shape = n * trigamma_positive(fit%shape)
        info_rate = n * fit%shape / (fit%rate * fit%rate)
        info_cross = -n / fit%rate
        determinant = info_shape * info_rate - info_cross * info_cross
        if (determinant <= 0.0_real64) return

        cov_shape = info_rate / determinant
        cov_rate = info_shape / determinant
        cov_cross = -info_cross / determinant
        fit%shape_se = sqrt(max(0.0_real64, cov_shape))
        fit%rate_se = sqrt(max(0.0_real64, cov_rate))

        dmean_dshape = 1.0_real64 / fit%rate
        dmean_drate = -fit%shape / (fit%rate * fit%rate)
        mean_variance = dmean_dshape**2 * cov_shape + dmean_drate**2 * cov_rate + &
            2.0_real64 * dmean_dshape * dmean_drate * cov_cross
        fit%mean_se = sqrt(max(0.0_real64, mean_variance))
    end subroutine asymptotic_standard_errors

end module generation_interval_fit_m
