module growth_inference_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: generation_interval_distribution
    public :: generation_interval_is_valid
    public :: generation_interval_mean, generation_interval_sd
    public :: generation_interval_laplace
    public :: reproduction_number_from_growth_rate
    public :: growth_rate_from_reproduction_number
    public :: doubling_time_from_growth_rate, halving_time_from_growth_rate
    public :: estimate_exponential_growth_rate
    public :: seir_generation_interval_density
    public :: build_seir_generation_interval
    public :: seir_early_growth_rate
    public :: seir_reproduction_number_from_growth_rate

    type :: generation_interval_distribution
        real(real64), allocatable :: time(:)
        real(real64), allocatable :: probability(:)
    end type generation_interval_distribution

contains

    logical function generation_interval_is_valid(distribution) result(valid)
        type(generation_interval_distribution), intent(in) :: distribution

        valid = .false.
        if (.not. allocated(distribution%time) .or. .not. allocated(distribution%probability)) return
        if (size(distribution%time) < 1) return
        if (size(distribution%time) /= size(distribution%probability)) return
        if (any(distribution%time <= 0.0_real64)) return
        if (any(distribution%probability < 0.0_real64)) return
        if (sum(distribution%probability) <= 0.0_real64) return
        valid = .true.
    end function generation_interval_is_valid

    real(real64) function generation_interval_mean(distribution) result(mean_time)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64) :: total_probability

        call validate_generation_interval(distribution)
        total_probability = sum(distribution%probability)
        mean_time = sum(distribution%time * distribution%probability) / total_probability
    end function generation_interval_mean

    real(real64) function generation_interval_sd(distribution) result(sd_time)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64) :: mean_time, total_probability, variance

        call validate_generation_interval(distribution)
        total_probability = sum(distribution%probability)
        mean_time = generation_interval_mean(distribution)
        variance = sum(distribution%probability * (distribution%time - mean_time)**2) / total_probability
        sd_time = sqrt(max(variance, 0.0_real64))
    end function generation_interval_sd

    real(real64) function generation_interval_laplace(distribution, growth_rate) result(transform)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64), intent(in) :: growth_rate
        real(real64) :: log_transform, log_huge, log_tiny

        call validate_generation_interval(distribution)
        log_transform = log_laplace_transform(distribution, growth_rate)
        log_huge = log(huge(1.0_real64))
        log_tiny = log(tiny(1.0_real64))

        if (log_transform >= log_huge) then
            transform = huge(1.0_real64)
        else if (log_transform <= log_tiny) then
            transform = 0.0_real64
        else
            transform = exp(log_transform)
        end if
    end function generation_interval_laplace

    real(real64) function reproduction_number_from_growth_rate(distribution, growth_rate) result(reproduction_number)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64), intent(in) :: growth_rate
        real(real64) :: log_reproduction_number, log_huge, log_tiny

        call validate_generation_interval(distribution)
        log_reproduction_number = -log_laplace_transform(distribution, growth_rate)
        log_huge = log(huge(1.0_real64))
        log_tiny = log(tiny(1.0_real64))

        if (log_reproduction_number >= log_huge) then
            reproduction_number = huge(1.0_real64)
        else if (log_reproduction_number <= log_tiny) then
            reproduction_number = 0.0_real64
        else
            reproduction_number = exp(log_reproduction_number)
        end if
    end function reproduction_number_from_growth_rate

    real(real64) function growth_rate_from_reproduction_number(distribution, reproduction_number, &
            tolerance, max_iterations) result(growth_rate)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64), intent(in) :: reproduction_number
        real(real64), intent(in), optional :: tolerance
        integer, intent(in), optional :: max_iterations
        real(real64) :: tol, log_r, lower, upper, midpoint, f_lower, f_upper, f_mid, scale
        integer :: iterations, idx

        call validate_generation_interval(distribution)
        if (reproduction_number <= 0.0_real64) error stop "reproduction number must be positive"

        tol = 1.0e-12_real64
        if (present(tolerance)) tol = tolerance
        if (tol <= 0.0_real64) error stop "tolerance must be positive"

        iterations = 200
        if (present(max_iterations)) iterations = max_iterations
        if (iterations < 1) error stop "max_iterations must be positive"

        log_r = log(reproduction_number)
        if (abs(log_r) <= tol) then
            growth_rate = 0.0_real64
            return
        end if

        scale = 1.0_real64 / generation_interval_mean(distribution)
        if (reproduction_number > 1.0_real64) then
            lower = 0.0_real64
            upper = scale
            f_lower = log_r
            f_upper = euler_lotka_log_residual(distribution, log_r, upper)
            do idx = 1, 128
                if (f_upper <= 0.0_real64) exit
                upper = 2.0_real64 * upper
                f_upper = euler_lotka_log_residual(distribution, log_r, upper)
            end do
        else
            upper = 0.0_real64
            lower = -scale
            f_upper = log_r
            f_lower = euler_lotka_log_residual(distribution, log_r, lower)
            do idx = 1, 128
                if (f_lower >= 0.0_real64) exit
                lower = 2.0_real64 * lower
                f_lower = euler_lotka_log_residual(distribution, log_r, lower)
            end do
        end if

        if (f_lower < 0.0_real64 .or. f_upper > 0.0_real64) then
            error stop "failed to bracket Euler-Lotka growth-rate root"
        end if

        do idx = 1, iterations
            midpoint = 0.5_real64 * (lower + upper)
            f_mid = euler_lotka_log_residual(distribution, log_r, midpoint)

            if (abs(f_mid) <= tol .or. abs(upper - lower) <= tol * max(1.0_real64, abs(midpoint))) then
                growth_rate = midpoint
                return
            end if

            if (f_mid > 0.0_real64) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do

        growth_rate = 0.5_real64 * (lower + upper)
    end function growth_rate_from_reproduction_number

    real(real64) function doubling_time_from_growth_rate(growth_rate) result(doubling_time)
        real(real64), intent(in) :: growth_rate

        if (growth_rate <= 0.0_real64) error stop "doubling time requires positive growth rate"
        doubling_time = log(2.0_real64) / growth_rate
    end function doubling_time_from_growth_rate

    real(real64) function halving_time_from_growth_rate(growth_rate) result(halving_time)
        real(real64), intent(in) :: growth_rate

        if (growth_rate >= 0.0_real64) error stop "halving time requires negative growth rate"
        halving_time = log(2.0_real64) / abs(growth_rate)
    end function halving_time_from_growth_rate

    real(real64) function estimate_exponential_growth_rate(times, incidence) result(growth_rate)
        real(real64), intent(in) :: times(:), incidence(:)
        real(real64), allocatable :: log_incidence(:)
        real(real64) :: mean_time, mean_log, denominator

        if (size(times) /= size(incidence)) error stop "times and incidence must have equal length"
        if (size(times) < 2) error stop "at least two observations are required"
        if (any(incidence <= 0.0_real64)) error stop "incidence must be positive for log-linear fitting"

        allocate(log_incidence(size(incidence)))
        log_incidence = log(incidence)
        mean_time = sum(times) / real(size(times), real64)
        mean_log = sum(log_incidence) / real(size(log_incidence), real64)
        denominator = sum((times - mean_time)**2)
        if (denominator <= 0.0_real64) error stop "times must not all be identical"

        growth_rate = sum((times - mean_time) * (log_incidence - mean_log)) / denominator
    end function estimate_exponential_growth_rate

    real(real64) function seir_generation_interval_density(time, sigma, gamma_rate) result(density)
        real(real64), intent(in) :: time, sigma, gamma_rate
        real(real64) :: rate_scale

        if (time < 0.0_real64) error stop "generation-interval time must be non-negative"
        if (sigma <= 0.0_real64 .or. gamma_rate <= 0.0_real64) then
            error stop "sigma and gamma must be positive"
        end if

        rate_scale = max(sigma, gamma_rate)
        if (abs(sigma - gamma_rate) <= 100.0_real64 * epsilon(1.0_real64) * rate_scale) then
            density = sigma * sigma * time * exp(-sigma * time)
        else
            density = sigma * gamma_rate * (exp(-sigma * time) - exp(-gamma_rate * time)) / &
                (gamma_rate - sigma)
        end if
    end function seir_generation_interval_density

    subroutine build_seir_generation_interval(sigma, gamma_rate, dt, max_time, distribution)
        real(real64), intent(in) :: sigma, gamma_rate, dt, max_time
        type(generation_interval_distribution), intent(out) :: distribution
        integer :: idx, n_bins
        real(real64) :: total_probability

        if (sigma <= 0.0_real64 .or. gamma_rate <= 0.0_real64) then
            error stop "sigma and gamma must be positive"
        end if
        if (dt <= 0.0_real64 .or. max_time <= 0.0_real64) then
            error stop "dt and max_time must be positive"
        end if

        n_bins = int(max_time / dt)
        if (n_bins < 1) error stop "max_time must contain at least one interval"
        allocate(distribution%time(n_bins), distribution%probability(n_bins))

        do idx = 1, n_bins
            distribution%time(idx) = (real(idx, real64) - 0.5_real64) * dt
            distribution%probability(idx) = seir_generation_interval_density( &
                distribution%time(idx), sigma, gamma_rate) * dt
        end do

        total_probability = sum(distribution%probability)
        if (total_probability <= 0.0_real64) error stop "SEIR generation-interval discretization has zero mass"
        distribution%probability = distribution%probability / total_probability
    end subroutine build_seir_generation_interval

    real(real64) function seir_early_growth_rate(beta, sigma, gamma_rate) result(growth_rate)
        real(real64), intent(in) :: beta, sigma, gamma_rate
        real(real64) :: discriminant

        if (beta < 0.0_real64) error stop "beta must be non-negative"
        if (sigma <= 0.0_real64 .or. gamma_rate <= 0.0_real64) then
            error stop "sigma and gamma must be positive"
        end if

        discriminant = (sigma - gamma_rate)**2 + 4.0_real64 * sigma * beta
        growth_rate = 0.5_real64 * (-(sigma + gamma_rate) + sqrt(discriminant))
    end function seir_early_growth_rate

    real(real64) function seir_reproduction_number_from_growth_rate(growth_rate, sigma, gamma_rate) &
            result(reproduction_number)
        real(real64), intent(in) :: growth_rate, sigma, gamma_rate

        if (sigma <= 0.0_real64 .or. gamma_rate <= 0.0_real64) then
            error stop "sigma and gamma must be positive"
        end if
        if (growth_rate <= -min(sigma, gamma_rate)) then
            error stop "growth rate implies a non-positive SEIR transmission rate"
        end if

        reproduction_number = (growth_rate + sigma) * (growth_rate + gamma_rate) / &
            (sigma * gamma_rate)
    end function seir_reproduction_number_from_growth_rate

    real(real64) function log_laplace_transform(distribution, growth_rate) result(log_transform)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64), intent(in) :: growth_rate
        real(real64), allocatable :: exponent(:)
        real(real64) :: shift, scaled_sum, total_probability

        allocate(exponent(size(distribution%time)))
        exponent = -growth_rate * distribution%time
        shift = maxval(exponent)
        total_probability = sum(distribution%probability)
        scaled_sum = sum(distribution%probability * exp(exponent - shift))
        if (scaled_sum <= 0.0_real64) error stop "generation-interval Laplace transform underflowed"
        log_transform = shift + log(scaled_sum / total_probability)
    end function log_laplace_transform

    real(real64) function euler_lotka_log_residual(distribution, log_reproduction_number, growth_rate) &
            result(residual)
        type(generation_interval_distribution), intent(in) :: distribution
        real(real64), intent(in) :: log_reproduction_number, growth_rate

        residual = log_reproduction_number + log_laplace_transform(distribution, growth_rate)
    end function euler_lotka_log_residual

    subroutine validate_generation_interval(distribution)
        type(generation_interval_distribution), intent(in) :: distribution

        if (.not. generation_interval_is_valid(distribution)) then
            error stop "invalid generation-interval distribution"
        end if
    end subroutine validate_generation_interval

end module growth_inference_m
