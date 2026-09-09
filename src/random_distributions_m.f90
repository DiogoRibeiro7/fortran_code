module random_distributions_m
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, uniform_open
    implicit none
    private

    public :: gamma_random, poisson_random

contains

    recursive real(real64) function gamma_random(rng, shape, rate) result(value)
        type(rng_state), intent(inout) :: rng
        real(real64), intent(in) :: shape, rate
        real(real64) :: d, c, x, v, u

        if (shape <= 0.0_real64 .or. rate <= 0.0_real64) then
            error stop "Gamma shape and rate must be positive"
        end if

        if (shape < 1.0_real64) then
            value = gamma_random(rng, shape + 1.0_real64, rate) * &
                uniform_open(rng) ** (1.0_real64 / shape)
            return
        end if

        d = shape - 1.0_real64 / 3.0_real64
        c = 1.0_real64 / sqrt(9.0_real64 * d)
        do
            x = standard_normal_random(rng)
            v = 1.0_real64 + c * x
            if (v <= 0.0_real64) cycle
            v = v * v * v
            u = uniform_open(rng)
            if (u < 1.0_real64 - 0.0331_real64 * x**4) exit
            if (log(u) < 0.5_real64 * x*x + d * (1.0_real64 - v + log(v))) exit
        end do
        value = d * v / rate
    end function gamma_random

    integer(int64) function poisson_random(rng, mean) result(draw)
        type(rng_state), intent(inout) :: rng
        real(real64), intent(in) :: mean
        real(real64) :: threshold, product
        real(real64) :: a, b, inv_alpha, v_r, u, v, us, lhs, rhs
        integer(int64) :: candidate

        if (mean < 0.0_real64) error stop "Poisson mean must be non-negative"
        if (mean == 0.0_real64) then
            draw = 0_int64
            return
        end if

        if (mean < 10.0_real64) then
            threshold = exp(-mean)
            product = 1.0_real64
            draw = -1_int64
            do
                draw = draw + 1_int64
                product = product * uniform_open(rng)
                if (product <= threshold) exit
            end do
            return
        end if

        b = 0.931_real64 + 2.53_real64 * sqrt(mean)
        a = -0.059_real64 + 0.02483_real64 * b
        inv_alpha = 1.1239_real64 + 1.1328_real64 / (b - 3.4_real64)
        v_r = 0.9277_real64 - 3.6224_real64 / (b - 2.0_real64)

        do
            u = uniform_open(rng) - 0.5_real64
            v = uniform_open(rng)
            us = 0.5_real64 - abs(u)
            if (us <= 0.0_real64) cycle

            candidate = int(floor((2.0_real64 * a / us + b) * u + mean + 0.43_real64), int64)
            if (candidate < 0_int64) cycle
            if (us >= 0.07_real64 .and. v <= v_r) then
                draw = candidate
                return
            end if
            if (us < 0.013_real64 .and. v > us) cycle

            lhs = log(v * inv_alpha / (a / (us * us) + b))
            rhs = -mean + real(candidate, real64) * log(mean) - &
                log_gamma(real(candidate + 1_int64, real64))
            if (lhs <= rhs) then
                draw = candidate
                return
            end if
        end do
    end function poisson_random

    real(real64) function standard_normal_random(rng) result(value)
        type(rng_state), intent(inout) :: rng
        real(real64) :: u1, u2
        real(real64), parameter :: two_pi = 2.0_real64 * acos(-1.0_real64)

        u1 = uniform_open(rng)
        u2 = uniform_open(rng)
        value = sqrt(-2.0_real64 * log(u1)) * cos(two_pi * u2)
    end function standard_normal_random

end module random_distributions_m
