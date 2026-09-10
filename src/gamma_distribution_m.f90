module gamma_distribution_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: gamma_quantile, regularized_gamma_p, regularized_gamma_q

contains

    real(real64) function gamma_quantile(probability, shape, rate) result(quantile)
        real(real64), intent(in) :: probability, shape, rate
        real(real64) :: lower, upper, midpoint
        integer :: iteration

        if (probability <= 0.0_real64) then
            quantile = 0.0_real64
            return
        end if
        if (probability >= 1.0_real64) then
            quantile = huge(1.0_real64)
            return
        end if
        if (shape <= 0.0_real64 .or. rate <= 0.0_real64) then
            error stop "invalid Gamma parameters"
        end if

        lower = 0.0_real64
        upper = max(1.0_real64 / rate, shape / rate + 10.0_real64 * sqrt(shape) / rate)
        do while (regularized_gamma_p(shape, rate * upper) < probability)
            upper = 2.0_real64 * upper
        end do

        do iteration = 1, 100
            midpoint = 0.5_real64 * (lower + upper)
            if (regularized_gamma_p(shape, rate * midpoint) < probability) then
                lower = midpoint
            else
                upper = midpoint
            end if
        end do
        quantile = 0.5_real64 * (lower + upper)
    end function gamma_quantile

    real(real64) function regularized_gamma_p(shape, x) result(value)
        real(real64), intent(in) :: shape, x

        if (shape <= 0.0_real64) error stop "Gamma shape must be positive"
        if (x <= 0.0_real64) then
            value = 0.0_real64
            return
        end if

        if (x < shape + 1.0_real64) then
            value = lower_gamma_series(shape, x)
        else
            value = 1.0_real64 - upper_gamma_continued_fraction(shape, x)
        end if
        value = min(1.0_real64, max(0.0_real64, value))
    end function regularized_gamma_p

    real(real64) function regularized_gamma_q(shape, x) result(value)
        real(real64), intent(in) :: shape, x

        if (shape <= 0.0_real64) error stop "Gamma shape must be positive"
        if (x <= 0.0_real64) then
            value = 1.0_real64
            return
        end if

        if (x < shape + 1.0_real64) then
            value = 1.0_real64 - lower_gamma_series(shape, x)
        else
            value = upper_gamma_continued_fraction(shape, x)
        end if
        value = min(1.0_real64, max(0.0_real64, value))
    end function regularized_gamma_q

    real(real64) function lower_gamma_series(shape, x) result(value)
        real(real64), intent(in) :: shape, x
        integer, parameter :: max_iterations = 10000
        real(real64), parameter :: tolerance = 1.0e-14_real64
        real(real64) :: sum_term, term, ap
        integer :: iteration

        ap = shape
        term = 1.0_real64 / shape
        sum_term = term
        do iteration = 1, max_iterations
            ap = ap + 1.0_real64
            term = term * x / ap
            sum_term = sum_term + term
            if (abs(term) <= abs(sum_term) * tolerance) exit
        end do
        value = sum_term * exp(-x + shape * log(x) - log_gamma(shape))
    end function lower_gamma_series

    real(real64) function upper_gamma_continued_fraction(shape, x) result(value)
        real(real64), intent(in) :: shape, x
        integer, parameter :: max_iterations = 10000
        real(real64), parameter :: tolerance = 1.0e-14_real64
        real(real64), parameter :: fpmin = tiny(1.0_real64) / tolerance
        real(real64) :: term, b, c, d, h, delta
        integer :: iteration

        b = x + 1.0_real64 - shape
        c = 1.0_real64 / fpmin
        d = 1.0_real64 / b
        h = d
        do iteration = 1, max_iterations
            term = -real(iteration, real64) * (real(iteration, real64) - shape)
            b = b + 2.0_real64
            d = term * d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + term / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_real64 / d
            delta = d * c
            h = h * delta
            if (abs(delta - 1.0_real64) <= tolerance) exit
        end do
        value = exp(-x + shape * log(x) - log_gamma(shape)) * h
    end function upper_gamma_continued_fraction

end module gamma_distribution_m
