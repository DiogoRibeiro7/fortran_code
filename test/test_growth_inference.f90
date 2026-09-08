program test_growth_inference
    use iso_fortran_env, only: real64
    use growth_inference_m, only: generation_interval_distribution, generation_interval_mean, &
        generation_interval_sd, generation_interval_laplace, reproduction_number_from_growth_rate, &
        growth_rate_from_reproduction_number, doubling_time_from_growth_rate, &
        halving_time_from_growth_rate, estimate_exponential_growth_rate, &
        build_seir_generation_interval, seir_early_growth_rate, &
        seir_reproduction_number_from_growth_rate
    implicit none

    call test_point_mass_euler_lotka()
    call test_shape_matters()
    call test_round_trip_growth_and_reproduction()
    call test_declining_epidemic()
    call test_log_linear_growth_fit()
    call test_seir_generation_interval()
    call test_seir_mechanistic_growth_relation()

    print *, "Growth inference tests passed"

contains

    subroutine test_point_mass_euler_lotka()
        type(generation_interval_distribution) :: distribution
        real(real64), parameter :: r = 0.10_real64
        real(real64) :: expected_r

        allocate(distribution%time(1), distribution%probability(1))
        distribution%time = [5.0_real64]
        distribution%probability = [1.0_real64]
        expected_r = exp(r * 5.0_real64)

        call assert_close(generation_interval_mean(distribution), 5.0_real64, 1.0e-12_real64, &
            "point-mass mean")
        call assert_close(generation_interval_sd(distribution), 0.0_real64, 1.0e-12_real64, &
            "point-mass standard deviation")
        call assert_close(generation_interval_laplace(distribution, r), exp(-0.5_real64), &
            1.0e-12_real64, "point-mass Laplace transform")
        call assert_close(reproduction_number_from_growth_rate(distribution, r), expected_r, &
            1.0e-12_real64, "point-mass Euler-Lotka R")
    end subroutine test_point_mass_euler_lotka

    subroutine test_shape_matters()
        type(generation_interval_distribution) :: concentrated, dispersed
        real(real64), parameter :: r = 0.10_real64
        real(real64) :: r_concentrated, r_dispersed

        allocate(concentrated%time(1), concentrated%probability(1))
        concentrated%time = [5.0_real64]
        concentrated%probability = [1.0_real64]

        allocate(dispersed%time(2), dispersed%probability(2))
        dispersed%time = [2.0_real64, 8.0_real64]
        dispersed%probability = [0.5_real64, 0.5_real64]

        call assert_close(generation_interval_mean(concentrated), generation_interval_mean(dispersed), &
            1.0e-12_real64, "equal generation-interval means")
        r_concentrated = reproduction_number_from_growth_rate(concentrated, r)
        r_dispersed = reproduction_number_from_growth_rate(dispersed, r)
        if (r_dispersed >= r_concentrated) then
            error stop "generation-interval shape should affect R at fixed growth rate"
        end if
    end subroutine test_shape_matters

    subroutine test_round_trip_growth_and_reproduction()
        type(generation_interval_distribution) :: distribution
        real(real64), parameter :: target_r = 0.07_real64
        real(real64) :: reproduction_number, recovered_r

        allocate(distribution%time(3), distribution%probability(3))
        distribution%time = [2.0_real64, 5.0_real64, 9.0_real64]
        distribution%probability = [0.2_real64, 0.5_real64, 0.3_real64]

        reproduction_number = reproduction_number_from_growth_rate(distribution, target_r)
        recovered_r = growth_rate_from_reproduction_number(distribution, reproduction_number)
        call assert_close(recovered_r, target_r, 1.0e-10_real64, "Euler-Lotka round trip")
        call assert_close(growth_rate_from_reproduction_number(distribution, 1.0_real64), &
            0.0_real64, 1.0e-12_real64, "R=1 threshold")
    end subroutine test_round_trip_growth_and_reproduction

    subroutine test_declining_epidemic()
        type(generation_interval_distribution) :: distribution
        real(real64) :: r

        allocate(distribution%time(2), distribution%probability(2))
        distribution%time = [3.0_real64, 7.0_real64]
        distribution%probability = [0.4_real64, 0.6_real64]

        r = growth_rate_from_reproduction_number(distribution, 0.8_real64)
        if (r >= 0.0_real64) error stop "R below one should imply negative growth"
        call assert_close(reproduction_number_from_growth_rate(distribution, r), 0.8_real64, &
            1.0e-10_real64, "declining epidemic round trip")
        if (halving_time_from_growth_rate(r) <= 0.0_real64) error stop "halving time must be positive"
    end subroutine test_declining_epidemic

    subroutine test_log_linear_growth_fit()
        real(real64) :: times(8), incidence(8), fitted_r
        real(real64), parameter :: true_r = 0.12_real64
        integer :: idx

        do idx = 1, size(times)
            times(idx) = real(idx - 1, real64)
            incidence(idx) = 7.0_real64 * exp(true_r * times(idx))
        end do

        fitted_r = estimate_exponential_growth_rate(times, incidence)
        call assert_close(fitted_r, true_r, 1.0e-12_real64, "log-linear growth-rate fit")
        call assert_close(doubling_time_from_growth_rate(fitted_r), log(2.0_real64) / true_r, &
            1.0e-12_real64, "doubling time")
    end subroutine test_log_linear_growth_fit

    subroutine test_seir_generation_interval()
        type(generation_interval_distribution) :: distribution
        real(real64), parameter :: sigma = 1.0_real64 / 3.0_real64
        real(real64), parameter :: gamma_rate = 1.0_real64 / 5.0_real64
        real(real64), parameter :: target_r = 0.10_real64
        real(real64) :: numerical_r0, analytic_r0

        call build_seir_generation_interval(sigma, gamma_rate, 0.01_real64, 100.0_real64, distribution)
        call assert_close(generation_interval_mean(distribution), 8.0_real64, 2.0e-3_real64, &
            "SEIR generation-interval mean")
        call assert_close(generation_interval_sd(distribution), sqrt(34.0_real64), 3.0e-3_real64, &
            "SEIR generation-interval standard deviation")

        numerical_r0 = reproduction_number_from_growth_rate(distribution, target_r)
        analytic_r0 = seir_reproduction_number_from_growth_rate(target_r, sigma, gamma_rate)
        call assert_close(numerical_r0, analytic_r0, 3.0e-5_real64, &
            "SEIR Euler-Lotka versus analytic growth relation")
    end subroutine test_seir_generation_interval

    subroutine test_seir_mechanistic_growth_relation()
        real(real64), parameter :: sigma = 0.25_real64
        real(real64), parameter :: gamma_rate = 0.10_real64
        real(real64), parameter :: beta = 0.28_real64
        real(real64) :: r, r0

        r = seir_early_growth_rate(beta, sigma, gamma_rate)
        r0 = seir_reproduction_number_from_growth_rate(r, sigma, gamma_rate)
        call assert_close(r0, beta / gamma_rate, 1.0e-12_real64, &
            "SEIR early-growth inversion")
    end subroutine test_seir_mechanistic_growth_relation

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, "Assertion failed: ", trim(message)
            print *, "  actual:   ", actual
            print *, "  expected: ", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_growth_inference
