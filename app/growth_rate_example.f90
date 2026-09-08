program growth_rate_example
    use iso_fortran_env, only: real64
    use growth_inference_m, only: generation_interval_distribution, generation_interval_mean, &
        generation_interval_sd, reproduction_number_from_growth_rate, build_seir_generation_interval, &
        seir_reproduction_number_from_growth_rate
    implicit none

    real(real64), parameter :: sigma = 1.0_real64 / 3.0_real64
    real(real64), parameter :: gamma_rate = 1.0_real64 / 5.0_real64
    real(real64), parameter :: observed_doubling_time = 5.0_real64
    type(generation_interval_distribution) :: distribution
    real(real64) :: growth_rate, mean_generation, sd_generation
    real(real64) :: euler_lotka_r0, mechanistic_r0, mean_only_r0

    growth_rate = log(2.0_real64) / observed_doubling_time
    call build_seir_generation_interval(sigma, gamma_rate, 0.01_real64, 100.0_real64, distribution)

    mean_generation = generation_interval_mean(distribution)
    sd_generation = generation_interval_sd(distribution)
    euler_lotka_r0 = reproduction_number_from_growth_rate(distribution, growth_rate)
    mechanistic_r0 = seir_reproduction_number_from_growth_rate(growth_rate, sigma, gamma_rate)
    mean_only_r0 = exp(growth_rate * mean_generation)

    print '(a, f10.6)', 'growth rate r (per day):          ', growth_rate
    print '(a, f10.4)', 'doubling time (days):            ', observed_doubling_time
    print '(a, f10.4)', 'generation-interval mean (days): ', mean_generation
    print '(a, f10.4)', 'generation-interval SD (days):   ', sd_generation
    print '(a, f10.4)', 'R0 from Euler-Lotka:             ', euler_lotka_r0
    print '(a, f10.4)', 'R0 from analytic SEIR relation:  ', mechanistic_r0
    print '(a, f10.4)', 'R0 using exp(r * mean GI):       ', mean_only_r0
end program growth_rate_example
