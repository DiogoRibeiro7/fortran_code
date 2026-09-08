program stochastic_seir_ensemble_example
    use iso_fortran_env, only: int64, real64
    use stochastic_seir_m, only: stochastic_seir_parameters, stochastic_seir_state, &
        seir_ensemble_summary, basic_reproduction_number_stochastic_seir, simulate_seir_ensemble
    implicit none

    type(stochastic_seir_parameters) :: params
    type(stochastic_seir_state) :: initial
    type(seir_ensemble_summary) :: summary

    params%population = 1000_int64
    params%beta = 0.30_real64
    params%sigma = 0.20_real64
    params%gamma = 0.10_real64

    initial%s = 999_int64
    initial%e = 1_int64
    initial%i = 0_int64
    initial%r = 0_int64
    initial%time = 0.0_real64

    call simulate_seir_ensemble(initial, params, 20260908_int64, 1000, 365.0_real64, &
        4000, 100_int64, summary)

    print '(a,f8.4)', 'R0 = ', basic_reproduction_number_stochastic_seir(params)
    print '(a,i0)', 'replicates = ', summary%n_replicates
    print '(a,f8.4)', 'extinction_probability = ', summary%extinction_probability
    print '(a,f8.4)', 'fadeout_probability = ', summary%fadeout_probability
    print '(a,f8.4)', 'major_outbreak_probability = ', summary%major_outbreak_probability
    print '(a,f10.3)', 'mean_final_size = ', summary%mean_final_size
    print '(a,f10.3)', 'sd_final_size = ', summary%sd_final_size
    print '(a,f10.3)', 'mean_peak_infectious = ', summary%mean_peak_infectious
end program stochastic_seir_ensemble_example
