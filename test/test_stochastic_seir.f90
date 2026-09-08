program test_stochastic_seir
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng
    use stochastic_seir_m, only: stochastic_seir_parameters, stochastic_seir_state, &
        stochastic_seir_rates, seir_ensemble_summary, event_none, event_infection, &
        event_progression, event_recovery, basic_reproduction_number_stochastic_seir, &
        seir_event_rates, simulate_seir_gillespie, simulate_seir_ensemble
    use seir_m, only: seir_parameters, seir_state, simulate_seir
    implicit none

    call test_rates_and_r0()
    call test_progression_recovery_chain()
    call test_absorbing_state()
    call test_seeded_event_prefix()
    call test_ensemble_reproducibility()
    call test_fadeout_probability()
    call test_deterministic_ensemble_agreement()

    print *, "Stochastic SEIR tests passed"

contains

    subroutine test_rates_and_r0()
        type(stochastic_seir_parameters) :: params
        type(stochastic_seir_state) :: state
        type(stochastic_seir_rates) :: rates

        params%population = 1000_int64
        params%beta = 0.30_real64
        params%sigma = 0.20_real64
        params%gamma = 0.10_real64
        state%s = 995_int64
        state%e = 3_int64
        state%i = 2_int64
        state%r = 0_int64

        call assert_close(basic_reproduction_number_stochastic_seir(params), 3.0_real64, &
            1.0e-12_real64, "stochastic SEIR R0")
        rates = seir_event_rates(state, params)
        call assert_close(rates%infection, 0.597_real64, 1.0e-12_real64, "infection hazard")
        call assert_close(rates%progression, 0.6_real64, 1.0e-12_real64, "progression hazard")
        call assert_close(rates%recovery, 0.2_real64, 1.0e-12_real64, "recovery hazard")
    end subroutine test_rates_and_r0

    subroutine test_progression_recovery_chain()
        type(stochastic_seir_parameters) :: params
        type(stochastic_seir_state) :: initial
        type(stochastic_seir_state) :: states(4)
        integer :: events(4), n_states
        type(rng_state) :: rng
        logical :: truncated

        params%population = 10_int64
        params%beta = 0.0_real64
        params%sigma = 1.0_real64
        params%gamma = 1.0_real64
        initial%s = 9_int64
        initial%e = 1_int64
        initial%i = 0_int64
        initial%r = 0_int64
        initial%time = 0.0_real64

        call seed_rng(rng, 17_int64)
        call simulate_seir_gillespie(initial, params, rng, 100.0_real64, states, events, n_states, truncated)

        if (truncated) error stop "progression/recovery chain was unexpectedly truncated"
        if (n_states /= 3) error stop "one exposed individual should produce exactly two events when beta=0"
        if (events(1) /= event_none .or. events(2) /= event_progression .or. &
                events(3) /= event_recovery) then
            error stop "expected progression followed by recovery"
        end if
        if (states(3)%e /= 0_int64 .or. states(3)%i /= 0_int64 .or. states(3)%r /= 1_int64) then
            error stop "progression/recovery chain ended in the wrong state"
        end if
    end subroutine test_progression_recovery_chain

    subroutine test_absorbing_state()
        type(stochastic_seir_parameters) :: params
        type(stochastic_seir_state) :: initial
        type(stochastic_seir_state) :: states(2)
        integer :: events(2), n_states
        type(rng_state) :: rng
        logical :: truncated

        params%population = 100_int64
        initial%s = 80_int64
        initial%e = 0_int64
        initial%i = 0_int64
        initial%r = 20_int64
        call seed_rng(rng, 99_int64)
        call simulate_seir_gillespie(initial, params, rng, 50.0_real64, states, events, n_states, truncated)

        if (truncated) error stop "absorbing state cannot truncate"
        if (n_states /= 1) error stop "disease-free SEIR state should be absorbing"
        if (events(1) /= event_none) error stop "absorbing state should contain no event"
    end subroutine test_absorbing_state

    subroutine test_seeded_event_prefix()
        type(stochastic_seir_parameters) :: params
        type(stochastic_seir_state) :: initial
        type(stochastic_seir_state), allocatable :: states(:)
        integer, allocatable :: events(:)
        integer :: n_states, idx
        integer, parameter :: expected(8) = [event_none, event_recovery, event_infection, &
            event_recovery, event_progression, event_progression, event_infection, event_recovery]
        type(rng_state) :: rng
        logical :: truncated
        integer(int64) :: population

        params%population = 1000_int64
        params%beta = 0.30_real64
        params%sigma = 0.20_real64
        params%gamma = 0.10_real64
        initial%s = 995_int64
        initial%e = 3_int64
        initial%i = 2_int64
        initial%r = 0_int64
        population = params%population

        allocate(states(10001), events(10001))
        call seed_rng(rng, 12345_int64)
        call simulate_seir_gillespie(initial, params, rng, 365.0_real64, states, events, n_states, truncated)
        if (truncated) error stop "seeded trajectory buffer was unexpectedly exhausted"
        if (n_states < size(expected)) error stop "seeded trajectory was shorter than expected"
        if (any(events(:size(expected)) /= expected)) error stop "seeded event prefix changed"

        do idx = 1, n_states
            if (states(idx)%s + states(idx)%e + states(idx)%i + states(idx)%r /= population) then
                error stop "stochastic SEIR population was not conserved"
            end if
            if (min(states(idx)%s, states(idx)%e, states(idx)%i, states(idx)%r) < 0_int64) then
                error stop "stochastic SEIR produced a negative compartment"
            end if
            if (idx > 1) then
                if (states(idx)%time <= states(idx - 1)%time) then
                    error stop "stochastic SEIR event times must increase strictly"
                end if
            end if
        end do
    end subroutine test_seeded_event_prefix

    subroutine test_ensemble_reproducibility()
        type(stochastic_seir_parameters) :: params
        type(stochastic_seir_state) :: initial
        type(seir_ensemble_summary) :: first, second

        params%population = 500_int64
        params%beta = 0.25_real64
        params%sigma = 0.20_real64
        params%gamma = 0.10_real64
        initial%s = 495_int64
        initial%e = 3_int64
        initial%i = 2_int64
        initial%r = 0_int64

        call simulate_seir_ensemble(initial, params, 777_int64, 100, 180.0_real64, &
            2000, 50_int64, first)
        call simulate_seir_ensemble(initial, params, 777_int64, 100, 180.0_real64, &
            2000, 50_int64, second)

        if (first%major_outbreak_count /= second%major_outbreak_count) then
            error stop "ensemble major-outbreak count is not reproducible"
        end if
        call assert_close(first%mean_final_size, second%mean_final_size, 0.0_real64, &
            "reproducible ensemble mean final size")
        call assert_close(first%sd_final_size, second%sd_final_size, 0.0_real64, &
            "reproducible ensemble final-size standard deviation")
    end subroutine test_ensemble_reproducibility

    subroutine test_fadeout_probability()
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

        call simulate_seir_ensemble(initial, params, 20260908_int64, 1000, 365.0_real64, &
            4000, 100_int64, summary)

        if (summary%fadeout_probability < 0.25_real64 .or. summary%fadeout_probability > 0.45_real64) then
            error stop "fadeout probability is inconsistent with the R0=3 branching regime"
        end if
        call assert_close(summary%fadeout_probability + summary%major_outbreak_probability, &
            1.0_real64, 1.0e-12_real64, "fadeout/major-outbreak partition")
    end subroutine test_fadeout_probability

    subroutine test_deterministic_ensemble_agreement()
        integer, parameter :: n_deterministic = 401
        real(real64), parameter :: dt = 0.05_real64
        real(real64), parameter :: population = 5000.0_real64
        type(stochastic_seir_parameters) :: stochastic_params
        type(stochastic_seir_state) :: stochastic_initial
        type(seir_ensemble_summary) :: summary
        type(seir_parameters) :: deterministic_params
        type(seir_state) :: deterministic_initial
        type(seir_state) :: deterministic_states(n_deterministic)

        stochastic_params%population = 5000_int64
        stochastic_params%beta = 0.24_real64
        stochastic_params%sigma = 0.25_real64
        stochastic_params%gamma = 0.12_real64
        stochastic_initial%s = 4900_int64
        stochastic_initial%e = 50_int64
        stochastic_initial%i = 50_int64
        stochastic_initial%r = 0_int64

        deterministic_params%population = population
        deterministic_params%beta = stochastic_params%beta
        deterministic_params%sigma = stochastic_params%sigma
        deterministic_params%gamma = stochastic_params%gamma
        deterministic_initial%s = 4900.0_real64
        deterministic_initial%e = 50.0_real64
        deterministic_initial%i = 50.0_real64
        deterministic_initial%r = 0.0_real64

        call simulate_seir(deterministic_initial, deterministic_params, dt, deterministic_states)
        call simulate_seir_ensemble(stochastic_initial, stochastic_params, 123456_int64, 300, &
            20.0_real64, 12000, 1000_int64, summary)

        call assert_fraction_close(summary%mean_s_at_horizon, deterministic_states(n_deterministic)%s, &
            population, "ensemble mean S versus deterministic SEIR")
        call assert_fraction_close(summary%mean_e_at_horizon, deterministic_states(n_deterministic)%e, &
            population, "ensemble mean E versus deterministic SEIR")
        call assert_fraction_close(summary%mean_i_at_horizon, deterministic_states(n_deterministic)%i, &
            population, "ensemble mean I versus deterministic SEIR")
        call assert_fraction_close(summary%mean_r_at_horizon, deterministic_states(n_deterministic)%r, &
            population, "ensemble mean R versus deterministic SEIR")
    end subroutine test_deterministic_ensemble_agreement

    subroutine assert_fraction_close(actual, expected, population, message)
        real(real64), intent(in) :: actual, expected, population
        character(len=*), intent(in) :: message

        if (abs(actual - expected) / population > 0.01_real64) then
            print *, "Agreement failure: ", trim(message)
            print *, "  stochastic mean: ", actual
            print *, "  deterministic:   ", expected
            error stop 1
        end if
    end subroutine assert_fraction_close

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

end program test_stochastic_seir
