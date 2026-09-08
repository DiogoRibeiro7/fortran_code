program test_stochastic_sir
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng, uniform_open
    use stochastic_sir_m, only: stochastic_sir_parameters, stochastic_sir_state, stochastic_sir_rates, &
        event_none, event_infection, event_recovery, basic_reproduction_number_stochastic_sir, &
        sir_event_rates, simulate_sir_gillespie
    implicit none

    call test_rng_reproducibility()
    call test_event_rates()
    call test_recovery_only_path()
    call test_seeded_trajectory()
    call test_absorbing_state()

    print *, "Stochastic SIR tests passed"

contains

    subroutine test_rng_reproducibility()
        type(rng_state) :: first, second
        real(real64) :: x1, x2
        integer :: idx

        call seed_rng(first, 20260908_int64)
        call seed_rng(second, 20260908_int64)
        do idx = 1, 20
            x1 = uniform_open(first)
            x2 = uniform_open(second)
            if (transfer(x1, 0_int64) /= transfer(x2, 0_int64)) then
                error stop "equal RNG seeds must reproduce the same stream"
            end if
            if (x1 <= 0.0_real64 .or. x1 >= 1.0_real64) error stop "RNG draw must lie in (0,1)"
        end do
    end subroutine test_rng_reproducibility

    subroutine test_event_rates()
        type(stochastic_sir_parameters) :: params
        type(stochastic_sir_state) :: state
        type(stochastic_sir_rates) :: rates

        params%population = 1000_int64
        params%beta = 0.30_real64
        params%gamma = 0.10_real64
        state%s = 990_int64
        state%i = 10_int64
        state%r = 0_int64

        call assert_close(basic_reproduction_number_stochastic_sir(params), 3.0_real64, &
            1.0e-12_real64, "R0")
        rates = sir_event_rates(state, params)
        call assert_close(rates%infection, 2.97_real64, 1.0e-12_real64, "infection rate")
        call assert_close(rates%recovery, 1.0_real64, 1.0e-12_real64, "recovery rate")
    end subroutine test_event_rates

    subroutine test_recovery_only_path()
        integer, parameter :: capacity = 16
        type(stochastic_sir_parameters) :: params
        type(stochastic_sir_state) :: initial, states(capacity)
        type(rng_state) :: rng
        integer :: events(capacity), n_states, idx
        logical :: truncated

        params%population = 100_int64
        params%beta = 0.0_real64
        params%gamma = 0.25_real64
        initial%s = 95_int64
        initial%i = 5_int64
        initial%r = 0_int64
        initial%time = 0.0_real64

        call seed_rng(rng, 7_int64)
        call simulate_sir_gillespie(initial, params, rng, 1000.0_real64, states, events, n_states, truncated)

        if (truncated) error stop "recovery-only trajectory should not truncate"
        if (n_states /= 6) error stop "five infectious individuals require five recovery events"
        do idx = 2, n_states
            if (events(idx) /= event_recovery) error stop "beta=0 trajectory must contain only recoveries"
            if (states(idx)%s /= initial%s) error stop "recoveries must not change susceptible count"
            if (states(idx)%time <= states(idx - 1)%time) error stop "event times must increase"
            call assert_population(states(idx), params%population)
        end do
        if (states(n_states)%i /= 0_int64 .or. states(n_states)%r /= 5_int64) then
            error stop "recovery-only trajectory must end in extinction"
        end if
    end subroutine test_recovery_only_path

    subroutine test_seeded_trajectory()
        integer, parameter :: capacity = 128
        type(stochastic_sir_parameters) :: params
        type(stochastic_sir_state) :: initial, states(capacity)
        type(rng_state) :: rng
        integer :: events(capacity), n_states, idx
        logical :: truncated
        integer, parameter :: expected_events(8) = [ &
            event_infection, event_infection, event_infection, event_infection, &
            event_infection, event_infection, event_infection, event_infection ]

        params%population = 200_int64
        params%beta = 0.35_real64
        params%gamma = 0.10_real64
        initial%s = 197_int64
        initial%i = 3_int64
        initial%r = 0_int64
        initial%time = 0.0_real64

        call seed_rng(rng, 123456_int64)
        call simulate_sir_gillespie(initial, params, rng, 20.0_real64, states, events, n_states, truncated)

        if (n_states < 9) error stop "seeded trajectory ended before regression prefix"
        do idx = 1, 8
            if (events(idx + 1) /= expected_events(idx)) error stop "seeded event regression changed"
        end do
        do idx = 1, n_states
            call assert_population(states(idx), params%population)
            if (min(states(idx)%s, states(idx)%i, states(idx)%r) < 0_int64) then
                error stop "stochastic trajectory produced a negative compartment"
            end if
        end do
    end subroutine test_seeded_trajectory

    subroutine test_absorbing_state()
        integer, parameter :: capacity = 4
        type(stochastic_sir_parameters) :: params
        type(stochastic_sir_state) :: initial, states(capacity)
        type(rng_state) :: rng
        integer :: events(capacity), n_states
        logical :: truncated

        params%population = 100_int64
        initial%s = 80_int64
        initial%i = 0_int64
        initial%r = 20_int64

        call seed_rng(rng, 1_int64)
        call simulate_sir_gillespie(initial, params, rng, 50.0_real64, states, events, n_states, truncated)
        if (n_states /= 1) error stop "disease-free state should be absorbing"
        if (events(1) /= event_none) error stop "initial state must not be labelled as an event"
        if (truncated) error stop "absorbing trajectory cannot be truncated"
    end subroutine test_absorbing_state

    subroutine assert_population(state, expected)
        type(stochastic_sir_state), intent(in) :: state
        integer(int64), intent(in) :: expected

        if (state%s + state%i + state%r /= expected) error stop "population not conserved"
    end subroutine assert_population

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_stochastic_sir
