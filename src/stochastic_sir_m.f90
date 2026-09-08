module stochastic_sir_m
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, uniform_open
    implicit none
    private

    integer, parameter, public :: event_none = 0
    integer, parameter, public :: event_infection = 1
    integer, parameter, public :: event_recovery = 2

    public :: stochastic_sir_parameters, stochastic_sir_state, stochastic_sir_rates
    public :: basic_reproduction_number_stochastic_sir, sir_event_rates
    public :: gillespie_step_sir, simulate_sir_gillespie

    type :: stochastic_sir_parameters
        real(real64) :: beta = 0.30_real64
        real(real64) :: gamma = 0.10_real64
        integer(int64) :: population = 1_int64
    end type stochastic_sir_parameters

    type :: stochastic_sir_state
        integer(int64) :: s = 0_int64
        integer(int64) :: i = 0_int64
        integer(int64) :: r = 0_int64
        real(real64) :: time = 0.0_real64
    end type stochastic_sir_state

    type :: stochastic_sir_rates
        real(real64) :: infection = 0.0_real64
        real(real64) :: recovery = 0.0_real64
    end type stochastic_sir_rates

contains

    pure real(real64) function basic_reproduction_number_stochastic_sir(params) result(r0)
        type(stochastic_sir_parameters), intent(in) :: params

        if (params%gamma <= 0.0_real64) then
            r0 = huge(1.0_real64)
        else
            r0 = params%beta / params%gamma
        end if
    end function basic_reproduction_number_stochastic_sir

    pure function sir_event_rates(state, params) result(rates)
        type(stochastic_sir_state), intent(in) :: state
        type(stochastic_sir_parameters), intent(in) :: params
        type(stochastic_sir_rates) :: rates

        if (params%population <= 0_int64) return

        rates%infection = params%beta * real(state%s, real64) * real(state%i, real64) / &
            real(params%population, real64)
        rates%recovery = params%gamma * real(state%i, real64)
    end function sir_event_rates

    subroutine gillespie_step_sir(state, params, rng, next_state, event_type, occurred)
        type(stochastic_sir_state), intent(in) :: state
        type(stochastic_sir_parameters), intent(in) :: params
        type(rng_state), intent(inout) :: rng
        type(stochastic_sir_state), intent(out) :: next_state
        integer, intent(out) :: event_type
        logical, intent(out) :: occurred
        type(stochastic_sir_rates) :: rates
        real(real64) :: total_rate, waiting_time, event_draw

        call validate_state_and_parameters(state, params)

        next_state = state
        event_type = event_none
        occurred = .false.

        rates = sir_event_rates(state, params)
        total_rate = rates%infection + rates%recovery
        if (state%i == 0_int64 .or. total_rate <= 0.0_real64) return

        waiting_time = -log(uniform_open(rng)) / total_rate
        event_draw = uniform_open(rng) * total_rate
        next_state%time = state%time + waiting_time

        if (event_draw < rates%infection) then
            if (state%s <= 0_int64) error stop "infection event selected without susceptible individuals"
            next_state%s = state%s - 1_int64
            next_state%i = state%i + 1_int64
            event_type = event_infection
        else
            if (state%i <= 0_int64) error stop "recovery event selected without infectious individuals"
            next_state%i = state%i - 1_int64
            next_state%r = state%r + 1_int64
            event_type = event_recovery
        end if

        occurred = .true.
    end subroutine gillespie_step_sir

    subroutine simulate_sir_gillespie(initial_state, params, rng, max_time, states, event_types, n_states, truncated)
        type(stochastic_sir_state), intent(in) :: initial_state
        type(stochastic_sir_parameters), intent(in) :: params
        type(rng_state), intent(inout) :: rng
        real(real64), intent(in) :: max_time
        type(stochastic_sir_state), intent(out) :: states(:)
        integer, intent(out) :: event_types(:)
        integer, intent(out) :: n_states
        logical, intent(out) :: truncated
        type(stochastic_sir_state) :: candidate, current
        integer :: event_type
        logical :: occurred

        if (size(states) < 1) error stop "states must contain at least one element"
        if (size(event_types) /= size(states)) error stop "event_types must match states size"
        if (max_time < initial_state%time) error stop "max_time precedes initial state time"
        call validate_state_and_parameters(initial_state, params)

        states(1) = initial_state
        event_types(1) = event_none
        n_states = 1
        truncated = .false.
        current = initial_state

        do while (current%i > 0_int64 .and. current%time < max_time)
            call gillespie_step_sir(current, params, rng, candidate, event_type, occurred)
            if (.not. occurred) exit
            if (candidate%time > max_time) exit

            if (n_states == size(states)) then
                truncated = .true.
                exit
            end if

            n_states = n_states + 1
            states(n_states) = candidate
            event_types(n_states) = event_type
            current = candidate
        end do
    end subroutine simulate_sir_gillespie

    subroutine validate_state_and_parameters(state, params)
        type(stochastic_sir_state), intent(in) :: state
        type(stochastic_sir_parameters), intent(in) :: params
        integer(int64) :: total_population

        if (params%population <= 0_int64) error stop "population must be positive"
        if (params%beta < 0.0_real64) error stop "beta must be non-negative"
        if (params%gamma < 0.0_real64) error stop "gamma must be non-negative"
        if (state%s < 0_int64 .or. state%i < 0_int64 .or. state%r < 0_int64) then
            error stop "compartment counts must be non-negative"
        end if

        total_population = state%s + state%i + state%r
        if (total_population /= params%population) error stop "state population does not match parameters"
    end subroutine validate_state_and_parameters

end module stochastic_sir_m
