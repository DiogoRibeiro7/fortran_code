module stochastic_seir_m
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng, uniform_open
    implicit none
    private

    integer, parameter, public :: event_none = 0
    integer, parameter, public :: event_infection = 1
    integer, parameter, public :: event_progression = 2
    integer, parameter, public :: event_recovery = 3

    public :: stochastic_seir_parameters, stochastic_seir_state, stochastic_seir_rates
    public :: seir_ensemble_summary
    public :: basic_reproduction_number_stochastic_seir, seir_event_rates
    public :: gillespie_step_seir, simulate_seir_gillespie, simulate_seir_ensemble

    type :: stochastic_seir_parameters
        real(real64) :: beta = 0.30_real64
        real(real64) :: sigma = 0.20_real64
        real(real64) :: gamma = 0.10_real64
        integer(int64) :: population = 1_int64
    end type stochastic_seir_parameters

    type :: stochastic_seir_state
        integer(int64) :: s = 0_int64
        integer(int64) :: e = 0_int64
        integer(int64) :: i = 0_int64
        integer(int64) :: r = 0_int64
        real(real64) :: time = 0.0_real64
    end type stochastic_seir_state

    type :: stochastic_seir_rates
        real(real64) :: infection = 0.0_real64
        real(real64) :: progression = 0.0_real64
        real(real64) :: recovery = 0.0_real64
    end type stochastic_seir_rates

    type :: seir_ensemble_summary
        integer :: n_replicates = 0
        integer :: extinction_count = 0
        integer :: major_outbreak_count = 0
        integer :: fadeout_count = 0
        real(real64) :: extinction_probability = 0.0_real64
        real(real64) :: major_outbreak_probability = 0.0_real64
        real(real64) :: fadeout_probability = 0.0_real64
        real(real64) :: mean_final_size = 0.0_real64
        real(real64) :: sd_final_size = 0.0_real64
        real(real64) :: mean_peak_infectious = 0.0_real64
        real(real64) :: mean_s_at_horizon = 0.0_real64
        real(real64) :: mean_e_at_horizon = 0.0_real64
        real(real64) :: mean_i_at_horizon = 0.0_real64
        real(real64) :: mean_r_at_horizon = 0.0_real64
    end type seir_ensemble_summary

contains

    pure real(real64) function basic_reproduction_number_stochastic_seir(params) result(r0)
        type(stochastic_seir_parameters), intent(in) :: params

        if (params%gamma <= 0.0_real64) then
            r0 = huge(1.0_real64)
        else
            r0 = params%beta / params%gamma
        end if
    end function basic_reproduction_number_stochastic_seir

    pure function seir_event_rates(state, params) result(rates)
        type(stochastic_seir_state), intent(in) :: state
        type(stochastic_seir_parameters), intent(in) :: params
        type(stochastic_seir_rates) :: rates

        if (params%population <= 0_int64) return

        rates%infection = params%beta * real(state%s, real64) * real(state%i, real64) / &
            real(params%population, real64)
        rates%progression = params%sigma * real(state%e, real64)
        rates%recovery = params%gamma * real(state%i, real64)
    end function seir_event_rates

    subroutine gillespie_step_seir(state, params, rng, next_state, event_type, occurred)
        type(stochastic_seir_state), intent(in) :: state
        type(stochastic_seir_parameters), intent(in) :: params
        type(rng_state), intent(inout) :: rng
        type(stochastic_seir_state), intent(out) :: next_state
        integer, intent(out) :: event_type
        logical, intent(out) :: occurred
        type(stochastic_seir_rates) :: rates
        real(real64) :: total_rate, waiting_time, event_draw

        call validate_state_and_parameters(state, params)

        next_state = state
        event_type = event_none
        occurred = .false.

        rates = seir_event_rates(state, params)
        total_rate = rates%infection + rates%progression + rates%recovery
        if (state%e + state%i == 0_int64 .or. total_rate <= 0.0_real64) return

        waiting_time = -log(uniform_open(rng)) / total_rate
        event_draw = uniform_open(rng) * total_rate
        next_state%time = state%time + waiting_time

        if (event_draw < rates%infection) then
            if (state%s <= 0_int64) error stop "infection event selected without susceptible individuals"
            next_state%s = state%s - 1_int64
            next_state%e = state%e + 1_int64
            event_type = event_infection
        else if (event_draw < rates%infection + rates%progression) then
            if (state%e <= 0_int64) error stop "progression event selected without exposed individuals"
            next_state%e = state%e - 1_int64
            next_state%i = state%i + 1_int64
            event_type = event_progression
        else
            if (state%i <= 0_int64) error stop "recovery event selected without infectious individuals"
            next_state%i = state%i - 1_int64
            next_state%r = state%r + 1_int64
            event_type = event_recovery
        end if

        occurred = .true.
    end subroutine gillespie_step_seir

    subroutine simulate_seir_gillespie(initial_state, params, rng, max_time, states, event_types, n_states, truncated)
        type(stochastic_seir_state), intent(in) :: initial_state
        type(stochastic_seir_parameters), intent(in) :: params
        type(rng_state), intent(inout) :: rng
        real(real64), intent(in) :: max_time
        type(stochastic_seir_state), intent(out) :: states(:)
        integer, intent(out) :: event_types(:)
        integer, intent(out) :: n_states
        logical, intent(out) :: truncated
        type(stochastic_seir_state) :: candidate, current
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

        do while (current%e + current%i > 0_int64 .and. current%time < max_time)
            call gillespie_step_seir(current, params, rng, candidate, event_type, occurred)
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
    end subroutine simulate_seir_gillespie

    subroutine simulate_seir_ensemble(initial_state, params, base_seed, n_replicates, max_time, &
            max_events, major_outbreak_threshold, summary)
        type(stochastic_seir_state), intent(in) :: initial_state
        type(stochastic_seir_parameters), intent(in) :: params
        integer(int64), intent(in) :: base_seed
        integer, intent(in) :: n_replicates, max_events
        real(real64), intent(in) :: max_time
        integer(int64), intent(in) :: major_outbreak_threshold
        type(seir_ensemble_summary), intent(out) :: summary
        type(stochastic_seir_state), allocatable :: states(:)
        integer, allocatable :: event_types(:)
        type(rng_state) :: rng
        integer(int64) :: final_size, peak_i
        integer :: replicate, n_states
        logical :: truncated
        real(real64) :: value, delta, mean_final, m2_final

        if (n_replicates <= 0) error stop "n_replicates must be positive"
        if (max_events < 1) error stop "max_events must be positive"
        if (major_outbreak_threshold < 0_int64 .or. major_outbreak_threshold > params%population) then
            error stop "major_outbreak_threshold must lie within the population"
        end if
        call validate_state_and_parameters(initial_state, params)

        allocate(states(max_events + 1), event_types(max_events + 1))
        summary = seir_ensemble_summary()
        summary%n_replicates = n_replicates
        mean_final = 0.0_real64
        m2_final = 0.0_real64

        do replicate = 1, n_replicates
            call seed_rng(rng, base_seed + 104729_int64 * int(replicate - 1, int64))
            call simulate_seir_gillespie(initial_state, params, rng, max_time, states, event_types, n_states, truncated)
            if (truncated) error stop "ensemble event buffer exhausted; increase max_events"

            final_size = params%population - states(n_states)%s
            peak_i = maxval(states(:n_states)%i)
            if (states(n_states)%e + states(n_states)%i == 0_int64) then
                summary%extinction_count = summary%extinction_count + 1
            end if
            if (final_size >= major_outbreak_threshold) then
                summary%major_outbreak_count = summary%major_outbreak_count + 1
            else
                summary%fadeout_count = summary%fadeout_count + 1
            end if

            value = real(final_size, real64)
            delta = value - mean_final
            mean_final = mean_final + delta / real(replicate, real64)
            m2_final = m2_final + delta * (value - mean_final)

            summary%mean_peak_infectious = summary%mean_peak_infectious + real(peak_i, real64)
            summary%mean_s_at_horizon = summary%mean_s_at_horizon + real(states(n_states)%s, real64)
            summary%mean_e_at_horizon = summary%mean_e_at_horizon + real(states(n_states)%e, real64)
            summary%mean_i_at_horizon = summary%mean_i_at_horizon + real(states(n_states)%i, real64)
            summary%mean_r_at_horizon = summary%mean_r_at_horizon + real(states(n_states)%r, real64)
        end do

        summary%extinction_probability = real(summary%extinction_count, real64) / real(n_replicates, real64)
        summary%major_outbreak_probability = real(summary%major_outbreak_count, real64) / real(n_replicates, real64)
        summary%fadeout_probability = real(summary%fadeout_count, real64) / real(n_replicates, real64)
        summary%mean_final_size = mean_final
        if (n_replicates > 1) summary%sd_final_size = sqrt(m2_final / real(n_replicates - 1, real64))
        summary%mean_peak_infectious = summary%mean_peak_infectious / real(n_replicates, real64)
        summary%mean_s_at_horizon = summary%mean_s_at_horizon / real(n_replicates, real64)
        summary%mean_e_at_horizon = summary%mean_e_at_horizon / real(n_replicates, real64)
        summary%mean_i_at_horizon = summary%mean_i_at_horizon / real(n_replicates, real64)
        summary%mean_r_at_horizon = summary%mean_r_at_horizon / real(n_replicates, real64)
    end subroutine simulate_seir_ensemble

    subroutine validate_state_and_parameters(state, params)
        type(stochastic_seir_state), intent(in) :: state
        type(stochastic_seir_parameters), intent(in) :: params
        integer(int64) :: total_population

        if (params%population <= 0_int64) error stop "population must be positive"
        if (params%beta < 0.0_real64) error stop "beta must be non-negative"
        if (params%sigma < 0.0_real64) error stop "sigma must be non-negative"
        if (params%gamma < 0.0_real64) error stop "gamma must be non-negative"
        if (state%s < 0_int64 .or. state%e < 0_int64 .or. state%i < 0_int64 .or. state%r < 0_int64) then
            error stop "compartment counts must be non-negative"
        end if

        total_population = state%s + state%e + state%i + state%r
        if (total_population /= params%population) error stop "state population does not match parameters"
    end subroutine validate_state_and_parameters

end module stochastic_seir_m
