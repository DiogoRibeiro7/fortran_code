program stochastic_sir_example
    use iso_fortran_env, only: int64, real64
    use rng_m, only: rng_state, seed_rng
    use stochastic_sir_m, only: stochastic_sir_parameters, stochastic_sir_state, &
        basic_reproduction_number_stochastic_sir, event_none, event_infection, event_recovery, &
        simulate_sir_gillespie
    implicit none

    integer, parameter :: capacity = 10000
    type(stochastic_sir_parameters) :: params
    type(stochastic_sir_state) :: initial, states(capacity)
    type(rng_state) :: rng
    integer :: events(capacity), n_states, idx
    logical :: truncated

    params%population = 1000_int64
    params%beta = 0.30_real64
    params%gamma = 0.10_real64

    initial%s = 995_int64
    initial%i = 5_int64
    initial%r = 0_int64
    initial%time = 0.0_real64

    call seed_rng(rng, 20260908_int64)
    call simulate_sir_gillespie(initial, params, rng, 160.0_real64, states, events, n_states, truncated)

    print '(a,f8.4)', '# R0 = ', basic_reproduction_number_stochastic_sir(params)
    print '(a)', 'event,time,susceptible,infectious,recovered'
    do idx = 1, n_states
        print '(a,",",f12.6,",",i0,",",i0,",",i0)', &
            event_label(events(idx)), states(idx)%time, states(idx)%s, states(idx)%i, states(idx)%r
    end do

    if (truncated) error stop "event buffer exhausted before simulation completed"

contains

    pure function event_label(event_type) result(label)
        integer, intent(in) :: event_type
        character(len=9) :: label

        select case (event_type)
        case (event_none)
            label = 'initial'
        case (event_infection)
            label = 'infection'
        case (event_recovery)
            label = 'recovery'
        case default
            label = 'unknown'
        end select
    end function event_label

end program stochastic_sir_example
