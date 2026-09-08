program seir_intervention_example
    use iso_fortran_env, only: real64
    use seir_m, only: seir_parameters, seir_state, simulate_seir_scheduled
    use transmission_m, only: transmission_schedule, beta_at
    implicit none

    integer, parameter :: n_steps = 241
    real(real64), parameter :: dt = 0.5_real64
    type(seir_parameters) :: params
    type(seir_state) :: initial_state
    type(seir_state) :: states(n_steps)
    type(transmission_schedule) :: schedule
    real(real64) :: time
    integer :: idx

    params = seir_parameters(beta=0.30_real64, sigma=0.20_real64, &
        gamma=0.10_real64, population=100000.0_real64)
    initial_state = seir_state(s=99990.0_real64, e=0.0_real64, &
        i=10.0_real64, r=0.0_real64)

    schedule%baseline_beta = params%beta
    schedule%change_times = [20.0_real64, 60.0_real64]
    schedule%multipliers = [0.35_real64, 0.70_real64]

    call simulate_seir_scheduled(initial_state, params, schedule, 0.0_real64, dt, states)

    print '(a)', 'day,beta,susceptible,exposed,infectious,recovered'
    do idx = 1, n_steps
        time = real(idx - 1, real64) * dt
        print '(f7.2,a,f8.5,a,f12.3,a,f12.3,a,f12.3,a,f12.3)', &
            time, ',', beta_at(schedule, time), ',', states(idx)%s, ',', &
            states(idx)%e, ',', states(idx)%i, ',', states(idx)%r
    end do
end program seir_intervention_example
