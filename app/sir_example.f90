program sir_example
    use iso_fortran_env, only: real64
    use epidemiology_m, only: sir_parameters, sir_state, simulate_sir, basic_reproduction_number
    implicit none

    integer, parameter :: n_steps = 161
    real(real64), parameter :: dt = 0.5_real64
    type(sir_parameters) :: params
    type(sir_state) :: initial_state
    type(sir_state) :: states(n_steps)
    integer :: idx

    params = sir_parameters(beta=0.30_real64, gamma=0.10_real64, population=100000.0_real64)
    initial_state = sir_state(s=99990.0_real64, i=10.0_real64, r=0.0_real64)

    call simulate_sir(initial_state, params, dt, states)

    print '(a,f6.3)', '# R0 = ', basic_reproduction_number(params)
    print '(a)', 'day,susceptible,infectious,recovered'
    do idx = 1, n_steps
        print '(f8.2,a,f12.3,a,f12.3,a,f12.3)', &
            real(idx - 1, real64) * dt, ',', states(idx)%s, ',', states(idx)%i, ',', states(idx)%r
    end do
end program sir_example
