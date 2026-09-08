program test_sir
    use iso_fortran_env, only: real64
    use epidemiology_m, only: sir_parameters, sir_state, simulate_sir, basic_reproduction_number
    implicit none

    integer, parameter :: n_steps = 401
    real(real64), parameter :: dt = 0.25_real64
    real(real64), parameter :: tol = 1.0e-8_real64
    type(sir_parameters) :: params
    type(sir_state) :: initial_state
    type(sir_state) :: states(n_steps)
    real(real64) :: total_population
    integer :: idx

    params = sir_parameters(beta=0.30_real64, gamma=0.10_real64, population=1000.0_real64)
    initial_state = sir_state(s=999.0_real64, i=1.0_real64, r=0.0_real64)

    if (abs(basic_reproduction_number(params) - 3.0_real64) > tol) then
        error stop "incorrect basic reproduction number"
    end if

    call simulate_sir(initial_state, params, dt, states)

    do idx = 1, n_steps
        total_population = states(idx)%s + states(idx)%i + states(idx)%r
        if (abs(total_population - params%population) > 1.0e-6_real64) then
            error stop "SIR simulation does not conserve population"
        end if
        if (min(states(idx)%s, states(idx)%i, states(idx)%r) < -1.0e-9_real64) then
            error stop "SIR compartment became negative"
        end if
    end do

    if (states(2)%i <= states(1)%i) then
        error stop "infectious population should initially grow when R0 > 1"
    end if

    if (states(n_steps)%s >= states(1)%s) then
        error stop "susceptible population should decline during the epidemic"
    end if
end program test_sir
