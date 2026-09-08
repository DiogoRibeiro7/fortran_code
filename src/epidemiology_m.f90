module epidemiology_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: sir_parameters, sir_state, basic_reproduction_number
    public :: sir_derivative, rk4_step_sir, simulate_sir, incidence_rate

    type :: sir_parameters
        real(real64) :: beta = 0.30_real64
        real(real64) :: gamma = 0.10_real64
        real(real64) :: population = 1.0_real64
    end type sir_parameters

    type :: sir_state
        real(real64) :: s = 0.0_real64
        real(real64) :: i = 0.0_real64
        real(real64) :: r = 0.0_real64
    end type sir_state

contains

    pure real(real64) function basic_reproduction_number(params) result(r0)
        type(sir_parameters), intent(in) :: params

        if (params%gamma <= 0.0_real64) then
            r0 = huge(1.0_real64)
        else
            r0 = params%beta / params%gamma
        end if
    end function basic_reproduction_number

    pure real(real64) function incidence_rate(state, params) result(incidence)
        type(sir_state), intent(in) :: state
        type(sir_parameters), intent(in) :: params

        if (params%population <= 0.0_real64) then
            incidence = 0.0_real64
        else
            incidence = params%beta * state%s * state%i / params%population
        end if
    end function incidence_rate

    pure function sir_derivative(state, params) result(dstatedt)
        type(sir_state), intent(in) :: state
        type(sir_parameters), intent(in) :: params
        type(sir_state) :: dstatedt
        real(real64) :: incidence

        incidence = incidence_rate(state, params)
        dstatedt%s = -incidence
        dstatedt%i = incidence - params%gamma * state%i
        dstatedt%r = params%gamma * state%i
    end function sir_derivative

    pure function rk4_step_sir(state, params, dt) result(next_state)
        type(sir_state), intent(in) :: state
        type(sir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(sir_state) :: next_state
        type(sir_state) :: k1, k2, k3, k4, tmp

        k1 = sir_derivative(state, params)
        tmp = add_scaled(state, k1, 0.5_real64 * dt)
        k2 = sir_derivative(tmp, params)
        tmp = add_scaled(state, k2, 0.5_real64 * dt)
        k3 = sir_derivative(tmp, params)
        tmp = add_scaled(state, k3, dt)
        k4 = sir_derivative(tmp, params)

        next_state%s = state%s + dt * (k1%s + 2.0_real64*k2%s + 2.0_real64*k3%s + k4%s) / 6.0_real64
        next_state%i = state%i + dt * (k1%i + 2.0_real64*k2%i + 2.0_real64*k3%i + k4%i) / 6.0_real64
        next_state%r = state%r + dt * (k1%r + 2.0_real64*k2%r + 2.0_real64*k3%r + k4%r) / 6.0_real64
    end function rk4_step_sir

    subroutine simulate_sir(initial_state, params, dt, states)
        type(sir_state), intent(in) :: initial_state
        type(sir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(sir_state), intent(out) :: states(:)
        integer :: idx

        if (size(states) < 1) error stop "states must contain at least one element"
        if (dt <= 0.0_real64) error stop "dt must be positive"
        if (params%population <= 0.0_real64) error stop "population must be positive"

        states(1) = initial_state
        do idx = 2, size(states)
            states(idx) = rk4_step_sir(states(idx - 1), params, dt)
        end do
    end subroutine simulate_sir

    pure function add_scaled(state, derivative, scale) result(updated)
        type(sir_state), intent(in) :: state, derivative
        real(real64), intent(in) :: scale
        type(sir_state) :: updated

        updated%s = state%s + scale * derivative%s
        updated%i = state%i + scale * derivative%i
        updated%r = state%r + scale * derivative%r
    end function add_scaled

end module epidemiology_m
