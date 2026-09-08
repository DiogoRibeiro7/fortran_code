module seirs_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: seirs_parameters, seirs_state
    public :: basic_reproduction_number_seirs, disease_free_susceptible_fraction
    public :: control_reproduction_number, seirs_incidence_rate
    public :: seirs_derivative, rk4_step_seirs, simulate_seirs

    type :: seirs_parameters
        real(real64) :: beta = 0.30_real64
        real(real64) :: sigma = 0.20_real64
        real(real64) :: gamma = 0.10_real64
        real(real64) :: vaccination_rate = 0.0_real64
        real(real64) :: waning_rate = 0.0_real64
        real(real64) :: mortality_rate = 0.0_real64
        real(real64) :: population = 1.0_real64
    end type seirs_parameters

    type :: seirs_state
        real(real64) :: s = 0.0_real64
        real(real64) :: e = 0.0_real64
        real(real64) :: i = 0.0_real64
        real(real64) :: r = 0.0_real64
    end type seirs_state

contains

    pure real(real64) function basic_reproduction_number_seirs(params) result(r0)
        type(seirs_parameters), intent(in) :: params
        real(real64) :: latent_exit, infectious_exit

        latent_exit = params%sigma + params%mortality_rate
        infectious_exit = params%gamma + params%mortality_rate

        if (latent_exit <= 0.0_real64 .or. infectious_exit <= 0.0_real64) then
            r0 = huge(1.0_real64)
        else
            r0 = params%beta * params%sigma / (latent_exit * infectious_exit)
        end if
    end function basic_reproduction_number_seirs

    pure real(real64) function disease_free_susceptible_fraction(params) result(fraction)
        type(seirs_parameters), intent(in) :: params
        real(real64) :: immune_exit, denominator

        immune_exit = params%waning_rate + params%mortality_rate
        denominator = immune_exit + params%vaccination_rate

        if (denominator <= 0.0_real64) then
            fraction = 1.0_real64
        else
            fraction = immune_exit / denominator
        end if
    end function disease_free_susceptible_fraction

    pure real(real64) function control_reproduction_number(params) result(rc)
        type(seirs_parameters), intent(in) :: params

        rc = basic_reproduction_number_seirs(params) * disease_free_susceptible_fraction(params)
    end function control_reproduction_number

    pure real(real64) function seirs_incidence_rate(state, params) result(incidence)
        type(seirs_state), intent(in) :: state
        type(seirs_parameters), intent(in) :: params

        if (params%population <= 0.0_real64) then
            incidence = 0.0_real64
        else
            incidence = params%beta * state%s * state%i / params%population
        end if
    end function seirs_incidence_rate

    pure function seirs_derivative(state, params) result(dstatedt)
        type(seirs_state), intent(in) :: state
        type(seirs_parameters), intent(in) :: params
        type(seirs_state) :: dstatedt
        real(real64) :: incidence, births

        incidence = seirs_incidence_rate(state, params)
        births = params%mortality_rate * params%population

        dstatedt%s = births - incidence - params%vaccination_rate * state%s + &
            params%waning_rate * state%r - params%mortality_rate * state%s
        dstatedt%e = incidence - (params%sigma + params%mortality_rate) * state%e
        dstatedt%i = params%sigma * state%e - (params%gamma + params%mortality_rate) * state%i
        dstatedt%r = params%gamma * state%i + params%vaccination_rate * state%s - &
            (params%waning_rate + params%mortality_rate) * state%r
    end function seirs_derivative

    pure function rk4_step_seirs(state, params, dt) result(next_state)
        type(seirs_state), intent(in) :: state
        type(seirs_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(seirs_state) :: next_state
        type(seirs_state) :: k1, k2, k3, k4, tmp

        k1 = seirs_derivative(state, params)
        tmp = add_scaled(state, k1, 0.5_real64 * dt)
        k2 = seirs_derivative(tmp, params)
        tmp = add_scaled(state, k2, 0.5_real64 * dt)
        k3 = seirs_derivative(tmp, params)
        tmp = add_scaled(state, k3, dt)
        k4 = seirs_derivative(tmp, params)

        next_state = combine_rk4(state, k1, k2, k3, k4, dt)
    end function rk4_step_seirs

    subroutine simulate_seirs(initial_state, params, dt, states)
        type(seirs_state), intent(in) :: initial_state
        type(seirs_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(seirs_state), intent(out) :: states(:)
        integer :: idx

        call validate_simulation(params, dt, size(states))
        states(1) = initial_state
        do idx = 2, size(states)
            states(idx) = rk4_step_seirs(states(idx - 1), params, dt)
        end do
    end subroutine simulate_seirs

    subroutine validate_simulation(params, dt, n_states)
        type(seirs_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        integer, intent(in) :: n_states

        if (n_states < 1) error stop "states must contain at least one element"
        if (dt <= 0.0_real64) error stop "dt must be positive"
        if (params%population <= 0.0_real64) error stop "population must be positive"
        if (params%beta < 0.0_real64) error stop "beta must be non-negative"
        if (params%sigma <= 0.0_real64) error stop "sigma must be positive"
        if (params%gamma <= 0.0_real64) error stop "gamma must be positive"
        if (params%vaccination_rate < 0.0_real64) error stop "vaccination_rate must be non-negative"
        if (params%waning_rate < 0.0_real64) error stop "waning_rate must be non-negative"
        if (params%mortality_rate < 0.0_real64) error stop "mortality_rate must be non-negative"
    end subroutine validate_simulation

    pure function add_scaled(state, derivative, scale) result(updated)
        type(seirs_state), intent(in) :: state, derivative
        real(real64), intent(in) :: scale
        type(seirs_state) :: updated

        updated%s = state%s + scale * derivative%s
        updated%e = state%e + scale * derivative%e
        updated%i = state%i + scale * derivative%i
        updated%r = state%r + scale * derivative%r
    end function add_scaled

    pure function combine_rk4(state, k1, k2, k3, k4, dt) result(updated)
        type(seirs_state), intent(in) :: state, k1, k2, k3, k4
        real(real64), intent(in) :: dt
        type(seirs_state) :: updated

        updated%s = state%s + dt * (k1%s + 2.0_real64*k2%s + 2.0_real64*k3%s + k4%s) / 6.0_real64
        updated%e = state%e + dt * (k1%e + 2.0_real64*k2%e + 2.0_real64*k3%e + k4%e) / 6.0_real64
        updated%i = state%i + dt * (k1%i + 2.0_real64*k2%i + 2.0_real64*k3%i + k4%i) / 6.0_real64
        updated%r = state%r + dt * (k1%r + 2.0_real64*k2%r + 2.0_real64*k3%r + k4%r) / 6.0_real64
    end function combine_rk4

end module seirs_m
