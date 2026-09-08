module seir_m
    use iso_fortran_env, only: real64
    use transmission_m, only: transmission_schedule, beta_at, schedule_is_valid
    implicit none
    private

    public :: seir_parameters, seir_state
    public :: basic_reproduction_number_seir, seir_incidence_rate
    public :: seir_derivative, rk4_step_seir, simulate_seir
    public :: rk4_step_seir_scheduled, simulate_seir_scheduled

    type :: seir_parameters
        real(real64) :: beta = 0.30_real64
        real(real64) :: sigma = 0.20_real64
        real(real64) :: gamma = 0.10_real64
        real(real64) :: population = 1.0_real64
    end type seir_parameters

    type :: seir_state
        real(real64) :: s = 0.0_real64
        real(real64) :: e = 0.0_real64
        real(real64) :: i = 0.0_real64
        real(real64) :: r = 0.0_real64
    end type seir_state

contains

    pure real(real64) function basic_reproduction_number_seir(params) result(r0)
        type(seir_parameters), intent(in) :: params

        if (params%gamma <= 0.0_real64) then
            r0 = huge(1.0_real64)
        else
            r0 = params%beta / params%gamma
        end if
    end function basic_reproduction_number_seir

    pure real(real64) function seir_incidence_rate(state, params) result(incidence)
        type(seir_state), intent(in) :: state
        type(seir_parameters), intent(in) :: params

        if (params%population <= 0.0_real64) then
            incidence = 0.0_real64
        else
            incidence = params%beta * state%s * state%i / params%population
        end if
    end function seir_incidence_rate

    pure function seir_derivative(state, params) result(dstatedt)
        type(seir_state), intent(in) :: state
        type(seir_parameters), intent(in) :: params
        type(seir_state) :: dstatedt
        real(real64) :: incidence

        incidence = seir_incidence_rate(state, params)
        dstatedt%s = -incidence
        dstatedt%e = incidence - params%sigma * state%e
        dstatedt%i = params%sigma * state%e - params%gamma * state%i
        dstatedt%r = params%gamma * state%i
    end function seir_derivative

    pure function rk4_step_seir(state, params, dt) result(next_state)
        type(seir_state), intent(in) :: state
        type(seir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(seir_state) :: next_state
        type(seir_state) :: k1, k2, k3, k4, tmp

        k1 = seir_derivative(state, params)
        tmp = add_scaled(state, k1, 0.5_real64 * dt)
        k2 = seir_derivative(tmp, params)
        tmp = add_scaled(state, k2, 0.5_real64 * dt)
        k3 = seir_derivative(tmp, params)
        tmp = add_scaled(state, k3, dt)
        k4 = seir_derivative(tmp, params)

        next_state = combine_rk4(state, k1, k2, k3, k4, dt)
    end function rk4_step_seir

    pure function rk4_step_seir_scheduled(state, params, schedule, time, dt) result(next_state)
        type(seir_state), intent(in) :: state
        type(seir_parameters), intent(in) :: params
        type(transmission_schedule), intent(in) :: schedule
        real(real64), intent(in) :: time, dt
        type(seir_state) :: next_state
        type(seir_state) :: k1, k2, k3, k4, tmp
        type(seir_parameters) :: local_params

        local_params = params
        local_params%beta = beta_at(schedule, time)
        k1 = seir_derivative(state, local_params)

        local_params%beta = beta_at(schedule, time + 0.5_real64 * dt)
        tmp = add_scaled(state, k1, 0.5_real64 * dt)
        k2 = seir_derivative(tmp, local_params)
        tmp = add_scaled(state, k2, 0.5_real64 * dt)
        k3 = seir_derivative(tmp, local_params)

        local_params%beta = beta_at(schedule, time + dt)
        tmp = add_scaled(state, k3, dt)
        k4 = seir_derivative(tmp, local_params)

        next_state = combine_rk4(state, k1, k2, k3, k4, dt)
    end function rk4_step_seir_scheduled

    subroutine simulate_seir(initial_state, params, dt, states)
        type(seir_state), intent(in) :: initial_state
        type(seir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(seir_state), intent(out) :: states(:)
        integer :: idx

        call validate_simulation(params, dt, size(states))
        states(1) = initial_state
        do idx = 2, size(states)
            states(idx) = rk4_step_seir(states(idx - 1), params, dt)
        end do
    end subroutine simulate_seir

    subroutine simulate_seir_scheduled(initial_state, params, schedule, start_time, dt, states)
        type(seir_state), intent(in) :: initial_state
        type(seir_parameters), intent(in) :: params
        type(transmission_schedule), intent(in) :: schedule
        real(real64), intent(in) :: start_time, dt
        type(seir_state), intent(out) :: states(:)
        integer :: idx
        real(real64) :: time

        call validate_simulation(params, dt, size(states))
        if (.not. schedule_is_valid(schedule)) error stop "invalid transmission schedule"

        states(1) = initial_state
        time = start_time
        do idx = 2, size(states)
            states(idx) = rk4_step_seir_scheduled(states(idx - 1), params, schedule, time, dt)
            time = time + dt
        end do
    end subroutine simulate_seir_scheduled

    subroutine validate_simulation(params, dt, n_states)
        type(seir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        integer, intent(in) :: n_states

        if (n_states < 1) error stop "states must contain at least one element"
        if (dt <= 0.0_real64) error stop "dt must be positive"
        if (params%population <= 0.0_real64) error stop "population must be positive"
        if (params%beta < 0.0_real64) error stop "beta must be non-negative"
        if (params%sigma <= 0.0_real64) error stop "sigma must be positive"
        if (params%gamma <= 0.0_real64) error stop "gamma must be positive"
    end subroutine validate_simulation

    pure function add_scaled(state, derivative, scale) result(updated)
        type(seir_state), intent(in) :: state, derivative
        real(real64), intent(in) :: scale
        type(seir_state) :: updated

        updated%s = state%s + scale * derivative%s
        updated%e = state%e + scale * derivative%e
        updated%i = state%i + scale * derivative%i
        updated%r = state%r + scale * derivative%r
    end function add_scaled

    pure function combine_rk4(state, k1, k2, k3, k4, dt) result(updated)
        type(seir_state), intent(in) :: state, k1, k2, k3, k4
        real(real64), intent(in) :: dt
        type(seir_state) :: updated

        updated%s = state%s + dt * (k1%s + 2.0_real64*k2%s + 2.0_real64*k3%s + k4%s) / 6.0_real64
        updated%e = state%e + dt * (k1%e + 2.0_real64*k2%e + 2.0_real64*k3%e + k4%e) / 6.0_real64
        updated%i = state%i + dt * (k1%i + 2.0_real64*k2%i + 2.0_real64*k3%i + k4%i) / 6.0_real64
        updated%r = state%r + dt * (k1%r + 2.0_real64*k2%r + 2.0_real64*k3%r + k4%r) / 6.0_real64
    end function combine_rk4

end module seir_m
