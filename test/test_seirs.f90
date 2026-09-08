program test_seirs
    use iso_fortran_env, only: real64
    use seirs_m, only: seirs_parameters, seirs_state, basic_reproduction_number_seirs, &
        disease_free_susceptible_fraction, control_reproduction_number, seirs_derivative, simulate_seirs
    implicit none

    integer, parameter :: n_steps = 3651
    real(real64), parameter :: dt = 1.0_real64
    real(real64), parameter :: tolerance = 1.0e-6_real64
    type(seirs_parameters) :: base, vaccinated
    type(seirs_state) :: initial, deriv
    type(seirs_state), allocatable :: baseline_states(:), vaccinated_states(:)
    real(real64) :: r0, expected_r0, susceptible_fraction, expected_fraction
    real(real64) :: total, max_error
    integer :: idx

    base%population = 100000.0_real64
    base%beta = 0.30_real64
    base%sigma = 0.20_real64
    base%gamma = 0.10_real64
    base%mortality_rate = 1.0_real64 / (75.0_real64 * 365.0_real64)
    base%waning_rate = 1.0_real64 / (3.0_real64 * 365.0_real64)
    base%vaccination_rate = 0.0_real64

    expected_r0 = base%beta * base%sigma / &
        ((base%sigma + base%mortality_rate) * (base%gamma + base%mortality_rate))
    r0 = basic_reproduction_number_seirs(base)
    call assert_close(r0, expected_r0, 1.0e-12_real64, "R0 formula")

    vaccinated = base
    vaccinated%vaccination_rate = 0.0025_real64
    expected_fraction = (vaccinated%waning_rate + vaccinated%mortality_rate) / &
        (vaccinated%waning_rate + vaccinated%mortality_rate + vaccinated%vaccination_rate)
    susceptible_fraction = disease_free_susceptible_fraction(vaccinated)
    call assert_close(susceptible_fraction, expected_fraction, 1.0e-12_real64, &
        "disease-free susceptible fraction")
    if (control_reproduction_number(vaccinated) >= 1.0_real64) then
        error stop "vaccination scenario should have control reproduction number below one"
    end if
    if (control_reproduction_number(vaccinated) >= r0) then
        error stop "vaccination should reduce the reproduction number"
    end if

    initial%s = base%population - 30.0_real64
    initial%e = 20.0_real64
    initial%i = 10.0_real64
    initial%r = 0.0_real64

    deriv = seirs_derivative(initial, base)
    call assert_close(deriv%s + deriv%e + deriv%i + deriv%r, 0.0_real64, &
        1.0e-10_real64, "population derivative")

    allocate(baseline_states(n_steps), vaccinated_states(n_steps))
    call simulate_seirs(initial, base, dt, baseline_states)
    call simulate_seirs(initial, vaccinated, dt, vaccinated_states)

    max_error = 0.0_real64
    do idx = 1, n_steps
        total = baseline_states(idx)%s + baseline_states(idx)%e + &
            baseline_states(idx)%i + baseline_states(idx)%r
        max_error = max(max_error, abs(total - base%population))
        if (min(baseline_states(idx)%s, baseline_states(idx)%e, &
                baseline_states(idx)%i, baseline_states(idx)%r) < -tolerance) then
            error stop "baseline trajectory produced a negative compartment"
        end if
    end do
    if (max_error > 1.0e-5_real64) error stop "population was not conserved"

    if (vaccinated_states(n_steps)%i >= baseline_states(n_steps)%i) then
        error stop "vaccination should reduce long-run infectious prevalence"
    end if
    if (vaccinated_states(n_steps)%i >= initial%i) then
        error stop "subcritical vaccination scenario should drive infection down"
    end if

    call test_waning_flow(base)

    print *, "SEIRS demographic/vaccination tests passed"

contains

    subroutine test_waning_flow(template)
        type(seirs_parameters), intent(in) :: template
        type(seirs_parameters) :: params
        type(seirs_state) :: state, dstate

        params = template
        params%beta = 0.0_real64
        params%vaccination_rate = 0.0_real64
        params%mortality_rate = 0.0_real64
        params%waning_rate = 0.01_real64

        state%s = 900.0_real64
        state%e = 0.0_real64
        state%i = 0.0_real64
        state%r = 100.0_real64
        params%population = 1000.0_real64

        dstate = seirs_derivative(state, params)
        if (dstate%s <= 0.0_real64 .or. dstate%r >= 0.0_real64) then
            error stop "waning immunity should transfer recovered individuals to susceptible"
        end if
        call assert_close(dstate%s + dstate%r, 0.0_real64, 1.0e-12_real64, "waning conservation")
    end subroutine test_waning_flow

    subroutine assert_close(actual, expected, tol, message)
        real(real64), intent(in) :: actual, expected, tol
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tol) then
            print *, "Assertion failed: ", trim(message)
            print *, "  actual:   ", actual
            print *, "  expected: ", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_seirs
