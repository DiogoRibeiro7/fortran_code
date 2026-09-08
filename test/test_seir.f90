program test_seir
    use iso_fortran_env, only: real64
    use seir_m, only: seir_parameters, seir_state, basic_reproduction_number_seir, &
        simulate_seir, simulate_seir_scheduled
    use transmission_m, only: transmission_schedule, beta_at, schedule_is_valid
    implicit none

    integer, parameter :: n_steps = 241
    real(real64), parameter :: tol = 1.0e-6_real64
    type(seir_parameters) :: params
    type(seir_state) :: initial_state
    type(seir_state) :: baseline(n_steps), intervention(n_steps)
    type(transmission_schedule) :: schedule
    real(real64) :: total, expected_population
    integer :: idx

    params = seir_parameters(beta=0.30_real64, sigma=0.20_real64, &
        gamma=0.10_real64, population=100000.0_real64)
    initial_state = seir_state(s=99990.0_real64, e=0.0_real64, &
        i=10.0_real64, r=0.0_real64)
    expected_population = params%population

    if (abs(basic_reproduction_number_seir(params) - 3.0_real64) > tol) &
        error stop "incorrect SEIR R0"

    call simulate_seir(initial_state, params, 0.5_real64, baseline)
    do idx = 1, n_steps
        total = baseline(idx)%s + baseline(idx)%e + baseline(idx)%i + baseline(idx)%r
        if (abs(total - expected_population) > 1.0e-5_real64) &
            error stop "SEIR population is not conserved"
        if (min(baseline(idx)%s, baseline(idx)%e, baseline(idx)%i, baseline(idx)%r) < -tol) &
            error stop "negative compartment"
    end do

    if (baseline(3)%e <= baseline(1)%e) &
        error stop "exposed compartment should initially grow"

    schedule%baseline_beta = params%beta
    schedule%change_times = [20.0_real64, 60.0_real64]
    schedule%multipliers = [0.35_real64, 0.70_real64]

    if (.not. schedule_is_valid(schedule)) error stop "valid schedule rejected"
    if (abs(beta_at(schedule, 10.0_real64) - 0.30_real64) > tol) &
        error stop "baseline beta mismatch"
    if (abs(beta_at(schedule, 30.0_real64) - 0.105_real64) > tol) &
        error stop "intervention beta mismatch"
    if (abs(beta_at(schedule, 80.0_real64) - 0.21_real64) > tol) &
        error stop "relaxed beta mismatch"

    call simulate_seir_scheduled(initial_state, params, schedule, 0.0_real64, &
        0.5_real64, intervention)

    if (intervention(n_steps)%r >= baseline(n_steps)%r) &
        error stop "intervention should reduce cumulative infections"
    if (maxval(intervention%i) >= maxval(baseline%i)) &
        error stop "intervention should reduce peak prevalence"

    print *, "SEIR tests passed"
end program test_seir
