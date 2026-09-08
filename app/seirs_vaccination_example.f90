program seirs_vaccination_example
    use iso_fortran_env, only: real64
    use seirs_m, only: seirs_parameters, seirs_state, basic_reproduction_number_seirs, &
        control_reproduction_number, simulate_seirs
    implicit none

    integer, parameter :: n_days = 3650
    integer, parameter :: n_states = n_days + 1
    real(real64), parameter :: dt = 1.0_real64
    type(seirs_parameters) :: baseline, vaccinated
    type(seirs_state) :: initial
    type(seirs_state), allocatable :: baseline_states(:), vaccinated_states(:)
    integer :: idx

    baseline%population = 100000.0_real64
    baseline%beta = 0.30_real64
    baseline%sigma = 0.20_real64
    baseline%gamma = 0.10_real64
    baseline%mortality_rate = 1.0_real64 / (75.0_real64 * 365.0_real64)
    baseline%waning_rate = 1.0_real64 / (3.0_real64 * 365.0_real64)
    baseline%vaccination_rate = 0.0_real64

    vaccinated = baseline
    vaccinated%vaccination_rate = 0.0025_real64

    initial%s = baseline%population - 30.0_real64
    initial%e = 20.0_real64
    initial%i = 10.0_real64
    initial%r = 0.0_real64

    allocate(baseline_states(n_states), vaccinated_states(n_states))
    call simulate_seirs(initial, baseline, dt, baseline_states)
    call simulate_seirs(initial, vaccinated, dt, vaccinated_states)

    print '(a,f8.4)', '# R0 = ', basic_reproduction_number_seirs(baseline)
    print '(a,f8.4)', '# Rc without vaccination = ', control_reproduction_number(baseline)
    print '(a,f8.4)', '# Rc with vaccination = ', control_reproduction_number(vaccinated)
    print '(a)', 'day,s_no_vax,e_no_vax,i_no_vax,r_no_vax,s_vax,e_vax,i_vax,r_vax'

    do idx = 1, n_states
        print '(f8.1,8(",",f14.6))', real(idx - 1, real64) * dt, &
            baseline_states(idx)%s, baseline_states(idx)%e, baseline_states(idx)%i, baseline_states(idx)%r, &
            vaccinated_states(idx)%s, vaccinated_states(idx)%e, vaccinated_states(idx)%i, vaccinated_states(idx)%r
    end do
end program seirs_vaccination_example
