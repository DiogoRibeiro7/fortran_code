program age_structured_sir_example
    use iso_fortran_env, only: real64
    use age_structured_sir_m, only: age_sir_parameters, age_sir_state, &
        basic_reproduction_number_age_sir, contact_reciprocity_residual, simulate_age_sir
    implicit none

    integer, parameter :: n_steps = 1801
    real(real64), parameter :: dt = 0.1_real64
    type(age_sir_parameters) :: baseline, intervention
    type(age_sir_state) :: initial
    type(age_sir_state), allocatable :: baseline_states(:), intervention_states(:)
    real(real64) :: baseline_peak(3), intervention_peak(3)
    integer :: idx

    call make_parameters(baseline)
    intervention = baseline

    ! Reduce contacts involving the child group while preserving reciprocity.
    intervention%contact(1,1) = 0.25_real64 * intervention%contact(1,1)
    intervention%contact(1,2) = 0.50_real64 * intervention%contact(1,2)
    intervention%contact(2,1) = 0.50_real64 * intervention%contact(2,1)
    intervention%contact(1,3) = 0.50_real64 * intervention%contact(1,3)
    intervention%contact(3,1) = 0.50_real64 * intervention%contact(3,1)

    initial%s = baseline%population - [10.0_real64, 20.0_real64, 5.0_real64]
    initial%i = [10.0_real64, 20.0_real64, 5.0_real64]
    initial%r = [0.0_real64, 0.0_real64, 0.0_real64]

    allocate(baseline_states(n_steps), intervention_states(n_steps))
    call simulate_age_sir(initial, baseline, dt, baseline_states)
    call simulate_age_sir(initial, intervention, dt, intervention_states)

    baseline_peak = 0.0_real64
    intervention_peak = 0.0_real64
    do idx = 1, n_steps
        baseline_peak = max(baseline_peak, baseline_states(idx)%i)
        intervention_peak = max(intervention_peak, intervention_states(idx)%i)
    end do

    print '(a,f8.4)', '# baseline R0 = ', basic_reproduction_number_age_sir(baseline)
    print '(a,f8.4)', '# intervention R0 = ', basic_reproduction_number_age_sir(intervention)
    print '(a,es12.4)', '# reciprocity residual = ', contact_reciprocity_residual(baseline)
    print '(a,3(f10.2,1x))', '# baseline peaks children/adults/older = ', baseline_peak
    print '(a,3(f10.2,1x))', '# intervention peaks children/adults/older = ', intervention_peak

contains

    subroutine make_parameters(params)
        type(age_sir_parameters), intent(out) :: params

        params%beta = 0.035_real64
        params%population = [20000.0_real64, 60000.0_real64, 20000.0_real64]
        params%gamma = [0.10_real64, 0.10_real64, 0.10_real64]
        params%susceptibility = [1.10_real64, 1.00_real64, 0.90_real64]
        params%infectiousness = [1.0_real64, 1.0_real64, 1.0_real64]
        allocate(params%contact(3,3))
        params%contact(1,:) = [8.0_real64, 6.0_real64, 0.5_real64]
        params%contact(2,:) = [2.0_real64, 5.0_real64, 1.0_real64]
        params%contact(3,:) = [0.5_real64, 3.0_real64, 3.0_real64]
    end subroutine make_parameters

end program age_structured_sir_example
