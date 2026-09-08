program test_age_structured_sir
    use iso_fortran_env, only: real64
    use age_structured_sir_m, only: age_sir_parameters, age_sir_state, &
        force_of_infection_age_sir, next_generation_matrix_age_sir, &
        basic_reproduction_number_age_sir, contact_reciprocity_residual, simulate_age_sir
    implicit none

    call test_one_group_reduction()
    call test_diagonal_next_generation_matrix()
    call test_reciprocity_and_intervention()
    call test_group_conservation()

    print *, "Age-structured SIR tests passed"

contains

    subroutine test_one_group_reduction()
        type(age_sir_parameters) :: params
        type(age_sir_state) :: state
        real(real64) :: lambda(1), ngm(1,1)

        params%beta = 0.03_real64
        params%population = [1000.0_real64]
        params%gamma = [0.10_real64]
        params%susceptibility = [1.0_real64]
        params%infectiousness = [1.0_real64]
        allocate(params%contact(1,1))
        params%contact(1,1) = 10.0_real64

        state%s = [990.0_real64]
        state%i = [10.0_real64]
        state%r = [0.0_real64]

        call force_of_infection_age_sir(state, params, lambda)
        call assert_close(lambda(1), 0.003_real64, 1.0e-12_real64, "one-group force of infection")
        call next_generation_matrix_age_sir(params, ngm)
        call assert_close(ngm(1,1), 3.0_real64, 1.0e-12_real64, "one-group NGM")
        call assert_close(basic_reproduction_number_age_sir(params), 3.0_real64, &
            1.0e-12_real64, "one-group R0")
    end subroutine test_one_group_reduction

    subroutine test_diagonal_next_generation_matrix()
        type(age_sir_parameters) :: params
        real(real64) :: ngm(2,2)

        params%beta = 0.02_real64
        params%population = [1000.0_real64, 2000.0_real64]
        params%gamma = [0.10_real64, 0.20_real64]
        params%susceptibility = [1.0_real64, 1.0_real64]
        params%infectiousness = [1.0_real64, 1.0_real64]
        allocate(params%contact(2,2))
        params%contact = 0.0_real64
        params%contact(1,1) = 10.0_real64
        params%contact(2,2) = 4.0_real64

        call next_generation_matrix_age_sir(params, ngm)
        call assert_close(ngm(1,1), 2.0_real64, 1.0e-12_real64, "group-one diagonal NGM")
        call assert_close(ngm(2,2), 0.4_real64, 1.0e-12_real64, "group-two diagonal NGM")
        call assert_close(ngm(1,2), 0.0_real64, 1.0e-12_real64, "off-diagonal NGM")
        call assert_close(basic_reproduction_number_age_sir(params), 2.0_real64, &
            1.0e-12_real64, "diagonal spectral radius")
    end subroutine test_diagonal_next_generation_matrix

    subroutine test_reciprocity_and_intervention()
        type(age_sir_parameters) :: base, reduced
        real(real64) :: base_r0, reduced_r0

        call make_three_group_parameters(base)
        call assert_close(contact_reciprocity_residual(base), 0.0_real64, &
            1.0e-12_real64, "reciprocal contact matrix")
        base_r0 = basic_reproduction_number_age_sir(base)

        reduced = base
        reduced%contact(1,1) = 0.25_real64 * reduced%contact(1,1)
        reduced%contact(1,2) = 0.50_real64 * reduced%contact(1,2)
        reduced%contact(2,1) = 0.50_real64 * reduced%contact(2,1)
        reduced%contact(1,3) = 0.50_real64 * reduced%contact(1,3)
        reduced%contact(3,1) = 0.50_real64 * reduced%contact(3,1)
        reduced_r0 = basic_reproduction_number_age_sir(reduced)

        if (base_r0 <= 1.0_real64) error stop "reference contact matrix should be supercritical"
        if (reduced_r0 >= base_r0) error stop "child-contact reduction should lower R0"
        call assert_close(contact_reciprocity_residual(reduced), 0.0_real64, &
            1.0e-12_real64, "reciprocity after symmetric intervention")
    end subroutine test_reciprocity_and_intervention

    subroutine test_group_conservation()
        integer, parameter :: n_steps = 1201
        real(real64), parameter :: dt = 0.1_real64
        type(age_sir_parameters) :: params
        type(age_sir_state) :: initial
        type(age_sir_state), allocatable :: states(:)
        integer :: idx, group
        real(real64) :: total

        call make_three_group_parameters(params)
        allocate(states(n_steps))
        initial%s = params%population - [10.0_real64, 20.0_real64, 5.0_real64]
        initial%i = [10.0_real64, 20.0_real64, 5.0_real64]
        initial%r = [0.0_real64, 0.0_real64, 0.0_real64]
        call simulate_age_sir(initial, params, dt, states)

        do idx = 1, n_steps
            do group = 1, 3
                total = states(idx)%s(group) + states(idx)%i(group) + states(idx)%r(group)
                call assert_close(total, params%population(group), 1.0e-6_real64, &
                    "within-group population conservation")
                if (min(states(idx)%s(group), states(idx)%i(group), states(idx)%r(group)) < -1.0e-7_real64) then
                    error stop "age-structured trajectory produced a negative compartment"
                end if
            end do
        end do
    end subroutine test_group_conservation

    subroutine make_three_group_parameters(params)
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
    end subroutine make_three_group_parameters

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message

        if (abs(actual - expected) > tolerance) then
            print *, "Assertion failed: ", trim(message)
            print *, "  actual:   ", actual
            print *, "  expected: ", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_age_structured_sir
