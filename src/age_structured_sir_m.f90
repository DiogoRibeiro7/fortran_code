module age_structured_sir_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: age_sir_parameters, age_sir_state
    public :: force_of_infection_age_sir, next_generation_matrix_age_sir
    public :: basic_reproduction_number_age_sir, contact_reciprocity_residual
    public :: age_sir_derivative, rk4_step_age_sir, simulate_age_sir

    type :: age_sir_parameters
        real(real64) :: beta = 0.0_real64
        real(real64), allocatable :: population(:)
        real(real64), allocatable :: gamma(:)
        real(real64), allocatable :: susceptibility(:)
        real(real64), allocatable :: infectiousness(:)
        real(real64), allocatable :: contact(:,:)
    end type age_sir_parameters

    type :: age_sir_state
        real(real64), allocatable :: s(:)
        real(real64), allocatable :: i(:)
        real(real64), allocatable :: r(:)
    end type age_sir_state

    interface
        subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
            import :: real64
            character(len=1), intent(in) :: jobvl, jobvr
            integer, intent(in) :: n, lda, ldvl, ldvr, lwork
            real(real64), intent(inout) :: a(lda, *)
            real(real64), intent(out) :: wr(*), wi(*), vl(ldvl, *), vr(ldvr, *), work(*)
            integer, intent(out) :: info
        end subroutine dgeev
    end interface

contains

    subroutine force_of_infection_age_sir(state, params, lambda)
        type(age_sir_state), intent(in) :: state
        type(age_sir_parameters), intent(in) :: params
        real(real64), intent(out) :: lambda(:)
        integer :: i_group, j_group, n_groups

        call validate_state_and_parameters(state, params)
        n_groups = size(params%population)
        if (size(lambda) /= n_groups) error stop "lambda size must match number of groups"

        lambda = 0.0_real64
        do i_group = 1, n_groups
            do j_group = 1, n_groups
                lambda(i_group) = lambda(i_group) + params%contact(i_group, j_group) * &
                    params%infectiousness(j_group) * state%i(j_group) / params%population(j_group)
            end do
            lambda(i_group) = params%beta * params%susceptibility(i_group) * lambda(i_group)
        end do
    end subroutine force_of_infection_age_sir

    subroutine next_generation_matrix_age_sir(params, matrix)
        type(age_sir_parameters), intent(in) :: params
        real(real64), intent(out) :: matrix(:,:)
        integer :: i_group, j_group, n_groups

        call validate_parameters(params)
        n_groups = size(params%population)
        if (size(matrix, 1) /= n_groups .or. size(matrix, 2) /= n_groups) then
            error stop "next-generation matrix has incorrect dimensions"
        end if

        do i_group = 1, n_groups
            do j_group = 1, n_groups
                matrix(i_group, j_group) = params%beta * params%susceptibility(i_group) * &
                    params%population(i_group) * params%contact(i_group, j_group) * &
                    params%infectiousness(j_group) / &
                    (params%population(j_group) * params%gamma(j_group))
            end do
        end do
    end subroutine next_generation_matrix_age_sir

    real(real64) function basic_reproduction_number_age_sir(params) result(r0)
        type(age_sir_parameters), intent(in) :: params
        real(real64), allocatable :: matrix(:,:)
        integer :: n_groups

        call validate_parameters(params)
        n_groups = size(params%population)
        allocate(matrix(n_groups, n_groups))
        call next_generation_matrix_age_sir(params, matrix)
        r0 = spectral_radius(matrix)
    end function basic_reproduction_number_age_sir

    real(real64) function contact_reciprocity_residual(params) result(residual)
        type(age_sir_parameters), intent(in) :: params
        real(real64) :: numerator, denominator
        integer :: i_group, j_group, n_groups

        call validate_parameters(params)
        residual = 0.0_real64
        n_groups = size(params%population)
        do i_group = 1, n_groups
            do j_group = i_group + 1, n_groups
                numerator = abs(params%population(i_group) * params%contact(i_group, j_group) - &
                    params%population(j_group) * params%contact(j_group, i_group))
                denominator = max(params%population(i_group) * params%contact(i_group, j_group), &
                    params%population(j_group) * params%contact(j_group, i_group), tiny(1.0_real64))
                residual = max(residual, numerator / denominator)
            end do
        end do
    end function contact_reciprocity_residual

    subroutine age_sir_derivative(state, params, dstatedt)
        type(age_sir_state), intent(in) :: state
        type(age_sir_parameters), intent(in) :: params
        type(age_sir_state), intent(out) :: dstatedt
        real(real64), allocatable :: lambda(:), incidence(:)
        integer :: n_groups

        call validate_state_and_parameters(state, params)
        n_groups = size(params%population)
        allocate(lambda(n_groups), incidence(n_groups))
        allocate(dstatedt%s(n_groups), dstatedt%i(n_groups), dstatedt%r(n_groups))

        call force_of_infection_age_sir(state, params, lambda)
        incidence = lambda * state%s
        dstatedt%s = -incidence
        dstatedt%i = incidence - params%gamma * state%i
        dstatedt%r = params%gamma * state%i
    end subroutine age_sir_derivative

    function rk4_step_age_sir(state, params, dt) result(next_state)
        type(age_sir_state), intent(in) :: state
        type(age_sir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(age_sir_state) :: next_state
        type(age_sir_state) :: k1, k2, k3, k4, tmp

        if (dt <= 0.0_real64) error stop "dt must be positive"
        call age_sir_derivative(state, params, k1)
        tmp = add_scaled(state, k1, 0.5_real64 * dt)
        call age_sir_derivative(tmp, params, k2)
        tmp = add_scaled(state, k2, 0.5_real64 * dt)
        call age_sir_derivative(tmp, params, k3)
        tmp = add_scaled(state, k3, dt)
        call age_sir_derivative(tmp, params, k4)
        next_state = combine_rk4(state, k1, k2, k3, k4, dt)
    end function rk4_step_age_sir

    subroutine simulate_age_sir(initial_state, params, dt, states)
        type(age_sir_state), intent(in) :: initial_state
        type(age_sir_parameters), intent(in) :: params
        real(real64), intent(in) :: dt
        type(age_sir_state), intent(out) :: states(:)
        integer :: idx

        if (size(states) < 1) error stop "states must contain at least one element"
        if (dt <= 0.0_real64) error stop "dt must be positive"
        call validate_state_and_parameters(initial_state, params)

        states(1) = initial_state
        do idx = 2, size(states)
            states(idx) = rk4_step_age_sir(states(idx - 1), params, dt)
        end do
    end subroutine simulate_age_sir

    real(real64) function spectral_radius(matrix) result(radius)
        real(real64), intent(in) :: matrix(:,:)
        real(real64), allocatable :: a(:,:), wr(:), wi(:), work(:)
        real(real64) :: vl(1,1), vr(1,1)
        integer :: n, info, lwork

        if (size(matrix, 1) /= size(matrix, 2)) error stop "spectral radius requires a square matrix"
        n = size(matrix, 1)
        if (n < 1) error stop "spectral radius requires a non-empty matrix"

        lwork = max(1, 4 * n)
        allocate(a(n,n), wr(n), wi(n), work(lwork))
        a = matrix
        call dgeev('N', 'N', n, a, n, wr, wi, vl, 1, vr, 1, work, lwork, info)
        if (info /= 0) error stop "LAPACK DGEEV failed while computing R0"
        radius = maxval(sqrt(wr * wr + wi * wi))
    end function spectral_radius

    function add_scaled(state, derivative, scale) result(updated)
        type(age_sir_state), intent(in) :: state, derivative
        real(real64), intent(in) :: scale
        type(age_sir_state) :: updated
        integer :: n_groups

        n_groups = size(state%s)
        allocate(updated%s(n_groups), updated%i(n_groups), updated%r(n_groups))
        updated%s = state%s + scale * derivative%s
        updated%i = state%i + scale * derivative%i
        updated%r = state%r + scale * derivative%r
    end function add_scaled

    function combine_rk4(state, k1, k2, k3, k4, dt) result(updated)
        type(age_sir_state), intent(in) :: state, k1, k2, k3, k4
        real(real64), intent(in) :: dt
        type(age_sir_state) :: updated
        integer :: n_groups

        n_groups = size(state%s)
        allocate(updated%s(n_groups), updated%i(n_groups), updated%r(n_groups))
        updated%s = state%s + dt * (k1%s + 2.0_real64*k2%s + 2.0_real64*k3%s + k4%s) / 6.0_real64
        updated%i = state%i + dt * (k1%i + 2.0_real64*k2%i + 2.0_real64*k3%i + k4%i) / 6.0_real64
        updated%r = state%r + dt * (k1%r + 2.0_real64*k2%r + 2.0_real64*k3%r + k4%r) / 6.0_real64
    end function combine_rk4

    subroutine validate_parameters(params)
        type(age_sir_parameters), intent(in) :: params
        integer :: n_groups

        if (params%beta < 0.0_real64) error stop "beta must be non-negative"
        if (.not. allocated(params%population) .or. .not. allocated(params%gamma) .or. &
                .not. allocated(params%susceptibility) .or. .not. allocated(params%infectiousness) .or. &
                .not. allocated(params%contact)) then
            error stop "all age-structured parameter arrays must be allocated"
        end if

        n_groups = size(params%population)
        if (n_groups < 1) error stop "at least one group is required"
        if (size(params%gamma) /= n_groups .or. size(params%susceptibility) /= n_groups .or. &
                size(params%infectiousness) /= n_groups .or. size(params%contact, 1) /= n_groups .or. &
                size(params%contact, 2) /= n_groups) then
            error stop "age-structured parameter dimensions are inconsistent"
        end if
        if (any(params%population <= 0.0_real64)) error stop "group populations must be positive"
        if (any(params%gamma <= 0.0_real64)) error stop "recovery rates must be positive"
        if (any(params%susceptibility < 0.0_real64)) error stop "susceptibility must be non-negative"
        if (any(params%infectiousness < 0.0_real64)) error stop "infectiousness must be non-negative"
        if (any(params%contact < 0.0_real64)) error stop "contact rates must be non-negative"
    end subroutine validate_parameters

    subroutine validate_state_and_parameters(state, params)
        type(age_sir_state), intent(in) :: state
        type(age_sir_parameters), intent(in) :: params
        integer :: n_groups
        real(real64), parameter :: tolerance = 1.0e-8_real64

        call validate_parameters(params)
        n_groups = size(params%population)
        if (.not. allocated(state%s) .or. .not. allocated(state%i) .or. .not. allocated(state%r)) then
            error stop "all age-structured state arrays must be allocated"
        end if
        if (size(state%s) /= n_groups .or. size(state%i) /= n_groups .or. size(state%r) /= n_groups) then
            error stop "age-structured state dimensions are inconsistent"
        end if
        if (any(state%s < -tolerance) .or. any(state%i < -tolerance) .or. any(state%r < -tolerance)) then
            error stop "age-structured compartments must be non-negative"
        end if
        if (any(abs(state%s + state%i + state%r - params%population) > tolerance * params%population)) then
            error stop "group compartment totals must match group populations"
        end if
    end subroutine validate_state_and_parameters

end module age_structured_sir_m
