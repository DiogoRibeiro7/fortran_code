module lapack_solvers_m
    use, intrinsic :: iso_fortran_env, only : real64
    implicit none (type, external)
    private

    integer, parameter, public :: solve_success = 0
    integer, parameter, public :: solve_invalid_shape = -1000
    public :: solve_linear_system

    interface
        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            import :: real64
            integer, intent(in) :: n
            integer, intent(in) :: nrhs
            integer, intent(in) :: lda
            integer, intent(out) :: ipiv(*)
            real(real64), intent(inout) :: a(lda, *)
            real(real64), intent(inout) :: b(ldb, *)
            integer, intent(in) :: ldb
            integer, intent(out) :: info
        end subroutine dgesv
    end interface

contains

    subroutine solve_linear_system(a, b, x, info)
        real(real64), intent(in) :: a(:, :)
        real(real64), intent(in) :: b(:)
        real(real64), intent(out) :: x(:)
        integer, intent(out) :: info

        real(real64), allocatable :: a_work(:, :)
        real(real64), allocatable :: rhs(:, :)
        integer, allocatable :: pivots(:)
        integer :: n

        n = size(a, dim=1)
        if (size(a, dim=2) /= n .or. size(b) /= n .or. size(x) /= n) then
            info = solve_invalid_shape
            return
        end if

        allocate(a_work(n, n), rhs(n, 1), pivots(n))
        a_work = a
        rhs(:, 1) = b

        call dgesv(n, 1, a_work, n, pivots, rhs, n, info)
        if (info == solve_success) x = rhs(:, 1)
    end subroutine solve_linear_system

end module lapack_solvers_m
