module rng_m
    use iso_fortran_env, only: int64, real64
    implicit none
    private

    integer(int64), parameter :: m1 = 4294967087_int64
    integer(int64), parameter :: m2 = 4294944443_int64
    integer(int64), parameter :: a12 = 1403580_int64
    integer(int64), parameter :: a13n = 810728_int64
    integer(int64), parameter :: a21 = 527612_int64
    integer(int64), parameter :: a23n = 1370589_int64
    integer(int64), parameter :: seed_modulus = 2147483647_int64
    integer(int64), parameter :: seed_multiplier = 48271_int64

    public :: rng_state, seed_rng, uniform_open

    type :: rng_state
        private
        integer(int64) :: s10 = 1_int64
        integer(int64) :: s11 = 1_int64
        integer(int64) :: s12 = 1_int64
        integer(int64) :: s20 = 1_int64
        integer(int64) :: s21 = 1_int64
        integer(int64) :: s22 = 1_int64
    end type rng_state

contains

    subroutine seed_rng(rng, seed)
        type(rng_state), intent(out) :: rng
        integer(int64), intent(in) :: seed
        integer(int64) :: x

        x = modulo(seed, seed_modulus - 1_int64) + 1_int64
        call advance_seed(x)
        rng%s10 = 1_int64 + modulo(x, m1 - 1_int64)
        call advance_seed(x)
        rng%s11 = 1_int64 + modulo(x, m1 - 1_int64)
        call advance_seed(x)
        rng%s12 = 1_int64 + modulo(x, m1 - 1_int64)
        call advance_seed(x)
        rng%s20 = 1_int64 + modulo(x, m2 - 1_int64)
        call advance_seed(x)
        rng%s21 = 1_int64 + modulo(x, m2 - 1_int64)
        call advance_seed(x)
        rng%s22 = 1_int64 + modulo(x, m2 - 1_int64)
    end subroutine seed_rng

    real(real64) function uniform_open(rng) result(value)
        type(rng_state), intent(inout) :: rng
        integer(int64) :: p1, p2, z

        p1 = modulo(a12 * rng%s11 - a13n * rng%s10, m1)
        rng%s10 = rng%s11
        rng%s11 = rng%s12
        rng%s12 = p1

        p2 = modulo(a21 * rng%s22 - a23n * rng%s20, m2)
        rng%s20 = rng%s21
        rng%s21 = rng%s22
        rng%s22 = p2

        if (p1 > p2) then
            z = p1 - p2
        else
            z = p1 - p2 + m1
        end if

        value = real(z, real64) / real(m1 + 1_int64, real64)
    end function uniform_open

    subroutine advance_seed(x)
        integer(int64), intent(inout) :: x

        x = modulo(seed_multiplier * x, seed_modulus)
        if (x == 0_int64) x = 1_int64
    end subroutine advance_seed

end module rng_m
