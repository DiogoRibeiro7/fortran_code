module transmission_m
    use iso_fortran_env, only: real64
    implicit none
    private

    public :: transmission_schedule, beta_at, schedule_is_valid

    type :: transmission_schedule
        real(real64) :: baseline_beta = 0.30_real64
        real(real64), allocatable :: change_times(:)
        real(real64), allocatable :: multipliers(:)
    end type transmission_schedule

contains

    pure logical function schedule_is_valid(schedule) result(valid)
        type(transmission_schedule), intent(in) :: schedule
        integer :: idx

        valid = schedule%baseline_beta >= 0.0_real64
        if (.not. valid) return

        if (allocated(schedule%change_times) .neqv. allocated(schedule%multipliers)) then
            valid = .false.
            return
        end if

        if (.not. allocated(schedule%change_times)) return
        if (size(schedule%change_times) /= size(schedule%multipliers)) then
            valid = .false.
            return
        end if
        if (any(schedule%multipliers < 0.0_real64)) then
            valid = .false.
            return
        end if

        do idx = 2, size(schedule%change_times)
            if (schedule%change_times(idx) <= schedule%change_times(idx - 1)) then
                valid = .false.
                return
            end if
        end do
    end function schedule_is_valid

    pure real(real64) function beta_at(schedule, time) result(beta)
        type(transmission_schedule), intent(in) :: schedule
        real(real64), intent(in) :: time
        integer :: idx

        beta = schedule%baseline_beta
        if (.not. allocated(schedule%change_times)) return

        do idx = 1, size(schedule%change_times)
            if (time < schedule%change_times(idx)) exit
            beta = schedule%baseline_beta * schedule%multipliers(idx)
        end do
    end function beta_at

end module transmission_m
