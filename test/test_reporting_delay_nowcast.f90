program test_reporting_delay_nowcast
    use iso_fortran_env, only: real64
    use reporting_delay_nowcast_m, only: delay_distribution, nowcast_summary, &
        validate_delay_distribution, reporting_completeness, nowcast_event_day
    implicit none

    type(delay_distribution) :: delay
    type(nowcast_summary) :: complete, recent

    allocate(delay%probability(4))
    delay%probability = [0.40_real64, 0.30_real64, 0.20_real64, 0.10_real64]

    if (.not. validate_delay_distribution(delay)) error stop "valid delay distribution rejected"
    call assert_close(reporting_completeness(delay, 0), 0.40_real64, 1.0e-12_real64, &
        "day-0 completeness")
    call assert_close(reporting_completeness(delay, 3), 1.0_real64, 1.0e-12_real64, &
        "complete reporting")

    call nowcast_event_day(20.0_real64, 3, delay, 0.001_real64, 0.001_real64, &
        0.95_real64, complete)
    call nowcast_event_day(8.0_real64, 0, delay, 0.001_real64, 0.001_real64, &
        0.95_real64, recent)

    if (.not. complete%valid .or. .not. recent%valid) error stop "valid nowcast marked invalid"
    if (recent%mean <= 8.0_real64) error stop "incomplete reporting should increase the nowcast"
    if (recent%expected_future_reports <= 0.0_real64) then
        error stop "recent event day should have an expected reporting backlog"
    end if
    if ((recent%upper - recent%lower) <= 0.0_real64) error stop "nowcast interval must have positive width"
    if (complete%expected_future_reports > 1.0e-12_real64) then
        error stop "fully reported event day should have no expected backlog"
    end if
    if (abs(complete%mean - 20.0_real64) > 0.1_real64) then
        error stop "weak-prior complete nowcast should recover the reported count"
    end if

    print *, "Reporting-delay nowcast tests passed"

contains

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

end program test_reporting_delay_nowcast
