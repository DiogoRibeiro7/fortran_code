program test_nowcast_rt_uncertainty
    use iso_fortran_env, only: int64, real64
    use reporting_delay_nowcast_m, only: delay_distribution
    use nowcast_rt_uncertainty_m, only: rt_monte_carlo_summary, propagate_nowcast_uncertainty_rt
    implicit none

    call test_reproducibility_and_complete_reporting()
    call test_incomplete_reporting_uncertainty()

    print *, "Nowcast/Rt uncertainty tests passed"

contains

    subroutine test_reproducibility_and_complete_reporting()
        integer, parameter :: n = 20
        real(real64) :: reported(n), generation_weights(1)
        integer :: ages(n)
        type(delay_distribution) :: delay
        type(rt_monte_carlo_summary) :: first(n), second(n)
        integer :: idx

        reported = 100.0_real64
        ages = 0
        generation_weights = [1.0_real64]
        allocate(delay%probability(1))
        delay%probability = [1.0_real64]

        call propagate_nowcast_uncertainty_rt(reported, ages, delay, generation_weights, &
            0.01_real64, 0.0001_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
            20260909_int64, 500, first)
        call propagate_nowcast_uncertainty_rt(reported, ages, delay, generation_weights, &
            0.01_real64, 0.0001_real64, 7, 1.0_real64, 1.0_real64, 0.95_real64, &
            20260909_int64, 500, second)

        do idx = 7, n
            if (.not. first(idx)%valid .or. .not. second(idx)%valid) then
                error stop "complete-reporting Rt summary should be valid"
            end if
            call assert_close(first(idx)%mean, second(idx)%mean, 0.0_real64, &
                "fixed-seed mean reproducibility")
            call assert_close(first(idx)%lower, second(idx)%lower, 0.0_real64, &
                "fixed-seed lower reproducibility")
            call assert_close(first(idx)%upper, second(idx)%upper, 0.0_real64, &
                "fixed-seed upper reproducibility")
        end do

        if (abs(first(n)%mean - 1.0_real64) > 0.08_real64) then
            error stop "complete constant incidence should imply Rt near one"
        end if
    end subroutine test_reproducibility_and_complete_reporting

    subroutine test_incomplete_reporting_uncertainty()
        integer, parameter :: n = 20
        real(real64) :: reported_complete(n), reported_incomplete(n), generation_weights(1)
        integer :: ages_complete(n), ages_incomplete(n), idx
        type(delay_distribution) :: complete_delay, incomplete_delay
        type(rt_monte_carlo_summary) :: complete_summary(n), incomplete_summary(n)
        real(real64) :: complete_width, incomplete_width

        reported_complete = 100.0_real64
        reported_incomplete = 100.0_real64
        ages_complete = 0
        ages_incomplete = 10
        generation_weights = [1.0_real64]

        allocate(complete_delay%probability(1))
        complete_delay%probability = [1.0_real64]
        allocate(incomplete_delay%probability(4))
        incomplete_delay%probability = [0.20_real64, 0.30_real64, 0.30_real64, 0.20_real64]

        ages_incomplete(n - 2) = 2
        ages_incomplete(n - 1) = 1
        ages_incomplete(n) = 0
        reported_incomplete(n - 2) = 80.0_real64
        reported_incomplete(n - 1) = 50.0_real64
        reported_incomplete(n) = 20.0_real64

        call propagate_nowcast_uncertainty_rt(reported_complete, ages_complete, complete_delay, &
            generation_weights, 0.01_real64, 0.0001_real64, 7, 1.0_real64, 1.0_real64, &
            0.95_real64, 777_int64, 1000, complete_summary)
        call propagate_nowcast_uncertainty_rt(reported_incomplete, ages_incomplete, incomplete_delay, &
            generation_weights, 0.01_real64, 0.0001_real64, 7, 1.0_real64, 1.0_real64, &
            0.95_real64, 777_int64, 1000, incomplete_summary)

        idx = n
        if (.not. complete_summary(idx)%valid .or. .not. incomplete_summary(idx)%valid) then
            error stop "final Rt summaries should be valid"
        end if

        complete_width = complete_summary(idx)%upper - complete_summary(idx)%lower
        incomplete_width = incomplete_summary(idx)%upper - incomplete_summary(idx)%lower
        if (incomplete_width <= complete_width) then
            error stop "incomplete reporting should widen propagated Rt uncertainty"
        end if
        if (abs(incomplete_summary(idx)%mean - 1.0_real64) > 0.20_real64) then
            error stop "nowcast uncertainty propagation should retain Rt near one for flat latent incidence"
        end if
    end subroutine test_incomplete_reporting_uncertainty

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

end program test_nowcast_rt_uncertainty
