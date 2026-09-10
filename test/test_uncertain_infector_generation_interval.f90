program test_uncertain_infector_generation_interval
    use iso_fortran_env, only: real64
    use uncertain_infector_generation_interval_m, only: fit_uncertain_infector_gamma, &
        uncertain_infector_gamma_log_likelihood
    use interval_censored_generation_interval_m, only: mixed_censored_gamma_log_likelihood, &
        interval_censored_gamma_fit
    implicit none

    integer, parameter :: start_single(4) = [1,2,3,4]
    integer, parameter :: end_single(4) = [1,2,3,4]
    real(real64), parameter :: lower_single(4) = [1.0_real64,2.0_real64,3.0_real64,4.0_real64]
    real(real64), parameter :: upper_single(4) = [1.8_real64,2.8_real64,3.8_real64,4.8_real64]
    logical, parameter :: right_single(4) = [.false.,.false.,.false.,.true.]
    real(real64), parameter :: weight_single(4) = [1.0_real64,1.0_real64,1.0_real64,1.0_real64]

    integer, parameter :: start_dup(2) = [1,3]
    integer, parameter :: end_dup(2) = [2,4]
    real(real64), parameter :: lower_dup(4) = [1.0_real64,1.0_real64,3.0_real64,3.0_real64]
    real(real64), parameter :: upper_dup(4) = [1.8_real64,1.8_real64,3.8_real64,3.8_real64]
    logical, parameter :: right_dup(4) = [.false.,.false.,.false.,.false.]
    real(real64), parameter :: weight_dup(4) = [0.5_real64,0.5_real64,0.2_real64,0.8_real64]

    integer, parameter :: start_hard(2) = [1,3]
    integer, parameter :: end_hard(2) = [2,4]
    real(real64), parameter :: lower_hard(4) = [1.0_real64,6.0_real64,2.0_real64,8.0_real64]
    real(real64), parameter :: upper_hard(4) = [1.8_real64,6.8_real64,2.8_real64,8.8_real64]
    logical, parameter :: right_hard(4) = [.false.,.false.,.false.,.false.]
    real(real64), parameter :: weight_hard(4) = [1.0_real64,1.0e-12_real64,1.0_real64,1.0e-12_real64]
    real(real64), parameter :: hard_lower(2) = [1.0_real64,2.0_real64]
    real(real64), parameter :: hard_upper(2) = [1.8_real64,2.8_real64]
    logical, parameter :: hard_right(2) = [.false.,.false.]

    type(interval_censored_gamma_fit) :: fit
    real(real64) :: ll_mixture, ll_reference, ll_duplicate, ll_hard

    ll_mixture = uncertain_infector_gamma_log_likelihood(start_single, end_single, lower_single, &
        upper_single, right_single, weight_single, 2.5_real64, 0.5_real64)
    ll_reference = mixed_censored_gamma_log_likelihood(lower_single, upper_single, right_single, &
        2.5_real64, 0.5_real64)
    call assert_close(ll_mixture, ll_reference, 1.0e-12_real64, "single-candidate reduction")

    ll_duplicate = uncertain_infector_gamma_log_likelihood(start_dup, end_dup, lower_dup, upper_dup, &
        right_dup, weight_dup, 2.5_real64, 0.5_real64)
    ll_reference = mixed_censored_gamma_log_likelihood([1.0_real64,3.0_real64], &
        [1.8_real64,3.8_real64], [.false.,.false.], 2.5_real64, 0.5_real64)
    call assert_close(ll_duplicate, ll_reference, 1.0e-12_real64, "duplicate-candidate invariance")

    ll_hard = uncertain_infector_gamma_log_likelihood(start_hard, end_hard, lower_hard, upper_hard, &
        right_hard, weight_hard, 2.5_real64, 0.5_real64)
    ll_reference = mixed_censored_gamma_log_likelihood(hard_lower, hard_upper, hard_right, &
        2.5_real64, 0.5_real64)
    call assert_close(ll_hard, ll_reference, 1.0e-9_real64, "hard-assignment limit")

    call fit_uncertain_infector_gamma(start_hard, end_hard, lower_hard, upper_hard, right_hard, &
        [0.85_real64,0.15_real64,0.75_real64,0.25_real64], fit)
    if (.not. fit%converged) error stop "uncertain-infector Gamma fit did not converge"
    if (fit%shape <= 0.0_real64 .or. fit%rate <= 0.0_real64 .or. fit%mean <= 0.0_real64) then
        error stop "uncertain-infector fit returned invalid parameters"
    end if

    print *, "Uncertain-infector generation-interval tests passed"

contains

    subroutine assert_close(actual, expected, tolerance, message)
        real(real64), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: message
        if (abs(actual-expected) > tolerance) then
            print *, trim(message), actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_uncertain_infector_generation_interval
