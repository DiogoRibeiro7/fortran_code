cmake_minimum_required(VERSION 3.25)

project(fortran_epidemiology VERSION 0.16.0 LANGUAGES Fortran)

include(CTest)
find_package(LAPACK REQUIRED)

function(enable_strict_fortran target_name)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target_name} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wconversion-extra
            -fcheck=all
            -fbacktrace
        )
    endif()
endfunction()

add_library(fortran_numerics
    src/factorial_m.f90
    src/lapack_solvers_m.f90
)
target_link_libraries(fortran_numerics PUBLIC LAPACK::LAPACK)
enable_strict_fortran(fortran_numerics)

add_library(epidemiology
    src/epidemiology_m.f90
    src/transmission_m.f90
    src/seir_m.f90
    src/seirs_m.f90
    src/rng_m.f90
    src/random_distributions_m.f90
    src/stochastic_sir_m.f90
    src/stochastic_seir_m.f90
    src/age_structured_sir_m.f90
    src/growth_inference_m.f90
    src/gamma_distribution_m.f90
    src/generation_interval_fit_m.f90
    src/renewal_rt_m.f90
    src/negative_binomial_rt_m.f90
    src/reporting_delay_nowcast_m.f90
    src/nowcasted_rt_pipeline_m.f90
    src/nowcast_rt_uncertainty_m.f90
    src/generation_interval_uncertainty_m.f90
    src/fitted_generation_interval_uncertainty_m.f90
)
target_link_libraries(epidemiology PUBLIC LAPACK::LAPACK)
enable_strict_fortran(epidemiology)

add_executable(factorial_example app/factorial_example.f90)
target_link_libraries(factorial_example PRIVATE fortran_numerics)
enable_strict_fortran(factorial_example)

add_executable(for_loop_example app/for_loop_example.f90)
enable_strict_fortran(for_loop_example)

add_executable(lapack_example app/lapack_example.f90)
target_link_libraries(lapack_example PRIVATE fortran_numerics)
enable_strict_fortran(lapack_example)

add_executable(sir_example app/sir_example.f90)
target_link_libraries(sir_example PRIVATE epidemiology)
enable_strict_fortran(sir_example)

add_executable(seir_intervention_example app/seir_intervention_example.f90)
target_link_libraries(seir_intervention_example PRIVATE epidemiology)
enable_strict_fortran(seir_intervention_example)

add_executable(seirs_vaccination_example app/seirs_vaccination_example.f90)
target_link_libraries(seirs_vaccination_example PRIVATE epidemiology)
enable_strict_fortran(seirs_vaccination_example)

add_executable(stochastic_sir_example app/stochastic_sir_example.f90)
target_link_libraries(stochastic_sir_example PRIVATE epidemiology)
enable_strict_fortran(stochastic_sir_example)

add_executable(stochastic_seir_ensemble_example app/stochastic_seir_ensemble_example.f90)
target_link_libraries(stochastic_seir_ensemble_example PRIVATE epidemiology)
enable_strict_fortran(stochastic_seir_ensemble_example)

add_executable(age_structured_sir_example app/age_structured_sir_example.f90)
target_link_libraries(age_structured_sir_example PRIVATE epidemiology)
enable_strict_fortran(age_structured_sir_example)

add_executable(growth_rate_example app/growth_rate_example.f90)
target_link_libraries(growth_rate_example PRIVATE epidemiology)
enable_strict_fortran(growth_rate_example)

add_executable(renewal_rt_example app/renewal_rt_example.f90)
target_link_libraries(renewal_rt_example PRIVATE epidemiology)
enable_strict_fortran(renewal_rt_example)

add_executable(negative_binomial_rt_example app/negative_binomial_rt_example.f90)
target_link_libraries(negative_binomial_rt_example PRIVATE epidemiology)
enable_strict_fortran(negative_binomial_rt_example)

add_executable(reporting_delay_nowcast_example app/reporting_delay_nowcast_example.f90)
target_link_libraries(reporting_delay_nowcast_example PRIVATE epidemiology)
enable_strict_fortran(reporting_delay_nowcast_example)

add_executable(nowcasted_rt_pipeline_example app/nowcasted_rt_pipeline_example.f90)
target_link_libraries(nowcasted_rt_pipeline_example PRIVATE epidemiology)
enable_strict_fortran(nowcasted_rt_pipeline_example)

add_executable(nowcast_rt_uncertainty_example app/nowcast_rt_uncertainty_example.f90)
target_link_libraries(nowcast_rt_uncertainty_example PRIVATE epidemiology)
enable_strict_fortran(nowcast_rt_uncertainty_example)

add_executable(generation_interval_uncertainty_example app/generation_interval_uncertainty_example.f90)
target_link_libraries(generation_interval_uncertainty_example PRIVATE epidemiology)
enable_strict_fortran(generation_interval_uncertainty_example)

add_executable(generation_interval_fit_example app/generation_interval_fit_example.f90)
target_link_libraries(generation_interval_fit_example PRIVATE epidemiology)
enable_strict_fortran(generation_interval_fit_example)

add_executable(fitted_generation_interval_uncertainty_example app/fitted_generation_interval_uncertainty_example.f90)
target_link_libraries(fitted_generation_interval_uncertainty_example PRIVATE epidemiology)
enable_strict_fortran(fitted_generation_interval_uncertainty_example)

if(BUILD_TESTING)
    add_executable(test_factorial test/test_factorial.f90)
    target_link_libraries(test_factorial PRIVATE fortran_numerics)
    enable_strict_fortran(test_factorial)
    add_test(NAME factorial COMMAND test_factorial)

    add_executable(test_linear_solve test/test_linear_solve.f90)
    target_link_libraries(test_linear_solve PRIVATE fortran_numerics)
    enable_strict_fortran(test_linear_solve)
    add_test(NAME linear_solve COMMAND test_linear_solve)

    add_executable(test_sir test/test_sir.f90)
    target_link_libraries(test_sir PRIVATE epidemiology)
    enable_strict_fortran(test_sir)
    add_test(NAME sir COMMAND test_sir)

    add_executable(test_seir test/test_seir.f90)
    target_link_libraries(test_seir PRIVATE epidemiology)
    enable_strict_fortran(test_seir)
    add_test(NAME seir COMMAND test_seir)

    add_executable(test_seirs test/test_seirs.f90)
    target_link_libraries(test_seirs PRIVATE epidemiology)
    enable_strict_fortran(test_seirs)
    add_test(NAME seirs COMMAND test_seirs)

    add_executable(test_stochastic_sir test/test_stochastic_sir.f90)
    target_link_libraries(test_stochastic_sir PRIVATE epidemiology)
    enable_strict_fortran(test_stochastic_sir)
    add_test(NAME stochastic_sir COMMAND test_stochastic_sir)

    add_executable(test_stochastic_seir test/test_stochastic_seir.f90)
    target_link_libraries(test_stochastic_seir PRIVATE epidemiology)
    enable_strict_fortran(test_stochastic_seir)
    add_test(NAME stochastic_seir COMMAND test_stochastic_seir)

    add_executable(test_age_structured_sir test/test_age_structured_sir.f90)
    target_link_libraries(test_age_structured_sir PRIVATE epidemiology)
    enable_strict_fortran(test_age_structured_sir)
    add_test(NAME age_structured_sir COMMAND test_age_structured_sir)

    add_executable(test_growth_inference test/test_growth_inference.f90)
    target_link_libraries(test_growth_inference PRIVATE epidemiology)
    enable_strict_fortran(test_growth_inference)
    add_test(NAME growth_inference COMMAND test_growth_inference)

    add_executable(test_renewal_rt test/test_renewal_rt.f90)
    target_link_libraries(test_renewal_rt PRIVATE epidemiology)
    enable_strict_fortran(test_renewal_rt)
    add_test(NAME renewal_rt COMMAND test_renewal_rt)

    add_executable(test_negative_binomial_rt test/test_negative_binomial_rt.f90)
    target_link_libraries(test_negative_binomial_rt PRIVATE epidemiology)
    enable_strict_fortran(test_negative_binomial_rt)
    add_test(NAME negative_binomial_rt COMMAND test_negative_binomial_rt)

    add_executable(test_reporting_delay_nowcast test/test_reporting_delay_nowcast.f90)
    target_link_libraries(test_reporting_delay_nowcast PRIVATE epidemiology)
    enable_strict_fortran(test_reporting_delay_nowcast)
    add_test(NAME reporting_delay_nowcast COMMAND test_reporting_delay_nowcast)

    add_executable(test_nowcasted_rt_pipeline test/test_nowcasted_rt_pipeline.f90)
    target_link_libraries(test_nowcasted_rt_pipeline PRIVATE epidemiology)
    enable_strict_fortran(test_nowcasted_rt_pipeline)
    add_test(NAME nowcasted_rt_pipeline COMMAND test_nowcasted_rt_pipeline)

    add_executable(test_nowcast_rt_uncertainty test/test_nowcast_rt_uncertainty.f90)
    target_link_libraries(test_nowcast_rt_uncertainty PRIVATE epidemiology)
    enable_strict_fortran(test_nowcast_rt_uncertainty)
    add_test(NAME nowcast_rt_uncertainty COMMAND test_nowcast_rt_uncertainty)

    add_executable(test_generation_interval_uncertainty test/test_generation_interval_uncertainty.f90)
    target_link_libraries(test_generation_interval_uncertainty PRIVATE epidemiology)
    enable_strict_fortran(test_generation_interval_uncertainty)
    add_test(NAME generation_interval_uncertainty COMMAND test_generation_interval_uncertainty)

    add_executable(test_generation_interval_fit test/test_generation_interval_fit.f90)
    target_link_libraries(test_generation_interval_fit PRIVATE epidemiology)
    enable_strict_fortran(test_generation_interval_fit)
    add_test(NAME generation_interval_fit COMMAND test_generation_interval_fit)

    add_executable(test_fitted_generation_interval_uncertainty test/test_fitted_generation_interval_uncertainty.f90)
    target_link_libraries(test_fitted_generation_interval_uncertainty PRIVATE epidemiology)
    enable_strict_fortran(test_fitted_generation_interval_uncertainty)
    add_test(NAME fitted_generation_interval_uncertainty COMMAND test_fitted_generation_interval_uncertainty)
endif()
