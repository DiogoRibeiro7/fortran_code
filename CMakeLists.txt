cmake_minimum_required(VERSION 3.25)

project(fortran_epidemiology VERSION 0.19.0 LANGUAGES Fortran)

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

function(add_epidemiology_example target_name)
    add_executable(${target_name} app/${target_name}.f90)
    target_link_libraries(${target_name} PRIVATE epidemiology)
    enable_strict_fortran(${target_name})
endfunction()

function(add_epidemiology_test test_name)
    add_executable(test_${test_name} test/test_${test_name}.f90)
    target_link_libraries(test_${test_name} PRIVATE epidemiology)
    enable_strict_fortran(test_${test_name})
    add_test(NAME ${test_name} COMMAND test_${test_name})
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
    src/interval_censored_generation_interval_m.f90
    src/renewal_rt_m.f90
    src/negative_binomial_rt_m.f90
    src/reporting_delay_nowcast_m.f90
    src/nowcasted_rt_pipeline_m.f90
    src/nowcast_rt_uncertainty_m.f90
    src/generation_interval_uncertainty_m.f90
    src/fitted_generation_interval_uncertainty_m.f90
    src/interval_censored_generation_interval_uncertainty_m.f90
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

set(EPIDEMIOLOGY_EXAMPLES
    sir_example
    seir_intervention_example
    seirs_vaccination_example
    stochastic_sir_example
    stochastic_seir_ensemble_example
    age_structured_sir_example
    growth_rate_example
    renewal_rt_example
    negative_binomial_rt_example
    reporting_delay_nowcast_example
    nowcasted_rt_pipeline_example
    nowcast_rt_uncertainty_example
    generation_interval_uncertainty_example
    generation_interval_fit_example
    fitted_generation_interval_uncertainty_example
    interval_censored_generation_interval_example
    interval_censored_generation_interval_uncertainty_example
    mixed_censored_generation_interval_example
)
foreach(example_name IN LISTS EPIDEMIOLOGY_EXAMPLES)
    add_epidemiology_example(${example_name})
endforeach()

if(BUILD_TESTING)
    add_executable(test_factorial test/test_factorial.f90)
    target_link_libraries(test_factorial PRIVATE fortran_numerics)
    enable_strict_fortran(test_factorial)
    add_test(NAME factorial COMMAND test_factorial)

    add_executable(test_linear_solve test/test_linear_solve.f90)
    target_link_libraries(test_linear_solve PRIVATE fortran_numerics)
    enable_strict_fortran(test_linear_solve)
    add_test(NAME linear_solve COMMAND test_linear_solve)

    set(EPIDEMIOLOGY_TESTS
        sir
        seir
        seirs
        stochastic_sir
        stochastic_seir
        age_structured_sir
        growth_inference
        renewal_rt
        negative_binomial_rt
        reporting_delay_nowcast
        nowcasted_rt_pipeline
        nowcast_rt_uncertainty
        generation_interval_uncertainty
        generation_interval_fit
        fitted_generation_interval_uncertainty
        interval_censored_generation_interval
        interval_censored_generation_interval_uncertainty
        mixed_censored_generation_interval
    )
    foreach(test_name IN LISTS EPIDEMIOLOGY_TESTS)
        add_epidemiology_test(${test_name})
    endforeach()
endif()
