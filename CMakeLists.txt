cmake_minimum_required(VERSION 3.25)

project(fortran_epidemiology VERSION 0.6.0 LANGUAGES Fortran)

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
    src/stochastic_sir_m.f90
    src/stochastic_seir_m.f90
)
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
endif()
