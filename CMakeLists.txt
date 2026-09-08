cmake_minimum_required(VERSION 3.25)

project(fortran_epidemiology VERSION 0.2.0 LANGUAGES Fortran)

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
endif()
