# Modern Fortran numerical examples

A small collection of modern Fortran examples focused on numerical computing, typed modules, explicit interfaces, and reproducible builds.

The repository currently contains:

- a reusable factorial module using `iso_fortran_env` kinds;
- a typed wrapper around LAPACK `DGESV` for dense linear systems;
- small executable examples under `app/`;
- numerical regression tests under `test/`;
- CMake/CTest build and test automation;
- GitHub Actions CI using GNU Fortran and system BLAS/LAPACK.

## Requirements

- CMake 3.25 or newer
- a Fortran compiler with Fortran 2018 support, such as GNU Fortran
- BLAS and LAPACK development libraries

On Debian or Ubuntu:

```bash
sudo apt-get install gfortran cmake libblas-dev liblapack-dev
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
```

Debug builds enable strict GNU Fortran diagnostics and runtime checks when `gfortran` is used.

## Test

```bash
ctest --test-dir build --output-on-failure
```

The current tests check the factorial implementation and verify that the LAPACK example solves

\[
A x = b
\]

with the expected solution

\[
x = (2, 3, -1)^T.
\]

## Run the examples

```bash
./build/factorial_example
./build/for_loop_example
./build/lapack_example
```

## Layout

```text
app/   executable examples
src/   reusable Fortran modules
test/  numerical tests
```

## Numerical conventions

Floating-point code uses the standard `real64` kind from `iso_fortran_env` rather than compiler-specific declarations such as `real*8`. External LAPACK calls use an explicit interface, and inputs are copied before calling `DGESV` because LAPACK overwrites the matrix and right-hand side in place.

## Next steps

Future work can add a maintained Python interface, broader LAPACK coverage, property-based numerical tests, and additional numerical-method examples without coupling those concerns to the core Fortran build.
