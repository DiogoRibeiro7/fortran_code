# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The current scientific core includes deterministic SIR and SEIR models together with piecewise time-varying transmission for intervention scenarios.

## Current epidemiology core

For a closed population of size \(N\), the SIR model is

\[
\frac{dS}{dt}=-\beta\frac{SI}{N},
\qquad
\frac{dI}{dt}=\beta\frac{SI}{N}-\gamma I,
\qquad
\frac{dR}{dt}=\gamma I.
\]

The SEIR extension introduces a latent compartment \(E\):

\[
\frac{dS}{dt}=-\beta\frac{SI}{N},
\qquad
\frac{dE}{dt}=\beta\frac{SI}{N}-\sigma E,
\]

\[
\frac{dI}{dt}=\sigma E-\gamma I,
\qquad
\frac{dR}{dt}=\gamma I.
\]

For both models, under the current assumptions,

\[
R_0=\frac{\beta}{\gamma}.
\]

The implementation currently provides:

- typed SIR and SEIR model parameters and compartment states;
- incidence \(\beta SI/N\);
- \(R_0\) calculation;
- fourth-order Runge-Kutta integration;
- complete deterministic epidemic simulation;
- piecewise transmission schedules for interventions through \(\beta(t)\);
- CSV-style output suitable for plotting or downstream analysis;
- tests for population conservation, non-negative compartments, early epidemic growth, exposed-compartment dynamics, intervention effects, and reproduction-number calculations.

## Build

Requirements:

- CMake 3.25 or newer;
- a modern Fortran compiler such as GNU Fortran;
- BLAS/LAPACK for the retained numerical examples.

On Debian or Ubuntu:

```bash
sudo apt-get install gfortran cmake libblas-dev liblapack-dev
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

## Run epidemic simulations

Baseline SIR:

```bash
./build/sir_example
```

SEIR with a transmission intervention:

```bash
./build/seir_intervention_example
```

The SEIR scenario starts with \(\beta=0.30\), reduces transmission to 35% of baseline on day 20, and relaxes it to 70% of baseline on day 60. Output contains

```text
day,beta,susceptible,exposed,infectious,recovered
```

so trajectories can be redirected directly to a data file.

## Repository layout

```text
src/   reusable epidemiological and numerical modules
app/   executable simulations and examples
test/  scientific and numerical regression tests
```

## Scientific roadmap

The intended development sequence is:

1. deterministic SIR and SEIR models;
2. time-varying transmission and intervention functions;
3. vaccination, births, deaths, waning immunity, and endemic models;
4. stochastic SIR/SEIR simulation using Gillespie-style event dynamics;
5. age-structured and contact-matrix models;
6. estimation of growth rates, generation intervals, and reproduction numbers;
7. likelihood-based parameter estimation and uncertainty quantification;
8. ensemble simulation and sensitivity analysis;
9. import/export interfaces for real outbreak data;
10. optional Python interoperability for visualization and analysis while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`. Model state and parameters are represented by derived types. Numerical integration and epidemiological equations are separated from executable programs so model behaviour can be tested independently.

Time-varying interventions are represented as transmission schedules rather than embedded directly in the differential equations. This keeps the biological model separate from policy or scenario assumptions and makes alternative intervention histories easy to compare.

The repository still contains a few generic numerical examples from its earlier life. They are retained temporarily because they provide tested numerical infrastructure; future PRs will either repurpose or remove them as the epidemiology library becomes self-contained.
