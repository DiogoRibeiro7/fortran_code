# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The current scientific core includes deterministic SIR, SEIR, and SEIRS models, time-varying transmission interventions, vaccination, waning immunity, and demographic turnover.

## Current epidemiology core

For a closed population of size \(N\), the SIR model is

\[
\frac{dS}{dt}=-\beta\frac{SI}{N},
\qquad
\frac{dI}{dt}=\beta\frac{SI}{N}-\gamma I,
\qquad
\frac{dR}{dt}=\gamma I.
\]

The SEIR extension introduces a latent compartment \(E\). The SEIRS model adds demographic turnover, vaccination, and waning immunity:

\[
\frac{dS}{dt}=\mu N-\beta\frac{SI}{N}-\nu S+\omega R-\mu S,
\]

\[
\frac{dE}{dt}=\beta\frac{SI}{N}-(\sigma+\mu)E,
\]

\[
\frac{dI}{dt}=\sigma E-(\gamma+\mu)I,
\]

\[
\frac{dR}{dt}=\gamma I+\nu S-(\omega+\mu)R.
\]

Here \(\mu\) is the natural mortality rate, with births entering at rate \(\mu N\), \(\nu\) is the susceptible vaccination rate, and \(\omega\) is the rate of waning immunity.

For the demographic SEIRS model,

\[
R_0=\frac{\beta\sigma}{(\sigma+\mu)(\gamma+\mu)}.
\]

At the vaccinated disease-free equilibrium, the susceptible fraction is

\[
\frac{S^*}{N}=\frac{\omega+\mu}{\omega+\mu+\nu},
\]

so the vaccination-adjusted control reproduction number is

\[
R_c=R_0\frac{S^*}{N}.
\]

The implementation currently provides:

- typed SIR, SEIR, and SEIRS model parameters and compartment states;
- standard-incidence transmission \(\beta SI/N\);
- basic and vaccination-adjusted reproduction-number calculations;
- fourth-order Runge-Kutta integration;
- complete deterministic epidemic and endemic simulations;
- piecewise transmission schedules for interventions through \(\beta(t)\);
- vaccination, waning immunity, births, and natural deaths;
- CSV-style output suitable for plotting or downstream analysis;
- tests for population conservation, non-negative compartments, early epidemic growth, intervention effects, waning flows, demographic balance, and reproduction-number calculations.

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

SEIRS with demographic turnover, waning immunity, and vaccination:

```bash
./build/seirs_vaccination_example
```

The SEIRS example compares a no-vaccination endemic scenario with a vaccinated scenario over ten years and reports \(R_0\) and \(R_c\) before emitting the paired trajectories.

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

Time-varying interventions are represented as transmission schedules rather than embedded directly in the differential equations. Demographic SEIRS dynamics use balanced births and natural deaths so that total population remains constant when initialized at the declared population size.

The repository still contains a few generic numerical examples from its earlier life. They are retained temporarily because they provide tested numerical infrastructure; future PRs will either repurpose or remove them as the epidemiology library becomes self-contained.
