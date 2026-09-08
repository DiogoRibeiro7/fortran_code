# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The current scientific core includes deterministic SIR, SEIR, and SEIRS models, time-varying transmission interventions, vaccination, waning immunity, demographic turnover, and stochastic SIR simulation with Gillespie's direct method.

## Current epidemiology core

For a closed population of size \(N\), the deterministic SIR model is

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
\qquad
\frac{dR}{dt}=\gamma I+\nu S-(\omega+\mu)R.
\]

For the demographic SEIRS model,

\[
R_0=\frac{\beta\sigma}{(\sigma+\mu)(\gamma+\mu)},
\qquad
R_c=R_0\frac{\omega+\mu}{\omega+\mu+\nu}.
\]

## Stochastic SIR

The stochastic SIR implementation treats infection and recovery as continuous-time Markov events with hazards

\[
a_{\mathrm{inf}}=\beta\frac{SI}{N},
\qquad
 a_{\mathrm{rec}}=\gamma I.
\]

Gillespie's direct method samples the next waiting time from an exponential distribution with rate \(a_0=a_{\mathrm{inf}}+a_{\mathrm{rec}}\), then selects infection or recovery proportionally to its hazard. Compartments remain integer-valued and total population is conserved exactly.

The stochastic layer uses an explicit reproducible pseudo-random number generator, so fixed seeds produce stable event trajectories across test runs.

The implementation currently provides:

- deterministic SIR, SEIR, and SEIRS models;
- standard-incidence transmission \(\beta SI/N\);
- basic and vaccination-adjusted reproduction numbers;
- fourth-order Runge-Kutta integration;
- time-varying transmission schedules;
- vaccination, waning immunity, births, and natural deaths;
- stochastic SIR simulation with Gillespie's direct method;
- reproducible seeded random-number streams;
- CSV/event-log output for plotting and downstream analysis;
- tests for deterministic invariants, event hazards, absorbing states, extinction paths, and seeded stochastic trajectories.

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

## Run simulations

```bash
./build/sir_example
./build/seir_intervention_example
./build/seirs_vaccination_example
./build/stochastic_sir_example
```

The stochastic executable emits

```text
event,time,susceptible,infectious,recovered
```

with one row per infection or recovery event.

## Repository layout

```text
src/   reusable epidemiological and numerical modules
app/   executable simulations and examples
test/  scientific and numerical regression tests
```

## Scientific roadmap

1. deterministic SIR and SEIR models;
2. time-varying transmission and intervention functions;
3. vaccination, births, deaths, waning immunity, and endemic models;
4. stochastic SIR/SEIR event dynamics;
5. age-structured and contact-matrix models;
6. estimation of growth rates, generation intervals, and reproduction numbers;
7. likelihood-based parameter estimation and uncertainty quantification;
8. ensemble simulation and sensitivity analysis;
9. import/export interfaces for real outbreak data;
10. optional Python interoperability for visualization and analysis while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`, while stochastic compartment counts use `int64`. Model state and parameters are represented by derived types. Numerical integration, stochastic event dynamics, and epidemiological equations are separated from executable programs so model behaviour can be tested independently.

Time-varying interventions are represented as transmission schedules rather than embedded directly in the differential equations. Demographic SEIRS dynamics use balanced births and natural deaths so that total population remains constant when initialized at the declared population size.
