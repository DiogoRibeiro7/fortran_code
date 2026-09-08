# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The current scientific core includes deterministic SIR, SEIR, and SEIRS models, time-varying interventions, vaccination, waning immunity, demographic turnover, stochastic SIR/SEIR simulation with Gillespie's direct method, reproducible ensembles, and age-structured transmission through contact matrices.

## Current epidemiology core

For a closed population of size \(N\), the deterministic SIR model is

\[
\frac{dS}{dt}=-\beta\frac{SI}{N},
\qquad
\frac{dI}{dt}=\beta\frac{SI}{N}-\gamma I,
\qquad
\frac{dR}{dt}=\gamma I.
\]

The deterministic SEIR extension adds a latent compartment \(E\), while the SEIRS model adds demographic turnover, vaccination, and waning immunity. For the demographic SEIRS model,

\[
R_0=\frac{\beta\sigma}{(\sigma+\mu)(\gamma+\mu)},
\qquad
R_c=R_0\frac{\omega+\mu}{\omega+\mu+\nu}.
\]

## Stochastic epidemic models

The stochastic SIR and SEIR implementations are continuous-time Markov chains simulated with Gillespie's direct method. For stochastic SEIR,

\[
a_{\mathrm{inf}}=\beta\frac{SI}{N},
\qquad
a_{\mathrm{prog}}=\sigma E,
\qquad
a_{\mathrm{rec}}=\gamma I.
\]

The stochastic layer uses a reproducible explicit pseudo-random number generator and supports seeded ensemble summaries including extinction, fadeout, major-outbreak probability, final-size moments, peak prevalence, and mean compartment values at a time horizon.

## Age-structured transmission

For population groups \(i=1,\ldots,m\), the force of infection is

\[
\lambda_i
=
\beta q_i
\sum_j
C_{ij} z_j\frac{I_j}{N_j},
\]

where \(C_{ij}\) is the per-person contact rate reported by group \(i\) with group \(j\), \(q_i\) is group-specific susceptibility, and \(z_j\) is group-specific infectiousness.

The group dynamics are

\[
\dot S_i=-\lambda_i S_i,
\qquad
\dot I_i=\lambda_i S_i-\gamma_i I_i,
\qquad
\dot R_i=\gamma_i I_i.
\]

At the disease-free equilibrium the next-generation matrix is

\[
K_{ij}
=
\frac{
\beta q_i N_i C_{ij} z_j
}{
N_j\gamma_j
},
\]

and the basic reproduction number is

\[
R_0=\rho(K),
\]

where \(\rho(K)\) is the spectral radius. The implementation computes this through LAPACK and exposes a reciprocity diagnostic based on

\[
N_iC_{ij}=N_jC_{ji}.
\]

The reference three-group example uses children, adults, and older adults. Its reciprocal baseline contact matrix gives approximately

```text
baseline R0     = 3.9184
intervention R0 = 2.2936
```

when contacts involving the child group are reduced while preserving reciprocity.

## Build

Requirements:

- CMake 3.25 or newer;
- a modern Fortran compiler such as GNU Fortran;
- BLAS and LAPACK.

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
./build/stochastic_seir_ensemble_example
./build/age_structured_sir_example
```

## Repository layout

```text
src/   reusable epidemiological and numerical modules
app/   executable simulations and examples
test/  scientific and numerical regression tests
```

## Scientific roadmap

1. deterministic SIR, SEIR, and SEIRS models;
2. time-varying transmission and intervention functions;
3. vaccination, births, deaths, waning immunity, and endemic dynamics;
4. stochastic SIR/SEIR event dynamics and ensembles;
5. age-structured and contact-matrix transmission;
6. generation intervals, growth rates, and reproduction-number estimation;
7. likelihood-based parameter estimation and uncertainty quantification;
8. metapopulation/spatial models and sensitivity analysis;
9. import/export interfaces for real outbreak data;
10. optional Python interoperability for visualization and analysis while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`, while stochastic compartment counts use `int64`. Model state and parameters are represented by derived types. Numerical integration, stochastic event dynamics, and epidemiological equations are separated from executable programs so model behaviour can be tested independently.

Time-varying interventions are represented as transmission schedules rather than embedded directly in the disease equations. Age-structured contact matrices use an explicit contact-rate convention and expose reciprocity diagnostics so heterogeneous transmission assumptions remain inspectable.
