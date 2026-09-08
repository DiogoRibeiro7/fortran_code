# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The scientific core now includes deterministic and stochastic compartmental models, intervention schedules, demographic turnover, age-structured transmission, contact-matrix reproduction numbers, early-epidemic growth inference, and renewal-equation estimation of time-varying reproduction numbers.

## Models

The repository includes deterministic SIR, SEIR, and SEIRS models together with stochastic SIR and SEIR continuous-time Markov chains simulated using Gillespie's direct method.

It also supports age-structured transmission through

\[
\lambda_i=\beta q_i\sum_j C_{ij}z_j\frac{I_j}{N_j},
\]

with basic reproduction number

\[
R_0=\rho(K),
\]

where \(K\) is the next-generation matrix.

## Early epidemic growth and generation intervals

For a generation-interval distribution \(w(\tau)\), exponential growth rate \(r\), and reproduction number \(R\), the code evaluates

\[
1=R\int_0^\infty e^{-r\tau}w(\tau)\,d\tau.
\]

The `growth_inference_m` module supports generation-interval moments, Euler-Lotka conversion between \(R\) and \(r\), doubling/halving times, log-linear growth fitting, and mechanistic SEIR generation intervals.

## Renewal-equation Rt inference

For observed incidence \(I_t\), the renewal model uses total infectiousness

\[
\Lambda_t=\sum_{s\ge1}I_{t-s}w_s
\]

and assumes

\[
I_t\mid R_t,\Lambda_t\sim\operatorname{Poisson}(R_t\Lambda_t).
\]

Within a sliding window where \(R_t\) is treated as constant, a Gamma prior

\[
R\sim\operatorname{Gamma}(a,b)
\]

with shape \(a\) and rate \(b\) gives the conjugate posterior

\[
R\mid I\sim\operatorname{Gamma}\left(a+\sum I_t,\;b+\sum\Lambda_t\right).
\]

The `renewal_rt_m` module provides:

- infectiousness convolution from incidence and discrete generation-interval weights;
- sliding-window posterior estimation;
- posterior mean and standard deviation;
- numerically evaluated Gamma credible intervals;
- explicit validity flags for windows without sufficient infectiousness.

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

## Run examples

```bash
./build/sir_example
./build/seir_intervention_example
./build/seirs_vaccination_example
./build/stochastic_sir_example
./build/stochastic_seir_ensemble_example
./build/age_structured_sir_example
./build/growth_rate_example
./build/renewal_rt_example
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
7. renewal-equation Rt inference and likelihood-based estimation;
8. uncertainty quantification, sensitivity analysis, and spatial/metapopulation models;
9. import/export interfaces for real outbreak data;
10. optional Python interoperability while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`, while stochastic compartment counts use `int64`. Model state, parameters, generation-interval distributions, numerical solvers, stochastic event dynamics, and statistical inference objects are represented explicitly and tested independently.
