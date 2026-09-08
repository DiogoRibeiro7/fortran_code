# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The scientific core includes deterministic and stochastic compartmental models, intervention schedules, demographic turnover, age-structured transmission, generation-interval inference, renewal-equation Rt estimation, and overdispersed surveillance likelihoods.

## Models

The repository includes deterministic SIR, SEIR, and SEIRS models together with stochastic SIR and SEIR continuous-time Markov chains simulated using Gillespie's direct method. It also supports age-structured transmission through contact matrices and next-generation-matrix reproduction numbers.

## Early epidemic growth and generation intervals

For a generation-interval distribution \(w(\tau)\), exponential growth rate \(r\), and reproduction number \(R\), the code evaluates

\[
1=R\int_0^\infty e^{-r\tau}w(\tau)\,d\tau.
\]

The `growth_inference_m` module supports generation-interval moments, Euler-Lotka conversion between \(R\) and \(r\), doubling/halving times, log-linear growth fitting, and mechanistic SEIR generation intervals.

## Renewal-equation Rt inference

For observed incidence \(I_t\), total infectiousness is

\[
\Lambda_t=\sum_{s\ge1}I_{t-s}w_s.
\]

The Poisson renewal model assumes

\[
I_t\mid R_t,\Lambda_t\sim\operatorname{Poisson}(R_t\Lambda_t),
\]

and uses Gamma-Poisson conjugacy to obtain posterior means, standard deviations, and Gamma credible intervals for sliding-window estimates of \(R_t\).

## Negative-binomial observation model

Real surveillance counts are often more variable than Poisson. The `negative_binomial_rt_m` module adds an NB2 renewal likelihood with

\[
\mu_t=R_t\Lambda_t,
\qquad
\operatorname{Var}(I_t\mid\mu_t)=\mu_t+\frac{\mu_t^2}{k},
\]

where \(k>0\) is the dispersion parameter. Large \(k\) approaches the Poisson model, while smaller \(k\) represents stronger overdispersion.

Because the negative-binomial likelihood is no longer Gamma-conjugate in \(R_t\), the module estimates \(R_t\) by maximum likelihood on each window and reports profile-likelihood confidence intervals. This keeps the uncertainty calculation consistent with the chosen observation model instead of reusing the Poisson posterior formula.

The module provides:

- NB2 variance and log-probability calculations;
- fixed-dispersion window and sliding-window Rt estimation;
- exact likelihood maximization for \(R_t\) via a bracketed score solve;
- profile-likelihood confidence intervals;
- explicit invalid-window handling when incidence has no renewal support;
- the Poisson limit as dispersion \(k\to\infty\).

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
./build/negative_binomial_rt_example
```

The negative-binomial example compares strongly overdispersed and near-Poisson profile intervals for the same incidence series.

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
7. Poisson and negative-binomial renewal-equation Rt inference;
8. uncertainty quantification, sensitivity analysis, and spatial/metapopulation models;
9. import/export interfaces for real outbreak data;
10. optional Python interoperability while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`, while stochastic compartment counts use `int64`. Model state, parameters, numerical solvers, stochastic event dynamics, and statistical inference objects are represented explicitly and tested independently.
