# Epidemiological modelling in modern Fortran

A compact epidemiological modelling codebase written in modern Fortran.

The project is developed as a numerically explicit, testable implementation of classical and computational epidemiology. The scientific core now includes deterministic and stochastic compartmental models, intervention schedules, demographic turnover, age-structured transmission, contact-matrix reproduction numbers, and early-epidemic growth inference from generation intervals.

## Models

The repository includes deterministic SIR, SEIR, and SEIRS models together with stochastic SIR and SEIR continuous-time Markov chains simulated using Gillespie's direct method.

It also supports age-structured transmission through

\[
\lambda_i
=
\beta q_i
\sum_j C_{ij}z_j\frac{I_j}{N_j},
\]

with basic reproduction number

\[
R_0=\rho(K),
\]

where \(K\) is the next-generation matrix and \(\rho(K)\) its spectral radius.

## Early epidemic growth and generation intervals

For a generation-interval distribution \(w(\tau)\), exponential growth rate \(r\), and reproduction number \(R\), the code evaluates the Euler-Lotka relation

\[
1
=
R\int_0^\infty e^{-r\tau}w(\tau)\,d\tau.
\]

The `growth_inference_m` module provides:

- arbitrary discrete generation-interval distributions;
- generation-interval mean and standard deviation;
- the Laplace transform of the generation interval;
- conversion from observed growth rate \(r\) to reproduction number \(R\);
- numerical inversion from \(R\) back to \(r\);
- doubling and halving times;
- log-linear fitting of an early exponential incidence curve;
- the mechanistic SEIR generation-interval density;
- analytic SEIR relations between \(r\), \(R_0\), \(\sigma\), and \(\gamma\).

For the standard latent-plus-infectious SEIR process, the generation interval has mean

\[
E[G]
=
\frac{1}{\sigma}+rac{1}{\gamma}.
\]

The early growth rate satisfies

\[
(r+\sigma)(r+\gamma)=\sigma\beta,
\]

so

\[
R_0
=
\frac{(r+\sigma)(r+\gamma)}{\sigma\gamma}.
\]

This is also recovered numerically by Euler-Lotka when the full SEIR generation-interval distribution is used.

A documented example uses a 3-day mean latent period, 5-day mean infectious period, and observed doubling time of 5 days. It gives approximately

```text
growth rate r                 = 0.138629 / day
generation-interval mean      = 8.0000 days
generation-interval SD        = 5.8309 days
R0 from Euler-Lotka           = 2.3973
R0 from analytic SEIR         = 2.3973
R0 from exp(r * mean GI)      = 3.0314
```

The difference in the last line is deliberate: using only the mean generation interval discards distributional information and can materially distort the inferred reproduction number.

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
10. optional Python interoperability while keeping the modelling kernel in Fortran.

## Numerical conventions

Floating-point calculations use `real64` from `iso_fortran_env`, while stochastic compartment counts use `int64`. Model state, parameters, generation-interval distributions, numerical solvers, and stochastic event dynamics are represented explicitly and tested independently.
