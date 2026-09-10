# Roadmap

This roadmap describes the planned scientific and engineering evolution of `fortran-epidemiology` beyond the current `0.19.0` baseline.

The project has already grown from a small modern-Fortran numerical codebase into a reproducible epidemiological modelling and surveillance-inference library. Future work should therefore prioritize methodological coherence, uncertainty propagation, real-data usability, and numerical reliability over adding isolated model variants.

## Guiding principles

Development should follow a few explicit rules:

1. **Statistical assumptions must be visible.** New inference routines should expose the observation model, parameterization, uncertainty semantics, and validity conditions rather than hiding them behind generic APIs.
2. **Uncertainty should be propagated, not silently conditioned away.** Whenever a downstream quantity depends on an estimated upstream object, the roadmap should prefer principled propagation over plug-in shortcuts unless the latter are clearly labelled.
3. **Numerical stability is part of scientific correctness.** Tail probabilities, likelihood differences, root solvers, matrix calculations, and random generators should be tested against analytic or independent references.
4. **Deterministic and stochastic implementations should cross-check one another whenever possible.** Large-population stochastic means, analytic reproduction numbers, and limiting cases provide important scientific regression targets.
5. **One coherent capability per pull request.** Feature branches should remain narrow and merge only after strict CI and scientific regression tests pass.
6. **Real-data interfaces should not weaken the modelling kernel.** The Fortran library should remain usable independently of Python, R, notebooks, or external services.

## Current baseline: 0.19.0

The current codebase already includes:

- deterministic SIR, SEIR, and SEIRS models;
- intervention schedules and vaccination/demography;
- stochastic SIR and SEIR simulation with Gillespie's direct method;
- reproducible stochastic ensembles and fadeout analysis;
- age-structured transmission and next-generation-matrix `R0`;
- generation-interval and early-growth inference through Euler-Lotka;
- Poisson renewal-equation `Rt` inference;
- negative-binomial `Rt` inference with profile likelihood;
- reporting-delay nowcasting;
- Monte Carlo propagation of reporting-delay uncertainty into `Rt`;
- sensitivity analysis for uncertain generation-interval PMFs;
- exact-pair Gamma generation-interval fitting;
- bootstrap propagation of fitted generation-interval uncertainty;
- interval-censored generation-interval fitting;
- bootstrap uncertainty propagation for interval-censored pair records;
- mixed finite-interval and right-censored generation-interval likelihoods;
- stable Gamma CDF, survival, and quantile infrastructure;
- portable seeded random-number generation;
- strict CMake/CTest/GitHub Actions validation.

The current strongest end-to-end surveillance pipeline can therefore be viewed schematically as

\[
\text{transmission-pair information}
\rightarrow
\text{generation interval}
\rightarrow
\text{reporting-delay model}
\rightarrow
\text{latent incidence}
\rightarrow
\text{renewal infectiousness}
\rightarrow
R_t.
\]

The next releases should make every uncertain link in that chain estimable from realistic data.

---

# Near-term roadmap

## 0.20.x — Complete mixed-censoring uncertainty propagation

### Priority: very high

The point-estimate generation-interval model now handles finite intervals and right-censored observations, but the bootstrap uncertainty pipeline still assumes finite censoring windows only.

### Planned work

- Extend pair-record bootstrapping to a mixed-censoring record type.
- Preserve censoring type when resampling.
- Refit the mixed Gamma likelihood inside every bootstrap replicate.
- Propagate the resulting renewal-weight uncertainty into the existing reporting-delay/`Rt` Monte Carlo pipeline.
- Report convergence rates separately for finite-only and mixed-censoring bootstrap regimes.
- Add regression tests showing that heavier right censoring increases uncertainty in the fitted generation-time tail.

### Acceptance criteria

- Exact fixed-seed reproducibility.
- All-finite input reproduces the existing interval-censored bootstrap path.
- Mixed-censoring bootstrap diagnostics expose shape, rate, mean generation interval, and tail-mass uncertainty.
- `Rt` summaries are based only on successful fitted draws, with explicit effective draw counts.

---

## 0.21.x — Uncertain infector identity

### Motivation

Real transmission-pair datasets often contain several plausible infectors for one infectee. Treating the selected infector as known can bias the generation-interval distribution.

### Planned model

For infectee `i` with candidate infectors `j`, use candidate probabilities or likelihood weights

\[
\pi_{ij},\qquad \sum_j\pi_{ij}=1,
\]

and marginalize rather than hard-assign:

\[
L_i(k,\beta)
=
\sum_j
\pi_{ij}
P(L_{ij}<G\le U_{ij}\mid k,\beta).
\]

Right-censored candidate links should use survival contributions analogously.

### Planned work

- Candidate-transmission record structure.
- Mixture likelihood with numerically stable log-sum-exp evaluation.
- Optional user-supplied prior weights on candidate infectors.
- Diagnostics for posterior/normalized candidate-link weights.
- Bootstrap resampling at the infectee level, not the candidate-link level.

### Acceptance criteria

- Single-candidate records reduce exactly to current censoring likelihoods.
- Equal duplicate candidates do not change the likelihood.
- Extreme candidate probabilities recover the corresponding hard-assignment limit.

---

## 0.22.x — Observation-window and truncation corrections for transmission-pair studies

### Motivation

Transmission pairs are generally sampled inside finite study windows. Long generation intervals can be preferentially excluded near the study boundaries, producing backward/forward truncation bias.

### Planned work

- Explicit study-start and study-end dates.
- Conditional likelihood correcting for pair eligibility under the observation window.
- Tests demonstrating bias in naive fitting under strong truncation.
- Comparison of naive, interval-censored, and truncation-corrected fits on synthetic data.
- Bootstrap uncertainty respecting study-window eligibility.

### Scientific target

Make generation-interval inference depend on the actual ascertainment process rather than only on observed pair intervals.

---

## 0.23.x — Estimate reporting delays from reporting triangles

### Motivation

The current nowcasting layer treats the reporting-delay PMF as known. In operational surveillance this distribution should usually be estimated from historical revisions.

### Planned data structure

A reporting triangle containing cumulative or incremental reports by event date and report date.

### Planned models

Initial implementation should support a transparent discrete delay model:

\[
D\sim p_0,p_1,\ldots,p_K,
\]

with likelihood derived from historical reporting increments.

Later extensions may include:

- weekday effects;
- calendar-time-varying delay distributions;
- hierarchical smoothing;
- right truncation of recent reporting cohorts.

### Acceptance criteria

- Recovery of known delay probabilities in simulated triangles.
- Correct normalization and truncation handling.
- Plug-compatible output with `reporting_delay_nowcast_m`.
- Bootstrap or posterior uncertainty that can feed the `Rt` Monte Carlo layer.

---

## 0.24.x — Joint reporting-delay uncertainty

Once reporting-delay estimation exists, the current assumption of a fixed delay PMF should be removed from the full pipeline.

Each Monte Carlo iteration should be able to draw

\[
p^{(m)}(D),
\]

then use that draw to compute reporting completeness, latent incidence, and `Rt`.

The target uncertainty chain becomes

\[
\text{pair study uncertainty}
+
\text{generation interval}
+
\text{reporting-delay distribution}
+
\text{latent incidence}
+
\text{renewal }R_t.
\]

---

# Medium-term roadmap

## 0.25.x — Observation-model integration

The project currently has Poisson and negative-binomial renewal models, while the full nowcasting uncertainty pipeline is primarily coupled to the conjugate Poisson model.

### Planned work

- Joint nowcast-plus-NB2 `Rt` inference.
- Propagation of reporting and generation-time uncertainty under overdispersed observations.
- Estimation or profiling of NB dispersion rather than requiring it to be fixed in all workflows.
- Comparison between Poisson and NB2 model adequacy on simulated overdispersed data.

### Longer-term option

Consider Poisson-lognormal or state-space observation models only if they add a clearly justified capability rather than model proliferation.

---

## 0.26.x — Time-varying `Rt` smoothing and state-space inference

Sliding windows are transparent and robust, but they impose piecewise-constant transmission over a chosen window.

### Planned work

- Random-walk model on `log Rt`.
- Penalized/state-space estimation.
- Separation of process evolution and surveillance observation noise.
- Simulation-based calibration against known trajectories.
- Comparison with existing sliding-window estimates.

### Principle

Keep the current renewal likelihood visible. Smoothing should be an explicit prior/process layer, not a black-box filter replacing the epidemiological model.

---

## 0.27.x — Susceptible depletion and renewal models

The current renewal formulation is primarily suited to settings where susceptible depletion is negligible or implicitly absorbed into effective transmission.

### Planned work

- Optional susceptible fraction correction.
- Relationship between renewal `Rt`, effective reproduction number, and compartmental susceptible depletion.
- Consistency tests against deterministic SIR/SEIR trajectories.

---

## 0.28.x — Importation and local transmission

### Planned work

Distinguish imported incidence from locally generated incidence:

\[
I_t=I_t^{\mathrm{local}}+I_t^{\mathrm{imported}}.
\]

Only locally infectious pressure should be attributed to local reproduction dynamics according to an explicitly chosen model.

This is important for low-incidence settings and regional surveillance.

---

## 0.29.x — Spatial and metapopulation transmission

### Planned scientific layer

For regions `i,j`, introduce mobility/contact coupling

\[
\lambda_i(t)
=
\sum_j M_{ij}(t)\,\mathcal I_j(t),
\]

with explicit definitions for mobility and infectiousness.

### Planned capabilities

- Deterministic metapopulation SIR/SEIR models.
- Region-specific transmission rates.
- Mobility matrices.
- Imported/local force-of-infection decomposition.
- Spectral reproduction-number calculations for structured systems.
- Stochastic metapopulation simulation where computationally reasonable.

---

# Advanced inference roadmap

## Uncertain infection times from exposure histories

Move beyond precomputed generation-interval censoring windows and infer infection times from exposure intervals directly.

Potential hierarchy:

\[
T_i^{\mathrm{infection}}
\sim
p(T_i\mid\text{exposure history}),
\]

followed by generation intervals derived from latent infection times.

This is substantially more complex than interval censoring and should be added only after candidate-infector and study-truncation models are stable.

## Serial interval versus generation interval

Add explicit symptom-onset models so the package can distinguish:

- generation interval;
- serial interval;
- incubation period;
- onset-to-transmission timing.

The library should avoid treating these quantities as interchangeable.

## Time-varying generation intervals

Interventions and behavior changes can alter generation intervals as well as transmission intensity.

Potential model:

\[
w_t(s)
\]

rather than a single fixed `w(s)`.

This requires careful identifiability work because simultaneous changes in `Rt` and `w_t` can be difficult to distinguish from incidence alone.

## Hierarchical generation-interval models

Potential grouping variables include:

- epidemic phase;
- variant/pathogen strain;
- age group;
- vaccination status;
- setting;
- calendar period.

A hierarchical model should be introduced only with strong regression coverage and clear shrinkage semantics.

---

# Deterministic and stochastic modelling roadmap

## Additional intervention mechanics

Possible additions include:

- vaccination with imperfect efficacy;
- waning vaccine protection;
- treatment reducing infectious duration;
- isolation/quarantine compartments;
- importation seeding;
- seasonally varying transmission.

Each extension should have a clear epidemiological quantity or validation target rather than simply increasing compartment count.

## Adaptive and higher-order ODE integration

Current deterministic models use explicit RK4. Potential improvements:

- adaptive Runge-Kutta integration;
- user-defined tolerances;
- event handling at intervention boundaries;
- convergence comparisons against RK4 reference trajectories.

This should remain optional because fixed-step RK4 is simple, reproducible, and useful for testing.

## Faster stochastic simulation

For large populations or high event rates:

- tau-leaping;
- hybrid stochastic/deterministic methods;
- bounded error diagnostics against exact Gillespie simulation.

Exact SSA remains the scientific reference implementation.

---

# Real-data and reproducibility roadmap

## Stable data input layer

Add explicit parsers/interfaces for:

- daily incidence CSV files;
- reporting triangles;
- transmission-pair tables;
- censoring metadata;
- contact matrices;
- mobility matrices.

Inputs should be validated aggressively with informative failure modes.

## Reproducible command-line workflows

Provide small CLI-style executables that can run complete analyses without editing source code.

Potential workflows:

```text
fit-generation-interval
estimate-rt
nowcast-incidence
estimate-nowcasted-rt
simulate-seir
simulate-stochastic-seir
```

Configuration should eventually be file-driven rather than compiled into examples.

## Reference datasets

Add small, redistributable synthetic/reference datasets for integration tests and examples.

Never place confidential or restricted public-health data in the repository.

## End-to-end synthetic validation

Create a simulation-recovery suite:

1. simulate an epidemic with known transmission parameters;
2. generate transmission-pair observations;
3. impose censoring and ascertainment;
4. generate reporting delays;
5. infer generation intervals;
6. nowcast latent incidence;
7. estimate `Rt`;
8. compare estimates with known truth.

This should become one of the strongest scientific validation assets in the repository.

---

# Numerical roadmap

## Distribution functions

Continue strengthening numerical distribution infrastructure with explicit tests against analytic cases and independent references.

Candidate additions:

- log-CDF and log-survival functions;
- stable log-difference helpers;
- Beta distribution functions;
- normal quantiles;
- negative-binomial random variates;
- robust multinomial sampling.

## Optimization

The project currently uses problem-specific root finding and Nelder-Mead where appropriate.

Future work should consider a small reusable optimization layer providing:

- bracketed scalar roots;
- bounded one-dimensional minimization;
- Nelder-Mead diagnostics;
- convergence codes rather than only logical flags;
- iteration and objective histories where useful.

Do not introduce a general optimizer abstraction unless it simplifies multiple existing inference modules.

## Random-number generation

Planned checks:

- documented RNG algorithm and period;
- reproducibility contract across supported compilers;
- basic moment tests for Gamma, Poisson, and Dirichlet draws;
- optional stream splitting for parallel Monte Carlo.

---

# Software-engineering roadmap

## API stabilization toward 1.0

Before declaring a stable `1.0.0` API:

- review naming conventions across all modules;
- identify duplicated summary types;
- standardize validity/convergence reporting;
- standardize shape/rate terminology;
- separate model parameters, fit objects, and Monte Carlo summaries consistently;
- document which APIs are intended as public/stable.

## Error handling

The current library uses `error stop` for invalid caller inputs in many places.

Potential improvement:

- optional status-return APIs for applications that cannot terminate the process;
- structured error/status codes;
- retain strict fail-fast wrappers for examples/tests.

## Compiler matrix

Expand CI beyond GNU Fortran when practical:

- GNU Fortran;
- Intel ifx if accessible;
- LLVM/flang when sufficiently mature for the required features.

Scientific tests must remain compiler-independent.

## Sanitizers and runtime diagnostics

Where supported:

- bounds checks;
- floating-point exception checks;
- uninitialized-variable diagnostics;
- address/undefined-behavior sanitizers for mixed-language layers.

## Packaging

Investigate packaging through the Fortran Package Manager (`fpm`) while retaining CMake support if it continues to serve LAPACK and multi-language integration well.

The goal should be dual usability, not replacing one build system merely for fashion.

---

# Interoperability roadmap

## Python interface

A Python layer is useful for data handling and visualization, but it should remain a wrapper around the Fortran scientific kernel.

Preferred direction:

- stable C ABI or modern F2PY-compatible wrapper layer;
- NumPy arrays for trajectories and inference outputs;
- Python packaging only after the Fortran API is sufficiently stable.

Avoid restoring the previous repository's stale platform-specific Python extension artifact.

## R interface

Potentially useful for epidemiological and statistical workflows, especially generation-interval fitting and `Rt` estimation.

Possible implementation routes:

- C ABI plus `.Call`/`.Fortran` wrappers;
- generated shared library consumed by an R package.

This is lower priority than Python unless a concrete user workflow requires it.

---

# Documentation roadmap

## Methods documentation

Each major statistical method should eventually have a standalone note containing:

- model definition;
- parameterization;
- likelihood/posterior;
- numerical method;
- uncertainty interpretation;
- assumptions;
- known limitations;
- references where appropriate;
- regression/reference values.

## Architecture documentation

Add an architecture document showing the dependency structure among:

```text
compartmental models
stochastic simulation
random distributions
generation-interval inference
reporting-delay inference
renewal inference
joint uncertainty propagation
```

## Examples as scientific demonstrations

Examples should remain small and inspectable. They should demonstrate one scientifically meaningful concept each rather than becoming undocumented mini-applications.

---

# Release roadmap

## Toward 0.5

A `0.5.x` maturity milestone in the new roadmap sense should require:

- mixed-censoring bootstrap uncertainty;
- uncertain infector handling;
- reporting-delay estimation from data;
- joint delay-distribution uncertainty propagation;
- improved real-data input contracts;
- end-to-end synthetic recovery tests.

The project version is already numerically beyond `0.5`; this milestone therefore describes maturity, not a literal future semantic version number.

## Toward 1.0.0

`1.0.0` should mean a stable scientific API, not simply a large feature count.

Minimum criteria:

- documented public modules and stable naming conventions;
- deterministic, stochastic, surveillance, and generation-interval components covered by reference/regression tests;
- explicit uncertainty propagation for all major fitted upstream quantities;
- real-data interfaces with strong validation;
- stable compiler/build support;
- no known silent numerical-failure modes;
- documented statistical assumptions and limitations;
- changelog and roadmap kept current for every release.

---

# Explicit non-goals

The project should avoid several tempting but low-value directions unless a concrete scientific need appears:

- adding arbitrary compartment models solely to increase model count;
- opaque machine-learning prediction models without epidemiological structure;
- hidden automatic model selection;
- silently substituting serial intervals for generation intervals;
- reporting `Rt` point estimates without corresponding uncertainty semantics;
- platform-specific compiled binaries committed to source control;
- Python-first architecture that makes the Fortran core secondary.

---

# Immediate next sequence

The recommended development order from the current `0.19.0` baseline is:

1. mixed finite/right-censored generation-interval bootstrap into `Rt`;
2. uncertain infector identity;
3. study-window truncation correction;
4. reporting-triangle delay estimation;
5. reporting-delay distribution uncertainty propagation;
6. negative-binomial integration with the full nowcasting pipeline;
7. end-to-end synthetic simulation/recovery benchmark;
8. real-data interfaces;
9. spatial/metapopulation modelling;
10. API consolidation and interoperability work.

This order keeps the project centered on a single coherent goal: **transparent epidemiological inference in which the data-generating process and every important source of uncertainty remain visible from raw surveillance data to reproduction-number estimates.**
