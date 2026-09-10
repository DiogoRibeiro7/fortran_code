# Roadmap

This roadmap describes the planned scientific and engineering evolution of `fortran-epidemiology` beyond the current `0.21.0` baseline.

The project has grown from a small modern-Fortran numerical codebase into a reproducible epidemiological modelling and surveillance-inference library. Future work should prioritize methodological coherence, uncertainty propagation, realistic observation processes, real-data usability, and numerical reliability over adding isolated model variants.

## Guiding principles

Development should follow a few explicit rules:

1. **Statistical assumptions must be visible.** New inference routines should expose the observation model, parameterization, uncertainty semantics, and validity conditions rather than hiding them behind generic APIs.
2. **Uncertainty should be propagated, not silently conditioned away.** Whenever a downstream quantity depends on an estimated upstream object, prefer principled propagation over plug-in shortcuts unless the latter are clearly labelled.
3. **Numerical stability is part of scientific correctness.** Tail probabilities, likelihood differences, root solvers, matrix calculations, optimization, and random generators should be tested against analytic or independent references.
4. **Deterministic and stochastic implementations should cross-check one another whenever possible.** Large-population stochastic means, analytic reproduction numbers, and limiting cases provide strong scientific regression targets.
5. **One coherent capability per pull request.** Feature branches should remain narrow and merge only after strict CI and scientific regression tests pass.
6. **Real-data interfaces should not weaken the modelling kernel.** The Fortran library should remain usable independently of Python, R, notebooks, or external services.
7. **Observation processes matter.** Censoring, truncation, reporting delay, ascertainment, and uncertain transmission links should be modelled explicitly when they affect inference.

## Current baseline: 0.21.0

The current codebase includes:

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
- mixed-censoring bootstrap propagation into `Rt`;
- uncertain-infector generation-interval inference through candidate-mixture likelihoods;
- stable Gamma CDF, survival, and quantile infrastructure;
- portable seeded random-number generation;
- strict CMake/CTest/GitHub Actions validation;
- repository-level `CHANGELOG.md` and `ROADMAP.md`.

The strongest current surveillance/inference chain can be viewed schematically as

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

The next releases should make the observation and ascertainment mechanisms in that chain increasingly estimable from realistic data.

---

# Completed near-term milestones

## 0.20.0 — Mixed-censoring uncertainty propagation

**Status: completed.**

The point-estimate generation-interval model already handled finite intervals and right-censored observations. Version `0.20.0` extended uncertainty propagation so complete mixed-censoring records are resampled together.

Implemented:

- pair-record bootstrapping of `(L_i,U_i,delta_i)`;
- preservation of censoring type under resampling;
- mixed Gamma refitting inside each bootstrap replicate;
- propagation of fitted renewal-weight uncertainty through reporting-delay completion and renewal `Rt` inference;
- explicit requested-versus-converged bootstrap counts;
- fixed-seed regression coverage;
- dedicated example and methods documentation.

Scientific rule retained: bootstrap samples that do not contain enough finite information for a stable fit are skipped rather than repaired with invented data.

## 0.21.0 — Uncertain infector identity

**Status: completed at the point-estimate likelihood level.**

For infectee `i` with candidate infectors `j`, candidate identity is marginalized rather than hard-assigned:

\[
L_i(k,\beta)=\sum_j \tilde\pi_{ij}L_{ij}(k,\beta),
\qquad
\tilde\pi_{ij}=\frac{\pi_{ij}}{\sum_h\pi_{ih}}.
\]

Each candidate can contribute either a finite censoring probability or a right-censored survival probability.

Implemented:

- flat candidate-transmission representation with contiguous per-infectee ranges;
- externally supplied non-negative candidate/link weights;
- numerically stable log-sum-exp mixture likelihood;
- log-scale Gamma shape/rate optimization;
- single-candidate reduction to the existing censoring likelihood;
- duplicate-candidate invariance;
- hard-assignment limiting behavior under extreme candidate weights;
- example and methods documentation.

Remaining extensions for uncertain infectors are tracked below: posterior/normalized candidate diagnostics, infectee-level bootstrap propagation, and eventual joint models for link evidence.

---

# Near-term roadmap

## 0.22.x — Observation-window and truncation corrections for transmission-pair studies

### Priority: very high

Transmission pairs are generally observed inside finite study windows. Long generation intervals can be preferentially excluded near study boundaries, creating backward or forward truncation bias even when censoring itself is modelled correctly.

### Planned statistical model

Introduce explicit pair-eligibility conditioning. If a latent generation interval `G` would only enter the observed sample when it satisfies a study-specific eligibility event `A_i`, use a conditional likelihood of the form

\[
L_i(\theta\mid A_i)
=
\frac{P_\theta(\text{observed pair information}_i)}
     {P_\theta(A_i)}.
\]

The exact denominator must reflect the study design rather than a generic truncation correction.

### Planned work

- explicit study-start and study-end times;
- candidate-specific infector/infectee observation windows where needed;
- conditional likelihood correcting for pair eligibility;
- support for finite, right-censored, and uncertain-infector records where mathematically identifiable;
- synthetic studies demonstrating bias in naive fitting under strong truncation;
- comparison of naive, censoring-aware, and truncation-corrected fits;
- bootstrap uncertainty that respects study-window eligibility rather than resampling impossible records.

### Acceptance criteria

- no correction when the study window is effectively unbounded;
- recovery of known Gamma parameters in synthetic truncated samples;
- demonstrable reduction of long-interval bias near study boundaries;
- stable likelihood evaluation when eligibility probabilities are small;
- explicit failure when the supplied study design does not identify the requested correction.

---

## 0.23.x — Bootstrap uncertainty with uncertain infector identity

### Motivation

Version `0.21.0` marginalizes candidate infectors in the likelihood, but uncertainty propagation should also preserve the clustered structure of candidates belonging to one infectee.

### Planned work

- bootstrap at the infectee level, not the candidate-link level;
- preserve all candidate links and their weights when an infectee is resampled;
- refit the uncertain-infector mixture likelihood in each bootstrap replicate;
- propagate fitted generation-interval uncertainty into renewal weights and `Rt`;
- report effective bootstrap sample size and convergence diagnostics;
- expose candidate-weight sensitivity diagnostics.

### Acceptance criteria

- one-candidate-per-infectee data reduce to the current mixed-censoring bootstrap pipeline;
- duplicate identical candidates remain invariant under the bootstrap likelihood;
- fixed-seed reproducibility;
- increasing the number of independent infectees reduces fitted-parameter uncertainty in controlled synthetic experiments.

---

## 0.24.x — Estimate reporting delays from reporting triangles

### Motivation

The current nowcasting layer treats the reporting-delay PMF as known. In operational surveillance the delay distribution should generally be estimated from historical revisions.

### Planned data structure

A reporting triangle containing cumulative or incremental reports by event date and report date.

### Initial model

Start with a transparent discrete delay distribution

\[
D\sim p_0,p_1,\ldots,p_K,
\qquad
\sum_{d=0}^Kp_d=1.
\]

The likelihood should be derived from historical reporting increments with explicit treatment of incomplete recent cohorts.

### Planned work

- validated reporting-triangle data structure;
- conversion between cumulative and incremental triangle representations;
- estimation of a stationary discrete delay PMF;
- right-truncation correction for recent reporting cohorts;
- weekday/report-day effects as a later extension;
- calendar-time-varying delay distributions only after the stationary model is validated.

### Acceptance criteria

- recovery of known delay probabilities in simulated triangles;
- correct normalization and support handling;
- no leakage from future reports into earlier analysis dates;
- plug-compatible output with `reporting_delay_nowcast_m`;
- bootstrap or posterior uncertainty suitable for the `Rt` Monte Carlo layer.

---

## 0.25.x — Joint reporting-delay distribution uncertainty

Once reporting-delay estimation exists, remove the full-pipeline assumption that the delay PMF is fixed.

Each Monte Carlo iteration should be able to draw

\[
p^{(m)}(D),
\]

then use that draw to compute reporting completeness, latent incidence, and `Rt`.

The target uncertainty chain becomes

\[
\text{pair-study uncertainty}
+
\text{generation interval}
+
\text{reporting-delay distribution}
+
\text{latent incidence}
+
\text{renewal }R_t.
\]

Acceptance should include decomposition diagnostics showing how much uncertainty enters from the delay distribution versus latent-incidence completion and renewal inference.

---

# Medium-term roadmap

## 0.26.x — Observation-model integration

The project has both Poisson and negative-binomial renewal models, while the full uncertainty pipeline remains primarily coupled to the conjugate Poisson formulation.

Planned work:

- joint nowcast-plus-NB2 `Rt` inference;
- propagation of reporting and generation-time uncertainty under overdispersed observations;
- estimation or profiling of NB dispersion rather than requiring it to be fixed in every workflow;
- model-adequacy comparisons on simulated overdispersed data;
- limiting checks showing convergence to the Poisson pipeline as dispersion tends to infinity.

Poisson-lognormal or other observation models should be considered only if they add a clearly justified capability rather than model proliferation.

## 0.27.x — Time-varying `Rt` smoothing and state-space inference

Sliding windows are transparent but impose piecewise-constant transmission over a chosen window.

Planned work:

- random-walk model on `log Rt`;
- explicit process variance;
- separation of process evolution from surveillance observation noise;
- penalized/state-space estimation;
- simulation-based calibration against known `Rt` trajectories;
- comparisons with current sliding-window estimates.

Smoothing must remain an explicit prior/process layer rather than a black-box replacement for the renewal likelihood.

## 0.28.x — Susceptible depletion and renewal models

Add an optional susceptible-fraction correction and document the distinction among basic, effective, and instantaneous reproduction measures.

Key validation target: consistency between renewal estimates and deterministic SIR/SEIR trajectories when susceptible depletion becomes material.

## 0.29.x — Importation and local transmission

Distinguish

\[
I_t=I_t^{\mathrm{local}}+I_t^{\mathrm{imported}}.
\]

Imported cases should enter infectious pressure according to an explicit model while avoiding attribution of imported incidence to local reproduction dynamics.

This is particularly important for low-incidence and regional surveillance.

## 0.30.x — Spatial and metapopulation transmission

For regions `i,j`, introduce mobility/contact coupling such as

\[
\lambda_i(t)=\sum_jM_{ij}(t)\mathcal I_j(t),
\]

with explicit definitions for mobility and infectiousness.

Planned capabilities:

- deterministic metapopulation SIR/SEIR models;
- region-specific transmission rates;
- mobility matrices;
- imported/local force-of-infection decomposition;
- structured next-generation matrices and spectral reproduction numbers;
- stochastic metapopulation simulation where computationally reasonable.

---

# Advanced generation-interval and transmission-pair inference

## Candidate-link diagnostics

The uncertain-infector likelihood currently treats candidate weights as external evidence. Useful next diagnostics include

\[
P(j\mid i,\hat\theta)
\propto
\pi_{ij}L_{ij}(\hat\theta),
\]

reported as normalized candidate responsibilities without pretending they are causal probabilities unless the link-weight model justifies that interpretation.

## Uncertain infection times from exposure histories

Move beyond precomputed generation-interval censoring windows and infer infection times from exposure histories directly:

\[
T_i^{\mathrm{infection}}
\sim
p(T_i\mid\text{exposure history}).
\]

Generation intervals would then be derived from latent infection times. This should follow candidate-infector and study-truncation corrections because the joint model otherwise becomes difficult to identify and validate.

## Serial interval versus generation interval

Add explicit symptom-onset models so the package distinguishes:

- generation interval;
- serial interval;
- incubation period;
- onset-to-transmission timing.

The library must never silently substitute one for another.

## Time-varying generation intervals

Interventions, isolation, behavior, or pathogen evolution can change generation intervals as well as transmission intensity. A possible future model uses

\[
w_t(s)
\]

instead of one fixed `w(s)`.

This requires careful identifiability analysis because simultaneous changes in `Rt` and `w_t` can be difficult to distinguish from incidence alone.

## Hierarchical generation-interval models

Potential grouping variables include epidemic phase, strain/variant, age, vaccination status, setting, or calendar period. Hierarchical models should only be introduced with clear shrinkage semantics and strong synthetic recovery tests.

---

# Deterministic and stochastic modelling roadmap

## Additional intervention mechanics

Potential extensions include:

- vaccination with imperfect efficacy;
- waning vaccine protection;
- treatment reducing infectious duration;
- isolation/quarantine compartments;
- importation seeding;
- seasonal transmission.

Each extension should have a clear epidemiological quantity or validation target rather than simply increasing compartment count.

## Adaptive ODE integration

Current deterministic models use explicit RK4. Potential optional improvements:

- adaptive Runge-Kutta integration;
- user-defined tolerances;
- exact event handling at intervention boundaries;
- convergence comparisons against fixed-step RK4 reference trajectories.

Fixed-step RK4 should remain available because it is simple, reproducible, and valuable for testing.

## Faster stochastic simulation

For large populations or high event rates consider:

- tau-leaping;
- hybrid stochastic/deterministic methods;
- bounded-error diagnostics against exact Gillespie simulation.

Exact SSA remains the scientific reference implementation.

---

# Real-data and reproducibility roadmap

## Stable data input layer

Add validated interfaces for:

- daily incidence CSV files;
- reporting triangles;
- transmission-pair tables;
- candidate-infector tables;
- censoring and study-window metadata;
- contact matrices;
- mobility matrices.

Inputs should fail loudly on ambiguous units, invalid date ordering, negative counts, impossible censoring bounds, or malformed candidate groups.

## Reproducible command-line workflows

Provide small CLI-style executables that can run complete analyses without editing source code, for example:

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

Add small, redistributable synthetic/reference datasets for integration tests and examples. Never place confidential or restricted public-health data in the repository.

## End-to-end synthetic validation

Create a simulation-recovery benchmark that:

1. simulates an epidemic with known transmission parameters;
2. generates transmission-pair candidates;
3. imposes censoring, uncertain links, and study-window ascertainment;
4. generates reporting delays and surveillance counts;
5. re-estimates the generation interval;
6. estimates reporting delays;
7. nowcasts latent incidence;
8. estimates `Rt`;
9. compares estimates and interval coverage with known truth.

This should become one of the strongest scientific validation assets in the repository.

---

# Numerical roadmap

## Distribution functions

Continue strengthening numerical distribution infrastructure with analytic and independent-reference tests.

Candidate additions:

- log-CDF and log-survival functions;
- stable log-difference helpers;
- Beta distribution functions;
- normal quantiles;
- negative-binomial random variates;
- robust multinomial/categorical sampling.

## Optimization

The project currently uses problem-specific root finding and Nelder-Mead where appropriate. A small reusable optimization layer may eventually provide:

- bracketed scalar roots;
- bounded one-dimensional minimization;
- Nelder-Mead diagnostics;
- convergence/status codes;
- iteration and objective histories;
- reusable numerical Hessian tools.

Do not introduce a general optimizer abstraction unless it simplifies several existing inference modules.

## Random-number generation

Planned checks:

- document the RNG algorithm and period;
- define a reproducibility contract across supported compilers;
- moment/distribution tests for Gamma, Poisson, and Dirichlet draws;
- categorical sampling for latent candidate-infector workflows;
- optional stream splitting for parallel Monte Carlo.

---

# Software-engineering roadmap

## API stabilization toward 1.0

Before declaring `1.0.0`:

- review naming conventions across all modules;
- identify duplicated summary types;
- standardize validity/convergence reporting;
- standardize shape/rate terminology;
- separate model parameters, fit objects, diagnostics, and Monte Carlo summaries consistently;
- document which APIs are intended as public/stable;
- establish deprecation rules rather than silently renaming public procedures.

## Error handling

The library currently uses `error stop` for many invalid caller inputs. Potential improvement:

- optional status-return APIs for embedded applications;
- structured error/status codes;
- strict fail-fast wrappers for examples and tests.

## Compiler matrix

Expand CI beyond GNU Fortran when practical:

- GNU Fortran;
- Intel ifx when accessible;
- LLVM/flang when sufficiently mature for required features.

Scientific tests must remain compiler-independent.

## Sanitizers and runtime diagnostics

Where supported, add floating-point exception checks, uninitialized-variable diagnostics, bounds checks, and sanitizers for future mixed-language layers.

## Packaging

Investigate Fortran Package Manager (`fpm`) support while retaining CMake if it continues to serve LAPACK and multi-language integration well. The goal is dual usability rather than replacing one build system for fashion.

---

# Interoperability roadmap

## Python interface

A Python layer is useful for data handling and visualization but should remain a wrapper around the Fortran scientific kernel.

Preferred direction:

- stable C ABI or modern F2PY-compatible wrapper layer;
- NumPy arrays for trajectories and inference outputs;
- Python packaging after the Fortran API is sufficiently stable.

Do not restore the previous repository's stale platform-specific Python extension artifact.

## R interface

Potentially useful for epidemiological/statistical workflows, especially generation-interval fitting and `Rt` estimation.

Possible implementation routes include a C ABI with `.Call`/`.Fortran` wrappers or a generated shared library consumed by an R package. This remains lower priority than the scientific inference roadmap unless a concrete workflow requires it.

---

# Documentation roadmap

Each major statistical method should eventually have a standalone note containing:

- model definition;
- parameterization;
- likelihood or posterior;
- numerical method;
- uncertainty interpretation;
- assumptions;
- known limitations;
- references where appropriate;
- regression/reference values.

An architecture document should show dependencies among compartmental models, stochastic simulation, random distributions, generation-interval inference, reporting-delay inference, renewal inference, and joint uncertainty propagation.

Examples should remain small scientific demonstrations rather than undocumented mini-applications.

---

# Release roadmap

## Maturity milestone before 1.0

Before API stabilization, the project should complete at least:

- study-window truncation correction;
- uncertain-infector bootstrap propagation;
- reporting-delay estimation from data;
- joint delay-distribution uncertainty propagation;
- overdispersed full-pipeline inference;
- real-data input contracts;
- end-to-end synthetic recovery tests.

## Toward 1.0.0

`1.0.0` should mean a stable scientific API, not simply a large feature count.

Minimum criteria:

- documented public modules and stable naming conventions;
- deterministic, stochastic, surveillance, and generation-interval components covered by reference/regression tests;
- explicit uncertainty propagation for major fitted upstream quantities;
- real-data interfaces with strong validation;
- stable compiler/build support;
- no known silent numerical-failure modes;
- documented statistical assumptions and limitations;
- changelog and roadmap kept current for every release.

---

# Explicit non-goals

Avoid unless a concrete scientific requirement appears:

- arbitrary compartment models solely to increase model count;
- opaque machine-learning prediction models without epidemiological structure;
- hidden automatic model selection;
- silently substituting serial intervals for generation intervals;
- `Rt` point estimates without explicit uncertainty semantics;
- platform-specific compiled binaries committed to source control;
- a Python-first architecture that makes the Fortran core secondary.

---

# Immediate next sequence

The recommended development order from the current `0.21.0` baseline is:

1. **study-window/truncation correction for transmission-pair inference**;
2. infectee-level bootstrap propagation for uncertain infector identity;
3. reporting-triangle delay estimation;
4. reporting-delay distribution uncertainty propagation;
5. negative-binomial integration with the full nowcasting pipeline;
6. end-to-end synthetic simulation/recovery benchmark;
7. real-data interfaces and file-driven workflows;
8. susceptible depletion/importation extensions;
9. spatial/metapopulation modelling;
10. API consolidation, packaging, and interoperability work.

This order keeps the project centered on one coherent goal: **transparent epidemiological inference in which the data-generating process and every important source of uncertainty remain visible from transmission-pair and surveillance data to reproduction-number estimates.**
