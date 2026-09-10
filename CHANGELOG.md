# Changelog

All notable changes to `fortran-epidemiology` are documented in this file.

The project evolved from the earlier `fortran_code` repository into a focused epidemiological modelling and inference library. The entries below follow the semantic versions used by the CMake project and summarize both scientific and engineering changes. Where relevant, the associated pull request is noted.

## [Unreleased]

### Planned

- Observation-window and truncation correction for transmission-pair studies.
- Estimation of reporting-delay distributions from reporting triangles rather than treating delay probabilities as fixed inputs.
- Propagation of fitted reporting-delay uncertainty through nowcasting and `Rt` inference.
- Broader real-data interfaces and reproducible surveillance workflows.

## [0.21.0] - 2026-09-10

### Added

- Generation-interval inference for infectees with multiple plausible candidate infectors.
- Candidate-level finite-interval and right-censored likelihood contributions.
- Marginalization over candidate infector identity using externally supplied non-negative link weights normalized within infectee.
- Numerically stable log-sum-exp evaluation of each infectee-level mixture likelihood.
- Flat candidate-array representation with contiguous start/end ranges per infectee.
- `uncertain_infector_generation_interval_example`.
- `docs/uncertain-infector-generation-interval.md`.

### Scientific notes

- The infectee-level likelihood is
  \[
  L_i(k,\beta)=\sum_j \tilde\pi_{ij}L_{ij}(k,\beta),
  \]
  where each candidate contribution may be either a finite censoring interval or a right-censored survival probability.
- Candidate weights are treated as external evidence rather than estimated jointly with the generation-interval model.
- Structural regression tests verify single-candidate reduction, duplicate-candidate invariance, and the hard-assignment limit under extreme candidate weights.

### Changed

- Project version advanced to `0.21.0`.

### Pull request

- PR #22: `Marginalize uncertain infector identity in generation-interval inference`.

## [0.20.0] - 2026-09-10

### Added

- Bootstrap propagation of mixed finite-interval and right-censored generation-interval uncertainty into renewal-equation `Rt` inference.
- Resampling of complete transmission-pair records `(L_i,U_i,delta_i)`, preserving censoring type within each bootstrap observation.
- Mixed-censoring Gamma refitting inside every successful bootstrap replicate.
- Propagation of fitted generation-interval uncertainty together with reporting-delay completion uncertainty and renewal posterior uncertainty.
- Explicit requested-versus-converged bootstrap counts when resamples become non-identifiable or numerical fits fail.
- `mixed_censored_generation_interval_uncertainty_example`.
- `docs/mixed-censored-generation-interval-uncertainty.md`.

### Scientific notes

- Bootstrap samples with fewer than two finite records are skipped rather than repaired with invented information.
- The finite-only uncertainty path remains available as a separate API, keeping censoring semantics explicit.

### Changed

- Project version advanced to `0.20.0`.

### Pull request

- PR #21: `Propagate mixed-censoring generation-interval uncertainty into Rt`.

## [0.19.0] - 2026-09-10

### Added

- Mixed finite-interval and right-censored Gamma generation-interval inference.
- Right-censored records with contribution
  \[
  P(G>L)=Q_\Gamma(L;k,\beta).
  \]
- `regularized_gamma_q` for numerically stable Gamma survival probabilities.
- A mixed-censoring regression fixture and executable example.
- `docs/right-censored-generation-interval.md`.

### Changed

- Finite interval probabilities can be evaluated using either CDF differences or survival-function differences, choosing the numerically safer positive representation.
- Renewal-grid tail probabilities are computed directly through the Gamma survival function instead of `1 - CDF`.
- Project version advanced to `0.19.0`.

### Scientific notes

- The mixed-censoring reference fit with 16 finite intervals and 4 right-censored observations yields approximately
  \[
  \hat k=2.26534655,\qquad \hat\beta=0.34202380,
  \]
  with fitted mean generation interval about `6.62336` days.
- Right censoring materially changes the fitted tail and therefore the inferred generation-time distribution.

### Pull request

- PR #19: `Support right-censored generation intervals`.

## [0.18.0] - 2026-09-10

### Added

- Bootstrap propagation of interval-censored generation-interval uncertainty into renewal-equation `Rt` estimates.
- Resampling of complete censoring records `(L_i,U_i)` rather than lower and upper endpoints independently.
- Bootstrap diagnostics for Gamma shape, rate, mean generation interval, convergence counts, and renewal-grid tail mass.
- `interval_censored_generation_interval_uncertainty_example`.
- `docs/interval-censored-generation-interval-uncertainty.md`.

### Changed

- CMake example and test target declarations were consolidated into reusable helper functions.
- Project version advanced to `0.18.0`.

### Scientific notes

- The Monte Carlo pipeline now combines uncertainty from interval-censored transmission-pair sampling, reporting-delay completion, and renewal-model `Rt` inference.

### Pull request

- PR #18: `Bootstrap interval-censored generation intervals into Rt uncertainty`.

## [0.17.0] - 2026-09-10

### Added

- Interval-censored Gamma generation-interval fitting.
- Likelihood contribution
  \[
  P(L_i<G_i\le U_i)=F_\Gamma(U_i;k,\beta)-F_\Gamma(L_i;k,\beta).
  \]
- Two-parameter optimization on the log-shape/log-rate scale using Nelder-Mead.
- Observed-information standard errors from a numerical Hessian.
- Daily renewal-weight discretization from the fitted continuous Gamma distribution.
- `interval_censored_generation_interval_example`.
- `docs/interval-censored-generation-interval.md`.

### Scientific notes

- The reference interval-censored fit gives approximately
  \[
  \hat k=2.99616785,\qquad \hat\beta=0.50320350,
  \]
  with fitted mean near `5.95419` days.

### Pull request

- PR #17: `Fit interval-censored generation intervals`.

## [0.16.0] - 2026-09-09

### Added

- Nonparametric bootstrap of exact transmission-pair generation intervals.
- Refit of the Gamma generation-interval model inside each bootstrap replicate.
- Propagation of fitted generation-interval uncertainty into the nowcast-plus-renewal `Rt` Monte Carlo pipeline.
- Bootstrap diagnostics for fitted shape, rate, mean generation interval, and omitted renewal-grid tail probability.
- `fitted_generation_interval_uncertainty_example`.
- `docs/fitted-generation-interval-uncertainty.md`.

### Scientific notes

- This version replaced an externally selected Dirichlet concentration with uncertainty induced by the observed finite transmission-pair sample.
- Degenerate bootstrap samples that do not admit a finite Gamma MLE are skipped explicitly.

### Pull request

- PR #16: `Bootstrap fitted generation intervals into Rt uncertainty`.

## [0.15.0] - 2026-09-09

### Added

- Direct maximum-likelihood fitting of a Gamma generation-interval distribution from exact positive transmission-pair intervals.
- Shape/rate parameterization
  \[
  G_i\sim\operatorname{Gamma}(k,\beta).
  \]
- Safeguarded root solve for the Gamma shape equation
  \[
  \log k-\psi(k)=\log \bar G-\overline{\log G}.
  \]
- Analytic rate estimate `beta = k / mean(G)`.
- Large-sample standard errors from the Fisher information matrix.
- Conversion of the continuous fitted distribution into discrete renewal weights
  \[
  w_j=P(j-1<G\le j).
  \]
- Explicit reporting of omitted probability mass beyond the configured maximum lag.
- `generation_interval_fit_example`.
- `docs/generation-interval-fit.md`.

### Scientific notes

- The 20-pair reference sample yields approximately shape `2.9875687421`, rate `0.5016908047`, mean `5.9550` days, and SD `3.4453` days.

### Pull request

- PR #15: `Fit Gamma generation intervals from transmission pairs`.

## [0.14.0] - 2026-09-09

### Added

- Generation-interval uncertainty in the Monte Carlo `Rt` pipeline using
  \[
  w^{(m)}\sim\operatorname{Dirichlet}(\kappa\bar w).
  \]
- Reusable Dirichlet random sampler built from Gamma variates.
- Joint propagation of reporting-delay uncertainty, generation-interval uncertainty, and renewal-estimation uncertainty.
- Diagnostics for sampled mean generation intervals.
- `generation_interval_uncertainty_example`.
- `docs/generation-interval-uncertainty.md`.

### Scientific notes

- The concentration parameter `kappa` is interpreted as a sensitivity parameter controlling uncertainty around a baseline PMF; it is not automatically equated with a sample size.

### Pull request

- PR #14: `Propagate generation-interval uncertainty into Rt`.

## [0.13.0] - 2026-09-09

### Added

- Full Monte Carlo propagation of reporting-delay nowcast uncertainty into renewal-equation `Rt` inference.
- Gamma and Poisson random variate generation using the repository RNG.
- Posterior predictive completion of incomplete event-day incidence:
  \[
  \lambda_t^{(m)}\sim p(\lambda_t\mid Y_{t,a}),
  \]
  \[
  U_t^{(m)}\sim\operatorname{Poisson}(\lambda_t^{(m)}[1-F(a)]),
  \]
  \[
  I_t^{(m)}=Y_{t,a}+U_t^{(m)}.
  \]
- Empirical Monte Carlo `Rt` summaries from conditional Gamma `Rt` draws.
- `nowcast_rt_uncertainty_example`.
- `docs/nowcast-rt-uncertainty.md`.

### Scientific notes

- Complete reporting collapses exactly to observed incidence because future reports are identically zero when `F(a)=1`.

### Pull request

- PR #13: `Propagate reporting-delay uncertainty into Rt`.

## [0.12.0] - 2026-09-09

### Added

- End-to-end plug-in surveillance pipeline from reported incidence through reporting-delay nowcasting to renewal `Rt` inference.
- Parallel output of raw-report and nowcast-corrected `Rt` series.
- Regression test showing that right truncation can manufacture a false decline in recent `Rt`.
- `nowcasted_rt_pipeline_example`.

### Scientific notes

- This version used posterior-mean nowcast incidence as a plug-in input and therefore corrected right truncation without yet propagating nowcast uncertainty through `Rt`.

### Pull request

- PR #12: `Connect reporting-delay nowcasts to Rt inference`.

## [0.11.0] - 2026-09-08

### Added

- Reporting-delay nowcasting for right-truncated surveillance counts.
- Discrete reporting-delay PMFs and event-day completeness `F(a)`.
- Poisson thinning model
  \[
  Y_{t,a}\mid\lambda_t\sim\operatorname{Poisson}(\lambda_tF(a)).
  \]
- Gamma posterior for latent incidence.
- Posterior mean, SD, credible intervals, and expected future reports.
- Shared `gamma_distribution_m` numerical utilities.
- `reporting_delay_nowcast_example`.

### Changed

- Gamma CDF/quantile code was extracted from renewal inference into shared infrastructure.

### Pull request

- PR #11: `Add reporting-delay nowcasting`.

## [0.10.0] - 2026-09-08

### Added

- Negative-binomial NB2 renewal observation model for overdispersed surveillance counts.
- Variance model
  \[
  \operatorname{Var}(I_t\mid\mu_t)=\mu_t+\mu_t^2/k.
  \]
- Maximum-likelihood `Rt` estimation for fixed dispersion.
- Profile-likelihood confidence intervals.
- Poisson limiting-regime tests.
- `negative_binomial_rt_example`.

### Scientific notes

- Gamma-Poisson conjugacy is intentionally not reused for the negative-binomial likelihood.

### Pull request

- PR #10: `Add negative-binomial Rt observation model`.

## [0.9.0] - 2026-09-08

### Added

- Renewal-equation inference for time-varying reproduction numbers.
- Total infectiousness
  \[
  \Lambda_t=\sum_{s\ge1}I_{t-s}w_s.
  \]
- Poisson renewal observation model
  \[
  I_t\mid R_t,\Lambda_t\sim\operatorname{Poisson}(R_t\Lambda_t).
  \]
- Sliding-window Gamma-Poisson posterior inference.
- Exact Gamma credible intervals by numerical CDF inversion.
- `renewal_rt_example`.

### Pull request

- PR #9: `Add renewal-equation Rt inference`.

## [0.8.0] - 2026-09-08

### Added

- Generation-interval distributions and Euler-Lotka inference.
- `R <-> r` conversion through
  \[
  1=R\int_0^\infty e^{-r\tau}w(\tau)\,d\tau.
  \]
- Doubling and halving times.
- Log-linear early-growth estimation.
- Mechanistic SEIR generation-interval distribution.
- Analytic SEIR growth relation
  \[
  (r+\sigma)(r+\gamma)=\sigma\beta.
  \]
- Tests showing that generation-interval shape matters, not only its mean.
- `growth_rate_example`.

### Pull request

- PR #8: `Add generation-interval and epidemic growth inference`.

## [0.7.0] - 2026-09-08

### Added

- Deterministic age-structured SIR transmission.
- Explicit contact matrices with heterogeneous susceptibility, infectiousness, and recovery rates.
- Force of infection
  \[
  \lambda_i=\beta q_i\sum_j C_{ij}z_jI_j/N_j.
  \]
- Next-generation matrix and `R0` as its spectral radius.
- LAPACK `DGEEV` eigenvalue calculation.
- Population-weighted contact reciprocity diagnostics.
- Three-group age-structured example.

### Pull request

- PR #7: `Add age-structured transmission and contact-matrix R0`.

## [0.6.0] - 2026-09-08

### Added

- Stochastic SEIR continuous-time Markov model using Gillespie's direct method.
- Infection, progression, and recovery event hazards.
- Reproducible stochastic ensembles.
- Extinction, fadeout, major-outbreak, final-size, peak-prevalence, and horizon-state summaries.
- Tests comparing large-population stochastic ensemble means with deterministic SEIR trajectories.
- `stochastic_seir_ensemble_example`.

### Pull request

- PR #6: `Add stochastic SEIR ensembles and fadeout analysis`.

## [0.5.0] - 2026-09-08

### Added

- Stochastic SIR continuous-time Markov model using Gillespie's direct method.
- Explicit infection and recovery hazards.
- Portable seeded pseudo-random number generator.
- Exact integer population conservation and natural extinction semantics.
- Fixed-seed stochastic regression tests.
- `stochastic_sir_example`.

### Pull request

- PR #5: `Add stochastic SIR simulation with Gillespie SSA`.

## [0.4.0] - 2026-09-08

### Added

- Endemic SEIRS model with births, deaths, waning immunity, and vaccination.
- Demographic reproduction number
  \[
  R_0=\frac{\beta\sigma}{(\sigma+\mu)(\gamma+\mu)}.
  \]
- Vaccinated disease-free susceptible fraction and control reproduction number.
- Long-horizon vaccination example and endemic-dynamics regression tests.

### Pull request

- PR #4: `Add SEIRS demography, waning immunity, and vaccination`.

## [0.3.0] - 2026-09-08

### Added

- Deterministic SEIR model.
- Reusable piecewise transmission schedules.
- RK4 integration with transmission rate evaluated at each RK4 stage.
- Intervention example demonstrating reduced peak and cumulative incidence.

### Pull request

- PR #3: `Add SEIR dynamics and intervention schedules`.

## [0.2.0] - 2026-09-08

### Added

- First epidemiological model: deterministic closed-population SIR.
- Typed SIR parameters and state.
- Incidence, derivative, RK4 stepping, complete trajectory simulation, and `R0=beta/gamma`.
- Scientific invariants for population conservation, non-negativity, infectious growth, and susceptible depletion.
- `sir_example`.

### Changed

- Repository scientific direction shifted from generic numerical demonstrations to epidemiological modelling.
- CMake project renamed to `fortran_epidemiology`.
- README rewritten around the new epidemiology scope.

### Pull request

- PR #2: `Build epidemiology modelling foundation`.

## [0.1.0] - 2026-09-08

### Added

- Modern CMake/CTest Fortran project structure.
- Reusable modules under `src/`, executables under `app/`, and tests under `test/`.
- Strict GNU Fortran warning and runtime-check flags.
- GitHub Actions CI with BLAS/LAPACK.
- Typed LAPACK `DGESV` wrapper with shape validation and explicit status handling.

### Changed

- Replaced legacy `real*8` declarations with `real(real64)` from `iso_fortran_env`.
- Removed committed executables, stale macOS/Python binary artifacts, and ad hoc build scripts.
- Replaced the original one-line README with reproducible build and numerical documentation.

### Historical note

- This release marks the modernization of the repository that was then known as `fortran_code`.

### Pull request

- PR #1: `Modernize Fortran numerical examples`.
