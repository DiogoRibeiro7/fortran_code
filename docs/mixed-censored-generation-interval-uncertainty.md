# Mixed-censored generation-interval uncertainty

This module extends the interval-censored generation-interval bootstrap to datasets containing both finite censoring windows and right-censored transmission-pair records.

## Observation records

Each transmission-pair record is represented by a lower bound `L_i`, an upper bound `U_i`, and a logical right-censoring flag `delta_i`.

For finite records (`delta_i = false`),

\[
L_i < G_i \le U_i,
\]

and the Gamma likelihood contribution is

\[
F_\Gamma(U_i;k,\beta)-F_\Gamma(L_i;k,\beta).
\]

For right-censored records (`delta_i = true`),

\[
G_i > L_i,
\]

and the likelihood contribution is

\[
Q_\Gamma(L_i;k,\beta).
\]

The upper bound of a right-censored record is ignored by the likelihood.

## Bootstrap semantics

The bootstrap resampling unit is the complete observed record

\[
(L_i,U_i,\delta_i).
\]

Lower bounds, upper bounds, and censoring indicators are never resampled independently. This preserves the empirical joint distribution of observed censoring patterns.

For bootstrap replicate `m`:

1. sample `n` complete pair records with replacement;
2. fit the mixed-censored Gamma generation-interval model;
3. discretize the fitted Gamma distribution into renewal weights;
4. sample reporting-delay completion uncertainty for recent incidence;
5. estimate the conditional renewal-equation `Rt` posterior;
6. draw `Rt` from that conditional posterior.

The resulting empirical `Rt` distribution therefore propagates uncertainty from the finite transmission-pair study, censoring structure, incomplete surveillance reporting, and renewal inference.

## Degenerate bootstrap samples

The current mixed-censoring fitter requires at least two finite interval records for a stable two-parameter Gamma fit. A bootstrap resample containing fewer than two finite records is skipped. Failed numerical fits are also skipped. The returned bootstrap summary reports both requested and converged draw counts so this loss of Monte Carlo support is visible.

A low convergence fraction should be treated as a warning about study design, sample size, or censoring severity rather than hidden by synthetic parameter substitutions.

## Renewal discretization

Each successful Gamma fit is converted into daily weights

\[
w_j=P(j-1<G\le j),\qquad j=1,\ldots,L.
\]

The probability mass beyond the configured maximum lag is recorded before the retained weights are renormalized. A large mean tail probability indicates that the renewal grid is too short for the fitted generation-time distribution.

## Interpretation

This procedure is a nonparametric bootstrap over observed transmission-pair records combined with parametric Gamma refitting. It quantifies finite-sample uncertainty conditional on the observed censoring mechanism and on the Gamma family.

It does not correct for transmission-pair ascertainment bias, epidemic-phase bias, preferential observation of short intervals, uncertain infector identity, dependence between pairs, or misspecification of the Gamma family. Those are separate modelling problems in the roadmap.

## Relationship to earlier uncertainty models

The repository now contains three progressively more data-grounded generation-interval uncertainty paths:

- `generation_interval_uncertainty_m`: Dirichlet sensitivity around a supplied discrete PMF;
- `fitted_generation_interval_uncertainty_m`: bootstrap of exact positive transmission-pair intervals;
- `interval_censored_generation_interval_uncertainty_m`: bootstrap of finite censoring windows;
- `mixed_censored_generation_interval_uncertainty_m`: bootstrap of finite and right-censored records.

The mixed-censoring path is the preferred API when a transmission-pair study contains open-ended right-censored observations.
