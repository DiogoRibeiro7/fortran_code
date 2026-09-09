# Fitting a Gamma generation-interval distribution

The `generation_interval_fit_m` module estimates a Gamma generation-interval distribution from strictly positive observed transmission-pair intervals.

The continuous model is

\[
G_i \sim \operatorname{Gamma}(k,\beta),
\]

with shape `k` and rate `beta`. The mean and standard deviation are

\[
E[G]=\frac{k}{\beta},
\qquad
\operatorname{SD}(G)=\frac{\sqrt{k}}{\beta}.
\]

For observed intervals \(G_1,\ldots,G_n\), the profile likelihood equation for the shape is

\[
\log k-\psi(k)
=
\log \bar G-\overline{\log G},
\]

where \(\psi\) is the digamma function. The implementation solves this equation with a bracketed bisection method and then sets

\[
\hat\beta=\frac{\hat k}{\bar G}.
\]

The fit object reports the MLEs, fitted mean and standard deviation, log likelihood, and large-sample standard errors from the two-parameter Fisher information matrix.

## Renewal discretization

The continuous fit is converted to daily renewal weights using

\[
w_j
=
P(j-1 < G \le j),
\qquad j=1,\ldots,L.
\]

For a finite maximum lag \(L\), the implementation reports the omitted tail probability

\[
P(G>L)
\]

before renormalizing the retained bins to sum to one. This makes truncation explicit rather than silently absorbing the tail.

## Scope

The current estimator assumes the observed transmission-pair intervals are independent, positive, and directly observed without interval censoring, truncation, or pair-assignment uncertainty. It is therefore a first parametric layer, not yet a full transmission-pair likelihood for real contact-tracing data.

Future extensions should handle interval-censored exposure windows, uncertain infectors, right truncation in pair ascertainment, and bootstrap or likelihood-based parameter uncertainty propagated into the renewal model.
