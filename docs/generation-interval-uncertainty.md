# Generation-interval uncertainty in Rt inference

The renewal equation depends on the discrete generation-interval probability mass function

\[
w=(w_1,\ldots,w_L),
\qquad
w_s\ge0,
\qquad
\sum_s w_s=1.
\]

Treating `w` as fixed can understate uncertainty in `Rt`, especially when the generation interval is itself estimated from limited epidemiological data.

This implementation uses a Dirichlet uncertainty model around a baseline PMF `w_bar`:

\[
w^{(m)}\sim\operatorname{Dirichlet}(\kappa\,\bar w).
\]

The concentration parameter `kappa` controls how tightly draws cluster around the baseline distribution:

- large `kappa`: generation interval is nearly fixed;
- small `kappa`: substantially more uncertainty in generation timing.

For every draw, the sampled PMF is non-negative and sums exactly to one, and

\[
E[w^{(m)}]=\bar w.
\]

The parameter `kappa` is an uncertainty parameter, not automatically the number of observed transmission pairs. Mapping study evidence to `kappa` requires an explicit statistical model for how the baseline generation interval was estimated.

## Joint Monte Carlo propagation

For each Monte Carlo draw:

1. sample a generation-interval PMF from the Dirichlet model;
2. sample/impute complete incidence under the reporting-delay model;
3. recompute renewal infectiousness using the sampled PMF;
4. estimate the conditional Gamma posterior for `Rt`;
5. sample one `Rt` value from that conditional posterior.

The empirical distribution therefore propagates three sources of uncertainty:

\[
\text{reporting delay}
+
\text{generation interval}
+
\text{renewal estimation}.
\]

The implementation also records the mean generation interval for every Monte Carlo draw, allowing the dispersion induced by `kappa` to be checked directly.

## Scope

This is a sensitivity/uncertainty model for a discrete generation-interval PMF. It does not estimate the PMF from transmission-pair data, correct serial intervals for incubation-period censoring, or model correlation between reporting delays and generation intervals.
