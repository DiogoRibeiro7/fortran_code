# Interval-censored generation-interval fitting

Transmission-pair infection times are often not observed exactly. Instead, contact tracing can identify an interval in which a secondary infection occurred. For pair `i`, let the generation interval satisfy

\[
L_i < G_i \le U_i.
\]

The interval-censored Gamma model assumes

\[
G_i\sim\operatorname{Gamma}(k,\beta),
\]

with shape `k` and rate `beta`. The likelihood contribution is therefore

\[
P(L_i<G_i\le U_i)
=
F_\Gamma(U_i;k,\beta)-F_\Gamma(L_i;k,\beta).
\]

The implementation maximizes the summed log likelihood on the unconstrained log-parameter scale

\[
(\log k,\log\beta),
\]

using a two-dimensional Nelder-Mead search. Positivity of the fitted parameters is therefore guaranteed by construction.

Observed-information standard errors are computed numerically from the Hessian of the negative log likelihood on the log-parameter scale and transformed back to shape, rate, and mean-generation-interval uncertainty.

For renewal inference, the fitted continuous Gamma distribution is discretized as

\[
w_j=P(j-1<G\le j),
\]

for integer lag `j`. The probability mass beyond the configured maximum lag is returned explicitly before the retained bins are renormalized.

## Current scope

This implementation supports finite interval censoring with `0 <= L_i < U_i < infinity`. It does not yet support right-censored open-ended intervals, uncertain infector identity, pair-selection bias, or joint inference over exposure windows for both infector and infectee.

The fixed regression fixture uses windows of width 0.8 days around the exact-pair reference sample. An independent optimization gives approximately

```text
shape = 2.996168
rate  = 0.503204
mean  = 5.954187 days
```

and the Fortran implementation is regression-tested against those values.
