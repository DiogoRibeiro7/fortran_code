# Right-censored generation intervals

Transmission-pair studies may provide a finite interval for some pairs and only a lower bound for others. The generation-interval fitter therefore supports two record types under the same Gamma(shape, rate) model.

For a finite interval,

\[
L_i < G_i \le U_i,
\]

the likelihood contribution is

\[
F_\Gamma(U_i;k,\beta)-F_\Gamma(L_i;k,\beta).
\]

For a right-censored record,

\[
G_i > L_i,
\]

the contribution is the Gamma survival probability

\[
Q_\Gamma(L_i;k,\beta)=1-F_\Gamma(L_i;k,\beta).
\]

The implementation evaluates the survival function directly with the upper incomplete-Gamma continued fraction in the upper tail instead of forming `1 - CDF`, which avoids cancellation when the survival probability is small.

The optimizer remains a two-dimensional Nelder-Mead search on `log(shape)` and `log(rate)`. Observed-information standard errors are calculated from the numerical Hessian on that log-parameter scale and transformed back to shape, rate, and the mean generation interval.

The existing finite-interval API is preserved. `fit_interval_censored_gamma` is now a special case of `fit_mixed_censored_gamma` with all censoring indicators set to false.

The fitted continuous Gamma distribution is converted to renewal weights using

\[
w_j=P(j-1<G\le j),
\]

and the omitted probability beyond the configured maximum lag is evaluated directly from the Gamma survival function before the retained bins are renormalized.

## Scope

This slice supports finite interval censoring and right censoring. Left censoring, exact zero-width observations, uncertain infector identity, dependent transmission pairs, and ascertainment correction remain outside the current model.
