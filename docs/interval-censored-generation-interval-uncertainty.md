# Interval-censored generation-interval uncertainty

This layer propagates uncertainty from interval-censored transmission-pair data into renewal-equation estimates of the time-varying reproduction number.

For pair `i`, the generation interval is only known to satisfy

\[
L_i < G_i \le U_i.
\]

A bootstrap replicate resamples the observed pair records `(L_i, U_i)` with replacement. Lower and upper bounds are always resampled together; they are never treated as independent observations.

For each bootstrap sample, the Gamma interval-censored likelihood is refitted,

\[
\ell(k,\beta)
=
\sum_i \log\left[F_\Gamma(U_i;k,\beta)-F_\Gamma(L_i;k,\beta)\right],
\]

and the fitted continuous distribution is discretized into daily renewal weights

\[
w_j=P(j-1<G\le j).
\]

The probability beyond the configured maximum lag is recorded before the retained bins are renormalized.

Each successful bootstrap fit is then combined with reporting-delay completion uncertainty and the conditional renewal `Rt` posterior. The resulting empirical `Rt` distribution therefore includes sampling uncertainty in the interval-censored transmission-pair study, recent-reporting uncertainty, and renewal-estimation uncertainty.

The implementation reports bootstrap means and standard deviations for Gamma shape, rate, and mean generation interval, together with the number of converged fits and the average omitted tail probability.

This remains a nonparametric case bootstrap over the observed pair records. It does not model dependence between transmission pairs, ascertainment bias, uncertain infector identity, or open-ended censoring. Those require a richer observation model rather than a different resampling scheme.
