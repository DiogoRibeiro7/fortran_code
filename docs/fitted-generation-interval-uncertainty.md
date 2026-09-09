# Fitted generation-interval uncertainty

This layer propagates uncertainty from observed transmission-pair generation intervals into renewal-equation estimates of the time-varying reproduction number.

Given observed positive pair intervals

\[
G_1,\ldots,G_n,
\]

each Monte Carlo replicate performs a nonparametric bootstrap of the pair sample:

\[
G_1^{(m)},\ldots,G_n^{(m)}
\sim
\widehat F_n.
\]

A Gamma generation-interval model is then refitted to that bootstrap sample,

\[
G^{(m)}\sim \operatorname{Gamma}(k^{(m)},\beta^{(m)}),
\]

and discretized to renewal weights

\[
w_j^{(m)}=P(j-1<G^{(m)}\le j).
\]

For a finite renewal lag window, the omitted Gamma tail probability is recorded before the retained weights are renormalized.

The same replicate also propagates reporting-delay uncertainty. Recent incomplete event-day counts are completed by drawing from the reporting-delay posterior predictive distribution. The resulting completed incidence path and bootstrap-fitted generation weights are passed to the renewal estimator, and one draw is then taken from the conditional Gamma posterior for \(R_t\).

Thus each Monte Carlo draw contains three uncertainty sources:

1. transmission-pair sampling uncertainty in the fitted generation interval;
2. reporting-delay uncertainty in recent incidence;
3. renewal-estimation uncertainty in \(R_t\).

The returned bootstrap diagnostics include the number of converged Gamma fits, mean and standard deviation of bootstrap shape/rate estimates, mean and standard deviation of the fitted mean generation interval, and mean truncated tail probability.

## Interpretation

This bootstrap answers a different question from the Dirichlet sensitivity model. The Dirichlet approach treats uncertainty around a supplied discrete generation-interval PMF through a chosen concentration parameter. The bootstrap approach instead derives generation-interval uncertainty from the finite transmission-pair dataset itself.

## Current scope

The procedure assumes observed transmission-pair intervals are independent, positive, and measured without error. It does not yet account for uncertain infector assignment, interval-censored infection times, pair ascertainment bias, or dependence between transmission pairs. Bootstrap replicates with degenerate intervals can fail to yield a finite Gamma MLE; those fits are skipped and the number of successful replicates is reported explicitly.
