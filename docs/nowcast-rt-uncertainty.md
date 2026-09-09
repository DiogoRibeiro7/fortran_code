# Propagating reporting-delay uncertainty into Rt

The plug-in nowcast pipeline uses the posterior mean of latent event-day incidence as the input to renewal-equation inference. This is useful for bias correction, but it does not propagate reporting-delay uncertainty into the final Rt interval.

The Monte Carlo layer in `nowcast_rt_uncertainty_m` propagates both uncertainty sources.

For event day `t`, cumulative reports observed at age `a` satisfy

\[
Y_{t,a}\mid\lambda_t\sim\operatorname{Poisson}(\lambda_t F(a)),
\]

with posterior

\[
\lambda_t\mid Y_{t,a}\sim\operatorname{Gamma}(a_t,b_t).
\]

For Monte Carlo replicate `m`, the implementation draws

\[
\lambda_t^{(m)}\sim\operatorname{Gamma}(a_t,b_t)
\]

and then imputes only the reports that remain outstanding:

\[
U_t^{(m)}\sim\operatorname{Poisson}\!\left(\lambda_t^{(m)}[1-F(a)]\right).
\]

The simulated completed incidence is

\[
I_t^{(m)}=Y_{t,a}+U_t^{(m)}.
\]

This construction preserves reports that are already observed. In particular, if reporting is complete, `F(a)=1` and therefore `U_t=0` exactly.

Each completed incidence path is passed through the renewal-equation estimator. Conditional on that path, the Poisson renewal model yields a Gamma posterior for Rt. One draw is then taken from that conditional posterior. Repeating this produces draws from the mixture induced by reporting-delay uncertainty and renewal uncertainty.

The reported Monte Carlo summary contains the empirical mean and central quantile interval of these Rt draws. Fixed seeds make the full propagation reproducible.

The implementation currently assumes a known reporting-delay distribution and known generation-interval distribution. Uncertainty in those distributions is not yet propagated.
