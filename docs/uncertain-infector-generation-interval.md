# Uncertain infector identity in generation-interval inference

Transmission-pair studies do not always identify one infector with certainty. An infectee may have several plausible candidate sources. Hard-assigning the most likely candidate before fitting the generation-interval distribution discards this uncertainty and can distort the fitted timing distribution.

This module keeps candidate identity uncertainty inside the likelihood.

## Candidate representation

For infectee `i`, candidate infectors `j=1,...,m_i` are represented by contiguous ranges in flat candidate arrays. Each candidate carries:

- a lower generation-interval bound `L_ij`;
- an upper bound `U_ij` when finitely interval-censored;
- a right-censoring indicator;
- a non-negative prior/link weight `pi_ij`.

Weights are normalized within each infectee, so only their relative values matter:

\[
\tilde\pi_{ij}=\frac{\pi_{ij}}{\sum_h\pi_{ih}}.
\]

The weights are treated as externally supplied evidence about candidate plausibility. They are not estimated by the current model.

## Marginal likelihood

For a finite candidate interval,

\[
L_{ij}(k,\beta)=F_\Gamma(U_{ij};k,\beta)-F_\Gamma(L_{ij};k,\beta).
\]

For a right-censored candidate,

\[
L_{ij}(k,\beta)=Q_\Gamma(L_{ij};k,\beta).
\]

Candidate identity is marginalized:

\[
L_i(k,\beta)=\sum_j\tilde\pi_{ij}L_{ij}(k,\beta),
\]

and the full log likelihood is

\[
\ell(k,\beta)=\sum_i\log L_i(k,\beta).
\]

The implementation evaluates the candidate mixture with a log-sum-exp calculation to avoid unnecessary underflow when some candidate links are very unlikely under a proposed Gamma distribution.

## Optimization

Shape and rate are optimized on the logarithmic scale,

\[
(\log k,\log\beta),
\]

using a two-dimensional Nelder-Mead search. Positive parameter support is therefore enforced by construction.

The fitted object uses the same `interval_censored_gamma_fit` type as the finite/mixed-censoring modules, so fitted shape, rate, mean, standard deviation, likelihood, iteration count, and convergence state remain consistent across generation-interval workflows.

## Reduction properties

The implementation is regression-tested against three structural identities:

1. **Single-candidate reduction.** One candidate per infectee gives exactly the existing mixed-censoring likelihood.
2. **Duplicate-candidate invariance.** Splitting one candidate into identical duplicates with weights summing to the original weight leaves the likelihood unchanged.
3. **Hard-assignment limit.** As one candidate weight approaches one and the others approach zero, the mixture likelihood approaches the likelihood obtained by selecting that candidate directly.

These are stronger checks than simply asserting optimizer convergence.

## Interpretation

Candidate weights should represent information available independently of the generation-interval density being fitted, such as contact-tracing evidence, genomic linkage scores, household membership, or another externally defined transmission-link model.

If candidate weights themselves depend on the same generation-interval parameters, the simple marginal likelihood above is no longer the appropriate joint model. That requires a coupled transmission-tree model and is outside the current scope.

## Current limitations

The first implementation deliberately does not yet:

- infer candidate-link weights;
- model dependence among infectees or transmission-tree constraints;
- enforce that one source can or cannot infect multiple recipients;
- model genomic sequence evolution;
- correct candidate ascertainment for finite study windows;
- bootstrap infectees and their complete candidate sets;
- propagate uncertain-infector generation-interval uncertainty into `Rt`.

Those are natural follow-up layers. The immediate next uncertainty extension should resample at the infectee level so that all candidate links belonging to an infectee remain coupled inside a bootstrap replicate.
