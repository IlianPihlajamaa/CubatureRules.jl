# Error estimates

## Embedded rules

An [`EmbeddedRule`](@ref)`(fine, coarse)` pairs a rule with a coarser rule whose nodes are a
subset of its own. When it is constructed, the coarse weights are placed on the matching
nodes of the fine rule (with zeros elsewhere), so that one pass over the nodes evaluates the
integrand once per node and accumulates both sums. Their difference estimates the error of
the fine result, at no extra cost in integrand evaluations.

Nodes are matched with a relative tolerance. Double-exponential rules on unbounded domains
have nodes spread over tens of orders of magnitude, and an absolute tolerance would match
every node near the origin with every other.

`embedded(family, domain; degree)` builds such a pair for families that nest.
Gauss–Kronrod contains its Gauss rule by construction; consecutive levels of the
double-exponential rules share nodes, and so do second-kind Fejér rules at degrees
`2^k - 1`.

`integrate(f, e; error = true)` returns an [`IntegrationResult`](@ref) with the value, the
error estimate, the number of evaluations, and the family and degree of the fine rule.

## Why there is no adaptive driver

Earlier versions had `integrate(f, domain; rtol)`, which built rules of increasing degree
from one family until two successive results agreed. It was removed in v0.5, for three
reasons:

- It assumed that one family forms a clean sequence of rules. Seeded families have gaps and
  an upper degree, and families round requested degrees up. Two requests could give the same
  rule, whose difference is exactly zero; this made the driver report convergence on a
  result that was wrong in the fifth digit.
- Each step builds a new rule, and for seeded families refines it. At high precision that
  can take much longer than the integration itself, and nothing in the call shows it.
- It was the one place where the package constructed rules out of sight, which the rest of
  the design avoids.

Adaptive integration is a separate problem, and good packages exist for it: QuadGK.jl on
intervals, HCubature.jl on boxes and HAdaptiveIntegration.jl on simplices. They subdivide
the domain, which also handles integrands that a single rule of any degree cannot.
