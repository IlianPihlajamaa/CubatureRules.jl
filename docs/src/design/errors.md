# Error estimates

## Embedded rules

An [`EmbeddedRule`](@ref)`(fine, coarse)` pairs a rule with a coarser rule whose nodes are a
subset of its own. When it is constructed, the coarse weights are placed on the matching
nodes of the fine rule (with zeros elsewhere), so that one pass over the nodes evaluates the
integrand once per node and accumulates both sums, at no extra cost in integrand
evaluations. Their difference is the error of the coarse result, up to the much smaller
error of the fine one: an estimate, usually pessimistic, of the error of the fine result.
In the Fejér example of [I want an error estimate](../tutorial/errors.md) the estimate is
several orders of magnitude above the finer rule's actual error.

Nodes are matched with a relative tolerance. Double-exponential rules on unbounded domains
have nodes spread over tens of orders of magnitude, and an absolute tolerance would match
every node near the origin with every other.

`embedded(family, domain; degree)` builds such a pair for the two families whose rules nest
by construction: Gauss–Kronrod, which contains its Gauss rule, and Gauss–Patterson, each of
whose rules contains the level below it. Other rules nest too — consecutive levels of the
double-exponential rules share nodes, and so do second-kind Fejér rules at degrees `2^k - 1`
— and for them the pair is built by hand, as `EmbeddedRule(fine, coarse)`.

`integrate(f, e; error = true)` returns an [`IntegrationResult`](@ref) with the value, the
error estimate, the number of evaluations, and the family and degree of the fine rule.

## Why there is no adaptive driver

Before the first release the package had `integrate(f, domain; rtol)`, which built rules of
increasing degree from one family until two successive results agreed. It was removed for
three reasons:

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
