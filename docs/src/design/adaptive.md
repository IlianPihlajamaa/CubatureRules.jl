# Sequences and adaptivity

## Sequences of rules

A [`RuleSequence`](@ref)`(family, domain; degrees)` is a lazy sequence of rules of increasing
degree. Iterating over it constructs the rules one by one. The default schedule is
`1, 3, 7, 15, 31, 63`, limited to the degrees the family offers.

A family often answers a requested degree with a higher one that it has. Requested degrees
that would give the same rule are removed when the sequence is created, so a sequence never
contains the same rule twice. (Before this was done, two identical rules in a row gave a
difference of exactly zero, and tolerance-based integration stopped too early.)

A [`LevelSequence`](@ref)`(family, domain; levels)` is the same for families that are
parameterised by a level instead of a degree, such as `TanhSinh`, `ExpSinh` and `SinhSinh`.

## Integration to a tolerance

`integrate(f, domain; rtol, atol)` walks a sequence and stops when two successive results
differ by less than `max(atol, rtol |value|)`. It returns an [`IntegrationResult`](@ref) with
the last value, the last difference as error estimate, the total number of function
evaluations, the family and degree of the last rule, and whether the tolerance was met.

The family is chosen by `adaptive_family`: among the candidates at a moderate degree it
prefers nested families, whose rules share nodes with the rules of lower degree (for example
Gauss–Kronrod), and otherwise takes the one with the fewest points. On domains where no
family has a polynomial degree, it uses a double-exponential family and walks its levels.

The difference between successive results estimates the error of the less accurate of the
two, so returning the more accurate one is conservative.

## Embedded rules

An [`EmbeddedRule`](@ref)`(fine, coarse)` pairs a rule with a coarser rule whose nodes are a
subset of its own. Both estimates are computed from the same function values, and their
difference estimates the error of the fine rule without extra evaluations of `f`. This is
how Gauss–Kronrod rules are normally used.

`embedded(family, domain; degree)` builds such a pair for families that are nested;
Gauss–Kronrod contains its Gauss rule by construction. `integrate(f, e; error = true)` returns
an `IntegrationResult`.

## What this is not

All of the above raises the order of a single rule on the whole domain. It never subdivides
the domain in response to the integrand. That works well when the integrand is smooth, where
the error decreases quickly with the degree, and badly when the integrand has a local feature
such as a peak or an interior singularity. For those, space-adaptive methods such as
QuadGK.jl (in one dimension) and HCubature.jl (in several) are the right tools.
