# Selection

This page describes how `rule(domain; degree)` finds and ranks families.

## Finding families

The registry has no list of families. Each time `rule`, `available` or `compare` is called,
it collects all concrete subtypes of `RuleFamily` that are currently loaded, using
`subtypes`. A family defined in another package is therefore found as soon as that package
is loaded, with no registration step.

This lookup deliberately happens at call time. Doing it once at load time (for example in a
`const`) would run during precompilation and miss every family defined in a package loaded
later. The test suite includes a separately precompiled package that defines a family, to
check that it is found.

## Candidates

For a request on domain `dom` with claim `PolynomialDegree(d)`, each family type `F` is asked
`candidates(F, reference(dom), claim)`. It returns the instances that can serve the request,
usually zero or one. Parameterised families return the instance that matches the domain;
`GaussJacobi` reads `α` and `β` from a `JacobiWeight`, for example.

Combinator families (`ConicalProduct`, `TensorProduct`, `BallProduct`, `GaussianProduct`)
return one candidate for each suitable inner family, so every combination is ranked.

For each candidate the registry records, without building anything: the number of points,
the achieved degree, whether it is derived or seeded, whether its weights are positive and
its nodes interior, its symmetry, and an estimated cost.

## Ranking

Candidates are sorted by

1. number of points (fewest first),
2. derived before seeded,
3. family name.

The third key only exists to make the order deterministic. Filters are applied after
sorting: `positive = true` and `interior = true` remove candidates without those properties,
and candidates that cannot produce the requested output type (`supports_type`) are removed.
Then, if the first remaining candidate has negative weights and a positive one needs at most
25% more points ([`POSITIVE_PREFERENCE`](@ref CubatureRules.POSITIVE_PREFERENCE)), the
positive one moves to the front. `rule` builds the first candidate.

The number of points is the main criterion because it is what determines the cost of using
the rule. The construction cost is reported by `compare` but not used for ranking, since a
rule is usually built once and applied many times.

Points alone would hand out rules that lose digits. Grundmann–Möller at degree 31 on the
tetrahedron has 3876 points against the conical product's 4096, but its weights alternate
in sign and amplify rounding errors in the integrand by `Σ|w|/Σw = 1.3e5`: `exp`
integrated to 7e-11 instead of 8e-16. `available` shows that factor in its `amplification`
column, from [`weight_amplification`](@ref CubatureRules.weight_amplification), without
building anything. When a rule with negative weights is still chosen because the positive
alternative is much larger — Grundmann–Möller in four dimensions, where the conical product
needs five times the points — `rule` says so and names the alternative; `positive = true`
takes it.

The ranking is recorded in `provenance(r).selection`.

## Families that are not selected automatically

A family can define `selectable(f) = false`. It is then listed by `available` but never
chosen by `rule(domain; degree)` unless the caller passes `copyleft = true` or names the
family. This is used for families whose rules come with licence terms the caller may not
want, such as rules refined from the GPL tables in Lebedev.jl; see
[Provenance and licensing](provenance.md).

When such a family would have given a rule with fewer points than the one chosen, `rule`
prints a warning, and so it does when a package that is not loaded would have: Lebedev.jl
on the sphere (see [`unloaded_alternatives`](@ref CubatureRules.unloaded_alternatives)). A
rule with negative weights is never suggested in place of a positive one. Each combination
of family and degree warns once; `CubatureRules.selection_warnings!(false)` turns the
warnings off.

## Explicit families

`rule(family, domain; ...)` skips the ranking. It still checks that the family supports the
domain, degree and output type, and throws a [`NoRuleError`](@ref) otherwise, naming the
nearest degree the family does offer.

## Error messages

When no candidate is left, the [`NoRuleError`](@ref) lists every candidate with the reason it
was rejected, and the nearest request that would succeed. When a family needs a package
that is not loaded, the message names the package.
