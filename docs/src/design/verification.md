# Verification

[`check`](@ref)`(r)` tests a rule against its claim and returns a [`Verification`](@ref CubatureRules.Verification).
[`verify`](@ref)`(r; degree, bits)` does the same with explicit settings.

```@repl ver
using CubatureRules
v = check(rule(Sphere{3}(); degree = 11))
```

For a [`PolynomialDegree`](@ref)`(d)` claim there are three parts.

## Exactness

Every function in a basis of the polynomials of degree `≤ d` is integrated with the rule and
compared with its exact integral. The basis depends on the domain:

| Domain | Basis |
|---|---|
| `Interval` | orthonormal Legendre polynomials |
| `Interval` with a Jacobi weight | orthonormal Jacobi polynomials |
| `LaguerreRay`, `HermiteLine` | the weight's own orthonormal polynomials, by their recurrence |
| `Simplex{2}` | orthonormal Dubiner polynomials |
| `Simplex{3}` | orthonormal tetrahedral Dubiner polynomials |
| `Simplex{D}`, `D ≥ 4` | barycentric monomials with exact Dirichlet moments |
| `Orthotope{D}` | tensor Legendre polynomials, graded by total degree |
| `Sphere{2}` | Fourier modes on the circle |
| `Sphere{3}` | real spherical harmonics |
| `Sphere{D}`, `D ≥ 4` | monomials with exact sphere moments |
| `Ball{D}` | monomials with exact ball moments |
| `GaussianSpace(D)` | monomials with exact Gaussian moments |
| moment-defined weight | the auxiliary polynomials, against the given moments |

Orthonormal bases are used where they exist because they keep the check well conditioned.
Where monomials are used, it is either because the moments are exact rationals (simplices,
balls) or because no orthonormal basis is implemented. On a sphere, monomials are linearly
dependent (`Σxᵢ² = 1`), so there they are a spanning set rather than a basis; this costs
redundant equations but does not affect the result.

The check runs at twice the rule's precision. Each residual is compared with the error that
rounding the delivered nodes and weights could cause, computed from the size of the basis
function and its gradient at the nodes. Rules in exact rational arithmetic must have a
residual of exactly zero.

## Sharpness

The rule must *not* be exact at degree `d + 1`. This keeps a rule from reporting a lower
degree than it has. The largest residual at degree `d + 1` is compared with the tolerance:
more than 100 times the tolerance counts as not exact (`sharp = true`), at most the
tolerance counts as exact (`sharp = false`, and verification fails), and anything in
between is reported as inconclusive (`sharp = nothing`).

This is also why families whose rules are automatically exact at the next odd degree, such
as the centrally symmetric Lebedev rules, only offer odd degrees: an even-degree request is
answered with the odd rule above it, which is then sharp.

## Structure

- The weights add up to `measure(domain)`.
- The nodes lie in the domain, and strictly inside it if the rule claims interior nodes.
- The weights are positive, if the rule claims positive weights.
- The rule has its claimed symmetry exactly: applying each group element maps the node set
  onto itself, with equal weights.

## Rules without a claim

A rule with [`NoClaim`](@ref) has nothing to integrate exactly. For such families,
`CubatureRules.verify_convergence(rules, f, reference)` integrates a function with known
integral using a sequence of rules and checks that the error decreases towards the
precision limit. The resulting `Verification` has `empirical = true`.

## Verification and construction are independent

The verification bases are not the equations the rules were solved from. The symmetric
triangle rules, for example, are solved in terms of `S₃`-invariant moments of orbit
parameters, and verified against the full Dubiner basis on the expanded node set. The
Lebedev rules are solved in terms of the invariants `p₄` and `p₆`, and verified against all
spherical harmonics. A mistake in the construction therefore shows up in the verification.

## In test suites

```julia
using Test
@test_exact r 10      # exact to degree 10 and not to 11
```
