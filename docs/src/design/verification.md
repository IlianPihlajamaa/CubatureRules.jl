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
| `Simplex{D}`, `D ≥ 3` | orthonormal Dubiner polynomials of that dimension |
| `Orthotope{D}` | tensor Legendre polynomials, graded by total degree |
| `Sphere{2}` | Fourier modes on the circle |
| `Sphere{3}` | real spherical harmonics |
| `Sphere{D}`, `D ≥ 4` | monomials with exact sphere moments |
| `Disk()`, `Ball{3}()` | orthonormal polynomials: Jacobi in `r²` times Fourier modes or spherical harmonics |
| `Ball{D}`, `D ≥ 4` | monomials with exact ball moments |
| `GaussianSpace(D)` | tensor orthonormal Hermite polynomials, graded by total degree |
| `ExponentialSpace(D)` | monomials with exact moments of `e^{-r}` |
| moment-defined weight | the auxiliary polynomials, against the given moments |

Orthonormal bases are used where they exist because they keep the check well conditioned:
with monomials, functions that are nearly dependent can leave a large error in an
orthonormal direction while every monomial residual sits at rounding level. Monomials remain
where no orthonormal basis is implemented here: above the 2-sphere and the 3-ball there are
no harmonics in this package, and `e^{-r}` has no classical orthogonal polynomials in `r²`.
In exact arithmetic the simplices and boxes use the same polynomials unnormalised. On a
sphere, monomials are linearly dependent (`Σxᵢ² = 1`), so there they are a spanning set
rather than a basis; this costs redundant equations but does not affect the result.

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

In one dimension, the function tested for sharpness matters. Every function of degree
`d + 1` has the same error up to its leading coefficient, so the choice only affects how the
error compares with rounding. The top orthonormal polynomial `p_{d+1}` is a poor choice: for
an `n`-point Gauss rule its error is the ratio of leading coefficients `k_{2n} / k_n²`, which
for Laguerre is `(n!)² / (2n)! ≈ 4⁻ⁿ` and drops below the rounding tolerance once `n` exceeds
about 1.7 times the number of digits. The check then reported a correct rule as exact one
degree too high. Sharpness is therefore tested on `p_a p_b` with `a + b = d + 1`, whose
integral is 1 or 0. For a Gauss rule that is `p_n²`, and its error is exactly 1 for any
weight and any `n`.

Some rules come closer to the next degree than any choice of test function can fix. The
127-point Gauss–Patterson rule misses degree 192 by `1e-20`, so in `Float64` it is exact
there to within rounding. That is not the same as being exact, and the claim of degree 191 is
correct. A family whose rules behave like this measures the miss at working precision during
construction and records it in the certificate (`next_error`). When the delivered rule then
looks exact at `d + 1`, verification leaves sharpness undecided (`sharp = nothing`), says why
in `sharp_note`, and passes. A rule without that record still fails, as before. At enough
digits the miss becomes visible and sharpness is tested normally.

## Structure

- The weights add up to `measure(domain)`.
- The nodes lie in the domain, and strictly inside it if the rule claims interior nodes.
- The weights are positive, if the rule claims positive weights.
- The rule has its claimed symmetry exactly: applying each group element maps the node set
  onto itself, with equal weights. The group is a type (`provenance(r).symmetry`, a
  `SymmetryGroup`), and every group is checked the same way, by orbits: two nodes lie in one
  orbit exactly when their canonical keys agree (the coordinates the group permutes, sorted;
  those it reflects, made positive), so grouping the nodes by key finds the candidate orbits
  in `O(N log N)`, and the rule is invariant when every group has distinct members, as many
  as the orbit of its key, and a single weight.

## Symmetric rules on the sphere

Checked the general way, a Lebedev rule of degree 125 would be tested against 16 129
spherical harmonics at each of its 5294 nodes, and its symmetry by comparing each of the 48
images of every node with every node. At 100 digits that takes on the order of an hour and
a half. For a rule that claims `O_h` symmetry, `check` does two cheaper things instead.

**Symmetry.** By orbits, as for every group: two points of `ℝ³` lie in the same `O_h` orbit
exactly when their sorted absolute coordinates agree, and the orbits found give the
representatives and total weights the exactness check below works on.

**Exactness.** For an invariant rule, a polynomial and its group average integrate the same
way, both under the rule and over the sphere, so only invariant test functions matter. The
group average of a monomial `xᵅ` is zero when an exponent is odd, and otherwise depends only
on the exponents up to order. On the sphere, the monomials of one even degree `m` span all
even polynomials of degree `≤ m`, and odd polynomials integrate to zero on both sides by
central symmetry. So exactness to degree `d` is checked with the sorted even exponents of the
largest even degree `≤ d`, and sharpness with those of the next even degree. Each test
function is evaluated once per orbit, at the orbit's representative, and weighted by the
orbit's total weight.

At degree 125 this is 352 test functions for exactness (the number of invariants, as it
must be) and 363 for sharpness, over 132 orbits. The test functions are monomials on the
delivered nodes, not the invariants `p₄`, `p₆` in orbit parameters that the rule was solved
from, so the check is independent of the construction. At degree 59 it takes 0.3 s where the
general check takes 120 s.

If the nodes do not group into complete orbits, the general check is used, and the rule is
reported as not symmetric. `check(r; use_symmetry = false)` forces the general check.

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
