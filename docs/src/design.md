# Design

## Seed → refine → certify

Every rule goes through one pipeline:

- **Seed**: a cheap Float64 approximation. It comes from a stored table, from the rule's
  orbit structure plus multistart search, or from Golub–Welsch for Gauss rules.
- **Refine**: Newton or Gauss–Newton on the rule's *defining equations*, in BigFloat at the
  requested precision plus guard digits. The result is rounded exactly once at the end.
  MPFR arithmetic is correctly rounded, which is why `rule(...)` gives bitwise-identical
  output on every platform.
- **Certify**: the defining-equation residual of the delivered rule is evaluated at twice
  the precision and attached as a [`Certificate`](@ref).

## Certificate versus verification

A `Certificate` answers *did Newton converge on the system it was given?* A
[`Verification`](@ref) answers *is this rule what it claims to be?* They are different
checks. A rule seeded from a mistranscribed table can have a tiny defining-equation
residual and still be the wrong rule, and only verification catches that.

## Exactness claims

A rule claims [`PolynomialDegree`](@ref)`(d)`, [`SpanOf`](@ref)`(basis)` or
[`NoClaim`](@ref)`()`. Transport operations state their effect on the claim:

| Operation | Effect on `PolynomialDegree(d)` |
|---|---|
| `map_to(r, dom)` (affine only) | preserved |
| `subdivide(r, dom, n)` | preserved |
| `transform(r, φ, Jφ)` | destroyed → `NoClaim` (unless asserted by the caller, which is recorded) |
| `duffy(r)` | destroyed → `NoClaim` |

## Guard digits

Guard digits are measured, not assumed. The refiner estimates the Jacobian's condition
number and carries `32 + log₂ κ` extra bits, rounded up to a multiple of 8 so that small
cross-platform differences in the estimate cannot change the working precision. Both the
condition number and the guard are recorded in the certificate.

## Symmetric rules and the orbit algebra

A fully symmetric rule on a `D`-simplex is a list of orbits of `S_{D+1}`, which permutes the
barycentric coordinates. Each orbit type is a pattern of equal coordinates. On the triangle
the types are the centroid `[3]`, the vertex type `[2,1]` and the general type `[1,1,1]`.

The moment equations are projected onto the invariant subspace. That subspace is the
range of the Reynolds operator on each orthonormal degree block, and its dimension is
checked against the Molien series. The system has exactly as many equations as there are
independent invariant moments. That basis is computed in Float64. An inexact basis changes
only the conditioning, never the converged rule, because a symmetric rule's residual
already lies in the invariant subspace.

## No caching

`rule(...)` constructs and returns; it does not memoise. Rules are immutable values with
content-based `==` and `hash`, so a user who wants memoisation can use a `Dict`.
