# The pipeline

Every rule is produced by the same three steps.

1. **Seed.** A starting approximation in `Float64`. Depending on the family it comes from a
   stored table, from a search that starts from the orbit structure alone, from the rule one
   degree lower, or from the Golub–Welsch eigenvalue method for Gauss rules. See
   [Seed strategies](seeds.md).
2. **Refine.** Newton's method or Gauss–Newton on the equations that define the rule, in
   `BigFloat` at the requested precision plus guard digits. The number of guard digits is
   derived from the condition number of the system. See
   [Precision and guard digits](precision.md).
3. **Certify.** The residual of the defining equations is computed for the final rule, at
   twice the working precision, and stored in a [`Certificate`](@ref
   CubatureRules.Certificate) together with the condition number, the guard and the number
   of iterations. See [Certificates](certificates.md).

The result is rounded to the output type once, at the end, and returned as an immutable
[`QuadratureRule`](@ref) with a [`Provenance`](@ref CubatureRules.Provenance) describing all
of the above.

There is one shortcut. A seeded family's table holds each rule as the correctly rounded
`Float64` values of the refined rule, with its residual recorded when the table was
checked. At 53 bits or fewer there is nothing for steps 2 and 3 to add, so a `Float64`
request returns the stored rule, and a `Float32` or `Float16` request rounds it. The
certificate then quotes the recorded residual. See
[Seed strategies](seeds.md#Float64-requests). Any higher precision goes through all three
steps.

Verification is separate from this pipeline. The certificate describes how the solver
converged; [`check`](@ref) tests the finished rule against an independent basis. See
[Verification](verification.md).

## Derived and seeded families

Families come in two kinds, recorded by `derivation(F)`:

- **Derived** families compute their rules from a formula or a well-conditioned algorithm:
  Gauss–Jacobi from the three-term recurrence, Grundmann–Möller from a closed formula,
  product rules from their factors. They exist at every degree.
- **Seeded** families start from stored data and, above `Float64`, refine it: the minimal
  symmetric rules on triangles, tetrahedra and spheres. They exist only at the degrees for
  which a seed is available.

When two candidates have the same number of points, the selector prefers the derived one.

## No caching

`rule(...)` constructs the rule and returns it. Nothing is memoised. Rules are immutable
values with content-based `==` and `hash`, so if you want to reuse rules across calls you
can store them in a `Dict` yourself.

The reason is predictability: a function that sometimes takes a second and sometimes takes
a microsecond, depending on what was called before, is hard to reason about in a hot loop.
It is also why there is no `integrate(f, domain; degree)`.

## Determinism

The same call gives bitwise identical output on every platform. This follows from three
choices:

- All arithmetic that determines the result is done either in MPFR (`BigFloat`), which is
  correctly rounded, or in plain `Float64` operations in a fixed order. Both are platform
  independent.
- The linear solves in the refinement use factorisations written out in Julia rather than
  LAPACK or BLAS, whose results can differ between machines. A well-conditioned Newton step
  is solved with a `Float64` QR and refined in `BigFloat`; an ill-conditioned one uses a QR
  in `BigFloat`.
- Guard digits are rounded up to a multiple of 8 bits, so that small platform differences in
  the `Float64` condition estimate cannot change the working precision.

The test suite checks this with pinned SHA hashes of a set of reference rules.

## Cancellation and progress

A long construction can be stopped from another task with a
[`CancellationToken`](@ref CubatureRules.CancellationToken): pass it as `cancel`, and call
`cancel!` on it. The refinement checks the token once per iteration and throws a
`CancelledError`.

`verbose = true` reports what is about to be computed before the expensive part starts, and
then one line per iteration. The messages are logged with `@info`, so the usual logging
tools can redirect or filter them. Both the token and the verbosity level travel to the
family's `build` method in the [`BuildContext`](@ref CubatureRules.BuildContext).
