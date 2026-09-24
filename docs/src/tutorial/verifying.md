# I want to verify a rule

Every rule states which polynomials it integrates exactly. [`check`](@ref) tests that
statement independently:

```@repl v
using CubatureRules
r = rule(Simplex{2}(); degree = 12);
check(r)
```

`passed(check(r))` gives the outcome as a `Bool`.

## What is checked

1. **Exactness.** Every function in an orthonormal basis up to the claimed degree is
   integrated and compared with its exact integral. The basis depends on the domain:
   Dubiner polynomials on triangles, Legendre or Jacobi polynomials on intervals, spherical
   harmonics on the sphere, and monomials with exactly known integrals elsewhere. The check
   runs at twice the precision of the rule, and each error is compared with what rounding the
   nodes and weights to the rule's precision could explain. Exact rational rules must give
   an error of zero.
2. **Sharpness.** The rule must *not* be exact at one degree higher. Without this check, a
   rule could claim a lower degree than it actually has.
3. **Structure.** The weights add up to the size of the domain, the nodes lie inside it, the
   weights are positive, and the rule has the symmetry it claims, where these properties are
   claimed.

## Verification and the certificate

Each rule also has a certificate:

```@repl v
certificate(r)
```

The two answer different questions. The [certificate](../design/certificates.md) says
whether the solver converged: how small the residual of the defining equations is, and how
well conditioned they were. The verification says whether the result is correct: whether
the rule actually integrates the polynomials it should.

A rule can pass the first and fail the second, for instance when the seed came from a table
with a typo and the solver converged to a different rule. That is why the package does both.

## In a test suite

```julia
using Test, CubatureRules

@test_exact r 12            # exact to degree 12, and not to degree 13
@test passed(check(r))
```

## Rules without a degree

Families like `TanhSinh` are not exact for any polynomial space, so `check` has nothing to
test. Instead, [`verify_convergence`](@ref CubatureRules.verify_convergence) checks a
sequence of rules on an integral whose value you know:

```@repl v
I = 2log(big(2)) - 2;                  # ∫ log(1 - x) dx over [-1, 1]
seq = [rule(TanhSinh(m), Interval(); digits = 50) for m in 1:6];
CubatureRules.verify_convergence(seq, x -> log(1 - x), I)
```

It computes the error of each rule and checks exactly two things:

1. **No growth.** The error never grows by more than half from one rule to the next, except
   once it has reached the floor: four times the smallest error seen, or the rounding error
   of the last rule's sum. Above, the last two errors are at that floor.
2. **Last error.** The last error is at most `rtol` times `max(|I|, 1)`. By default `rtol`
   is `16√ε` at the rules' precision `ε`, about `1e-24` here.

It does *not* check how fast the errors fall, or that they fall at all once they are below
the target. It catches divergence and a wrong answer, not slow convergence. For example,
`1/√(1 - x²)` passes in `Float64` although its error never improves:

```@repl v
seq = [rule(TanhSinh(m), Interval()) for m in 2:6];
CubatureRules.verify_convergence(seq, x -> 1 / sqrt(1 - x^2), Float64(π))
```

Here the integrand blows up at the endpoints, where `1 - x` has lost half of its digits, so
every level stops at about `√ε ≈ 1e-8`. That is the reason for the default `rtol`. Pass a
smaller `rtol` if you need more.

The result is marked `empirical`, since it is based on observed errors rather than an exact
property.

## High-precision rules

```@repl v
r50 = rule(Simplex{2}(); degree = 12, digits = 50);
v = check(r50);
passed(v), Float64(v.max_residual), Float64(v.tolerance)
```

The tolerance follows the precision of the rule, so a 50-digit rule is checked to about 50
digits.

## Next steps

- [I want to cite a rule](citing.md)
- [Verification](../design/verification.md) describes the bases and tolerances in detail.
