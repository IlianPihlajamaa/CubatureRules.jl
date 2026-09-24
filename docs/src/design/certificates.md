# Certificates

A [`Certificate`](@ref CubatureRules.Certificate) records how the refinement of a rule went.

```@repl cert
using CubatureRules
certificate(rule(Simplex{2}(); degree = 12, digits = 40))
```

| Field | Content |
|---|---|
| `equations` | which system of equations defines the rule |
| `residual` | the residual of those equations for the delivered rule |
| `residual_bits` | the precision at which the residual was computed |
| `digits` | the target number of decimal digits |
| `guard_digits` | the extra digits used during the computation |
| `cond` | the condition number of the system |
| `iterations` | the number of Newton or Gauss–Newton steps |

The residual is computed after rounding to the output type, at a precision well above it, so
it measures the rule you actually receive.

## Certificates and verification answer different questions

A certificate says whether the solver converged on the system it was given. A
[`Verification`](@ref CubatureRules.Verification) says whether the resulting rule has the properties it claims. Both
are needed, because each can pass while the other fails.

A certificate can look perfect for a wrong rule. If a seed table contained a typo, the
refinement could converge to a different solution of the same equations: small residual,
reasonable condition number, but not the intended rule. Verification tests the rule against
an independent basis and would catch this.

A small residual can also hide a large error when the equations are badly conditioned. The
clearest example in this package is a Gauss rule computed from ordinary moments: at 40
points, a rule whose nodes are wrong in the sixth decimal still satisfies all 80 moment
equations to about `10⁻⁵²`, because the condition number is around `10⁴⁶`. For this reason
the moment-based family does not use its residual to decide when to stop; see
[Measures given by moments](moments.md).

The condition number in the certificate is what connects the two: the error in the rule is
roughly the residual times the condition number.

## Families without refinement

Some families do not refine anything. A Grundmann–Möller rule, for example, is computed from
a closed formula. Its certificate says so, reports no iterations and `cond = 1`, and uses as
residual the relative rounding error of the weights:

```@repl cert
certificate(rule(GrundmannMöller(), Simplex{2}(); degree = 5))
```
