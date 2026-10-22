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
| `next_error` | how far the rule misses the next degree, measured at working precision, or `nothing` |

For the Gauss-type families and the seeded ones above `Float64`, the residual is computed
after rounding to the output type, at `2·target + guard` bits, so it measures the rule you
actually receive. Families that find their working precision by repeating the computation
(modified Chebyshev, Gauss–Patterson) report the residual at the precision they reached. The
discretised Stieltjes family and the Duffy families for singular kernels have no system of
equations that defines the rule: their `residual` is the change between the last two
resolutions of the discretisation, an estimate of the error rather than a residual.

`next_error` is recorded by families whose rules come so close to the next degree that a rule
delivered at ordinary precision cannot show the difference: the stored sphere rules,
Gauss–Patterson and Filon. A Gauss–Patterson rule misses degree `d + 1` by about `10⁻²⁰` at
127 points. When [`check`](@ref) finds such a rule exact at `d + 1`, it reports sharpness as
not resolvable at that precision, quoting this number, instead of failing.

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

## Rules returned from a table

A `Float64` request for a seeded family returns the stored rule without refining it (see
[Seed strategies](seeds.md#Float64-requests)). Its certificate reports no iterations and no
guard digits, and quotes the residual recorded when the table was checked, at 256 bits.
That is the residual of the rule you receive, since the stored and delivered parameters are
the same numbers:

```@repl cert
certificate(rule(Simplex{2}(); degree = 12))
```

For `Float32` and `Float16`, the residual of the rounded rule is evaluated in `Float64`.

## Families without refinement

Some families do not refine anything. A Grundmann–Möller rule, for example, is computed from
a closed formula. Its certificate says so, reports no iterations and `cond = 1`, and uses as
residual the relative rounding error of the weights:

```@repl cert
certificate(rule(GrundmannMöller(), Simplex{2}(); degree = 5))
```
