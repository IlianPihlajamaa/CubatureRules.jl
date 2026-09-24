# Precision and guard digits

## Output precision

The output type is chosen with `T` or `digits` (see
[I want arbitrary precision](../tutorial/precision.md)). Internally this becomes a
`BuildContext{T}` with a number of target bits: 53 for `Float64`, `ceil(d log₂ 10) + 1` for
`digits = d`, and the ambient precision for `T = BigFloat`.

## Working precision

Families compute in `BigFloat` at the target precision plus a number of *guard bits*. For
the families built on Newton or Gauss–Newton refinement, the guard is measured rather than
fixed: the refiner estimates the condition number `κ` of the Jacobian of the defining
equations and uses

```
guard = 32 + log₂ κ   (rounded up to a multiple of 8)
```

The rounding to multiples of 8 keeps small platform differences in the `Float64` estimate of
`κ` from changing the working precision.

If the condition number met during the refinement is larger than the one estimated at the
seed, the refinement is repeated once with a correspondingly larger guard. Both the
condition number and the guard are stored in the certificate.

For the Gauss families built on a three-term recurrence, the equations are well conditioned
and the guard is a function of the number of points only (`24 + 2 log₂(n+1)` bits).

For rules given by moments the condition number cannot be measured cheaply before the
computation, so the working precision is found by repeating the computation at increasing
precision until two successive results agree; see [Measures given by moments](moments.md).

## Estimating the condition number

`κ` is estimated from a `Float64` copy of the Jacobian, which is cheap. If that estimate is
above `10¹²`, where `Float64` cannot resolve it, the singular values are recomputed in
256-bit arithmetic. That fallback is expensive for large systems (a 352 × 352 matrix takes
longer than the solve itself), so the Gauss–Newton loop computes the accurate estimate only
at the first and last iterate and uses the free ratio `|R₁₁| / |R_rr|` from the pivoted QR
factorisation in between.

## Rounding once

All values are kept in `BigFloat` at working precision until the rule is complete, and then
rounded to the output type once, by `finalize_number`. Rounding intermediate results would
make the output depend on the order of operations, and so on the platform.

Because MPFR arithmetic is correctly rounded and the linear algebra is pure Julia, the same
call produces bitwise identical rules everywhere.

## Exact output

For `T = Rational{BigInt}` the family computes in exact rational arithmetic, and there is
no rounding. Only families with rational nodes and weights support this (`GrundmannMöller`,
`NewtonCotes`, tensor products of those); the others are excluded by `supports_type`.

## Other float types

Any `AbstractFloat` can be the output type. The rule is computed in `BigFloat` with enough
bits for that type and converted at the end, so `Float32`, `Double64` and `Float128` rules
are as accurate as their type allows.
