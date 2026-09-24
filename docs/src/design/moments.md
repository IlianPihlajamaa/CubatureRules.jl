# Measures given by moments

A [`MomentWeight`](@ref) describes a measure `w(x) dx` on an interval by its *modified
moments*

```math
m_k = \int \pi_k(x)\, w(x)\, dx, \qquad k = 0, 1, 2, \dots
```

against an auxiliary family of monic polynomials `πₖ` with a known three-term recurrence.
From `2n` moments, Wheeler's algorithm computes the first `n` recurrence coefficients of the
measure's own orthogonal polynomials, and the Gauss rule follows from those in the usual way.
This is implemented by the [`ModifiedChebyshev`](@ref) family.

## Conditioning

With the monomials as auxiliary family, the moments are the ordinary moments `∫ xᵏ w dx`, and
the problem is badly conditioned. This is the classical Chebyshev algorithm, and the reason
most libraries (PolyChaos.jl, for example) use the Stieltjes or Lanczos procedures instead.

Recovering Gauss–Legendre on `[0, 1]` from `mₖ = 1/(k+1)` and comparing with the rule
computed directly gives these errors in the nodes:

| points | `Float64` | 200 digits |
|---|---|---|
| 4  | 8.2e-14 | 1.1e-198 |
| 8  | 1.4e-08 | 1.9e-193 |
| 12 | 7.1e-04 | 4.6e-188 |
| 16 | fails (β₁₃ < 0) | 1.2e-181 |
| 25 | fails | 1.9e-168 |
| 40 | fails | 5.7e-146 |

About 1.4 decimal digits are lost per point in both columns. In `Float64` that leaves nothing
after 11 or 12 points, and from 16 points on a recurrence coefficient becomes negative, so
the result no longer describes a measure. At 200 digits the same algorithm is usable up to
roughly 140 points.

With an auxiliary family that is orthogonal on the same interval, the problem is well
conditioned. For `w(x) = log(1/x)` on `[0, 1]` with monic shifted Legendre polynomials, the
error against the exact values `∫₀¹ xʲ log(1/x) dx = 1/(j+1)²` is:

| points | `Float64` | 200 digits |
|---|---|---|
| 4  | 5.2e-16 | 1.1e-200 |
| 16 | 2.8e-15 | 2.9e-200 |
| 50 | 6.3e-15 | 2.1e-199 |
| 80 | 7.9e-15 | 3.6e-199 |

So a good auxiliary family should be used whenever one is available. Ordinary moments are
supported for the cases where nothing else is known, and there the extra precision is what
makes the computation possible.

## Choosing the working precision

The condition number cannot be estimated cheaply beforehand, so `ModifiedChebyshev` finds the
working precision by experiment. It computes the rule at one precision and again 64 bits
higher, and doubles the difference until two successive results agree to the requested
accuracy. The more accurate of the two is returned.

The residual of the moment equations is not used as a stopping criterion, because it does
not reflect the accuracy of the nodes when the problem is badly conditioned. At 40 points
from ordinary moments, a rule with nodes wrong in the sixth decimal still satisfies all 80
moment equations to about `10⁻⁵²`. An early version of the code stopped on this residual
and returned such a rule with a certificate claiming 50 digits.

The certificate reports the measured condition number: the difference between the last two
results divided by the rounding unit of the less precise one. For 40 points from ordinary
moments this is about `4.7 × 10⁵⁵`, which matches the 54 digits lost in the table above. For
the `log(1/x)` example it is about 2.4.

A breakdown of the recursion (`βₖ ≤ 0`) is treated as a sign of too little precision and
triggers a retry at higher precision. If the precision limit is reached, a
`RefinementError` is thrown, noting whether the recursion kept breaking down or the results
kept disagreeing.

## Requirements on the moments

The moment function is called with the number type to use, and at whatever precision the
search reaches. It must therefore compute the moments (from a formula, or as exact
fractions); stored `Float64` values cannot be refined, and a measure known only to `Float64`
cannot give a rule accurate to more digits than that.

If the moments do not come from a positive measure on the interval, the recursion breaks
down at every precision and the construction fails.

## Exact change of basis

Many weights with an endpoint singularity have ordinary moments that are exact rationals.
[`LogWeight`](@ref)`(α; power)`, the weight `x^α log(1/x)^m` on `[0, 1]`, is one:

```math
\int_0^1 x^{j+\alpha} \log(1/x)^m \, dx = \frac{m!}{(j + \alpha + 1)^{m+1}} ,
```

and a floating-point `α` is itself an exact rational. Fed to Wheeler's algorithm directly,
ordinary moments are the badly conditioned case above. But the change of basis to modified
moments against the shifted Legendre polynomials can be done *exactly*, in rational
arithmetic, before anything is rounded. The only rounding then happens after the change of
basis, and the problem is well conditioned: at 20 points the measured condition numbers are
between 1 and about 10⁴, and the rules reproduce all `2n` moments to the requested accuracy.

A moment weight can record the interval its moments were computed on with
`support = (a, b)`. `LogWeight` does; pairing it with any other interval then throws an
`ArgumentError` instead of producing a rule for the wrong measure.

## Christoffel modification

`CubatureRules.christoffel(domain, z; power)` gives the Gauss rules of `|x − z| dμ`
(`power = 1`) or `dμ / |x − z|` (`power = -1`) for a measure `μ` whose monic orthogonal
polynomials `πₖ` are known (Legendre, Jacobi, Laguerre) and a point `z` outside its support.
Both are modified-moment problems against the `πₖ` themselves, which is the well-conditioned
case.

Multiplying needs only two moments. With `s` the sign of `x − z` on the support, and
`x − z = π₁ + (α₀ − z)`,

```math
\int \pi_k \, |x - z| \, d\mu = s \bigl( (\alpha_0 - z)\beta_0 \, \delta_{k0} + \beta_0\beta_1 \, \delta_{k1} \bigr).
```

Dividing gives `∫ πₖ dμ / |x − z| = −s ρₖ(z)`, where `ρₖ(z) = ∫ πₖ(x) dμ(x) / (z − x)` are the
functions of the second kind. They satisfy the same three-term recurrence as the `πₖ`, with
`ρ₋₁ = 1`, but as its minimal solution, so forward recursion loses them to rounding. They are
computed by backward recursion instead, which is a continued fraction, starting deeper until
the values no longer change. The convergence is geometric on a bounded interval, slower on
the half line, and slows as `z` approaches the support.

Checked against references that do not use the moments: multiplying Legendre by `1 + x`
reproduces Gauss–Jacobi(0, 1) to the last bit, and the other cases match integrals computed
with a 200-point rule of the base weight, or exact values on the half line, to the requested
precision. The measured condition numbers are between 1 and about 10⁵.

## Helpers

- [`LogWeight`](@ref) gives `x^α log(1/x)^power` on `[0, 1]`, as above.
- [`OrdinaryMoments`](@ref) builds a `MomentWeight` with the monomials as auxiliary family.
- `CubatureRules.shifted_legendre_recurrence()` is the monic shifted Legendre family on
  `[0, 1]` with exact rational coefficients.
- `CubatureRules.monic(rec)` converts one of the package's orthonormal recurrences
  (`jacobi_recurrence`, `laguerre_recurrence`, `hermite_recurrence`) into monic form.
- `CubatureRules.shift(aux, lo, hi)` moves a monic family from `[-1, 1]` to `[lo, hi]`.
- `CubatureRules.wheeler(w, n, T)` runs Wheeler's algorithm on its own and returns the
  recurrence coefficients.

## References

- J. C. Wheeler, Modified moments and Gaussian quadratures, *Rocky Mountain J. Math.* 4
  (1974) 287–296.
- W. Gautschi, On the construction of Gaussian quadrature rules from modified moments,
  *Math. Comp.* 24 (1970) 245–260.
- R. A. Sack and A. F. Donovan, An algorithm for Gaussian quadrature given modified
  moments, *Numer. Math.* 18 (1972) 465–478.
