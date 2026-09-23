# Weights given by moments

Every other family in this package knows its weight in closed form. `GaussJacobi` knows that
the orthogonal polynomials of `(1-x)^α (1+x)^β` are the Jacobi polynomials and knows their
three-term recurrence; the rule follows from the recurrence in a few lines.

Most weights are not like that. If you need a Gauss rule for `w(x) = log(1/x)`, or for
`exp(-x) / (1 + x²)`, or for a weight that only exists as a table of moments from an
experiment, no classical family applies. What you usually *can* get is moments — integrals of
known polynomials against the weight, often in closed form.

Wheeler's algorithm turns `2n` of those into the `n`-point Gauss rule. This page is about why
that is normally considered a bad idea, and why arbitrary precision changes the answer.

## The short version

```julia
using CubatureRules

# Lebesgue measure on [0,1], described only by ∫₀¹ xᵏ dx = 1/(k+1)
w = OrdinaryMoments((k, T) -> one(T) / (k + 1); label = "Lebesgue on [0,1]")
r = rule(WeightedDomain(Interval(0, 1), w); degree = 79, digits = 50)
```

That returns the 40-point Gauss–Legendre rule on `[0,1]`, correct to 50 digits. The same
computation in `Float64` does not return a poor answer — it returns no answer at all, because
a recurrence coefficient goes negative and the output stops describing a measure.

## Why this is supposed to be impossible

The map from moments to recurrence coefficients is the Hilbert matrix in disguise. Its
inverse amplifies error by a factor that grows geometrically in the number of points, so the
classical Chebyshev algorithm is a standard example of a numerically hopeless procedure.
PolyChaos.jl, which solves the same problem, says so explicitly and implements Stieltjes and
Lanczos instead.

That reputation is deserved, and nothing here contradicts it. Recovering Gauss–Legendre on
`[0,1]` from `mₖ = 1/(k+1)`, comparing against the rule built directly:

| points | `Float64` | 200 digits |
|---|---|---|
| 4  | 8.2e-14 | 1.1e-198 |
| 8  | 1.4e-08 | 1.9e-193 |
| 12 | 7.1e-04 | 4.6e-188 |
| 16 | **β₁₃ < 0** | 1.2e-181 |
| 25 | **β₁₃ < 0** | 1.9e-168 |
| 40 | **β₁₃ < 0** | 5.7e-146 |

The conditioning is identical in both columns — roughly 1.4 decimal digits lost per point, in
a straight line. `Float64` has about 16 digits to spend, so it runs out at 11 or 12 points,
which is exactly where the table shows it failing. Two hundred digits buys about 140 points
before the same thing happens.

This is the package's central claim in one table. The algorithm is not rehabilitated; the
constraint that made it useless was the arithmetic, and the arithmetic is a parameter.

## A weight with no classical family

Ordinary moments are the worst case, chosen above to make a point. In practice you pick a
better auxiliary family and the conditioning largely goes away.

The modified moments of `w(x) = log(1/x)` against the monic shifted Legendre polynomials on
`[0,1]` are known in closed form, so:

```julia
using CubatureRules
import CubatureRules: monic, shift, jacobi_recurrence

# monic shifted Legendre on [0,1], from the Legendre recurrence the package already has
aux = shift(monic(jacobi_recurrence(0, 0)), 0, 1)

logw = MomentWeight(aux,
                    (k, T) -> k == 0 ? one(T) :
                              T((-1)^k) / (T(k) * (k + 1) * T(binomial(big(2k), big(k))));
                    label = "log(1/x) on [0,1]")

r = rule(WeightedDomain(Interval(0, 1), logw); degree = 79, digits = 50)
integrate(f, r)     # ∫₀¹ f(x) log(1/x) dx
```

At 40 points this reports a measured condition number of about 2.4 and 30 guard digits —
which is just the fixed cost of the stability check described below — against 4.7e55 and 164
guard digits for the same size of rule from ordinary moments. Choosing the auxiliary family well is worth more than any amount of arithmetic — the
point of the previous section is that you are no longer *required* to choose it well.

## How the precision is chosen

The driver is given a target accuracy, not a working precision, and has to find the second
from the first. The obvious way to check its work is wrong, and worth spelling out because it
is a trap that the certificate machinery would otherwise walk straight into.

The natural test is whether the rule reproduces the moments it was built from: does
`Σᵢ wᵢ πₖ(xᵢ) = mₖ` hold for `k < 2n`? It is the defining property, it is cheap, and it is
checked against the input data rather than against anything derived from it.

It is also nearly useless here. At 40 points from ordinary moments, a rule whose nodes are
wrong in the **sixth** decimal still reproduces all 80 moments to 10⁻⁵². That is not a bug in
the test; it is what a condition number of 10⁴⁶ means. Stopping when the residual looks small
would ship six correct digits under a certificate claiming fifty.

So [`ModifiedChebyshev`](@ref) does not stop on the residual. It computes the rule twice, at
working precisions 64 bits apart, and doubles the precision until the two agree to the
accuracy requested. Consecutive iterates bracket the error in the thing the caller actually
receives, and the conditioning is measured rather than assumed: the ratio of the error to the
rounding unit that produced it is reported as `certificate(r).cond`.

```julia
julia> c = certificate(rule(WeightedDomain(Interval(0, 1), w); degree = 79, digits = 50));

julia> c.cond, c.guard_digits
(4.699e55, 164)
```

Both failure modes are handled the same way. A breakdown in Wheeler's recursion — some
`βₖ ≤ 0`, which cannot happen for a positive measure — is treated as a symptom of too little
precision and retried, because that is what it usually is. If it survives to the precision
cap, or if successive precisions never agree, you get a `RefinementError` that says which.

## What you have to supply

The moment function is called at whatever precision the search reaches, so it must be able to
answer in `BigFloat` at the ambient precision. A closed form or an exact rational works; a
lookup of stored `Float64` values does not, and no amount of working precision will rescue it.
That is a real limitation rather than an implementation gap: if the measure is only known to
16 digits, a 50-digit rule for it does not exist to be computed.

If the moments do not come from a positive measure on the interval you named, the
construction fails rather than returning something meaningless. Verification then has nothing
to check, which is the intended outcome.

## Reference

```@docs
MomentWeight
OrdinaryMoments
ModifiedChebyshev
CubatureRules.MonicRecurrence
CubatureRules.monic
CubatureRules.shift
CubatureRules.monomial_recurrence
CubatureRules.wheeler
CubatureRules.MomentBreakdownError
```
