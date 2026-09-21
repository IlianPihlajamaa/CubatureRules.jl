# Rule families

| Family | Domain | Degrees | Derivation | Points | Output types |
|---|---|---|---|---|---|
| `XiaoGimbutas` | triangle | 1–40 | seeded | Xiao–Gimbutas counts, or fewer (degrees 27–33, 38–40) | floating |
| `FullySymmetric` | tetrahedron | 1–15 | seeded | smallest found by in-house search | floating |
| `GrundmannMöller` | any simplex | any (odd) | derived | ``\sum_{i=0}^{s} \binom{s-i+D}{D}`` | floating, `Rational{BigInt}` |
| `ConicalProduct` | any simplex | any | derived | ``\lceil (d+1)/2 \rceil^D`` | floating |
| `TensorProduct` | box (`Orthotope`) | any | derived | ``\prod_i m_i`` | floating, exact over exact factors |
| `GaussJacobi` | interval, Jacobi weight | any | derived | ``\lceil (d+1)/2 \rceil`` | floating |
| `NewtonCotes` | interval | any | derived | `d`, or `d+1` when even | floating, `Rational{BigInt}` |
| `Fejer` | interval | any | derived | `d`, or `d+1` when even | floating |
| `TanhSinh` | interval | none (`NoClaim`) | derived | set by the level | floating |

## XiaoGimbutas

```@docs
XiaoGimbutas
```

Fully symmetric, positive and interior triangle rules at the minimal point counts of
H. Xiao and Z. Gimbutas, *Comput. Math. Appl.* 59 (2010) 663–676,
doi:10.1016/j.camwa.2009.10.027.

The seeds are MIT-licensed. They were generated in-house from orbit structures alone, and
no published numbers were used. Above degree 20 they come from node elimination, which at
degree 27 reached a 139-point rule where the paper reports 141. Minimal rules are not unique; the shipped rule is the one
with the largest minimum barycentric coordinate among those found.

## FullySymmetric

```@docs
FullySymmetric
```

Fully symmetric (S₄), positive and interior tetrahedron rules. Unlike the triangle table,
nothing here comes from a paper. The point count at each degree is the smallest that the
in-house search reached, walking upwards from a count below which no fully symmetric rule
can exist. The counts are not proven minimal, and the family is not named after any
published one.

## GrundmannMöller

```@docs
GrundmannMöller
```

A. Grundmann and H. M. Möller, *SIAM J. Numer. Anal.* 15 (1978) 282–290,
doi:10.1137/0715019.

## ConicalProduct

```@docs
ConicalProduct
```

A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall, 1971.

## GaussJacobi

```@docs
GaussJacobi
GaussLegendre
```

Seeds come from Golub–Welsch (G. H. Golub and J. H. Welsch, *Math. Comp.* 23 (1969)
221–230). They are refined by Newton on the three-term recurrence, and the weights come
from the Christoffel function.

## TensorProduct

```@docs
TensorProduct
```

Tensor products of 1D rules on an [`Orthotope`](@ref), one family per axis, built from any
1D family. Two rules can also be combined directly with `r₁ ⊗ r₂`.

The claim is deliberately conservative. A tensor rule is exact on a tensor-product
polynomial space, which is larger than the total-degree space it reports — a 4×4 Gauss rule
integrates `x⁶y⁶` exactly but claims only total degree 7. The selector ranks candidates by
node count against a total degree, so that is what is reported, and `show` never prints
something ambiguous.

## TanhSinh

```@docs
TanhSinh
```

The family that forces the [`ExactnessClaim`](@ref) hierarchy to exist. Tanh-sinh is exact
on no polynomial space — its claim is a convergence *rate* — so it reports
[`NoClaim`](@ref), is never offered for a `degree` request, and is checked by a convergence
sweep rather than by exact integration:

```julia
seq = [rule(TanhSinh(m), Interval()) for m in 2:6]
CubatureRules.verify_convergence(seq, x -> 1/sqrt(1 - x^2), π)
```

It earns its place on integrands Gauss rules handle badly. On `1/√(1-x²)` at 201 points,
tanh-sinh reaches 5e-8 where Gauss–Legendre manages 9e-3.

Accuracy is limited by how well `1 - x` survives in floating point. The outermost node sits
a few units of roundoff from the endpoint, so an integrand evaluated at `x` near ±1 loses
about half the working digits — 3e-8 in `Float64`, 3e-26 at 50 digits. Ask for more digits,
or substitute so that the integrand is written in terms of the distance to the endpoint.
