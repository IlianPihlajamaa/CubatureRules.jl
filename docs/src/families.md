# Rule families

| Family | Domain | Degrees | Derivation | Points | Output types |
|---|---|---|---|---|---|
| `XiaoGimbutas` | triangle | 1–46 | seeded | Xiao–Gimbutas counts, or fewer (degrees 27–33 and 35–46) | floating |
| `FullySymmetric` | tetrahedron | 1–16 | seeded | smallest found by in-house search | floating |
| `GrundmannMöller` | any simplex | any (odd) | derived | ``\sum_{i=0}^{s} \binom{s-i+D}{D}`` | floating, `Rational{BigInt}` |
| `ConicalProduct` | any simplex | any | derived | ``\lceil (d+1)/2 \rceil^D`` | floating |
| `TensorProduct` | box (`Orthotope`) | any | derived | ``\prod_i m_i`` | floating, exact over exact factors |
| `GaussJacobi` | interval, Jacobi weight | any | derived | ``\lceil (d+1)/2 \rceil`` | floating |
| `GaussLaguerre` | `LaguerreRay(α)` | any | derived | ``\lceil (d+1)/2 \rceil`` | floating |
| `GaussHermite` | `HermiteLine()` | any | derived | ``\lceil (d+1)/2 \rceil`` | floating |
| `NewtonCotes` | interval | any | derived | `d`, or `d+1` when even | floating, `Rational{BigInt}` |
| `Fejer` | interval | any | derived | `d`, or `d+1` when even | floating |
| `TanhSinh` | interval | none (`NoClaim`) | derived | set by the level | floating |
| `ExpSinh` | `HalfLine()` | none (`NoClaim`) | derived | set by the level | floating |
| `SinhSinh` | `RealLine()` | none (`NoClaim`) | derived | set by the level | floating |
| `SphereProduct` | `Sphere{2}`, `Sphere{3}` | any | derived | `d+1` on the circle, ``\lceil (d+1)/2 \rceil (d+1)`` on the sphere | floating |

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

## GaussLaguerre and GaussHermite

```@docs
GaussLaguerre
GaussHermite
```

The two classical rules on unbounded domains, and the first families whose domain is
weighted by something other than a Jacobi weight:

| Family | Domain | Weight | ``\sum_i w_i`` |
|---|---|---|---|
| `GaussLaguerre(α)` | [`LaguerreRay`](@ref)`(α)` = ``[0, \infty)`` | ``x^\alpha e^{-x}``, ``\alpha > -1`` | ``\Gamma(\alpha+1)`` |
| `GaussHermite()` | [`HermiteLine`](@ref)`()` = ``\mathbb{R}`` | ``e^{-x^2}`` | ``\sqrt{\pi}`` |

The weight is part of the domain, not of the integrand, so `integrate` applies it for you:

```julia
r = rule(LaguerreRay(); degree = 11)     # 6 points
integrate(x -> x^5, r)                   # ∫₀^∞ x⁵ e^{-x} dx = 120
```

Both share the driver behind [`GaussJacobi`](@ref) — Golub–Welsch in `Float64` for the
seed, Newton on the three-term recurrence at working precision, Christoffel weights — so a
family only supplies ``a_k``, ``b_k`` and ``\mu_0``. Laguerre nodes spread out to about
``4n`` rather than staying in ``[-1, 1]``, so the guard allowance grows with ``n``; the
degree-15 rule still reproduces ``7! `` to 60 digits.

Verification cannot use a monomial basis here — ``\int x^k w`` grows factorially and the
conditioning with it — so these domains verify against the weight's *own* orthonormal
polynomials, evaluated by the same recurrence.

The unbounded base domains [`HalfLine`](@ref) and [`RealLine`](@ref) report
`measure(...) == Inf` and exist to be weighted; a rule's weights sum to the weight's mass,
which is finite.

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

## SphereProduct

```@docs
SphereProduct
```

The fallback on a sphere, as [`ConicalProduct`](@ref) is on a simplex: available at every
degree and every precision, positive throughout, and needing no table.

On ``S^2``, with ``t = \cos\theta`` the surface integral separates,

```math
\int_{S^2} f \, d\sigma = \int_{-1}^{1}\!\!\int_0^{2\pi} f \, d\varphi \, dt,
```

into Gauss–Legendre in ``t`` against the trapezoid rule in ``\varphi``. The trapezoid rule
is spectrally exact on a circle — ``m`` equally spaced points integrate ``e^{ik\varphi}``
exactly for ``|k| < m`` — so ``m = d+1`` and ``n = \lceil (d+1)/2 \rceil`` give degree
``d`` in ``\lceil (d+1)/2 \rceil (d+1)`` points, roughly twice a minimal rule.

On the circle the trapezoid rule is the whole story, and there it is optimal: `d+1` points
for degree `d`, which nothing can beat.

```julia
r = rule(Sphere{3}(); degree = 15)        # 128 points
integrate(x -> x[3]^2, r)                 # 4π/3
rule(Sphere((1.0, 2.0, 3.0), 2.5); degree = 9)    # a similarity preserves the claim
```

Rules on a sphere are verified against real spherical harmonics, not monomials: monomials
restricted to a sphere are linearly dependent (``\sum x_i^2 = 1``), while the harmonics of
degree ``\le d`` span exactly the polynomials of degree ``\le d`` restricted to it, which
is what a degree claim on a sphere means.

## ExpSinh and SinhSinh

```@docs
ExpSinh
SinhSinh
```

The same double-exponential idea on the unbounded domains: exp-sinh on
[`HalfLine`](@ref)`()` and sinh-sinh on [`RealLine`](@ref)`()`, both `NoClaim`, both nested
in the level.

```julia
r = rule(ExpSinh(5), HalfLine())          # 289 points
integrate(x -> exp(-x)/sqrt(x), r)        # √π, to 4e-16 — singular at 0, decaying at ∞
```

Truncation is the one real design choice, and it differs from tanh-sinh. On a finite
interval the sum can stop where the weight falls below the working resolution, because the
missing tail is then provably negligible. Here the weights *grow* double-exponentially, so
no such bound exists without knowing how fast the integrand decays — which is the premise
of the transformation, not something the rule can check. These rules therefore truncate at
a fixed dynamic range: nodes span ``[2^{-2(p-2)}, 2^{2(p-2)}]`` in magnitude at output
precision ``p``, about ``10^{\pm 31}`` in `Float64`.

The range is squared rather than plain because what truncation drops is a head and a tail,
not a weight. An integrand like ``x^{-1+\delta}`` at the origin leaves ``x_{\min}^\delta``
behind, and one decaying like ``x^{-1-\delta}`` leaves ``x_{\max}^{-\delta}``. The worst
case worth serving is ``\delta = 1/2``: over the plain range it caps the accuracy of
``\int_0^\infty x^{-1/2} e^{-x}\,dx`` at 5e-8, and over the squared range at roundoff, for
18% more nodes — ``t`` enters through `asinh`, so a squared range is a constant more of it.

What this does not do is rescue an integrand that does not decay, or one whose mass sits
outside the range. The certificate only says the emitted numbers satisfy the closed form;
the convergence sweep is what tells you the rule has resolved *your* integrand.
