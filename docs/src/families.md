# Rule families

| Family | Domain | Degrees | Derivation | Points | Output types |
|---|---|---|---|---|---|
| `XiaoGimbutas` | triangle | 1–47 | seeded | Xiao–Gimbutas counts, or fewer (degrees 27–33 and 35–47) | floating |
| `FullySymmetric` | tetrahedron | 1–20 | seeded | smallest found by in-house search | floating |
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
| `LebedevRule{LebedevJLSeeds}` | `Sphere{3}` | to 125, from Lebedev.jl | seeded elsewhere | Lebedev's counts | floating |
| `BallProduct` | `Ball{D}` | any | derived | radial × angular | floating |
| `Lebedev` | `Sphere{3}` | odd, tabulated | seeded | Lebedev counts (78 at degree 13) | floating |
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

## BallProduct

```@docs
BallProduct
```

A ball separates into a radius and a sphere,

```math
\int_{B^D} f \, dx = \int_0^1 r^{D-1} \!\! \int_{S^{D-1}} f(r\omega) \, d\sigma(\omega) \, dr,
```

and nothing couples the two factors: a monomial of degree `d` restricted to a ray is `r^d`
times a monomial on the sphere. So the radial factor is a Gauss–Jacobi rule for the weight
``r^{D-1}`` and the angular factor is any sphere family — which is why this is a combinator.

```julia
rule(Ball{3}(); degree = 9)                       # BallProduct(Lebedev), 190 points
rule(BallProduct(SphereProduct()), Ball{3}(); degree = 9)    # 250, no table needed
rule(Disk(); degree = 5)
```

The selector offers every combination and ranks them by node count, so `Ball{3}()` gets the
Lebedev-angular rule wherever a Lebedev seed exists and falls back to the product rule
elsewhere.

Rules on a ball are verified against monomials with exact moments: unlike on a sphere, the
coordinates of a ball satisfy no relation, so monomials really are a basis there. Their
moments come from the sphere's, divided by ``|\alpha| + D``.

## Lebedev

```@docs
LebedevRule
```

Octahedrally symmetric, positive-weight rules on the sphere: the `SphereProduct` costs about
twice the nodes of one of these, so the selector prefers Lebedev wherever a seed is shipped.

```julia
rule(Sphere{3}(); degree = 11)            # 50 points, not the product rule's 72
rule(LebedevRule(), Sphere{3}(); degree = 15, digits = 60)
```

The rules are invariant under the 48 signed permutations of the coordinates, whose orbits
are the classical six — 6, 12 and 8 fixed points, two 24-point families with one parameter
each, and a 48-point family with two. The moment system is written in the invariants
``p_4 = \sum x_i^4`` and ``p_6 = x^2y^2z^2`` rather than in spherical harmonics, because a
rule built from whole orbits annihilates every non-invariant harmonic identically: at degree
131 that would be 17424 equations of which a few hundred say anything.

Only odd degrees are tabulated. Every orbit is centrally symmetric, so the odd harmonics
integrate to zero whatever the parameters, and a rule of degree ``2k`` is automatically of
degree ``2k+1``; an even-degree request is answered by the odd rule above it.

The seeds are MIT-licensed and were generated in-house by searching point counts upward
from the orbit structures alone — no published table was used. The counts found agree with
Lebedev's at every degree searched except degree 13, where his 74-point rule has a negative
weight: the search requires positive weights, so it reports a 78-point rule instead. That
74-point structure was checked exhaustively; every solution found has a negative weight.

### Rules from an installed Lebedev.jl

The type parameter says where the seeds came from. `LebedevRule()` is `LebedevRule{InHouseSeeds}()`,
the table above. `UpstreamLebedev()` is `LebedevRule{LebedevJLSeeds}()`, the caller's own
[Lebedev.jl](https://github.com/stefabat/Lebedev.jl) — GPL-3, and tabulated to degree 125.
This package ships none of those numbers.

```julia
using CubatureRules
import Lebedev as LebedevJL       # aliased: this package exports a `Lebedev` type of its own

available(Sphere{3}(); degree = 19)
#  Lebedev (Lebedev.jl, GPL-3)   146 points
#  SphereProduct                 200 points

rule(Sphere{3}(); degree = 19)                   # SphereProduct, and a warning
rule(Sphere{3}(); degree = 19, copyleft = true)  # the 146-point rule, licence attached
rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 40)   # refined, still GPL-3
```

[`available`](@ref) lists both, with the licence carried in the family name. Only the
in-house one is chosen automatically: a rule whose terms this package cannot pass on should
not be the silent answer to a request that said only "degree 19". Ask for it by name, or
pass `copyleft = true`. When a cheaper rule is passed over for this reason, `rule` says so
once per case — see [`license_warnings!`](@ref CubatureRules.license_warnings!) to silence it.

Refinement is offered on both. A rule refined from Lebedev.jl's table is a derived work of
it, so the GPL-3 licence travels into the refined rule's `provenance` unchanged and the
derivation path records where it came from. What the seed policy forbids is the other thing:
absorbing someone else's numbers into `src/data` under this package's own licence. Nothing
from Lebedev.jl is ever written there.

The orbit structure is recovered from the raw points before refinement, by
`CubatureRules.classify_octahedral`, which decomposes any sphere rule into `O_h` orbits or
reports that it does not decompose. That is what makes an outside table refinable here, and
what will let rules from other sources be ingested later.

!!! note "The name `Lebedev`"
    This package exports a type called `Lebedev`, so `using Lebedev` collides with it. Write
    `import Lebedev as LebedevJL`; the extension activates either way, and the alias never
    has to be used afterwards.

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
