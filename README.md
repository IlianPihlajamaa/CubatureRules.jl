# CubatureRules.jl

Quadrature and cubature rules **generated on demand**, at arbitrary order and arbitrary
precision, instead of tabulated in double precision.

```julia
using CubatureRules

r = rule(Simplex{2}(); degree = 20, digits = 200)   # construct once
integrate(x -> exp(x[1] * x[2]), r)                 # apply
```

That first line returns a 79-point, fully symmetric, positive-weight, interior-node rule
on the triangle. It has the minimal point count of Xiao & Gimbutas (2010) and is correct to
200 digits. No such rule is published anywhere else.

```
julia> r
QuadratureRule{2,BigFloat} on Simplex{2}()
  family    : XiaoGimbutas (seeded, Newton-refined)
  exactness : polynomial degree 20 (claimed; `check(rule)` verifies)
  points    : 79, all interior
  weights   : all positive, Σw = 0.5
  precision : 200 digits (guard 12, est. cond(J) = 1.6e+02)
  residual  : 7.6e-201  (defining equations at 413-digit arithmetic: S₃-invariant moment system …)
  reference : Xiao & Gimbutas (2010), doi:10.1016/j.camwa.2009.10.027

julia> check(r)
Verification: PASSED
  basis     : orthonormal Dubiner, 1332-bit arithmetic
  exactness : degree 20 ✓  (max residual 1.00e-200, tolerance 4.26e-198)
  sharpness : not exact at degree 21 ✓  (residual 3.56e-01)
  structure : Σw = measure ✓, interior ✓, positive ✓, symmetric ✓
```

## Seed → refine → certify

Every rule goes through the same pipeline:

1. **Seed**: a cheap Float64 starting point, from a stored table, the rule's orbit structure
   plus multistart search, or Golub–Welsch.
2. **Refine**: Newton or Gauss–Newton on the rule's defining equations, in BigFloat. The
   solve is rank-revealing. Guard digits come from the measured cond(J), and a
   cancellation token is checked every iteration.
3. **Certify**: the defining-equation residual of the *delivered* rule is evaluated at
   twice the precision and attached as a `Certificate`.

**Verification is a separate check.** `check(rule)` tests exactness against an
orthonormal basis of the claimed space. It also tests sharpness (the rule is *not* exact
one degree higher), the weight sum, interiority, positivity and the claimed symmetry.

## Asking for what you need, not for a scheme name

```julia
julia> available(Simplex{2}(); degree = 17)
Candidates on Simplex{2}() for degree ≥ 17, ranked:
  family                       npoints  degree  derivation  positive  interior  symmetry
  XiaoGimbutas                 60       17      seeded      yes       yes       S3
  ConicalProduct(GaussJacobi)  81       17      derived     yes       yes       none
  GrundmannMöller              165      17      derived     no        yes       S3

julia> rule(Simplex{2}(); degree = 9, T = Rational{BigInt})     # exact rationals → Grundmann–Möller
julia> rule(Simplex{3}(); degree = 11, digits = 50)             # any dimension
julia> rule(Simplex{2}())                                       # errors, listing what is available
```

Families are found with `subtypes(RuleFamily)` at call time. A downstream package can add
a family by defining one type and a few methods, with no registration step.

## Status: v0.1 (the cornerstone release) in progress

| | |
|---|---|
| Domains | `Interval`, `Simplex{D}` (any D), `WeightedDomain` with Jacobi weights; the others are stubbed with the release that brings them |
| Families | `XiaoGimbutas` (triangle, degrees 1–27, any precision), `FullySymmetric` (tetrahedron, degrees 1–15, any precision), `GrundmannMöller` (any D, exact rational), `ConicalProduct` (any D, any degree), `GaussJacobi` |
| Claims | `PolynomialDegree`, `SpanOf`, `NoClaim`; claim preservation under `map_to`, `subdivide`, `transform`, `duffy` |
| Application | `integrate` (generic return types; allocation-free on `static(rule)`), over a domain, over a mesh (threaded on request; 1.14× a hand-written FEM kernel over 10⁶ triangles), batched |
| Verification | `check` / `verify` / `@test_exact`, dispatched on the claim |
| Tooling | `show`, `cite` (BibTeX/APA/plain), content `hash`/`==`, `rule_hash` for bitwise reproducibility, `benchmark_construction` |

See [ROADMAP.md](ROADMAP.md) for the stages and [PLAN.md](PLAN.md) for the design. The v0.0
feasibility spike that decided the go/no-go gates is written up in
[notes/v0.0-spike.md](notes/v0.0-spike.md), and the tetrahedron gate in
[notes/v0.2-tetrahedra.md](notes/v0.2-tetrahedra.md).

## Seed data and licensing

All shipped seed data is MIT-licensed and generated in-house from orbit structures alone
(`scripts/generate_triangle_seeds.jl`, `scripts/generate_tetrahedron_seeds.jl`). No numbers
were copied from any published table,
and none from quadpy. Every file in `src/data` has an entry in `src/data/PROVENANCE.toml`,
and CI fails the build if one is missing.

Point counts are the smallest our searches reached, not proven minima — at degree 27 on the
triangle, node elimination found a 139-point rule where Xiao & Gimbutas (2010) report 141
(also 148 against 150 at degree 28, and 157 against 159 at degree 29).
Minimal symmetric rules are also **not unique**. The shipped rule at each degree is the most
interior of the valid rules the search found, so it need not match the published table
node for node. To refine a published table whose licence allows it, pass
`seed = ExplicitSeed(θ; source = "...")`.

## Related packages

1D is largely covered elsewhere and is not reimplemented here beyond what the conical
product needs: [QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl),
[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) and
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).
Adjacent packages, none with arbitrary-precision generation: Lebedev.jl,
GrundmannMoeller.jl, HAdaptiveIntegration.jl, Cubature.jl, GaussQuadrature.jl.
