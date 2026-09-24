# CubatureRules.jl

[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://ilianpihlajamaa.github.io/CubatureRules.jl/dev/)
[![CI](https://github.com/IlianPihlajamaa/CubatureRules.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/IlianPihlajamaa/CubatureRules.jl/actions/workflows/CI.yml)

Quadrature and cubature rules computed on demand, for any polynomial degree and at any
precision, on intervals, simplices, boxes, spheres, balls and unbounded domains.

```julia
using CubatureRules

r = rule(Simplex{2}(); degree = 20, digits = 50)
integrate(x -> exp(x[1] * x[2]), r)
```

```
julia> r
QuadratureRule{2,BigFloat} on Simplex{2}()
  family    : XiaoGimbutas (seeded, Newton-refined)
  exactness : polynomial degree 20 (claimed; `check(rule)` verifies)
  points    : 79, all interior
  weights   : all positive, Σw = 0.5
  precision : 50 digits (guard 12, est. cond(J) = 1.6e+02)
  ...

julia> check(r)
Verification: PASSED
  basis     : orthonormal Dubiner, 336-bit arithmetic
  exactness : degree 20 ✓  (max residual 7.32e-51, tolerance 6.07e-48)
  sharpness : not exact at degree 21 ✓  (residual 3.56e-01)
  structure : Σw = measure ✓, interior ✓, positive ✓, symmetric ✓
```

## What it does differently

- **Rules are computed, not looked up.** The minimal symmetric rules on triangles,
  tetrahedra and spheres are refined by Newton's method to the precision you ask for, with
  guard digits chosen from the conditioning of the problem.
- **Each rule can be checked.** `check(r)` verifies exactness against an orthonormal basis
  at twice the rule's precision, and that the rule is *not* exact one degree higher.
- **Each rule has a record.** `provenance(r)` shows how it was made and why it was chosen;
  `cite(r)` gives the references; the licence of the underlying data travels with the rule.
- **You ask for a degree, not a scheme.** `rule(domain; degree)` ranks all available
  families by number of points; `available(domain; degree)` shows the ranking.
- **Unusual weights.** Gauss rules for a weight given only by its moments, including the
  ordinary moments that are too ill-conditioned to use in `Float64`.

## When to use something else

This package provides rules; it does not integrate adaptively to a tolerance.

| If you need | Use |
|---|---|
| Adaptive integration on an interval | [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) |
| Adaptive integration on a box | [HCubature.jl](https://github.com/JuliaMath/HCubature.jl) |
| Adaptive integration on a triangle or tetrahedron | [HAdaptiveIntegration.jl](https://github.com/zmoitier/HAdaptiveIntegration.jl) |
| Gauss rules with very many points in `Float64` | [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) |

Gauss–Kronrod rules come from QuadGK.jl, and Lobatto, Radau and Clenshaw–Curtis rules from
[QuadratureRules.jl](https://github.com/JuliaGNI/QuadratureRules.jl) when it is loaded.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/IlianPihlajamaa/CubatureRules.jl")
```

Julia 1.11 or later.

## Documentation

- [Tutorial](https://ilianpihlajamaa.github.io/CubatureRules.jl/dev/tutorial/first-rule/),
  starting from "I just want a rule"
- [Catalogue](https://ilianpihlajamaa.github.io/CubatureRules.jl/dev/catalogue/) of domains
  and rule families
- [Architecture & Design](https://ilianpihlajamaa.github.io/CubatureRules.jl/dev/design/pipeline/)
- [ROADMAP.md](ROADMAP.md) for planned work

## Licence

MIT. All tables shipped with the package were computed for it from orbit structures; no
published tables were copied. Rules refined from other packages' data, such as the GPL
tables of Lebedev.jl, keep that licence and are only used when you ask for them.
