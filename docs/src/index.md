# CubatureRules.jl

CubatureRules.jl constructs quadrature and cubature rules on demand, for any polynomial
degree and at any precision.

```@repl
using CubatureRules
r = rule(Simplex{2}(); degree = 20)
integrate(x -> exp(x[1] * x[2]), r)
```

Most quadrature libraries store tables of rules. This package computes the rule when you
ask for it, to the precision you ask for, and keeps a record of how it was computed.

- **Any precision.** A rule with 200 digits is requested the same way as one with 16. The
  nodes and weights are computed at that precision, with extra working digits chosen from
  the conditioning of the problem.
- **A record of each rule.** Every rule stores which family built it, what it started from,
  how it was refined, under which licence its data falls, and why it was chosen. `cite(r)`
  returns the references as BibTeX.
- **Verification.** `check(r)` tests a rule against an orthonormal basis at twice its
  precision: that it is exact to its degree, not exact to the next degree, and that its
  weights, nodes and symmetry are as stated.

## Contents

- **[Tutorial](tutorial/first-rule.md)**: short task-based guides, starting with
  [I just want a rule](tutorial/first-rule.md).
- **[Architecture & Design](design/pipeline.md)**: how the package works and why.
- **[Catalogue](catalogue/index.md)**: all domains and the rule families available on each.
- **[API](api.md)**: reference for all exported and public functions.

## Scope

The package provides rules. It does not integrate adaptively to a tolerance: you choose the
rule, and can estimate its error by comparing two degrees or with an embedded pair (see
[I want an error estimate](tutorial/errors.md)). For adaptive integration, use
[QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) on intervals,
[HCubature.jl](https://github.com/JuliaMath/HCubature.jl) on boxes, or
[HAdaptiveIntegration.jl](https://github.com/zmoitier/HAdaptiveIntegration.jl) on
simplices.

Rules that other packages already compute in generic arithmetic are taken from them rather
than reimplemented: Gauss–Kronrod from QuadGK.jl, and Lobatto, Radau and Clenshaw–Curtis from
QuadratureRules.jl when it is loaded. See [External providers](design/providers.md).
