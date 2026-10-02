# Box

```@setup box
using CubatureRules, Markdown
function family_table(dom)
    fmt(r) = last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"
    rows = ["| Family | Degrees | Derivation |", "|---|---|---|"]
    for row in available(dom)
        der = row.derivation isa CubatureRules.Derived ? "derived" : "seeded"
        push!(rows, "| `$(row.family)` | $(fmt(row.degrees)) | $der |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

## Domain

```@docs
Orthotope
```

On the cube:

```@example box
family_table(Orthotope{3}())
```

## Fully symmetric rules

On the square and the cube, [`FullySymmetric`](@ref) gives positive-weight rules with
interior nodes that are invariant under all 8 (square) or 48 (cube) symmetries of the
domain, the signed permutations of the coordinates. At the same total degree they need far
fewer points than the tensor Gauss rule:

```@example box
using Markdown
rows = ["| Degree | Square: symmetric | Square: tensor | Cube: symmetric | Cube: tensor |", "|---|---|---|---|---|"]
for d in 3:2:11
    sq = npoints(rule(Orthotope{2}(); degree = d))
    cu = d <= last(CubatureRules.degree_range(FullySymmetric(), Orthotope{3}())) ? npoints(rule(Orthotope{3}(); degree = d)) : "—"
    m = cld(d + 1, 2)
    push!(rows, "| $d | $sq | $(m^2) | $cu | $(m^3) |")
end
Markdown.parse(join(rows, "\n"))
```

The tensor rule integrates the larger space of tensor-product polynomials, `xᵃyᵇzᶜ` with
each exponent up to the degree. That is what the mass and stiffness matrices of
tensor-product (`Q_k`) elements need, and there the tensor rule is the one to use: a total
degree high enough to cover those integrands would need many more points. The symmetric
rules are for total-degree integrands — `P_k` or serendipity elements on quadrilaterals and
hexahedra, discontinuous Galerkin with total-degree bases, nonlinear terms.

The seeds were found in-house: random starts on every orbit structure of each point count at
low degree, then each degree grown from the one below and thinned by node elimination. The
point counts are compared with Witherden & Vincent (2015), who tabulate the same class to
degree 21 on the square and 11 on the cube. On the cube at degree 3 the table has 8 points
where they give 6: a 6-point rule of degree 3 has to put its nodes on the faces, at
`(±1, 0, 0)`, and these rules keep every node inside. Only odd degrees are tabulated: every
orbit is centrally symmetric, so a rule of degree `2k` is one of degree `2k + 1`.

## Tensor product

```@docs
TensorProduct
```

A tensor product of one-dimensional rules, one per axis, built from any one-dimensional
family. `TensorProduct(GaussLegendre(), 3)` uses Gauss–Legendre on all three axes, and
`TensorProduct((f₁, f₂))` uses different families per axis. Two existing rules can also be
combined with `r₁ ⊗ r₂` (`using CubatureRules: ⊗`).

The claim is the total degree, even though the rule is exact on the larger tensor-product
space: a 4 × 4 Gauss rule integrates `x⁶y⁶` exactly but claims degree 7. The selector
compares rules by total degree, so that is what is reported.

Exact output (`T = Rational{BigInt}`) is available when all factors support it, for example
with `NewtonCotes`.
