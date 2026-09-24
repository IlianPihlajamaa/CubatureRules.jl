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
