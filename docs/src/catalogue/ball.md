# Ball

```@setup ball
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

## Domains

```@docs
Ball
Disk
```

On the ball in `ℝ³`:

```@example ball
family_table(Ball{3}())
```

## Ball product

```@docs
BallProduct
```

A ball separates into a radius and a sphere:

```math
\int_{B^D} f \, dx = \int_0^1 r^{D-1} \int_{S^{D-1}} f(r\omega) \, d\sigma(\omega) \, dr .
```

A polynomial of degree `d` restricted to a ray is `r^d` times a polynomial on the sphere, so
the two factors can be treated separately: a Gauss–Jacobi rule for the weight `r^{D-1}` in
the radius, and any sphere family for the angles. The selector tries every sphere family as
the angular factor, so `Ball{3}()` gets a Lebedev angular rule where a Lebedev seed exists
and a product rule elsewhere.

```@example ball
npoints(rule(Ball{3}(); degree = 9)), npoints(rule(BallProduct(SphereProduct()), Ball{3}(); degree = 9))
```

Rules on a ball are verified against monomials with exact moments. Unlike on the sphere, the
coordinates in a ball satisfy no relation, so the monomials are a basis. Their moments are
those of the sphere divided by `|α| + D`.
