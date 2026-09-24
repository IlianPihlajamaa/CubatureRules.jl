# Catalogue

This section lists every domain and the rule families available on it. The tables are
generated from [`available`](@ref) when the documentation is built, so they show the degrees
the current version actually supports.

| Page | Domains |
|---|---|
| [Interval](interval.md) | `Interval`, Jacobi-weighted intervals, weights given by moments |
| [Simplex](simplex.md) | `Simplex{D}`: triangle, tetrahedron and higher |
| [Box](box.md) | `Orthotope{D}` |
| [Sphere](sphere.md) | `Sphere{D}`: circle, sphere and higher |
| [Ball](ball.md) | `Ball{D}`, `Disk` |
| [Unbounded](unbounded.md) | `LaguerreRay`, `HermiteLine`, `GaussianSpace`, `HalfLine`, `RealLine` |

## All families at a glance

```@setup overview
using CubatureRules, Markdown
function overview(doms)
    fmt(r) = last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"
    rows = ["| Domain | Family | Degrees | Derivation |", "|---|---|---|---|"]
    for (label, dom) in doms
        for row in available(dom)
            der = row.derivation isa CubatureRules.Derived ? "derived" : "seeded"
            push!(rows, "| $label | `$(row.family)` | $(fmt(row.degrees)) | $der |")
        end
    end
    Markdown.parse(join(rows, "\n"))
end
```

```@example overview
overview([
    "`Interval()`" => Interval(),
    "`Simplex{2}()`" => Simplex{2}(),
    "`Simplex{3}()`" => Simplex{3}(),
    "`Simplex{4}()`" => Simplex{4}(),
    "`Orthotope{3}()`" => Orthotope{3}(),
    "`Sphere{2}()`" => Sphere{2}(),
    "`Sphere{3}()`" => Sphere{3}(),
    "`Ball{3}()`" => Ball{3}(),
    "`LaguerreRay()`" => LaguerreRay(),
    "`HermiteLine()`" => HermiteLine(),
    "`GaussianSpace(3)`" => GaussianSpace(3),
])
```

A *derived* family computes its rules from a formula and exists at every degree. A *seeded*
family refines stored starting values and exists at the degrees for which a seed is shipped;
see [Seed strategies](../design/seeds.md).

Not listed above:

- Families without a polynomial degree (`TanhSinh`, `ExpSinh`, `SinhSinh`), which are
  requested by level rather than degree.
- `Lobatto`, `Radau` and `ClenshawCurtis`, which appear once QuadratureRules.jl is loaded.
- `UpstreamLebedev`, which appears once Lebedev.jl is loaded.
- `ModifiedChebyshev`, which applies to intervals with a [`MomentWeight`](@ref).

## Planned domains

`Polytope`, `Wedge` and `Pyramid` are defined so that requests for them give a clear error
naming the release in which they are planned.
