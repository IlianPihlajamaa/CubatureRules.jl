# Catalogue

This section lists every domain and the rule families available on it. The tables are
generated from [`available`](@ref) when the documentation is built, so they show the degrees
the current version actually supports.

| Page | Domains |
|---|---|
| [Interval](interval.md) | `Interval`, Jacobi-weighted intervals, weights given by moments or by functions |
| [Simplex](simplex.md) | `Simplex{D}`: triangle, tetrahedron and higher; `SurfaceTriangle`, a triangle in space; the inverse-distance kernel of boundary elements |
| [Box](box.md) | `Orthotope{D}` |
| [Sphere](sphere.md) | `Sphere{D}`: circle, sphere and higher |
| [Ball](ball.md) | `Ball{D}`, `Disk` |
| [Wedge and pyramid](wedge_pyramid.md) | `Wedge`, `Pyramid` |
| [Unbounded](unbounded.md) | `LaguerreRay`, `HermiteLine`, `GaussianSpace`, `ExponentialSpace`, `HalfLine`, `RealLine` |

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
    "`Wedge()`" => Wedge(),
    "`Pyramid()`" => Pyramid(),
    "`LaguerreRay()`" => LaguerreRay(),
    "`HermiteLine()`" => HermiteLine(),
    "`GaussianSpace(3)`" => GaussianSpace(3),
    "`ExponentialSpace(3)`" => ExponentialSpace(3),
])
```

A *derived* family computes its rules from a formula and exists at every degree. A *seeded*
family refines stored starting values and exists at the degrees for which a seed is shipped;
see [Seed strategies](../design/seeds.md).

## Seed tables

The seeded families are the fully symmetric rules: positive weights, interior nodes, and the
whole symmetry group of the domain. Each comes from a table of seeds found by the package's
own searches. Their point counts are the smallest those searches reached, not proven minima.
This table is generated when the documentation is built:

```@setup seedtables
using CubatureRules, Markdown
const CR = CubatureRules
function seed_tables()
    tables = [("Triangle", "XiaoGimbutas", "S₃", CR.xg_entries()),
              ("Tetrahedron", "FullySymmetric", "S₄", CR.simplex_entries(3)),
              ("4-simplex", "FullySymmetric", "S₅", CR.simplex_entries(4)),
              ("Square", "FullySymmetric", "D₄", CR.box_entries(2)),
              ("Cube", "FullySymmetric", "Oₕ", CR.box_entries(3)),
              ("4-cube", "FullySymmetric", "B₄", CR.box_entries(4)),
              ("Disk", "FullySymmetric", "D₄", CR.round_entries(2)),
              ("Ball", "FullySymmetric", "Oₕ", CR.round_entries(3)),
              ("Pyramid", "FullySymmetric", "C₄ᵥ", CR.pyramid_entries()),
              ("Wedge", "FullySymmetric", "D₃ₕ", CR.wedge_entries()),
              ("Sphere", "LebedevRule", "Oₕ", CR.lebedev_entries())]
    rows = ["| Domain | Family | Symmetry | Degrees | Entries | Largest rule |", "|---|---|---|---|---|---|"]
    for (dom, fam, group, es) in tables
        isempty(es) && continue
        ds = [e.degree for e in es]
        big = es[argmax([e.npoints for e in es])]
        push!(rows, "| $dom | `$fam` | $group | $(minimum(ds))–$(maximum(ds)) | $(length(es)) | " *
                    "$(big.npoints) points (degree $(big.degree)) |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

```@example seedtables
seed_tables()
```

A request between two stored degrees is answered by the smallest stored rule of at least
that degree. The sphere's table holds odd degrees only: every `Oₕ` orbit is centrally
symmetric, so a rule of degree `2k` is automatically of degree `2k + 1`.

Not listed above:

- Families without a polynomial degree (`TanhSinh`, `ExpSinh`, `SinhSinh`), which are
  requested by level rather than degree.
- `Lobatto`, `Radau` and `ClenshawCurtis`, which appear once QuadratureRules.jl is loaded.
- `UpstreamLebedev`, which appears once Lebedev.jl is loaded.
- `ModifiedChebyshev`, which applies to intervals with a [`MomentWeight`](@ref);
  `StieltjesDiscretization`, which applies to domains with a [`FunctionWeight`](@ref);
  `Filon`, which applies to intervals with an [`Oscillatory`](@ref) weight;
  `SingularGauss`, which applies to intervals with a [`PrincipalValue`](@ref) or
  [`FinitePart`](@ref) kernel; `DuffyGauss` and `DuffySinh`, which apply to triangles with an
  [`InverseDistance`](@ref) kernel; `DuffyFinitePart`, which applies to triangles with an
  [`InverseDistanceCubed`](@ref) or [`InverseDistanceGradient`](@ref) kernel; `DuffyCone`,
  which applies to tetrahedra with an [`InverseDistance`](@ref) kernel; and `DuffyCurved`, which
  applies to the reference triangle with a [`CurvedInverseDistance`](@ref) kernel.
