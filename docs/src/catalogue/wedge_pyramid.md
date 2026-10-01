# Wedge and pyramid

```@setup wp
using CubatureRules, Markdown
function family_table(dom; degree)
    rows = ["| Family | Points at degree $degree | Derivation |", "|---|---|---|"]
    for row in available(dom; degree)
        push!(rows, "| `$(row.family)` | $(row.npoints) | $(row.derivation isa CubatureRules.Derived ? "derived" : "seeded") |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

The two remaining finite-element shapes in three dimensions. Both are given by their
vertices and must be affine images of their reference shapes: a rule keeps its polynomial
degree under an affine map and under no other, so a wedge with a twisted top or a pyramid
with a non-parallelogram base is refused rather than integrated with a claim that no longer
holds. [`transform`](@ref) integrates over those, with `NoClaim`.

## Domains

```@docs
Wedge
Pyramid
```

A mesh of either is integrated like a mesh of simplices, with one reference rule mapped to
each element. The unit cube as two wedges over the triangles of its bottom face:

```@example wp
W1 = Wedge((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (0, 1, 1))
W2 = Wedge((1, 0, 0), (1, 1, 0), (0, 1, 0), (1, 0, 1), (1, 1, 1), (0, 1, 1))
r = rule(Wedge(); degree = 3)
integrate(x -> x[1]^2 * x[2] + x[3]^3, r, [W1, W2])     # 1/6 + 1/4
```

## Wedge: product rules

```@docs
WedgeProduct
```

A polynomial of total degree `d` on the wedge is a sum of terms `p(x, y) zʳ` with
`deg p + r ≤ d`, so a triangle rule of degree `d` times a Gauss–Legendre rule of degree `d`
integrates it exactly. The selector takes the triangle family with the fewest points: the
shipped [Xiao–Gimbutas](simplex.md) rules through their table, and the conical product
beyond it.

```@example wp
family_table(Wedge(); degree = 10)
```

## Pyramid: conical product

The collapsed coordinates `x = ξ(1 − ζ)`, `y = η(1 − ζ)`, `z = ζ` carry the cube
`[-1, 1]² × [0, 1]` onto the pyramid with Jacobian `(1 − ζ)²`. A monomial `xᵃ yᵇ zᶜ` becomes
`ξᵃ ηᵇ (1 − ζ)ᵃ⁺ᵇ ζᶜ`, of degree at most `a + b + c` in `ζ` against the weight `(1 − ζ)²`. So
[`ConicalProduct`](@ref), with `m` Gauss–Legendre points in `ξ` and `η` and `m` Gauss–Jacobi
`(2, 0)` points in `ζ`, is exact to total degree `2m − 1` with `m³` positive weights and
interior nodes, at any degree and precision.

```@example wp
family_table(Pyramid(); degree = 10)
```

Both families are verified against polynomials of order one on the reference shape. For
the wedge these are the orthonormal Dubiner polynomials of the triangle times Legendre
polynomials in `z`, which are orthonormal on the wedge. For the pyramid they are Legendre
polynomials on its bounding box, integrated over the pyramid slice by slice with the
closed form of `∫ P_a` over each slice and a Gauss–Legendre rule in `z`, which shares
nothing with the collapsed construction but one-dimensional Gauss nodes.

Rules with fewer points, symmetric under the wedge's or the pyramid's own symmetry group,
are not shipped yet; they would come from the same orbit search as the triangle and
tetrahedron tables.
