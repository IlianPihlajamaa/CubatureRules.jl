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

## Fully symmetric rules

[`FullySymmetric`](@ref) gives positive-weight rules with interior nodes that are invariant
under the shape's own symmetry group: on the wedge the 12 symmetries of the prism (`D₃ₕ`,
the permutations of the triangle's barycentric coordinates and `z ↦ −z`), on the pyramid the
8 symmetries of its square base (`C₄ᵥ`).

An orbit of the wedge is an orbit of the triangle, at `z = 0` or at the pair `±z`. The group
acts on `(x, y)` and on `z` separately, so its invariant polynomials are the triangle's
invariants times even polynomials in `z`; with orthonormal Legendre polynomials in `z` the
products are orthonormal on the wedge, and they are the moment equations.

An orbit of the pyramid is an orbit of the square in the collapsed coordinates
`(ξ, η) = (x, y)/(1 − z)`, at a free height `z`. The moment equations are in the invariant
polynomials `S_ij(ξ, η) (1 − z)^(i+j) q_k(z)`, with `S_ij` symmetrised products of even
Legendre polynomials and `q_k` orthonormal on `[0, 1]` for the weight `(1 − z)^(2(i+j)+2)`;
these are polynomials in `x, y, z` of degree `i + j + k`, and orthonormal on the pyramid.

The seeds were found in-house, by random starts on the orbit structures of each point count
at low degree and node elimination above (`scripts/generate_orbit_seeds.jl`). Point counts
against the product rules, and against Witherden & Vincent (2015), who tabulate the same
class to degree 10:

```@example wp
published = (wedge = [1, 5, 8, 11, 16, 28, 35, 46, 60, 85], pyramid = [1, 5, 6, 10, 15, 24, 31, 47, 62, 83])
rows = ["| Degree | Wedge: symmetric | Wedge: XG × Gauss | W&V | Pyramid: symmetric | Pyramid: conical | W&V |",
        "|---|---|---|---|---|---|---|"]
top(dom) = last(CubatureRules.degree_range(FullySymmetric(), dom))
for d in 1:max(top(Wedge()), top(Pyramid()))
    w = d <= top(Wedge()) ? npoints(rule(FullySymmetric(), Wedge(); degree = d)) : "—"
    p = d <= top(Pyramid()) ? npoints(rule(FullySymmetric(), Pyramid(); degree = d)) : "—"
    push!(rows, "| $d | $w | $(npoints(rule(WedgeProduct(), Wedge(); degree = d))) | " *
                "$(get(published.wedge, d, "—")) | $p | $(cld(d + 1, 2)^3) | $(get(published.pyramid, d, "—")) |")
end
Markdown.parse(join(rows, "\n"))
```

The counts are the smallest the searches reached, not proven minima.
