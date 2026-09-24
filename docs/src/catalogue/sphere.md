# Sphere

```@setup sphere
using CubatureRules, Markdown
function family_table(dom)
    fmt(r) = last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"
    rows = ["| Family | Degrees | Derivation | Symmetry |", "|---|---|---|---|"]
    for row in available(dom)
        der = row.derivation isa CubatureRules.Derived ? "derived" : "seeded"
        push!(rows, "| `$(row.family)` | $(fmt(row.degrees)) | $der | $(row.symmetry) |")
    end
    Markdown.parse(join(rows, "\n"))
end
```

## Domain

```@docs
Sphere
```

`Sphere{D}` is the unit sphere in `ℝ^D`, a surface of dimension `D - 1`: `Sphere{2}()` is
the circle and `Sphere{3}()` the ordinary sphere.

On the sphere in `ℝ³`:

```@example sphere
family_table(Sphere{3}())
```

## Lebedev

```@docs
LebedevRule
```

Rules invariant under the 48 signed permutations of the coordinates (the octahedral group
`O_h`), with positive weights. They need about half as many points as the product rules
below, so the selector prefers them where a seed is available.

The rules consist of whole `O_h` orbits, and the equations are written in terms of the
invariants `p₄ = Σxᵢ⁴` and `p₆ = x²y²z²` rather than in spherical harmonics, which reduces
their number considerably; see [Orbit algebra](../design/symmetry.md). Only odd degrees
exist: every orbit is centrally symmetric, so a rule of degree `2k` is automatically of
degree `2k + 1`, and an even-degree request gets the odd rule above it.

The seeds were generated for this package from the orbit structures alone, and are MIT
licensed. The point counts agree with Lebedev's at every degree searched except degree 13,
where Lebedev's 74-point rule has a negative weight. This package requires positive weights
and gives a 78-point rule there; the 74-point structure was searched exhaustively, and every
solution found has a negative weight.

Reference: V. I. Lebedev and D. N. Laikov, A quadrature formula for the sphere of the 131st
algebraic order of accuracy, *Doklady Mathematics* 59 (1999) 477–481.

### Rules from Lebedev.jl

`UpstreamLebedev()` is `LebedevRule{LebedevJLSeeds}()`. When
[Lebedev.jl](https://github.com/stefabat/Lebedev.jl) is loaded, it takes its seeds from that
package's tables, which reach degree 125, and refines them to the requested precision. The
orbit structure is recovered from the raw points by `CubatureRules.classify_octahedral`.

Lebedev.jl is GPL-3, and so are rules refined from its tables. They are listed by
[`available`](@ref) but not chosen automatically; request them by name or with
`copyleft = true`:

```julia
using CubatureRules, Lebedev

rule(Sphere{3}(); degree = 19, copyleft = true)                  # the 146-point rule
rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 40)   # refined, still GPL-3
```

See [Provenance and licensing](../design/provenance.md).

## Product rules

```@docs
SphereProduct
```

Available at every degree, dimension and precision, with positive weights and without a
table. On the sphere in `ℝ³`, with `t = cos θ`, the surface integral separates as

```math
\int_{S^2} f \, d\sigma = \int_{-1}^{1} \int_0^{2\pi} f \, d\varphi \, dt ,
```

and the rule is Gauss–Legendre in `t` combined with the trapezoid rule in `φ`. The
trapezoid rule with `m` equally spaced points integrates `e^{ikφ}` exactly for `|k| < m`, so
`m = d + 1` and `n = ⌈(d+1)/2⌉` give degree `d` with `⌈(d+1)/2⌉ (d+1)` points, about twice
as many as a Lebedev rule.

In higher dimensions the surface element in hyperspherical coordinates contains a power of
`sin θ` for each polar angle. The substitution `t = cos θ` turns each of these into a Jacobi
weight with `α = β`, so every polar angle gets a Gauss–Jacobi rule and the azimuth keeps the
trapezoid rule.

On the circle only the trapezoid rule is left, and it is optimal there: `d + 1` points for
degree `d`.

## Verification

Rules on `Sphere{3}` are verified against real spherical harmonics. Monomials are not a
basis on a sphere, because `Σxᵢ² = 1` makes them linearly dependent; the harmonics of degree
`≤ d` span exactly the polynomials of degree `≤ d` restricted to the sphere. In higher
dimensions the monomials are used as a spanning set.
