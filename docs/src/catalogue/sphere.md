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

The rules consist of whole `O_h` orbits, so there is one equation per invariant polynomial
rather than one per spherical harmonic, which reduces their number considerably. The
refinement solves them in the orthonormal `O_h`-invariant spherical harmonics, which keep the
system well conditioned at every degree; the certificate quotes the residual in the
invariants `p₄ᵃp₆ᵇ` (`p₄ = Σxᵢ⁴`, `p₆ = x²y²z²`), whose exact integrals are rationals. See
[Orbit algebra](../design/symmetry.md#On-the-sphere). Only odd degrees
exist: every orbit is centrally symmetric, so a rule of degree `2k` is automatically of
degree `2k + 1`, and an even-degree request gets the odd rule above it.

The seeds were generated for this package from the orbit structures alone, and are MIT
licensed. The point counts agree with Lebedev's at every degree searched except degree 13,
where Lebedev's 74-point rule has a negative weight. This package requires positive weights
and gives a 78-point rule there; the 74-point structure was searched exhaustively, and every
solution found has a negative weight.

The table also holds every odd degree from 133 to 201, from 6074 points at degree 133 to
13680 at degree 201. They come from a continuation search that grows a verified rule into
the next degree by adding orbits, then removes orbits until the system is square — as many
unknowns as invariant conditions — and solves it with positive weights. They are not claimed
minimal; their counts are 0.2% to 1.5% above the parameter-count estimate `(d+1)²/3`.
Every odd degree from 19 to 131 is shipped as well: in-house rules grown from the degree-17
rule by the pooled search (`scripts/search_campaign/`), at most three unknowns above square.
Their counts equal Lebedev's at degrees 19, 21, 23, 35 and 41 and are at most about 10% above
his elsewhere; they are not claimed minimal. Every rule from degree 18 on needs at least one
48-point orbit in general position, since `[xyz(x²−y²)(y²−z²)(z²−x²)]²` is an invariant of
degree 18 that vanishes, with its gradient, on every other orbit kind. The selector
(`rule(Sphere{3}(); degree)` without a family) takes these rules over the product rule at
every degree, on node count.
`Float64` rules are served from the table. More digits are refined in the orthonormal
invariant harmonics, whose condition number stays moderate (about 3·10³ at degree 133 and
3·10⁵ at 189): 160 bits took about 7 s at degree 133 and 30 s at degree 189 (October 2026). In `Float64` these rules come too close
to the next degree for `check` to resolve sharpness — degree 133 misses degree 134 by
3e-41, degree 201 misses 202 by 5e-61 — so the miss, measured at high precision, is
recorded in the table and reported instead.

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

rule(Sphere{3}(); degree = 29, copyleft = true)       # Lebedev.jl's 302 points, not the in-house 314
rule(UpstreamLebedev(), Sphere{3}(); degree = 29, digits = 40)   # refined, still GPL-3
```

```@setup sphere_counts
# the comment above names both counts; the build fails here if a campaign changes ours
using CubatureRules
@assert npoints(rule(LebedevRule(), Sphere{3}(); degree = 29)) == 314
@assert Dict(CubatureRules.UPSTREAM_LEBEDEV_POINTS)[29] == 302
```

Where the two have the same number of points, as at degree 19 (146), the in-house rule wins
the tie and `copyleft = true` changes nothing.

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

Lebedev rules are verified on their orbits instead, which is equivalent for a symmetric rule
and much faster: at degree 59 (about 1200 points) the check took 0.3 s instead of 120 s
(September 2026). See
[Verification](../design/verification.md#Symmetric-rules-on-the-sphere).
