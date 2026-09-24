# Unbounded domains

```@setup unb
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

An integral over an unbounded domain usually only converges with a weight, so most of these
domains carry one. The weight is part of the domain, and a rule on it integrates
`∫ f(x) w(x) dx`.

## Domains

```@docs
HalfLine
RealLine
RealSpace
LaguerreRay
HermiteLine
GaussianSpace
```

| Domain | Weight | `measure` |
|---|---|---|
| `LaguerreRay(α)` | `x^α e^{-x}` on `[0, ∞)` | `Γ(α + 1)` |
| `HermiteLine()` | `e^{-x²}` on `ℝ` | `√π` |
| `GaussianSpace(D)` | `e^{-|x|²}` on `ℝ^D` | `π^{D/2}` |
| `HalfLine()`, `RealLine()` | none | `Inf` |

## Gauss–Laguerre and Gauss–Hermite

```@docs
GaussLaguerre
GaussHermite
```

```@example unb
family_table(LaguerreRay())
```

Both use the same driver as [`GaussJacobi`](@ref): Golub–Welsch in `Float64` for the seed,
Newton on the three-term recurrence at working precision, and weights from the Christoffel
function. Laguerre nodes spread out to about `4n`, so the guard grows with `n`.

These domains are verified against the weight's own orthonormal polynomials rather than
monomials, because the moments `∫ xᵏ w` grow factorially and would make a monomial check
badly conditioned.

## Gaussian-weighted space

```@docs
GaussianProduct
```

```@example unb
family_table(GaussianSpace(3))
```

In polar coordinates the Gaussian weight separates into a radial and an angular part. With
`u = r²` the radial factor becomes a generalised Gauss–Laguerre problem with `α = (D - 2)/2`,
and odd degrees vanish through the angular factor, so only `⌈(d + 2)/4⌉` radial points are
needed.

## Double-exponential rules

```@docs
ExpSinh
SinhSinh
```

The double-exponential idea of [`TanhSinh`](@ref) on the half line and the real line, for
integrands without a weight. Both have no polynomial degree ([`NoClaim`](@ref)) and are
requested by level.

```@example unb
r = rule(ExpSinh(5), HalfLine());
npoints(r), integrate(x -> exp(-x) / sqrt(x), r) - sqrt(π)
```

On a finite interval a double-exponential sum can be cut off where the weights fall below
the working precision, because the rest is provably negligible. Here the weights grow instead,
so no such bound exists without knowing how fast the integrand decays. These rules therefore
use a fixed range: the nodes span `[2^{-2(p-2)}, 2^{2(p-2)}]` at output precision `p`, about
`10^{±31}` in `Float64`.

The range is squared because the truncation removes a head and a tail. An integrand behaving
like `x^{-1+δ}` near the origin leaves a remainder of about `x_min^δ`, and one decaying like
`x^{-1-δ}` leaves about `x_max^{-δ}`. For `δ = 1/2`, the unsquared range would limit the
accuracy of `∫₀^∞ x^{-1/2} e^{-x} dx` to `5e-8`; the squared range reaches rounding error
with 18% more nodes.

These rules cannot handle an integrand that does not decay, or one whose mass lies outside
the range. A convergence study, not the certificate, shows whether the rule has resolved a
particular integrand.
