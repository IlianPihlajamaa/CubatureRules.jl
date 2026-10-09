# Simplex

```@setup simplex
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
Simplex
```

On the triangle:

```@example simplex
family_table(Simplex{2}())
```

On the tetrahedron:

```@example simplex
family_table(Simplex{3}())
```

In four and more dimensions:

```@example simplex
family_table(Simplex{4}())
```

## Xiao–Gimbutas

```@docs
XiaoGimbutas
```

These are fully symmetric triangle rules with positive weights and interior nodes, at the
point counts of H. Xiao and Z. Gimbutas, *Comput. Math. Appl.* 59 (2010) 663–676,
doi:10.1016/j.camwa.2009.10.027, or fewer.

The seeds were generated for this package from the orbit structures alone; no published
numbers were used, and the table is MIT licensed. Above degree 20 they were found by growing
from the rule one degree lower and eliminating nodes (see [Seed strategies](../design/seeds.md)),
which found rules with fewer points than the published ones at most degrees from 27 on, for
example 139 instead of 141 at degree 27 and 412 instead of 423 at degree 48. Degrees 51 to
68 go beyond the published table. These are the smallest counts the searches reached, not
proven minima. Rules with the minimal number of points are not unique; the table contains
the one with the largest smallest barycentric coordinate among those found.

## Fully symmetric tetrahedron and 4-simplex rules

```@docs
FullySymmetric
```

Fully symmetric (`S₄`), positive-weight tetrahedron rules with interior nodes. Unlike the
triangle rules, these point counts are not taken from a paper: they are the smallest found by
the package's own search, starting from a number of points below which no fully symmetric
rule can exist, and they go up to degree 30. They are not proven to be minimal.

On the 4-simplex, the domain of space–time elements over tetrahedra, the same machinery
gives `S₅`-symmetric rules: orbits are patterns of equal barycentric coordinates, and the
moment equations are taken in the invariant subspace of an orthonormal basis on the
4-simplex. That basis is built the same way in every dimension (the Proriol–Koornwinder–
Dubiner construction, one homogenised Jacobi factor per coordinate), and the invariant
subspace is checked against the Molien series `Π_{j=2}^{5} 1/(1 − tʲ)`. The shipped table
covers low degrees only, as a check of the machinery:

```@example simplex
[(d, npoints(rule(FullySymmetric(), Simplex{4}(); degree = d)), npoints(rule(ConicalProduct(), Simplex{4}(); degree = d)))
 for d in 1:last(CubatureRules.degree_range(FullySymmetric(), Simplex{4}()))]
```

(degree, fully symmetric points, conical product points). The seeds are extended by
`scripts/generate_simplex4_seeds.jl`.

## Grundmann–Möller

```@docs
GrundmannMöller
```

A closed formula for rules of odd degree `2s + 1` on a simplex of any dimension `D`, with

```math
\sum_{i=0}^{s} \binom{s - i + D}{D}
```

points. All nodes and weights are rational, so these are the rules to use for exact
arithmetic. Their weights are negative from degree 3 on, and the number of points grows much
faster than for the other families.

Reference: A. Grundmann and H. M. Möller, *SIAM J. Numer. Anal.* 15 (1978) 282–290,
doi:10.1137/0715019.

## Conical product

```@docs
ConicalProduct
```

The simplex is mapped to a cube by a collapsed (Duffy-type) coordinate transformation, and a
Gauss–Jacobi rule is used in each direction, with Jacobi weights that absorb the Jacobian of
the transformation. With `m = ⌈(d+1)/2⌉` points per direction this gives `m^D` points,
available at every degree and precision.

Reference: A. H. Stroud, *Approximate Calculation of Multiple Integrals*, Prentice-Hall
(1971).

## The kernel 1/|y − x₀| on a triangle

```@docs
DuffyGauss
InverseDistance
```

The weakly singular kernel of boundary-element methods, `∫_T f(y) / |y − x₀| dy` with `x₀`
in the triangle, as the weight of a `WeightedDomain`. A rule for it carries the kernel in its
weights and is exact when `f` is a polynomial up to its degree. Here `x₀` is inside, so the
triangle is cut into three sub-triangles:

```@example simplex
T = Simplex((0, 0), (1, 0), (3 // 10, 8 // 10))
dom = WeightedDomain(T, InverseDistance((1 // 3, 1 // 4)))
r = rule(dom; degree = 9, digits = 30)
npoints(r), degree(r), integrate(y -> 1 + y[1] * y[2], r)
```

```@example simplex
v = check(r)
v.exact, v.sharp, v.positive, v.interior
```

The Duffy map from `x₀`, `y = x₀ + s (a + t (b − a))` on a sub-triangle with edges `a` and
`b`, has the Jacobian `s |a × b|` and gives `|y − x₀| = s √q(t)` with `q(t) = |a + t (b − a)|²`.
So `f(y) / |y − x₀| dy = f |a × b| / √q(t) ds dt`: the singularity cancels, and what remains is
a polynomial in `s` and `t` times the weight `1 / √q(t)`. The rule is Gauss–Legendre in `s`
times the Gauss rule of that weight in `t`. The weight depends on the shape of the
sub-triangle, so every rule is built for its triangle and its point, at any precision. Its
Gauss rule comes from the Stieltjes procedure on a discretisation of the weight. `q` has
complex roots `t* ± iη`, and the substitution `t = t* + η sinh u` turns the weight into the
constant `1/|b − a|`, so Gauss–Legendre in `u` discretises it with geometric convergence even
when the sub-triangle is thin.

An ordinary triangle rule converges only slowly on such an integrand, because the kernel is
not a polynomial: on the triangle above with `x₀` at a vertex, the 291-point Xiao–Gimbutas
rule of degree 40 integrates `1/|y − x₀|` itself with a relative error of 9e-4, where a rule
of this family integrates `f(y)/|y − x₀|` exactly for every polynomial `f` of its degree —
36 points for degree 11 (`notes/singular-bem.md`).

The Helmholtz kernel `e^{ikr}/r`, `r = |y − x₀|`, needs no other rule: integrate `f(y) e^{ikr}`
with the weight `1/r`. Its Taylor series in `r` has, besides the terms `p/r` the rule is exact
on, polynomial terms, which the Gauss rule of `1/√q` integrates with geometric convergence
since `√q` is analytic. With `x₀` at the right angle of the unit triangle the reference is
`∫₀^{π/2} sin(k R(θ))/k dθ`, `R(θ) = 1/(cos θ + sin θ)`, in polar coordinates:

```@example simplex
k = 5
dom = WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistance((0, 0)))
θr = rule(Interval(0, π / 2); degree = 99)
ref = sum(w * sin(k / (cos(θ[1]) + sin(θ[1]))) / k for (θ, w) in zip(nodes(θr), weights(θr)))
[(d, npoints(rule(dom; degree = d)), abs(integrate(y -> cos(k * hypot(y...)), rule(dom; degree = d)) - ref))
 for d in (5, 9, 13, 17, 21)]
```

A rule exact on the polynomial terms as well was built and measured, and did worse at equal
point count on every geometry tried (`notes/singular-bem.md`).

`check` verifies these rules against the orthonormal Dubiner polynomials of the triangle,
whose integrals with the kernel it computes independently, in polar coordinates about `x₀`.

The same works on a triangle in space, a [`SurfaceTriangle`](@ref), as a boundary element
of a surface mesh: the construction uses only distances and the area element, so it is the
same in three dimensions.

```@example simplex
E = SurfaceTriangle((0.1, 0.2, 0.3), (1.3, -0.4, 0.7), (0.2, 1.1, -0.5))
c = (E.vertices[1] + E.vertices[2] + E.vertices[3]) / 3        # the collocation point
r = rule(WeightedDomain(E, InverseDistance(c)); degree = 9)
npoints(r), check(r).exact
```

A point computed in floating point, like this centroid, is on the triangle only to rounding:
here it is 1.7e-17 off the plane. `x₀` within 64 units in the last place of the triangle's
plane or boundary is moved onto it and the rule is built for that point, as the provenance
says; exact coordinates are never moved. A point further off is near-singular, the next
section. The stronger singularities of double-layer and hypersingular kernels are not
covered yet.

Reference: M. G. Duffy, *SIAM J. Numer. Anal.* 19 (1982) 1260–1262, doi:10.1137/0719090.

## Near-singular integrals

```@docs
DuffySinh
```

With `x₀` off the triangle the integrand is smooth, but an ordinary rule needs ever more
points as `x₀` comes closer: the kernel then varies on the scale of the distance. Here
`x₀` is 1e-6 above a triangle in space:

```@example simplex
E = SurfaceTriangle((0, 0, 0), (1, 0, 0), (3 // 10, 8 // 10, 0))
r = rule(WeightedDomain(E, InverseDistance((1 // 3, 1 // 4, 1e-6))); degree = 9)
family(r), npoints(r), integrate(y -> 1 + y[1] * y[2], r)
```

The triangle is cut at its point `c` nearest to `x₀`, and mapped by the Duffy map from `c`.
Along each radial line the kernel is a weight on `[0, 1]`, `s / √(q s² + 2ℓ s + H²)` with `H`
the distance from `x₀` to `c`; its moments follow a three-term recurrence from a closed form,
and the line gets the Gauss rule of that weight, built from them at any precision. So the
radial direction is exact for polynomials however close `x₀` is, and the rule does not grow
as `x₀` approaches the triangle. Across the lines, what is left is analytic but not a
polynomial, and is integrated with Gauss–Legendre after a sinh substitution (P. R. Johnston
and D. Elliott, *Int. J. Numer. Meth. Engng* 62 (2005) 564–578, doi:10.1002/nme.1208), with
more points until the rule's integrals of every polynomial of its degree agree between two
resolutions to the precision asked for. The rule is therefore exact to working precision,
like every other rule here, but it grows with that precision, and `npoints(DuffySinh(), …)`
gives only its smallest size. `check` verifies it against integrals computed in polar
coordinates about `c`.

## Triangles in space

```@docs
SurfaceTriangle
```

A flat triangle in three dimensions takes every rule of the reference triangle, mapped onto
it with the area element, so the degree claim carries over to polynomials in the three
coordinates:

```@example simplex
E = SurfaceTriangle((0, 0, 0), (1, 0, 1), (0, 1, 1))
r = rule(E; degree = 6)
family(r), npoints(r), integrate(y -> y[3]^2, r), sqrt(3) / 4
```

To integrate over many triangles, build one rule on `Simplex{2}()` and pass the triangles to
`integrate`, which maps the rule onto each on the fly without building anything:

```@example simplex
r = rule(Simplex{2}(); degree = 6)
mesh = [SurfaceTriangle((0, 0, 0), (1, 0, 1), (0, 1, 1)), SurfaceTriangle((1, 0, 1), (1, 1, 2), (0, 1, 1))]
integrate(y -> y[3]^2, r, mesh)
```

## Strongly singular and hypersingular kernels

```@docs
DuffyFinitePart
InverseDistanceCubed
InverseDistanceGradient
```

Collocation at a point `x₀` of a panel meets, besides `1/r`, the kernels `1/r³` (the
hypersingular operator on a flat panel) and `(y − x₀)·e/r³` (the adjoint double layer at an
edge or vertex shared with a panel at an angle), `r = |y − x₀|`. Neither is integrable at `x₀`;
the rules compute the finite part with respect to `r`, the integral outside the disk `r < ε`
less its terms in `1/ε` and `ln ε`. That is additive over the panels around `x₀`, and for the
gradient kernel with `x₀` inside a panel it is the principal value. Here the principal value
of `(y − x₀)·e/r³` over a triangle with `x₀` inside, against the divergence theorem,
`−∮ (n·e)/r ds`, which is closed form along each edge:

```@example simplex
T = Simplex((0, 0), (1, 0), (3 // 10, 8 // 10))
x0, e = (2 // 5, 3 // 10), (3 // 5, 4 // 5)
r = rule(WeightedDomain(T, InverseDistanceGradient(x0, e)); degree = 6, digits = 40)
V = [big.(collect(v)) for v in CubatureRules.vertices(T)]
X0 = big.(collect(x0))
pv = -sum(zip(V, circshift(V, -1))) do (p, q)
    u = (q - p) / sqrt(sum(abs2, q - p))
    f = p + sum((X0 - p) .* u) * u
    h = sqrt(sum(abs2, X0 - f))
    (u[2] * e[1] - u[1] * e[2]) * (asinh(sum((q - f) .* u) / h) - asinh(sum((p - f) .* u) / h))
end
npoints(r), Float64(abs(integrate(y -> 1, r) - pv))
```

With the Duffy map `y = x₀ + s v(t)` on each piece of the triangle cut at `x₀`, `r = s √q(t)`,
so excluding `r < ε` excludes `s < ε/√q(t)`: the finite part in `r` is the finite part in `s`
plus terms in `ln √q(t)` times `φ(x₀)` or `∇φ(x₀)·v(t)` (Guiggiani's correction). The rule is
interpolatory in `s` on Gauss–Legendre nodes, weighted for the finite part; the logarithmic
terms are folded into its weights through the values along each ray, so no node sits at `x₀`;
the angular direction has the Gauss rule of `q(t)^{-3/2}`. The weights are signed and grow with
the degree — `Σ|w|/|Σw|` is about 2000 for `1/r³` at degree 9 — so a `Float64` rule loses up to
three digits to cancellation; ask for more digits where that matters.

For Helmholtz, both kernels arise multiplied by `e^{ikr}(1 − ikr) = 1 + k²r²/2 + ⋯` on a flat
panel: pass `f(y) e^{ikr}(1 − ikr)` as the integrand. The linear term in `r` cancels, which
is what keeps the finite part well defined (`notes/singular-bem.md`).

## The kernel 1/|y − x₀| on a tetrahedron

```@docs
DuffyCone
```

Volume integral equations need `∫_K f(y)/|y − x₀| dy` over a tetrahedron with `x₀` in it,
the same [`InverseDistance`](@ref) weight on a `Simplex` with three-dimensional vertices. The
tetrahedron is cut into the cones from `x₀` over its faces; along each cone the integrand is a
polynomial in the radial variable, and across it the near-singular kernel `1/|z − x₀|` on the
face, with `x₀` off the face's plane: the case of [`DuffySinh`](@ref). So the rules are exact
to working precision for every polynomial of their degree, wherever `x₀` lies in the
tetrahedron, and their size is set by that precision:

```@example simplex
K = Simplex((0, 0, 0), (1, 0, 0), (3 // 10, 9 // 10, 0), (1 // 5, 1 // 4, 4 // 5))
dom = WeightedDomain(K, InverseDistance((3 // 10, 3 // 10, 1 // 5)))
r = rule(dom; degree = 3)
v = check(r)
npoints(r), v.exact, v.sharp, npoints(rule(dom; degree = 3, T = Float32))
```

The points go mostly to the angular direction on the faces, where the integrand after the
radial integration is analytic but not polynomial: about twenty lines per piece of a face for
`Float64`, wherever `x₀` stands.

## Curved triangles

```@docs
QuadraticTriangle
CurvedInverseDistance
DuffyCurved
```

On a curved boundary element `Γ = χ(T̂)`, the six-node quadratic triangle of isoparametric
methods, the weakly singular integral is taken where the shape functions are polynomials, on
the reference triangle: `∫_Γ f/|y − x₀| dS = ∫_T̂ f(χ(ξ)) J(ξ)/|χ(ξ) − χ(ξ₀)| dξ`. The weight
`J/|χ − χ(ξ₀)|` is singular at `ξ₀` but not of the flat form, so `DuffyCurved` gives each ray
from `ξ₀` the Gauss rule of its own weight — smooth after the Duffy map, since
`χ(ξ₀ + s v) − χ(ξ₀) = s (Dχ(ξ₀) v + s ½ vᵀD²χ v)` exactly for a quadratic map — and refines the
angular direction until every polynomial of the degree is integrated to the precision asked
for. An octant of the unit sphere as one element, with the shape function of its first vertex:

```@example simplex
s = 1 / sqrt(2)
Γ = QuadraticTriangle((1, 0, 0), (0, 1, 0), (0, 0, 1), (s, s, 0), (0, s, s), (s, 0, s))
dom = WeightedDomain(Simplex{2}(), CurvedInverseDistance(Γ, (1 // 5, 3 // 10)))
r = rule(dom; degree = 5)
N₁(ξ) = (1 - ξ[1] - ξ[2]) * (1 - 2ξ[1] - 2ξ[2])
v = check(r)
npoints(r), integrate(N₁, r), v.exact, v.sharp
```

The nodes are in `ξ`; `Γ.(nodes(r))` are the points on the surface, for kernels such as
Helmholtz's `e^{ik|y − x₀|}` factor. A straight-sided element is a flat triangle, and there the
rules agree with [`DuffyGauss`](@ref) to rounding.

## Galerkin pairs

A Galerkin discretisation needs `∫_{T₁} ∫_{T₂} φ(x) ψ(y) G(x, y) dy dx` over pairs of
triangles, singular when they share a face, an edge or a vertex. The Sauter–Schwab
transformations regularise each configuration into an integral over `[0, 1]⁴` of an analytic
function, and [SauterSchwabQuadrature.jl](https://github.com/krcools/SauterSchwabQuadrature.jl)
implements them with the one-dimensional rule as an argument. It computes in the arithmetic
of that rule, so a Gauss–Legendre rule from this package takes it to any precision; this
package does not repeat it. The common face, against the closed form for a triangle with
sides `a`, `b`, `c` and area `A`, `(4A²/3) Σ (1/a) ln(((a + b)² − c²)/(b² − (c − a)²))`:

```@example simplex
using SauterSchwabQuadrature
setprecision(BigFloat, 256) do
    g = rule(GaussLegendre(), Interval(0, 1); degree = 2 * 16 - 1, digits = 60)
    qps = collect(zip(first.(nodes(g)), weights(g)))         # 16 points on [0, 1]
    v = (BigFloat[2, 0], BigFloat[3 // 10, 9 // 10], BigFloat[0, 0])
    len(p) = sqrt(sum(abs2, p))
    J = abs((v[1] - v[3])[1] * (v[2] - v[3])[2] - (v[1] - v[3])[2] * (v[2] - v[3])[1])
    x(u) = v[3] + u[1] * (v[1] - v[3]) + u[2] * (v[2] - v[3])  # the parametrisation it expects
    I = sauterschwab_parameterized((u, w) -> J^2 / len(x(u) - x(w)), CommonFace(qps))
    a, b, c = len(v[2] - v[3]), len(v[3] - v[1]), len(v[1] - v[2])
    t(a, b, c) = log(((a + b)^2 - c^2) / (b^2 - (c - a)^2)) / a
    Float64(I), Float64(abs(I - 4(J / 2)^2 / 3 * (t(a, b, c) + t(b, c, a) + t(c, a, b))))
end
```

The error falls geometrically with the number of points: 4e-5 at 5 points per axis, 1e-9 at
10, 1e-17 at 20 and 6e-33 at 40, where the work is `6n⁴` evaluations (90 s in `BigFloat` at
40). The common edge and the common vertex converge as fast, with the vertices ordered as the
package expects: a common vertex first in both triangles, a common edge at positions 1 and 3
of both, in the same order (other orderings converge only algebraically).
