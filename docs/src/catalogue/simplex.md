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

## Fully symmetric tetrahedron rules

```@docs
FullySymmetric
```

Fully symmetric (`S₄`), positive-weight tetrahedron rules with interior nodes. Unlike the
triangle rules, these point counts are not taken from a paper: they are the smallest found by
the package's own search, starting from a number of points below which no fully symmetric
rule can exist, and they go up to degree 30. They are not proven to be minimal.

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
