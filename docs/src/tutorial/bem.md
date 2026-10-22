# I want boundary-element integrals

A collocation or Nyström method on a triangulated surface needs, for every collocation point
`x₀` and every panel `T`, integrals like `∫_T f(y) G(x₀, y) dS(y)`, where `G` is singular at
`y = x₀`. This package builds rules that carry the kernel in their weights: for the right
kernel, a rule of degree `d` integrates `f · G` exactly for every polynomial `f` of degree `d`,
wherever `x₀` lies. This page says which kernel to use for which integral and where to go for
the details, which are in the [simplex catalogue](../catalogue/simplex.md).

## Which kernel

With `r = |y − x₀|`, `n` the normal of the panel, and `n₀` the normal at `x₀`:

| Integral | Kernel to use | Construction |
|---|---|---|
| single layer, `∫ f/r` | [`InverseDistance`](@ref)`(x₀)` | [`DuffyGauss`](@ref) on the panel, [`DuffySinh`](@ref) off it |
| adjoint double layer, `∫ f ∂/∂n₀ (1/r)` | [`InverseDistanceGradient`](@ref)`(x₀, n₀)` | [`DuffyFinitePart`](@ref) on the panel, [`DuffySinh`](@ref) off it |
| double layer, `∫ f ∂/∂n (1/r)` | minus [`InverseDistanceGradient`](@ref)`(x₀, n)` | the same |
| hypersingular, on a flat panel `∫ f/r³` | [`InverseDistanceCubed`](@ref)`(x₀)` | the same |

The gradient kernel is `(y − x₀)·e / r³`, and `∂/∂n (1/r) = −(y − x₀)·n / r³`, hence the sign.
On the panel that contains `x₀` both gradient integrals vanish, since `y − x₀` lies in the
panel's plane; they matter on the neighbouring panels, which meet it at an angle. With `x₀`
on the panel, `1/r³` and the gradient kernel are not integrable, and the rules compute the
finite part, which is additive over the panels around `x₀`; for the gradient kernel with `x₀`
inside a panel it is the principal value.

For the Helmholtz kernel `e^{ikr}/r`, use the `1/r` rule and put `e^{ikr}` in the integrand:
the rule is exact on the terms `p/r` and converges geometrically on the rest. For the
gradient and hypersingular Helmholtz kernels, pass `f(y) e^{ikr}(1 − ikr)` to the gradient or
`1/r³` rule. See [the catalogue](../catalogue/simplex.md#Strongly-singular-and-hypersingular-kernels)
for why that is enough.

## A panel and a collocation point

A panel is a [`SurfaceTriangle`](@ref), a triangle in space; the kernel is the weight of a
[`WeightedDomain`](@ref):

```@example bem
using CubatureRules
E = SurfaceTriangle((0.1, 0.2, 0.3), (1.3, -0.4, 0.7), (0.2, 1.1, -0.5))
c = (E.vertices[1] + E.vertices[2] + E.vertices[3]) / 3      # the centroid as collocation point
r = rule(WeightedDomain(E, InverseDistance(c)); degree = 9)
family(r), npoints(r), integrate(y -> 1 + y[1] * y[2], r)
```

The same call covers `x₀` anywhere: inside the panel, on an edge or at a vertex (where the
panel is cut into fewer pieces), and off the panel, near or far. A collocation point computed
in floating point, like this centroid, lies on the panel only to rounding; a point within 64
units in the last place of the panel's plane or boundary is moved onto it, and the provenance
says so. Exact coordinates are never moved.

```@example bem
using LinearAlgebra: cross, normalize
v = E.vertices
n = normalize(cross(v[2] - v[1], v[3] - v[1]))                  # the panel's unit normal
x₀ = Tuple(c + 1e-6 * n)                                       # 1e-6 off the panel
r = rule(WeightedDomain(E, InverseDistance(x₀)); degree = 9)
family(r), npoints(r)
```

## Many panels: the fast builder

The rules above are built in `BigFloat` and verified to the precision asked for, which takes
milliseconds to tenths of a second each. An assembly with hundreds of thousands of near and
singular pairs needs them faster. [`KernelRuleBuilder`](@ref) builds the same rules in
`Float64`, without allocating, in microseconds; it was validated against the certified rules
offline, and falls back to them outside the range it covers:

```@example bem
b = KernelRuleBuilder(InverseDistance, 9)
x, w = b(E, InverseDistance(c))                 # views into b's buffers, valid until the next call
integrate(y -> 1 + y[1] * y[2], b, E, InverseDistance(c))
```

Make one builder per kernel and degree, and one per thread. The
[catalogue](../catalogue/simplex.md#Fast-rules-for-boundary-element-assembly) has the accuracy
it was validated to.

## Curved panels

For isoparametric quadratic elements, give the six nodes of the element as a
[`QuadraticTriangle`](@ref) and use [`CurvedInverseDistance`](@ref) on the reference triangle,
with the collocation point either as a parameter point `ξ₀` of the element or as a point `x0`
in space, on or off the surface. The nodes of the rule are parameter points; map them with
the element, `Γ.(nodes(r))`, to evaluate a function on the surface. See
[curved triangles](../catalogue/simplex.md#Curved-triangles).

## Volume integrals

`∫_K f(y)/|y − x₀| dy` over a tetrahedron `K` with `x₀` inside uses the same
[`InverseDistance`](@ref) weight on a `Simplex` with three-dimensional vertices; the
construction is [`DuffyCone`](@ref).

## Galerkin methods

A Galerkin method needs double integrals over pairs of panels, singular when they share a
face, an edge or a vertex. [SauterSchwabQuadrature.jl](https://github.com/krcools/SauterSchwabQuadrature.jl)
regularises them and takes the one-dimensional rule as an argument, and a Gauss–Legendre rule
from this package takes it to any precision; see
[Galerkin pairs](../catalogue/simplex.md#Galerkin-pairs).

## Checking

Every certified rule here can be checked like any other: `check(r)` integrates the
orthonormal polynomials of the triangle against the kernel, computed independently in polar
coordinates about `x₀`, at twice the rule's precision.

```@example bem
v = check(rule(WeightedDomain(E, InverseDistance(c)); degree = 9))
v.exact, v.sharp
```
