# I just want a rule

Call [`rule`](@ref) with a domain and the polynomial degree you need:

```@repl first
using CubatureRules
r = rule(Simplex{2}(); degree = 20)
```

This is a rule on the reference triangle that integrates every polynomial of degree 20 or
less exactly. Apply it to a function with [`integrate`](@ref):

```@repl first
integrate(x -> exp(x[1] * x[2]), r)
```

In more than one dimension the integrand receives the point as an `SVector` of coordinates in
the rule's number type, so `x[1]` and `x[2]` are its coordinates. In one dimension it receives a
number.

## Build the rule once

Constructing a rule takes some work, and the package does not cache rules. Create the rule
once, outside any loop, and reuse it:

```julia
r = rule(Simplex{2}(); degree = 20)
total = integrate(f, r, mesh)              # mesh: a vector of triangles
```

`integrate(f, r, cell)` gives a single cell's integral, for a loop of your own.

For the same reason `integrate` always takes a rule, never a domain. To judge how accurate
the result is, see [I want an error estimate](errors.md).

## Domains

The most common domains, in their reference position:

```@repl first
rule(Interval(); degree = 9)          # the interval [-1, 1]
rule(Orthotope{3}(); degree = 9)      # the cube [-1, 1]³
rule(Sphere{3}(); degree = 9)         # the unit sphere in ℝ³
rule(Ball{3}(); degree = 9)           # the unit ball in ℝ³
```

`Simplex{D}()` works in any dimension, and `Disk()` is the unit disk.

Every bounded domain can also be given its own position and size. `rule` builds the rule on
the reference domain and maps it onto yours. The map is affine, so the rule keeps its
degree:

```@repl first
r = rule(Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 3.0)); degree = 9);
sum(weights(r))                                   # the triangle's area
r = rule(Orthotope((0, 0, 0), (1, 2, 3)); degree = 9);
sum(weights(r))                                   # the volume of [0, 1] × [0, 2] × [0, 3]
r = rule(Sphere((1.0, 0.0, 0.0), 2.0); degree = 9);
sum(weights(r))                                   # 4π r², with r = 2
```

| Domain | Reference | Your own |
|---|---|---|
| interval | `Interval()`, the interval `[-1, 1]` | `Interval(a, b)` |
| simplex | `Simplex{D}()`, vertices `0` and the unit vectors | `Simplex(v₀, v₁, …, v_D)` |
| box | `Orthotope{D}()`, the box `[-1, 1]^D` | `Orthotope(lo, hi)` |
| sphere | `Sphere{D}()`, the unit sphere | `Sphere(centre, radius)` |
| ball, disk | `Ball{D}()` and `Disk()`, the unit ball and disk | `Ball(centre, radius)` |
| wedge | `Wedge()`, the reference triangle times `[-1, 1]` | `Wedge(v₁, …, v₆)`, the top a translate of the bottom |
| pyramid | `Pyramid()`, base `[-1, 1]²`, apex `(0, 0, 1)` | `Pyramid(v₁, …, v₅)`, with a parallelogram base |
| triangle in space | none | `SurfaceTriangle(v₀, v₁, v₂)` |

Each of these is an affine image of its reference domain, and a shape that is not is
refused: a wedge with a twisted top, or a pyramid whose base is not a parallelogram. Such a
shape, or a curved element, needs a non-affine map, under which a rule no longer has a
degree; see [I want my own domain or mesh](domains.md), which also covers meshes.

On an unbounded domain an integral usually only converges with a weight, so these domains
come with one. The weight is part of the domain, and `integrate` includes it for you:

```@repl first
r = rule(LaguerreRay(); degree = 11);   # ∫₀^∞ f(x) e^{-x} dx
integrate(x -> x^5, r)                  # 5! = 120
```

`HermiteLine()` gives `∫ f(x) e^{-x²} dx` over the real line, and `GaussianSpace(D)` its
`D`-dimensional version. These come only in this standard form. For a shifted or scaled
weight, substitute: `∫ f(x) e^{-(x-μ)²/σ²} dx` is `σ ∫ f(μ + σt) e^{-t²} dt`, which is
`σ * integrate(t -> f(μ + σ * t), rule(HermiteLine(); degree = 20))`. The
[Catalogue](../catalogue/index.md) lists every domain and the rules available on it.

## What the rule contains

```@repl first
r = rule(Simplex{2}(); degree = 5)
npoints(r), degree(r), family(r)
nodes(r)[1], weights(r)[1]
```

The summary shows which family was used, the degree, the number of points, and whether the
weights are positive and the nodes inside the domain. It also reports how the rule was
computed: the working precision and the residual of the equations that define the rule.

The degree is listed as *claimed*. You can check it yourself with [`check`](@ref); see
[I want to verify a rule](verifying.md).

`degree(r)` can be higher than the degree you asked for. Some families only exist at certain
degrees, such as the Gauss rules on the interval, which have odd degrees only, and `rule` then
returns the next one up; see [About the degree](choosing.md#About-the-degree).

## Next steps

- [I want to choose a rule](choosing.md) if the default choice is not what you want.
- [I want arbitrary precision](precision.md) for more than 16 digits.
- [Troubleshooting](troubleshooting.md) if `rule` throws an error.
