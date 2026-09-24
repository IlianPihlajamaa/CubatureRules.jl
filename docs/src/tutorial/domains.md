# I want my own domain or mesh

Rules are constructed on a *reference* domain: the interval `[-1, 1]`, the triangle with
vertices `(0,0)`, `(1,0)` and `(0,1)`, the unit sphere, and so on. There are several ways to
use them on your own geometry.

## A rule on your domain

Pass your domain to `rule` directly:

```@repl dom
using CubatureRules
t = Simplex((0.0, 0.0), (2.0, 0.0), (0.0, 3.0))
r = rule(t; degree = 6);
domain(r) == t, sum(weights(r))     # the weights add up to the area
```

The rule is built on the reference triangle and mapped to yours. The map is affine, so the
rule keeps its degree.

The same works for intervals, boxes and spheres:

```@repl dom
rule(Interval(0, 1); degree = 9);
rule(Sphere((1.0, 2.0, 3.0), 2.5); degree = 9);
```

## A mesh

For many cells, build one reference rule and pass it to `integrate` together with the cells:

```@repl dom
cells = [Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
         Simplex((1.0, 0.0), (1.0, 1.0), (0.0, 1.0))]
r = rule(Simplex{2}(); degree = 6);
integrate(x -> x[1]^2 + x[2]^2, r, cells)       # sum over all cells
integrate(x -> x[1]^2 + x[2]^2, r, cells[1])    # a single cell
```

The nodes are mapped as they are used, so no new rule is built for each cell. See
[I want performance](performance.md) for a version that does not allocate.

If you need the mapped rule itself, use [`map_to`](@ref):

```@repl dom
mapped = map_to(r, cells[2]);
domain(mapped), degree(mapped)
```

(A triangle with the reference vertices, like `cells[1]`, prints as `Simplex{2}()`.)

## Subdividing

```@repl dom
r = rule(Simplex{2}(); degree = 4);
sub = subdivide(r, Simplex{2}(), 2);
npoints(sub), degree(sub)
```

[`subdivide`](@ref) splits the domain into smaller pieces of the same shape and puts a copy of
the rule on each. The degree stays the same and the number of points is multiplied. This can
help for integrands that are not smooth, where raising the degree does not help much.

## Non-affine maps

For an integrand with a singularity at a vertex of a simplex, [`duffy`](@ref) moves the
nodes of a rule towards the vertex at the origin:

```@repl dom
r = rule(Simplex{2}(); degree = 10);
d = duffy(r; power = 3);
exactness(d)
integrate(x -> 1 / sqrt(x[1]^2 + x[2]^2), d)     # singular at the origin
```

The map `x ↦ (x₁ + … + x_D)^(q-1) x` is not affine, so polynomials do not stay polynomials
and the new rule has no degree: its exactness is [`NoClaim`](@ref). The same applies to a
map you supply yourself with [`transform`](@ref). If you know that the degree is preserved
for your particular map, you can pass it as `claim`, and the provenance records that the
claim came from you.

The grading cancels a singularity of the form `r^{-(q-1)D}` at the vertex, where `q` is the
`power`. Apply `duffy` on the reference simplex and map the result to your triangle
afterwards.

## Summary

| Operation | Effect on the degree |
|---|---|
| `rule(domain; degree)` on a mapped domain | kept |
| `map_to(r, dom)` (affine) | kept |
| `subdivide(r, dom, n)` | kept |
| `transform(r, φ, Jφ)` | lost (`NoClaim`), unless you pass `claim` |
| `duffy(r; power)` | lost (`NoClaim`) |

## Next steps

- [I want performance](performance.md)
- [Application and transport](../design/application.md) explains how mapping is
  implemented.
