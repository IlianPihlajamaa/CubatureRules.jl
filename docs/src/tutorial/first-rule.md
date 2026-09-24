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
for cell in mesh
    total += integrate(f, r, cell)
end
```

For the same reason `integrate` always takes a rule, never a domain. To judge how accurate
the result is, see [I want an error estimate](errors.md).

## Domains

The most common domains:

```@repl first
rule(Interval(); degree = 9)          # the interval [-1, 1]
rule(Orthotope{3}(); degree = 9)      # the cube [-1, 1]³
rule(Sphere{3}(); degree = 9)         # the unit sphere in ℝ³
rule(Ball{3}(); degree = 9)           # the unit ball in ℝ³
```

`Simplex{D}()` works in any dimension, and `Disk()` is the unit disk.

On an unbounded domain an integral usually only converges with a weight, so these domains
come with one. The weight is part of the domain, and `integrate` includes it for you:

```@repl first
r = rule(LaguerreRay(); degree = 11);   # ∫₀^∞ f(x) e^{-x} dx
integrate(x -> x^5, r)                  # 5! = 120
```

`HermiteLine()` gives `∫ f(x) e^{-x²} dx` over the real line, and `GaussianSpace(D)` its
`D`-dimensional version. The [Catalogue](../catalogue/index.md) lists every domain and the
rules available on it.

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

`degree(r)` can be higher than the degree you asked for. Many families only exist at certain
degrees, and `rule` then returns the next one up.

## Next steps

- [I want to choose a rule](choosing.md) if the default choice is not what you want.
- [I want arbitrary precision](precision.md) for more than 16 digits.
- [Troubleshooting](troubleshooting.md) if `rule` throws an error.
