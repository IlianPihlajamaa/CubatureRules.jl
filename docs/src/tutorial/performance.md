# I want performance

The two most important points: build the rule outside your loop, and convert it to a static
rule.

## Build the rule once

`rule` computes the rule every time it is called; nothing is cached. Call it once and keep
the result:

```julia
r = rule(Simplex{2}(); degree = 6)      # once
for cell in cells
    total += integrate(f, r, cell)      # many times
end
```

## Static rules

```@repl perf
using CubatureRules
r = rule(Simplex{2}(); degree = 6);
s = static(r);
typeof(s)
```

[`static`](@ref) stores the nodes and weights as `SVector`s in a tuple. The rule then has a
fixed size known to the compiler, and integrating with it does not allocate:

```@repl perf
f(x) = exp(x[1] * x[2])
integrate(f, s)
@allocated integrate(f, s)
```

This is most useful when the integrand is cheap and you integrate many times. Each number of
points gives a different type, so code that uses many different static rules is compiled
once for each.

## Meshes

```@repl perf
cells = [Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
         Simplex((1.0, 0.0), (1.0, 1.0), (0.0, 1.0))]
integrate(f, s, cells)
```

The map to each cell is applied to the nodes as they are used, so nothing is allocated per
cell. For large meshes you can use threads:

```@repl perf
integrate(f, s, cells; threaded = true)
```

## Passing integrands through functions

If you write a function that takes an integrand and passes it on to `integrate`, add a type
parameter for it:

```julia
my_integral(f, r) = integrate(f, r)                   # slow
my_integral(f::F, r) where {F} = integrate(f, r)      # fast
```

Julia does not specialise on a function argument that is only passed on to another
function. Without the type parameter every evaluation of `f` goes through dynamic dispatch.

## Vectorised integrands

If your integrand is faster when evaluated on many points at once, you can work with the
nodes and weights directly:

```@repl perf
xs = nodes(r); ws = weights(r);
sum(ws .* map(f, xs))
```

`integrate(f, r; batch = true)` passes all nodes to `f` in a single call instead.

## Precision and speed

Rules are always constructed in `BigFloat`, even when you ask for `Float64`, so construction
is slower than a table lookup. Once built, a `Float64` rule is used with ordinary `Float64`
arithmetic.

## Next steps

- [Application and transport](../design/application.md)
- [The pipeline](../design/pipeline.md) explains what construction involves.
