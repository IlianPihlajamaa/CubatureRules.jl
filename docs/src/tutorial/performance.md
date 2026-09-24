# I want performance

The most important point is to build the rule once, outside your loop.

## Build the rule once

`rule` computes the rule every time it is called; nothing is cached. Call it once and keep
the result:

```julia
r = rule(Simplex{2}(); degree = 6)      # once
for cell in cells
    total += integrate(f, r, cell)      # many times
end
```

## Integration does not allocate

Once a rule is built, integrating with it allocates nothing, as long as the integrand does
not allocate either:

```@repl perf
using CubatureRules
r = rule(Simplex{2}(); degree = 6);
f(x) = exp(x[1] * x[2])
integrate(f, r)
@allocated integrate(f, r)
```

The cost is one evaluation of `f` per node, plus a multiply and an add.

## Meshes

```@repl perf
cells = [Simplex((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)),
         Simplex((1.0, 0.0), (1.0, 1.0), (0.0, 1.0))]
integrate(f, r, cells)
```

The map to each cell is applied to the nodes as they are used, so nothing is allocated per
cell. For large meshes you can use threads:

```@repl perf
integrate(f, r, cells; threaded = true)
```

## Passing integrands through your own functions

Julia compiles a function that only passes a function argument on, without calling it, once
for all functions rather than once per function. If you write such a wrapper and Julia does
not inline it, each call makes one dynamic dispatch into `integrate`:

```julia
my_integral(f, r) = integrate(f, r)
my_integral(f::F, r) where {F} = integrate(f, r)      # specialised on f
```

Evaluations of `f` inside `integrate` are not affected, so the cost is fixed per call. On
one machine it was about 50 ns and 176 bytes: 108 ns instead of 57 ns for a 12-point rule,
and 391 ns instead of 343 ns for a 79-point one. It matters only when a small rule is
applied many times, for example once per mesh cell. The type parameter removes it.

## Vectorised integrands

If your integrand is faster when evaluated on many points at once, pass `batch = true`. `f`
is then called once, with all nodes as a `D × N` matrix (a vector in one dimension), and
must return the `N` values:

```@repl perf
integrate(X -> exp.(X[1, :] .* X[2, :]), r; batch = true)
```

## Rules with a fixed size

Some code needs a rule of fixed size, for example inside a GPU kernel, where everything must
be `isbits`. Build one from the nodes and weights, as tuples or as `SVector`s from
StaticArrays.jl:

```@repl perf
xs = Tuple(nodes(r)); ws = Tuple(weights(r));
isbits(xs), isbits(ws)
```

For ordinary use on the CPU this is not faster than the rule itself; see
[Application and transport](../design/application.md) for the measurements.

## Construction time

Rules are always constructed in `BigFloat`, even when you ask for `Float64`, so construction
is slower than a table lookup. Once built, a `Float64` rule is used with ordinary `Float64`
arithmetic.

## Next steps

- [Application and transport](../design/application.md)
- [The pipeline](../design/pipeline.md) explains what construction involves.
