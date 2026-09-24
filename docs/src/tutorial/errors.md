# I want an error estimate

A rule of degree `d` integrates polynomials of degree `d` exactly, but your integrand is
usually not a polynomial. There are three ways to find out how accurate the result is.

## Compare two rules

The simplest check is to integrate with two rules of different degree and compare:

```@repl err
using CubatureRules
f(x) = exp(x[1] * x[2]) / (1 + x[1]^2)
a = integrate(f, rule(Simplex{2}(); degree = 10))
b = integrate(f, rule(Simplex{2}(); degree = 20))
abs(a - b)
```

The difference estimates the error of the lower-degree result. The higher-degree one is
usually much better.

## Use an embedded rule

An embedded pair is a rule together with a coarser rule whose nodes are a subset of its own.
Both results come from the same function values, so the error estimate costs no extra
evaluations of the integrand. Gauss–Kronrod rules are built this way:

```@repl err
e = embedded(GaussKronrod(), Interval(); degree = 21);
res = integrate(x -> exp(x) * cos(3x), e; error = true)
res.value, res.error_estimate, res.neval
```

You can also build a pair yourself from any two rules that share nodes, with
[`EmbeddedRule`](@ref). The second-kind Fejér rules nest, for example:

```@repl err
fe = EmbeddedRule(rule(Fejer(2), Interval(); degree = 7), rule(Fejer(2), Interval(); degree = 3));
integrate(x -> exp(x), fe; error = true)
```

If the coarse rule's nodes are not all nodes of the fine rule, `EmbeddedRule` throws an
error.

## Check the convergence

To see how quickly a family converges on your integrand, build a few rules and compare
their results with a known value, or with each other:

```@repl err
rules = [rule(TanhSinh(m), Interval()) for m in 2:6];
[integrate(x -> 1 / sqrt(1 - x^2), r) - π for r in rules]
```

`CubatureRules.verify_convergence(rules, f, reference)` does the same and checks that the
error decreases as it should; see [I want to verify a rule](verifying.md).

## Adaptive integration

This package does not integrate to a tolerance automatically: it provides rules and leaves
the choice of rule to you. If you need adaptive integration, other packages do it well:

| Domain | Package |
|---|---|
| Interval | [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) |
| Box | [HCubature.jl](https://github.com/JuliaMath/HCubature.jl) |
| Triangle, tetrahedron | [HAdaptiveIntegration.jl](https://github.com/zmoitier/HAdaptiveIntegration.jl) |

These subdivide the domain where the integrand needs it, which also handles peaks and
singularities that a single rule cannot resolve.

## Next steps

- [I want to verify a rule](verifying.md)
- [Error estimates](../design/errors.md) explains why the package has no adaptive driver.
