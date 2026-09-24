# I want to integrate to a tolerance

If you don't know which degree you need, give a tolerance instead:

```@repl tol
using CubatureRules
integrate(x -> exp(x[1] * x[2]), Simplex{2}(); rtol = 1e-12)
```

This form of `integrate` takes a domain rather than a rule. It builds rules of increasing
degree and stops when two successive results agree to within the tolerance.

## The result

The return value is an [`IntegrationResult`](@ref), not just a number:

```@repl tol
res = integrate(x -> exp(x[1] * x[2]), Simplex{2}(); rtol = 1e-12);
res.value, res.error_estimate, res.neval, res.converged
res.family, res.degree
```

It contains the value, an error estimate, the number of function evaluations, whether the
tolerance was reached, and which rule produced the final value. If `converged` is `false`,
the highest degree was reached before the tolerance was met, and `value` is the best
estimate available at that point.

## Unbounded domains

```@repl tol
integrate(x -> exp(-x^2), RealLine(); rtol = 1e-10)
integrate(x -> exp(-x) / sqrt(x), HalfLine(); rtol = 1e-10)
```

On unbounded domains without a weight there is no family with a polynomial degree, so the
package uses double-exponential rules (`SinhSinh`, `ExpSinh`) and increases their level
instead of the degree.

## Options

```@repl tol
integrate(x -> exp(x), Interval(); rtol = 1e-10, atol = 1e-14)
integrate(x -> exp(x), Interval(); rtol = 1e-10, family = GaussLegendre())
integrate(x -> exp(x), Interval(); rtol = 1e-30, digits = 40)
```

- `rtol` and `atol`: the result is accepted when the error estimate is below either one.
  `rtol` is required; there is no default tolerance.
- `family`: use a specific family instead of the one the package picks.
- `T` or `digits`: the number type, as for [`rule`](@ref).
- `maxdegree`: the highest degree to try (default 64).

## Limitations

The degree is raised for the whole domain at once; the domain is never subdivided. This
works well for smooth integrands, where the error drops quickly as the degree grows. It
works badly for integrands with a sharp peak or a singularity inside the domain. For those,
use [QuadGK.jl](https://github.com/JuliaMath/QuadGK.jl) in one dimension or
[HCubature.jl](https://github.com/JuliaMath/HCubature.jl) in several, which subdivide the
domain where needed.

If the singularity is at a known place, it is often better to move it into the weight of the
domain; see [I want an unusual weight](weights.md).

## Sequences of rules

The rules that `integrate` walks through are available directly, for example to study
convergence yourself:

```@repl tol
seq = RuleSequence(GaussLegendre(), Interval(); degrees = 1:2:15);
for r in seq
    println(npoints(r), " points: ", integrate(x -> exp(x), r))
end
```

[`LevelSequence`](@ref) does the same for families with levels instead of degrees, such as
`TanhSinh`.

## Next steps

- [I want performance](performance.md)
- [Sequences and adaptivity](../design/adaptive.md) explains how the error estimate is
  computed.
