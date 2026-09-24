# I want an unusual weight

If your integral has the form `∫ f(x) w(x) dx` and `w` is singular or rapidly decaying, it
is better to make `w` part of the domain than part of the integrand. The rule then places
its nodes according to `w`, and only `f` needs to be smooth.

## Built-in weights

```@repl w
using CubatureRules
integrate(x -> x^5, rule(LaguerreRay(); degree = 11))    # ∫₀^∞ x⁵ e^{-x} dx
integrate(x -> x^2, rule(HermiteLine(); degree = 9))     # ∫ x² e^{-x²} dx
```

On an interval, a Jacobi weight `(1-x)^α (1+x)^β` handles singularities at the endpoints:

```@repl w
dom = WeightedDomain(Interval(), JacobiWeight(-0.5, -0.5));   # 1/√(1-x²)
r = rule(dom; degree = 21);
integrate(x -> one(x), r)                                  # π
```

The weight is infinite at both endpoints, but the result is exact, because the singularity
is in the weight and not in the function that is being approximated.

## Other weights: using moments

For a weight the package does not know, you can describe it by its moments: the integrals of
a set of known polynomials against the weight. From `2n` moments the package computes the
`n`-point Gauss rule for the weight, using Wheeler's algorithm.

As an example, take `w(x) = log(1/x)` on `[0, 1]`. Its moments against the monic shifted
Legendre polynomials are known in closed form:

```@repl w
import CubatureRules: monic, shift, jacobi_recurrence
aux = shift(monic(jacobi_recurrence(0, 0)), 0, 1)    # monic shifted Legendre on [0, 1]
logw = MomentWeight(aux,
                    (k, T) -> k == 0 ? one(T) :
                              T((-1)^k) / (T(k) * (k + 1) * T(binomial(big(2k), big(k))));
                    label = "log(1/x) on [0,1]")
r = rule(WeightedDomain(Interval(0, 1), logw); degree = 41, digits = 30);
npoints(r)
integrate(x -> one(x), r)     # ∫₀¹ log(1/x) dx = 1
integrate(x -> x^3, r)        # 1/16
```

The moment function takes an index `k` and a number type `T`, and must return the `k`th
moment in type `T`. It is called at whatever precision the computation needs, so it should
compute the moments (for example from a formula or as exact fractions) rather than return
stored `Float64` values.

## Ordinary moments

If you only know the ordinary moments `∫ xᵏ w(x) dx`, use [`OrdinaryMoments`](@ref):

```@repl w
lebesgue = OrdinaryMoments((k, T) -> one(T) / (k + 1); label = "Lebesgue on [0,1]")
r = rule(WeightedDomain(Interval(0, 1), lebesgue); degree = 79, digits = 50);
npoints(r)
```

Computing a Gauss rule from ordinary moments is a badly conditioned problem: about 1.4
decimal digits are lost for every point. In `Float64` the computation breaks down before 16
points. The package compensates by raising the working precision until the result is stable,
and records how much was needed:

```@repl w
c = certificate(r);
c.cond, c.guard_digits
```

With a better choice of polynomials, such as the shifted Legendre polynomials above, the
problem is well conditioned and no extra precision is needed:

```@repl w
certificate(rule(WeightedDomain(Interval(0, 1), logw); degree = 79, digits = 50)).cond
```

So use modified moments against a suitable family when you can, and ordinary moments when
they are all you have. [Measures given by moments](../design/moments.md) has the details and
measurements.

## Next steps

- [I want to verify a rule](verifying.md)
- [Measures given by moments](../design/moments.md)
