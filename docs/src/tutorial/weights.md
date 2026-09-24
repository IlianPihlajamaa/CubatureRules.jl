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

## Logarithmic singularities

[`LogWeight`](@ref)`(α; power)` is the weight `x^α log(1/x)^power` on `[0, 1]`: a logarithmic
singularity at the origin, combined with an algebraic one when `α ≠ 0`.

```@repl w
r = rule(WeightedDomain(Interval(0, 1), LogWeight()); degree = 19, digits = 30);
integrate(x -> one(x), r)            # ∫₀¹ log(1/x) dx = 1
r = rule(WeightedDomain(Interval(0, 1), LogWeight(-1/2)); degree = 19);
integrate(x -> cos(x), r)            # ∫₀¹ cos(x) log(1/x) / √x dx
```

The weight is only defined on `[0, 1]`, and pairing it with another interval is an error. For
a logarithmic singularity at another point, change variables so that it sits at 0.

## Multiplying or dividing by a linear factor

`CubatureRules.christoffel(domain, z; power)` multiplies the weight of a domain by `|x − z|`
(`power = 1`) or divides it by `|x − z|` (`power = -1`), for a point `z` outside the domain:

```@repl w
using CubatureRules: christoffel
r = rule(christoffel(Interval(), -2; power = -1); degree = 19);   # 1/(x + 2) on [-1, 1]
integrate(x -> one(x), r)                                         # log 3
r = rule(christoffel(LaguerreRay(), -1; power = -1); degree = 19); # e^-x/(x + 1) on [0, ∞)
integrate(x -> one(x), r)                                         # e E₁(1) ≈ 0.596
```

The domain must be one whose orthogonal polynomials the package knows: an `Interval`, an
`Interval` with a `JacobiWeight`, or a `LaguerreRay`.

## Other weights: using moments

For a weight the package does not know, you can describe it by its moments: the integrals of
a set of known polynomials against the weight. From `2n` moments the package computes the
`n`-point Gauss rule for the weight, using Wheeler's algorithm.

`LogWeight` is built this way. Written out by hand for `w(x) = log(1/x)`, whose moments
against the monic shifted Legendre polynomials are known in closed form, it looks like this:

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

Moments only describe the weight on the interval they were computed for, and a domain with a
moment weight is never mapped to another interval. Pass `support = (0, 1)` to
`MomentWeight` to record the interval, so that pairing it with any other one is refused.

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
