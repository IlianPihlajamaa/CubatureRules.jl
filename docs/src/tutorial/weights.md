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

## Weights given by a function

If you can evaluate the weight but know nothing else about it, write it as a
[`FunctionWeight`](@ref). A weight in several pieces, or with point masses, is a sum:

```@repl w
w = FunctionWeight(x -> exp(-x), 0, 1; β = -1/2) + FunctionWeight(x -> 1 + x^2, 1, 2) + PointMass(2, 1/10);
r = rule(WeightedDomain(Interval(0, 2), w); degree = 19, digits = 30);
npoints(r), passed(verify(r))
```

The first piece is `e^{-x} / √x` on `[0, 1]`: the exponent `β = -1/2` puts the singularity
`(x − 0)^{-1/2}` in the weight's classical factor, and only the smooth `e^{-x}` is given as a
function. Each piece is integrated with the Gauss rule of its classical factor, so this
converges quickly, while a singularity inside the function converges slowly, and the
construction reports it:

```@repl w
rule(WeightedDomain(Interval(0, 1), FunctionWeight(sqrt, 0, 1)); degree = 9, digits = 30)
```

Written with the exponent instead, `FunctionWeight(x -> one(x), 0, 1; β = 1/2)`, the same
weight is no problem. Split a piece where the function has a kink.

The function is evaluated at `BigFloat` arguments, and must compute in the type of its
argument: a function that returns `Float64` describes the weight only to `Float64` and is
refused. Pieces can also lie on a half line, `FunctionWeight(g, a, Inf; β, rate)` for
`g(x) (x − a)^β e^{−rate (x − a)}`, or on the whole line, `FunctionWeight(g, -Inf, Inf)` for
`g(x) e^{−x²}`.

## Oscillatory integrands

For `∫ f(x) e^{iωx} dx` with large `ω`, put the oscillation in the weight with
[`Oscillatory`](@ref). The rule then needs only enough points for `f`:

```@repl w
r = rule(WeightedDomain(Interval(0, 1), Oscillatory(1000)); degree = 10, digits = 30);
npoints(r)
integrate(exp, r)                     # ∫₀¹ eˣ e^{1000ix} dx
(exp(big(1 + 1000im)) - 1) / (1 + 1000im)
```

The weights are complex. `Oscillatory(ω, cos)` and `Oscillatory(ω, sin)` give real rules
for `cos(ωx)` and `sin(ωx)`. The nodes include both ends of the interval.

## Principal values and finite parts

A Cauchy principal value `⨍ f(x) w(x) / (x − t) dx` is not an integral against a weight,
since the kernel changes sign at `t` and is not integrable there. It still has rules, and
[`PrincipalValue`](@ref) puts it in the same form as a weight:

```@repl w
dom = WeightedDomain(Interval(-1, 1), PrincipalValue(1//3; α = 1//2, β = 1//2));  # √(1-x²)/(x - 1/3)
r = rule(dom; degree = 20, digits = 30);
integrate(x -> one(x), r)             # -π/3
integrate(x -> x, r)                  # π/2 - π/9
```

[`FinitePart`](@ref) gives the Hadamard finite part `⨎ f(x) w(x) / (x − t)² dx` in the same
way:

```@repl w
r = rule(WeightedDomain(Interval(-1, 1), FinitePart(1//3)); degree = 20, digits = 30);
integrate(x -> one(x), r)             # -2/(1 - 1/9) = -9/4
```

Pass the smooth part `f` to `integrate`, not `f / (x − t)`. The point `t` is one of the
nodes, so `f` must be finite there, and the weights have both signs. `t` is taken exactly as
written: `0.3` is the binary number nearest 0.3, and `3//10` is 0.3.

## A singular kernel on a triangle

The weakly singular kernel of boundary-element methods, `∫_T f(y) / |y − x₀| dy` with `x₀` in
the triangle, is a weight on the triangle, [`InverseDistance`](@ref):

```@repl w
T = Simplex((0, 0), (1, 0), (0, 1));
r = rule(WeightedDomain(T, InverseDistance((0, 0))); degree = 9);
integrate(y -> 1.0, r)                # √2 log(1 + √2)
sqrt(2) * log(1 + sqrt(2))
```

Again pass only the smooth part `f`; the kernel is in the weights. `x₀` may be a vertex, on
an edge or inside, and the rule is built for that triangle and that point. See
[`DuffyGauss`](@ref) for how.

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
