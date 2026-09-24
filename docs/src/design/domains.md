# Domains and measures

A domain is a subtype of `Domain{D,T}`, where `D` is the dimension of the points and `T` the
type used to store the domain's own parameters (vertices, centre, radius).

## Reference and mapped domains

Rule families only construct rules on *reference* domains:

| Domain | Reference form |
|---|---|
| `Interval` | `[-1, 1]` |
| `Simplex{D}` | vertices `0, e₁, …, e_D` |
| `Orthotope{D}` | `[-1, 1]^D` |
| `Sphere{D}` | unit sphere in `ℝ^D`, centred at the origin |
| `Ball{D}` | unit ball in `ℝ^D`, centred at the origin |

`Simplex{2}()`, `Interval()` and so on construct the reference domain. A domain with other
vertices, endpoints, centre or radius is a mapped domain. `isreference(dom)` tells the two
apart, and `reference(dom)` returns the reference domain of the same kind.

When you call `rule` on a mapped domain, the rule is built on the reference domain and then
mapped with [`map_to`](@ref). All supported maps are affine (or similarities, for spheres),
so they preserve polynomial degree; see [Exactness claims](claims.md).

The advantage of this split is that families only ever need to handle one domain of each
kind, and the mapping code is written once.

## Weighted domains

A `WeightedDomain(base, weight)` is a domain with the measure `weight(x) dx`. The weight is
part of the domain, so a rule on a weighted domain integrates `∫ f(x) w(x) dx` when applied
to `f`. The built-in weights are:

| Weight | Domain | Shorthand |
|---|---|---|
| `JacobiWeight(α, β)` | `Interval()` | |
| `ExponentialWeight(α)`: `x^α e^{-x}` | `HalfLine()` | `LaguerreRay(α)` |
| `GaussianWeight()`: `e^{-|x|²}` | `RealLine()`, `RealSpace{D}()` | `HermiteLine()`, `GaussianSpace(D)` |
| `MomentWeight(aux, moments)` | an `Interval` | |

The unbounded base domains `HalfLine`, `RealLine` and `RealSpace{D}` have infinite measure;
they are mostly used as the base of a weighted domain. The double-exponential families
(`ExpSinh`, `SinhSinh`) do produce rules on them directly, but with [`NoClaim`](@ref).

Moving a singular or decaying factor from the integrand into the weight usually helps a lot:
the family can then place nodes according to the weight, and the remaining integrand is
smooth.

## Moment-defined measures

A [`MomentWeight`](@ref) describes a measure through its moments against a known family of
polynomials. Those moments are only valid on the interval they were computed for, so a
`WeightedDomain(Interval(a, b), MomentWeight(...))` is treated as its own reference domain
and is never mapped. See [Measures given by moments](moments.md).

## Measure

`measure(dom)` returns the total mass of the domain: its length, area or volume, or the
integral of its weight. Verification checks that the weights of a rule add up to it.
