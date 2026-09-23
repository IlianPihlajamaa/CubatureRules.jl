# A measure given by its moments rather than by a formula (PLAN §6 Tier 2).
#
# Every other weight in this package is a named function whose orthogonal polynomials are
# known in closed form. A weight can also arrive as a list of numbers: its integrals against
# a polynomial family that *is* known. Wheeler's algorithm turns 2n of those into the n-point
# Gauss rule, so the list is a complete description of the measure as far as the pipeline is
# concerned — no formula for w is ever evaluated.

"""
    MonicRecurrence(a, b)

The monic three-term recurrence of a polynomial family,

    π₋₁ = 0,  π₀ = 1,  π_{k+1}(x) = (x - a(k, T)) π_k(x) - b(k, T) π_{k-1}(x),

with coefficients returned in type `T`. This is the form Wheeler's algorithm takes and the
form the literature tabulates; [`monic`](@ref) converts the orthonormal
[`Recurrence`](@ref) used elsewhere in the package.
"""
struct MonicRecurrence{A,B}
    a::A
    b::B
end

"""
    monomial_recurrence()

The monomials `1, x, x², …` as a [`MonicRecurrence`](@ref): `a = b = 0`. Modified moments
against this family are the ordinary moments.
"""
monomial_recurrence() = MonicRecurrence((k, T) -> zero(T), (k, T) -> zero(T))

"""
    shift(aux::MonicRecurrence, lo, hi)

`aux`, monic-orthogonal on `[-1, 1]`, carried affinely to `[lo, hi]`. Monic families do not
survive a change of variable unchanged — `π̃ₖ(x) = hᵏ πₖ(t)` with `h = (hi - lo)/2` — but the
recurrence does, as `ãₖ = h aₖ + (hi + lo)/2` and `b̃ₖ = h² bₖ`.

Saves writing out a shifted classical family by hand:

```julia
shift(monic(jacobi_recurrence(0, 0)), 0, 1)     # monic shifted Legendre on [0,1]
```

The endpoints are converted at the working precision rather than stored, so this is exact at
any precision for rational `lo` and `hi`.
"""
shift(aux::MonicRecurrence, lo, hi) =
    MonicRecurrence((k, T) -> (T(hi) - T(lo)) / 2 * aux.a(k, T) + (T(hi) + T(lo)) / 2,
                    (k, T) -> ((T(hi) - T(lo)) / 2)^2 * aux.b(k, T))

"""
    MomentWeight(aux, moments; label = "MomentWeight")

A measure on an interval described by its **modified moments**

    mₖ = ∫ πₖ(x) w(x) dx,    k = 0, 1, 2, …

against the monic family `aux::MonicRecurrence`, where `moments(k, T)` returns `mₖ` in type
`T`. An `n`-point rule reads `mₖ` for `k < 2n + 2`. Pair it with the interval the moments
were computed on:

```julia
dom = WeightedDomain(Interval(0, 1), MomentWeight(aux, moments))
r = rule(dom; degree = 40, digits = 50)
```

The choice of `aux` is the whole game. Against the monomials ([`OrdinaryMoments`](@ref)) the
map from moments to rule is the classically ill-conditioned one, losing roughly 1.4 decimal
digits per node — which is why it is usually dismissed. Against a family orthogonal on the
same interval it loses essentially nothing. [`ModifiedChebyshev`](@ref) measures the loss
rather than assuming it, and raises working precision until the rule verifies against the
moments it came from, so the ill-conditioned case costs time instead of accuracy.

`moments` is called at whatever precision that search reaches, so it must be able to deliver
`mₖ` in `BigFloat` at the ambient precision — a closed form or an exact rational, not a
stored `Float64`. A measure known only to `Float64` cannot be rescued by precision, and the
error message says so.
"""
struct MomentWeight{A<:MonicRecurrence,M}
    aux::A
    moments::M
    label::String
end
MomentWeight(aux::MonicRecurrence, moments; label::AbstractString = "MomentWeight") =
    MomentWeight(aux, moments, String(label))

Base.show(io::IO, w::MomentWeight) = print(io, w.label)

"""
    OrdinaryMoments(moments; label = "OrdinaryMoments")

A measure described by its ordinary moments `∫ xᵏ w(x) dx`, where `moments(k, T)` returns the
`k`th in type `T`. Shorthand for a [`MomentWeight`](@ref) against
[`monomial_recurrence`](@ref).

This is the ill-conditioned case — the one the classical Chebyshev algorithm is named for and
the reason moment-based construction is usually avoided. It is supported because it is
sometimes all there is, and because arbitrary precision makes it work: the cost is about 1.4
extra digits of working precision per node.
"""
OrdinaryMoments(moments; label::AbstractString = "OrdinaryMoments") =
    MomentWeight(monomial_recurrence(), moments, String(label))

"""
    MomentInterval

An [`Interval`](@ref) carrying a measure given by its moments; see [`MomentWeight`](@ref).
"""
const MomentInterval = WeightedDomain{1,<:Any,<:Interval,<:MomentWeight}

# Modified moments are stated for the interval they were computed on, so there is no
# reference domain to map from and nothing to map: such a domain is its own reference.
# Without this the selector would build on [-1, 1] and map, silently pairing the moments
# with the wrong interval.
isreference(::MomentInterval) = true
reference(d::MomentInterval) = d

# π₀ = 1, so the zeroth modified moment is the total mass whatever the auxiliary family is.
measure(d::WeightedDomain{1,T,<:Interval,<:MomentWeight}) where {T} = d.weight.moments(0, float(T))
