# Named singular weights (PLAN §6 Tier 2, v0.5).
#
# A weight with an integrable endpoint singularity has no classical Gauss family, but its
# ordinary moments are often exact rationals. Fed to Wheeler's algorithm directly those would
# be the badly conditioned case (see `OrdinaryMoments`). Converting them to modified moments
# against the shifted Legendre polynomials *in exact arithmetic*, before any rounding, avoids
# the problem altogether: the only rounding happens after the change of basis, and Wheeler's
# algorithm on modified moments is well conditioned.

"""
    shifted_legendre_recurrence()

The monic shifted Legendre polynomials on `[0, 1]` as a [`MonicRecurrence`](@ref), with
exact rational coefficients `aₖ = 1/2` and `bₖ = k² / (4(4k² - 1))`.
"""
shifted_legendre_recurrence() =
    MonicRecurrence((k, T) -> one(T) / 2,
                    (k, T) -> k == 0 ? zero(T) : T(big(k)^2 // (4 * (4 * big(k)^2 - 1))))

"""
    log_weight_moment(k, α, power) -> Rational{BigInt}

The `k`th modified moment of `x^α log(1/x)^power` on `[0, 1]` against the monic shifted
Legendre polynomials, exactly. The monic polynomial has coefficients
`(-1)^(k+j) C(k,j) C(k+j,j) / C(2k,k)`, and the ordinary moments are
`∫₀¹ xʲ⁺ᵅ log(1/x)^m dx = m! / (j + α + 1)^(m+1)`.
"""
function log_weight_moment(k::Integer, α::Rational{BigInt}, power::Integer)
    fm = factorial(big(power))
    s = zero(Rational{BigInt})
    for j in 0:k
        c = (-1)^(k + j) * binomial(big(k), big(j)) * binomial(big(k + j), big(j))
        s += c * fm // (j + α + 1)^(power + 1)
    end
    return s // binomial(big(2k), big(k))
end

"""
    LogWeight(α = 0; power = 1)

The weight `x^α log(1/x)^power` on `[0, 1]`, for `α > -1` and a non-negative integer
`power`: a logarithmic singularity at the origin, together with an algebraic one when
`α ≠ 0`. Pair it with `Interval(0, 1)`; any other interval is refused.

```julia
dom = WeightedDomain(Interval(0, 1), LogWeight())          # log(1/x)
dom = WeightedDomain(Interval(0, 1), LogWeight(-1/2))      # log(1/x) / √x
r = rule(dom; degree = 39, digits = 50)
```

Rules come from [`ModifiedChebyshev`](@ref). The ordinary moments are exact rationals (a
floating-point `α` is itself an exact rational), so the modified moments against the shifted
Legendre polynomials are computed exactly before anything is rounded, and the construction
is well conditioned at every degree.
"""
function LogWeight(α::Real = 0; power::Integer = 1)
    α > -1 || throw(ArgumentError("LogWeight needs α > -1, got α = $α"))
    power >= 0 || throw(ArgumentError("LogWeight needs a non-negative integer power, got $power"))
    a = Rational{BigInt}(α)
    astr = α isa Rational ? string(numerator(α), "/", denominator(α)) : string(α)
    xpart = iszero(a) ? "" : isinteger(a) && a > 0 ? "x^$(astr) " : "x^($(astr)) "
    lpart = power == 0 ? "" : power == 1 ? "log(1/x)" : "log(1/x)^$power"
    label = strip(xpart * lpart) * " on [0,1]"
    return MomentWeight(shifted_legendre_recurrence(), (k, T) -> T(log_weight_moment(k, a, power));
                        label, support = (0, 1))
end
