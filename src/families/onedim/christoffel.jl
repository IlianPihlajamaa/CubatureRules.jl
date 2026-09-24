# Christoffel modification: a weight multiplied or divided by a linear factor
# (PLAN §6 Tier 2, v0.5).
#
# Given a measure μ whose monic orthogonal polynomials πₖ are known — Legendre, Jacobi,
# Laguerre — and a point z outside its support, the Gauss rules of |x - z| dμ and of
# dμ / |x - z| follow from their modified moments against those same πₖ, which Wheeler's
# algorithm turns into a recurrence. Against a measure's own orthogonal polynomials the map is
# well conditioned, so this is the easy case of the moment machinery.
#
# Multiplying needs almost nothing. With s the sign of x - z on the support and
# x - z = π₁ + (α₀ - z),
#
#     ∫ πₖ |x - z| dμ = s ((α₀ - z) β₀ δₖ₀ + β₀β₁ δₖ₁),
#
# two nonzero moments. Dividing gives ∫ πₖ dμ / |x - z| = -s ρₖ(z), where
#
#     ρₖ(z) = ∫ πₖ(x) dμ(x) / (z - x)
#
# are the functions of the second kind. They satisfy the same three-term recurrence as the πₖ
# (with ρ₋₁ = 1), but as its *minimal* solution, so forward recursion loses them to rounding.
# They are computed backwards instead, as a continued fraction, deepened until the result
# stops changing. Convergence is geometric on a bounded interval and slower on the half line,
# and slows as z approaches the support.

"""
    legendre_monic(a, b)

The monic Legendre polynomials on `[a, b]` as a [`MonicRecurrence`](@ref), with exact rational
coefficients when `a` and `b` are rational.
"""
legendre_monic(a, b) =
    MonicRecurrence((k, T) -> (T(a) + T(b)) / 2,
                    (k, T) -> k == 0 ? zero(T) : ((T(b) - T(a)) / 2)^2 * T(big(k)^2 // (4 * big(k)^2 - 1)))

# The monic recurrence, mass, base domain, support and a readable name, for the domains whose
# orthogonal polynomials the package knows.
christoffel_base(d::Interval) =
    (legendre_monic(d.a, d.b), T -> T(d.b) - T(d.a), d, d.a, d.b, "1 on [$(d.a), $(d.b)]")
function christoffel_base(d::WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight})
    isreference(d.base) ||
        throw(ArgumentError("a JacobiWeight is defined on [-1, 1], but was paired with $(d.base)"))
    α, β = d.weight.α, d.weight.β
    rec = jacobi_recurrence(α, β)
    return (monic(rec), rec.mass, d.base, -1, 1, "(1-x)^$α (1+x)^$β on [-1, 1]")
end
function christoffel_base(d::LaguerreDomain)
    α = d.weight.α
    rec = laguerre_recurrence(α)
    return (monic(rec), rec.mass, d.base, 0, Inf, iszero(α) ? "e^-x on [0, ∞)" : "x^$α e^-x on [0, ∞)")
end
christoffel_base(d) =
    throw(ArgumentError("christoffel needs a domain whose orthogonal polynomials are known — an " *
                        "Interval, an Interval with a JacobiWeight, or a LaguerreRay — but got $d"))

"""
    second_kind(aux, mass, z, K, T) -> Vector{T}

The functions of the second kind `ρₖ(z) = ∫ πₖ(x) dμ(x) / (z - x)`, for `k = 0, …, K`, of the
measure with monic recurrence `aux` and total mass `mass(T)`. One backward pass on
`rⱼ = ρⱼ / ρⱼ₋₁ = βⱼ / ((z - αⱼ) - rⱼ₊₁)`, with `β₀` the mass, gives all of them; the starting
depth doubles until every `ρₖ` agrees with the previous depth to the working precision.
"""
function second_kind(aux::MonicRecurrence, mass, z::Real, K::Integer, ::Type{T}) where {T}
    prev = nothing
    N = 2K + 64
    tol = 16 * eps(T)
    rs = Vector{T}(undef, K + 1)
    while true
        r = zero(T)
        for j in N:-1:0
            βj = j == 0 ? T(mass(T)) : aux.b(j, T)
            r = βj / ((T(z) - aux.a(j, T)) - r)
            j <= K && (rs[j + 1] = r)
        end
        ρ = cumprod(rs)
        prev !== nothing && all(i -> abs(ρ[i] - prev[i]) <= tol * abs(ρ[i]), eachindex(ρ)) && return ρ
        prev = ρ
        N *= 2
        N > 1 << 20 &&
            throw(ArgumentError("the functions of the second kind did not converge at z = $z within " *
                                "$(N ÷ 2) terms; z is too close to the support for this construction"))
    end
end

# The moments of a divided weight are read one at a time, but one backward pass gives all of
# them, and on the half line that pass is long. So they are kept per number type and working
# precision, inside this one weight, and extended when a higher index is asked for. This
# stores intermediate values of a single moment function; no rule is cached.
function second_kind_moments(aux, mass, z, s)
    store = Dict{Tuple{Type,Int},Vector}()
    guard = ReentrantLock()
    return function (k, ::Type{T}) where {T}
        key = (T, T === BigFloat ? precision(BigFloat) : 0)
        ρ = lock(guard) do
            v = get(store, key, nothing)
            if v === nothing || length(v) <= k
                v = second_kind(aux, mass, z, max(2k + 2, v === nothing ? 0 : 2 * length(v)), T)
                store[key] = v
            end
            v
        end
        return -s * T(ρ[k + 1])
    end
end

_factor_label(z) = iszero(z) ? "|x|" : z < 0 ? "|x + $(-z)|" : "|x - $z|"

"""
    christoffel(domain, z; power = 1)

`domain` with its weight multiplied by `|x - z|` (`power = 1`) or divided by it
(`power = -1`), for a point `z` outside the support. Returns a new weighted domain whose rules
are built by [`ModifiedChebyshev`](@ref).

```julia
christoffel(Interval(), -2)                        # (x + 2) on [-1, 1]
christoffel(Interval(), -2; power = -1)            # 1 / (x + 2) on [-1, 1]
christoffel(LaguerreRay(), -1; power = -1)         # e^-x / (x + 1) on [0, ∞)
```

`domain` must be one whose orthogonal polynomials are known: an [`Interval`](@ref), an
`Interval` with a [`JacobiWeight`](@ref), or a [`LaguerreRay`](@ref). For `power = 1`, `z` may
be an endpoint; for `power = -1` it must lie strictly outside.

The modified moments of the new weight against the old weight's orthogonal polynomials are
known in closed form (for `power = 1`) or through the functions of the second kind (for
`power = -1`), and the construction is well conditioned in both cases.
"""
function christoffel(dom, z::Real; power::Integer = 1)
    power in (1, -1) ||
        throw(ArgumentError("power must be 1 (multiply by |x - z|) or -1 (divide by it), got $power"))
    aux, mass, base, lo, hi, name = christoffel_base(dom)
    inside = power == 1 ? lo < z < hi : lo <= z <= hi
    inside && throw(ArgumentError(
        "z = $z lies " * (power == 1 ? "inside" : "in the closed hull of") * " the support [$lo, $hi]; " *
        (power == 1 ? "|x - z| would not be a polynomial factor there" :
                      "dividing by |x - z| would not leave an integrable positive weight")))
    s = z <= lo ? 1 : -1                     # the sign of x - z on the support
    factor = _factor_label(z)
    if power == 1
        moments = (k, T) -> k == 0 ? s * (aux.a(0, T) - T(z)) * T(mass(T)) :
                            k == 1 ? s * T(mass(T)) * aux.b(1, T) : zero(T)
        label = "$factor × $name"
    else
        moments = second_kind_moments(aux, mass, z, s)
        label = "$name / $factor"
    end
    return WeightedDomain(base, MomentWeight(aux, moments; label, support = (lo, hi)))
end
