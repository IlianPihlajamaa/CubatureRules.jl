# Cauchy principal-value and Hadamard finite-part integrals (PLAN §6 Tier 2, v0.5).
#
#     ⨍ₐᵇ f(x) w(x) / (x − t) dx        and        ⨎ₐᵇ f(x) w(x) / (x − t)² dx
#
# with a Jacobi weight w. Neither is an integral against a measure: the kernel changes sign
# at t and is not integrable there. Both are still linear functionals on smooth f, and a rule
# for them is exact on polynomials up to some degree like any other, with signed weights and
# with t among the nodes.
#
# Everything follows from one number, the Hilbert transform of the weight
# ρ₀(t) = ⨍ w/(x − t), and its derivative σ₀(t) = ρ₀'(t) = ⨎ w/(x − t)². For the orthogonal
# polynomials q_k of w, ρ_k(t) = ⨍ w q_k/(x − t) satisfies the recurrence of q_k with one
# extra term at k = 0, and σ_k = ρ_k' satisfies its derivative (families/onedim/singular_gauss.jl).

"""
    SingularKernel

The functional `f ↦ ⨍ f(x) w(x) / (x - t)^order dx` on an interval, with a Jacobi weight
`w`: a Cauchy principal value for `order = 1`, a Hadamard finite part for `order = 2`.
Made by [`PrincipalValue`](@ref) and [`FinitePart`](@ref).
"""
struct SingularKernel
    t::Real
    order::Int
    α::Real
    β::Real
    function SingularKernel(t::Real, order::Integer, α::Real, β::Real)
        isfinite(t) || throw(ArgumentError("the singular point must be finite, got t = $t"))
        (α > -1 && β > -1) || throw(ArgumentError("the weight (b - x)^α (x - a)^β needs α, β > -1, got ($α, $β)"))
        order in (1, 2) || throw(ArgumentError("order must be 1 (principal value) or 2 (finite part), got $order"))
        return new(t, Int(order), α, β)
    end
end

"""
    PrincipalValue(t; α = 0, β = 0)

The Cauchy principal value at `t`,

    ⨍ₐᵇ f(x) (b - x)^α (x - a)^β / (x - t) dx = lim_{ε→0} (∫ₐ^{t-ε} + ∫_{t+ε}^b) …,

as the weight of a [`WeightedDomain`](@ref) on an `Interval(a, b)` with `a < t < b`. `t` is
taken exactly as given: `0.3` is the binary number nearest 0.3, and `3//10` or `big"0.3"`
the decimal.

```julia
dom = WeightedDomain(Interval(-1, 1), PrincipalValue(3//10))
r = rule(dom; degree = 20, digits = 40)
integrate(exp, r)                        # ⨍₋₁¹ eˣ / (x - 0.3) dx
```

Rules come from [`SingularGauss`](@ref). Their weights are signed, and `t` is one of the
nodes, so the integrand must be finite there: the rule is for `f`, not for `f / (x - t)`.
"""
PrincipalValue(t::Real; α::Real = 0, β::Real = 0) = SingularKernel(t, 1, α, β)

"""
    FinitePart(t; α = 0, β = 0)

The Hadamard finite part at `t`,

    ⨎ₐᵇ f(x) (b - x)^α (x - a)^β / (x - t)² dx,

the derivative with respect to `t` of the principal value with the same weight. Use it as
the weight of a [`WeightedDomain`](@ref) on an `Interval(a, b)` with `a < t < b`; see
[`PrincipalValue`](@ref) for how `t` is read.

```julia
dom = WeightedDomain(Interval(-1, 1), FinitePart(1//2; α = 1//2, β = 1//2))
r = rule(dom; degree = 20, digits = 40)
```
"""
FinitePart(t::Real; α::Real = 0, β::Real = 0) = SingularKernel(t, 2, α, β)

function Base.show(io::IO, k::SingularKernel)
    print(io, k.order == 1 ? "PrincipalValue(" : "FinitePart(", k.t)
    (iszero(k.α) && iszero(k.β)) || print(io, "; α = ", k.α, ", β = ", k.β)
    print(io, ")")
end

"""
    KernelDomain

An [`Interval`](@ref) carrying a [`PrincipalValue`](@ref) or [`FinitePart`](@ref) kernel.
"""
const KernelDomain = WeightedDomain{1,<:Any,<:Interval,SingularKernel}

# The singular point is stated in the interval's own coordinates, and the finite part does
# not scale like a measure under a change of variable, so the rule is built where it is
# asked for: such a domain is its own reference.
isreference(::KernelDomain) = true
reference(d::KernelDomain) = d

function check_kernel(dom::KernelDomain)
    a, b, t = dom.base.a, dom.base.b, dom.weight.t
    (isfinite(a) && isfinite(b)) || throw(ArgumentError("$(dom.weight) needs a finite interval"))
    a < t < b || throw(ArgumentError("the singular point t = $t must lie inside ($a, $b); outside it the integral is not singular"))
    return nothing
end

# --- the Hilbert transform of a Jacobi weight -------------------------------------------------
#
# On [0, 1], with v = 1 − u,
#
#     H(α, β, u) = ⨍₀¹ sᵅ (1 − s)ᵝ / (s − u) ds,
#
# and H' = dH/du. For non-integer α (Gautschi & Wimp 1987)
#
#     H = −π cot(πα) uᵅ vᵝ + B(α, β + 1) ₂F₁(1, −α − β; 1 − α; u),
#
# a series in u, used for u ≤ 2/3. An integer exponent makes cot(πα) and Γ(α) infinite; a
# positive one is removed exactly, sᵐ = uᵐ + (s − u)(…), leaving exponent 0, for which
#
#     H(0, β, u) = vᵝ [log(v/u) − ψ(β + 1) − γ − Σ_{j≥1} C(β, j) rʲ / j],   r = u/v,
#
# used for u ≤ 1/3. Reflection, H(α, β, u) = −H(β, α, v), moves every other case into one of
# these, so each series converges at least like (2/3)ʲ.

# ₂F₁(1, b; c; z) and its z-derivative, 0 < z ≤ 2/3
function hyp_1bc(b, c, z)
    s, ds, term = one(z), zero(z), one(z)
    tol = eps(z) / 8
    n = 0
    while true
        term *= (b + n) / (c + n) * z
        n += 1
        s += term
        ds += n * term / z
        iszero(term) && break
        # past n ≈ |b| the terms shrink geometrically; the tail is below 3(n + 1)|term|
        n > abs(b) + 2 && 3 * (n + 1) * abs(term) <= tol * (abs(s) + z * abs(ds)) && break
        n > 100_000 && throw(RefinementError("SingularGauss", "the hypergeometric series did not converge"))
    end
    return s, ds
end

# H and H' for non-integer α, u ≤ 2/3
function hilbert_series(α::BigFloat, β::BigFloat, u, v)
    πb = BigFloat(π)
    K = πb * cot(πb * α)
    p = u^α * v^β
    F, dF = hyp_1bc(-α - β, 1 - α, u)
    s = α + β + 1
    # 1/Γ(α + β + 1) vanishes when α + β = −1
    C = isinteger(s) && s <= 0 ? zero(u) :
        SpecialFunctions.gamma(α) * SpecialFunctions.gamma(β + 1) / SpecialFunctions.gamma(s)
    return -K * p + C * F, -K * p * (α / u - β / v) + C * dF
end

# H(0, β, u) and H' for u ≤ 1/3
function hilbert_log(β::BigFloat, u, v)
    r = u / v
    S, dS = zero(u), zero(u)          # Σ C(β,j) rʲ/j and its r-derivative
    c = one(u)
    tol = eps(u) / 8
    j = 0
    while true
        c *= (β - j) / (j + 1) * r
        j += 1
        S += c / j
        dS += c / r
        iszero(c) && break
        # every term of the result is at least of order 1 (log(v/u) ≥ log 2, 1/u ≥ 3), and
        # the tails are below 2|c|
        j > abs(β) + 2 && abs(c) <= tol && break
        j > 100_000 && throw(RefinementError("SingularGauss", "the logarithmic series did not converge"))
    end
    E = log(v / u) - SpecialFunctions.digamma(β + 1) - BigFloat(Base.MathConstants.eulergamma) - S
    vb = v^β
    return vb * E, -β * vb / v * E + vb * (-1 / v - 1 / u - dS / v^2)
end

beta_fn(a, b) = SpecialFunctions.gamma(BigFloat(a)) * SpecialFunctions.gamma(BigFloat(b)) /
                SpecialFunctions.gamma(BigFloat(a) + BigFloat(b))

# H(α, β, u) and H'; the exponents are tested for being integers as given, before rounding
function hilbert01(α::Real, β::Real, u::BigFloat, v::BigFloat)
    if isinteger(α) && α > 0
        m = Int(α)
        h, dh = hilbert01(0, β, u, v)
        B = [beta_fn(i + 1, BigFloat(β) + 1) for i in 0:(m - 1)]
        H = u^m * h + sum(u^(m - 1 - i) * B[i + 1] for i in 0:(m - 1))
        dH = m * u^(m - 1) * h + u^m * dh + sum(((m - 1 - i) * u^(m - 2 - i) * B[i + 1] for i in 0:(m - 2)); init = zero(u))
        return H, dH
    end
    if isinteger(β) && β > 0
        n = Int(β)
        h, dh = hilbert01(α, 0, u, v)
        B = [beta_fn(BigFloat(α) + 1, j + 1) for j in 0:(n - 1)]
        H = v^n * h - sum(v^(n - 1 - j) * B[j + 1] for j in 0:(n - 1))
        dH = -n * v^(n - 1) * h + v^n * dh + sum(((n - 1 - j) * v^(n - 2 - j) * B[j + 1] for j in 0:(n - 2)); init = zero(u))
        return H, dH
    end
    ai, bi = isinteger(α), isinteger(β)              # integer exponents are 0 from here
    ai && bi && return log(v / u), -1 / v - 1 / u
    if u > v
        h, dh = hilbert01(β, α, v, u)
        return -h, dh
    end
    ai || return hilbert_series(BigFloat(α), BigFloat(β), u, v)
    3u <= 1 && return hilbert_log(BigFloat(β), u, v)
    h, dh = hilbert_series(BigFloat(β), BigFloat(α), v, u)   # 1/2 ≤ v < 2/3
    return -h, dh
end

"""
    jacobi_hilbert(α, β, t) -> (ρ₀, σ₀)

The Hilbert transform `ρ₀(t) = ⨍₋₁¹ (1-x)^α (1+x)^β / (x - t) dx` of a Jacobi weight and its
derivative `σ₀(t) = ⨎₋₁¹ (1-x)^α (1+x)^β / (x - t)² dx`, for `-1 < t < 1`, at the ambient
precision.
"""
function jacobi_hilbert(α::Real, β::Real, t::BigFloat)
    u, v = (1 + t) / 2, (1 - t) / 2
    h, dh = hilbert01(β, α, u, v)                      # s = (1 + x)/2 carries (1 + x)^β
    c = BigFloat(2)^(BigFloat(α) + BigFloat(β))
    return c * h, c * dh / 2
end

"""
    kernel_hilbert(k::SingularKernel, a, b, bits) -> (ρ₀, σ₀)

[`jacobi_hilbert`](@ref) at the singular point of `k` mapped from `[a, b]` to `[-1, 1]`,
correct to `bits`: computed at two precisions, with more guard bits until the two agree. An
exponent close to an integer makes the two terms of the series formula large and nearly
cancelling, which the guard bits alone would not reveal.
"""
function kernel_hilbert(k::SingularKernel, a, b, bits::Integer)
    μ = with_bits(() -> BigFloat(jacobi_mass(k.α, k.β)), bits + 16)
    prev = nothing
    for extra in (16, 48, 112, 240, 496)
        cur = with_bits(bits + extra) do
            h, c = (BigFloat(b) - BigFloat(a)) / 2, (BigFloat(b) + BigFloat(a)) / 2
            jacobi_hilbert(k.α, k.β, (BigFloat(k.t) - c) / h)
        end
        if prev !== nothing
            tol = ldexp(BigFloat(1), -bits)
            abs(cur[1] - prev[1]) <= tol * (abs(cur[1]) + μ) &&
                abs(cur[2] - prev[2]) <= tol * (abs(cur[2]) + μ) && return cur
        end
        prev = cur
    end
    throw(RefinementError("SingularGauss", "the Hilbert transform of the weight at t = $(k.t) did not settle to $bits bits"))
end

# h^(α + β) for the principal value, h^(α + β − 1) for the finite part: the factor a rule on
# [-1, 1] picks up on an interval of half-length h
kernel_scale(k::SingularKernel, h::BigFloat) = h^(BigFloat(k.α) + BigFloat(k.β) + 1 - k.order)

measure(d::KernelDomain) = kernel_moments(d, 0)[1]

"""
    kernel_moments(dom, K) -> Vector{BigFloat}

The functional of `dom` applied to `q_k((2x - a - b)/(b - a))` for `k = 0 … K`, at the
ambient precision, where `q_k` are the orthonormal Jacobi polynomials on `[-1, 1]` scaled to
`q₀ = 1`. These are the integrals a rule is verified against. Monic polynomials would do as
well in exact arithmetic, but on `[0, 1]` the monic one of degree 200 is of size 4⁻²⁰⁰, below
any sensible absolute tolerance.

They are computed without the recurrence the construction uses:

    ⨍ w q_k / (x - t) = ∫ w (q_k(x) - q_k(t)) / (x - t) + q_k(t) ρ₀(t),

where the divided difference is a polynomial, integrated by a Gauss–Jacobi rule with enough
points; the finite part subtracts one more Taylor term. Only `ρ₀` and `σ₀` are shared with
the construction.
"""
function kernel_moments(dom::KernelDomain, K::Integer)
    check_kernel(dom)
    k = dom.weight
    bits = precision(BigFloat)
    ρ0, σ0 = kernel_hilbert(k, dom.base.a, dom.base.b, bits)
    ρ0, σ0 = BigFloat(ρ0), BigFloat(σ0)               # back to the ambient precision
    h, c = (BigFloat(dom.base.b) - BigFloat(dom.base.a)) / 2, (BigFloat(dom.base.b) + BigFloat(dom.base.a)) / 2
    t = (BigFloat(k.t) - c) / h
    rec = jacobi_recurrence(k.α, k.β)
    a = [BigFloat(rec.a(j, BigFloat)) for j in 0:K]
    b = [j == 0 ? zero(t) : BigFloat(rec.b(j, BigFloat)) for j in 0:(K + 1)]   # b[j + 1] = b_j
    m = cld(K, 2) + 1                                   # exact to degree 2m − 1 ≥ K
    x, λ, _ = gauss_from_recurrence(m, rec, bits + 16)
    # b_{j+1} D_{j+1} = q_j(x) + (t − a_j) D_j − b_j D_{j−1} for the divided difference
    # D_j = (q_j(x) − q_j(t))/(x − t), and the same with D for q for the second one
    I1, I2 = zeros(BigFloat, K + 1), zeros(BigFloat, K + 1)
    for (xj, λj) in zip(x, λ)
        xj, λj = BigFloat(xj), BigFloat(λj)
        p0, p1 = zero(xj), one(xj)                      # q_{j−1}(x), q_j(x)
        d0, d1 = zero(xj), zero(xj)                     # first divided differences at t
        e0, e1 = zero(xj), zero(xj)                     # second divided differences
        for j in 0:K
            I1[j + 1] += λj * d1
            I2[j + 1] += λj * e1
            j == K && break
            p2 = ((xj - a[j + 1]) * p1 - b[j + 1] * p0) / b[j + 2]
            d2 = (p1 + (t - a[j + 1]) * d1 - b[j + 1] * d0) / b[j + 2]
            e2 = (d1 + (t - a[j + 1]) * e1 - b[j + 1] * e0) / b[j + 2]
            p0, p1, d0, d1, e0, e1 = p1, p2, d1, d2, e1, e2
        end
    end
    out = zeros(BigFloat, K + 1)
    s = kernel_scale(k, h)
    P0, P1, D0, D1 = zero(t), one(t), zero(t), zero(t)  # q_j(t) and q_j'(t)
    for j in 0:K
        out[j + 1] = s * (k.order == 1 ? I1[j + 1] + P1 * ρ0 : I2[j + 1] + P1 * σ0 + D1 * ρ0)
        j == K && break
        P2 = ((t - a[j + 1]) * P1 - b[j + 1] * P0) / b[j + 2]
        D2 = (P1 + (t - a[j + 1]) * D1 - b[j + 1] * D0) / b[j + 2]
        P0, P1, D0, D1 = P1, P2, D1, D2
    end
    return out
end
