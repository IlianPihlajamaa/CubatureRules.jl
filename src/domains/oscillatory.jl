# Oscillatory weights (PLAN §6 Tier 2, v0.5): e^{iωx}, cos ωx and sin ωx on an interval.
#
# On [-1, 1], with κ = ωh, the Legendre polynomials have the moments
#
#     ∫ P_k(x) e^{iκx} dx = 2 iᵏ j_k(κ),
#
# with j_k the spherical Bessel functions, and an interval c ± h only adds the factor
# h e^{iωc}. The construction computes j_k upwards from the closed forms of j₀ and j₁. Above
# k ≈ κ that loses digits, since j_k is the recessive solution there, so it runs at a
# precision raised by the loss. Verification computes them downwards by Miller's algorithm,
# which is stable for the recessive solution. The two share nothing but sin and cos.

"""
    OscillatoryWeight

The weight `e^{iωx}`, `cos(ωx)` or `sin(ωx)` on an interval. Made by
[`Oscillatory`](@ref).
"""
struct OscillatoryWeight
    ω::Real
    kind::Symbol                     # :exp, :cos or :sin
    function OscillatoryWeight(ω::Real, kind::Symbol)
        (isfinite(ω) && !iszero(ω)) ||
            throw(ArgumentError("the frequency must be finite and nonzero, got ω = $ω; for ω = 0 use the interval itself"))
        kind in (:exp, :cos, :sin) || throw(ArgumentError("kind must be :exp, :cos or :sin, got :$kind"))
        return new(ω, kind)
    end
end

"""
    Oscillatory(ω)
    Oscillatory(ω, cos)
    Oscillatory(ω, sin)

The weight `e^{iωx}`, `cos(ωx)` or `sin(ωx)`, as the weight of a [`WeightedDomain`](@ref) on
an `Interval`:

```julia
dom = WeightedDomain(Interval(0, 1), Oscillatory(200))
r = rule(dom; degree = 20, digits = 30)
integrate(x -> 1 / (1 + x), r)            # ∫₀¹ e^{200ix} / (1 + x) dx
```

Rules come from [`Filon`](@ref). For `e^{iωx}` the nodes are real and the weights complex;
for `cos` and `sin` both are real. `ω` is taken exactly as given, and must be nonzero.
"""
Oscillatory(ω::Real) = OscillatoryWeight(ω, :exp)
Oscillatory(ω::Real, ::typeof(cos)) = OscillatoryWeight(ω, :cos)
Oscillatory(ω::Real, ::typeof(sin)) = OscillatoryWeight(ω, :sin)

Base.show(io::IO, w::OscillatoryWeight) =
    print(io, "Oscillatory(", w.ω, w.kind === :exp ? "" : ", " * string(w.kind), ")")

"""
    OscillatoryDomain

An [`Interval`](@ref) carrying an [`Oscillatory`](@ref) weight.
"""
const OscillatoryDomain = WeightedDomain{1,<:Any,<:Interval,OscillatoryWeight}

# The weight is stated in the interval's own coordinates (e^{iωx} on [a, b] is not e^{iωx}
# on [-1, 1]), so the rule is built where it is asked for.
isreference(::OscillatoryDomain) = true
reference(d::OscillatoryDomain) = d

function check_oscillatory(dom::OscillatoryDomain)
    a, b = dom.base.a, dom.base.b
    (isfinite(a) && isfinite(b)) || throw(ArgumentError("$(dom.weight) needs a finite interval"))
    return nothing
end

# the extra bits the upward recurrence loses up to index N, from the growth of the dominant
# solution y_k over the recessive j_k above k ≈ |κ|, plus the cancellation in the closed form
# of j₁ when |κ| < 1
function bessel_up_loss(N::Integer, κ::BigFloat)
    lk = Float64(log2(abs(κ)))
    return 2max(0.0, -lk) + sum((max(0.0, 2(log2(2i + 1) - lk)) for i in 1:N); init = 0.0)
end

"""
    spherical_bessel_up(N, κ) -> Vector{BigFloat}

`j₀(κ) … j_N(κ)` at the ambient precision, by the upward recurrence
`j_{k+1} = (2k+1)/κ j_k − j_{k−1}` from the closed forms of `j₀` and `j₁`, run at a precision
raised by what the recurrence loses above `k ≈ |κ|`.
"""
function spherical_bessel_up(N::Integer, κ::BigFloat)
    bits = precision(BigFloat)
    j = with_bits(bits + 32 + ceil(Int, bessel_up_loss(N, κ))) do
        k = BigFloat(κ)
        s, c = sin(k), cos(k)
        out = Vector{BigFloat}(undef, N + 1)
        out[1] = s / k
        N >= 1 && (out[2] = s / k^2 - c / k)
        for n in 1:(N - 1)
            out[n + 2] = (2n + 1) / k * out[n + 1] - out[n]
        end
        out
    end
    return BigFloat.(j)
end

"""
    spherical_bessel_down(N, κ) -> Vector{BigFloat}

`j₀(κ) … j_N(κ)` at the ambient precision by Miller's algorithm: the recurrence run downwards
from an index far enough above `max(N, |κ|)`, normalised by the closed form of `j₀` or `j₁`,
whichever is larger. The starting index uses the decay rate `acosh(k/|κ|)` per step, which
holds near the turning point `k ≈ |κ|` as well as far above it.
"""
function spherical_bessel_down(N::Integer, κ::BigFloat)
    bits = precision(BigFloat)
    j = with_bits(bits + 32) do
        k = abs(BigFloat(κ))
        kf = max(Float64(k), floatmin(Float64))
        K, s = max(N, ceil(Int, kf)), 0.0
        while s < bits + 64
            K += 1
            s += acosh(max(1.0, (K + 0.5) / kf)) / log(2)
        end
        p = zeros(BigFloat, K + 2)                     # p[n + 1] ∝ j_n
        p[K + 1] = one(BigFloat)
        for n in K:-1:1
            p[n] = (2n + 1) / k * p[n + 1] - p[n + 2]
        end
        j0, j1 = sin(k) / k, sin(k) / k^2 - cos(k) / k
        scale = abs(j0) >= abs(j1) ? j0 / p[1] : j1 / p[2]
        [(κ < 0 && isodd(n) ? -scale : scale) * p[n + 1] for n in 0:N]
    end
    return BigFloat.(j)
end

# h, c and κ = ωh for the interval c ± h, at the ambient precision
function oscillatory_frame(dom::OscillatoryDomain)
    A, B = BigFloat(dom.base.a), BigFloat(dom.base.b)
    h, c = (B - A) / 2, (B + A) / 2
    return h, c, BigFloat(dom.weight.ω) * h
end

# what the weight makes of the complex moment h e^{iωc} 2 iᵏ j_k(κ)
oscillatory_part(kind::Symbol, z) = kind === :exp ? z : kind === :cos ? real(z) : imag(z)

"""
    oscillatory_moments(dom, K) -> Vector

`∫ P_k((2x - a - b)/(b - a)) w(x) dx` over `[a, b]` for `k = 0 … K`, with `P_k` the Legendre
polynomials and `w` the weight of `dom`, at the ambient precision: complex for `e^{iωx}`,
real for `cos` and `sin`. Computed with [`spherical_bessel_down`](@ref), not the recurrence
the construction uses.
"""
function oscillatory_moments(dom::OscillatoryDomain, K::Integer)
    check_oscillatory(dom)
    bits = precision(BigFloat)
    m = with_bits(bits + 32) do
        h, c, κ = oscillatory_frame(dom)
        j = spherical_bessel_down(K, κ)
        phase = h * cis(BigFloat(dom.weight.ω) * c)
        [oscillatory_part(dom.weight.kind, phase * 2 * im^k * j[k + 1]) for k in 0:K]
    end
    return [x isa Complex ? Complex{BigFloat}(x) : BigFloat(x) for x in m]
end

# ∫ w over the interval: h e^{iωc} 2 sin(κ)/κ, without the cancellation of the textbook
# (e^{iωb} − e^{iωa})/(iω) when ω(b − a) is small
function measure(d::OscillatoryDomain)
    check_oscillatory(d)
    m = with_bits(precision(BigFloat) + 32) do
        h, c, κ = oscillatory_frame(d)
        oscillatory_part(d.weight.kind, h * cis(BigFloat(d.weight.ω) * c) * 2 * sin(κ) / κ)
    end
    return m isa Complex ? Complex{BigFloat}(m) : BigFloat(m)
end
