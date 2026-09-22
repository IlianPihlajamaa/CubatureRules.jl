# Unbounded domains and their classical weights (PLAN §2.1).
#
# These are only useful weighted: nobody integrates against Lebesgue measure on the whole
# line. They therefore exist to be the `base` of a `WeightedDomain`, and the measure that
# matters — the one a rule's weights sum to — is the weight's mass ∫ w, which is finite.

"""
    HalfLine()

The ray `[0, ∞)`. Useful as the base of a [`WeightedDomain`](@ref); see [`LaguerreRay`](@ref).
"""
struct HalfLine{T} <: Domain{1,T} end
HalfLine() = HalfLine{Int}()

"""
    RealLine()

The whole line `(-∞, ∞)`. Useful as the base of a [`WeightedDomain`](@ref); see
[`HermiteLine`](@ref).
"""
struct RealLine{T} <: Domain{1,T} end
RealLine() = RealLine{Int}()

isreference(::HalfLine) = true
isreference(::RealLine) = true
reference(::HalfLine) = HalfLine()
reference(::RealLine) = RealLine()
convert_domain(::Type{S}, ::HalfLine) where {S} = HalfLine{S}()
convert_domain(::Type{S}, ::RealLine) where {S} = RealLine{S}()
measure(::HalfLine) = Inf
measure(::RealLine) = Inf
indomain(x, ::HalfLine; tol = 0) = only(x) >= -tol
indomain(x, ::RealLine; tol = 0) = isfinite(only(x))
isinterior(x, ::HalfLine; tol = 0) = only(x) > tol
isinterior(x, ::RealLine; tol = 0) = isfinite(only(x))
Base.show(io::IO, ::HalfLine) = print(io, "HalfLine()")
Base.show(io::IO, ::RealLine) = print(io, "RealLine()")

"""
    ExponentialWeight(α = 0)

The Laguerre weight `x^α e^{-x}` on the half line, with mass `Γ(α + 1)`.
"""
struct ExponentialWeight{T}
    α::T
    function ExponentialWeight(α::T = 0) where {T}
        α > -1 || throw(ArgumentError("the Laguerre weight needs α > -1, got $α"))
        return new{T}(α)
    end
end
(w::ExponentialWeight)(x) = x^w.α * exp(-x)

"""
    GaussianWeight()

The Hermite weight `e^{-x²}` on the whole line, with mass `√π`.
"""
struct GaussianWeight end
(::GaussianWeight)(x) = exp(-x^2)

Base.show(io::IO, w::ExponentialWeight) = print(io, iszero(w.α) ? "ExponentialWeight()" : "ExponentialWeight($(w.α))")
Base.show(io::IO, ::GaussianWeight) = print(io, "GaussianWeight()")

"""
    LaguerreRay(α = 0)

`[0, ∞)` with the weight `x^α e^{-x}`: the home of Gauss–Laguerre rules.
"""
LaguerreRay(α = 0) = WeightedDomain(HalfLine(), ExponentialWeight(α))

"""
    HermiteLine()

The whole line with the weight `e^{-x²}`: the home of Gauss–Hermite rules.
"""
HermiteLine() = WeightedDomain(RealLine(), GaussianWeight())

const LaguerreDomain = WeightedDomain{1,<:Any,<:HalfLine,<:ExponentialWeight}
const HermiteDomain = WeightedDomain{1,<:Any,<:RealLine,<:GaussianWeight}

# The measure of a weighted unbounded domain is the weight's mass — finite, and what a
# rule's weights sum to.
measure(d::LaguerreDomain) = laguerre_mass(d.weight.α)
measure(::HermiteDomain) = sqrt(big(π))

"`∫₀^∞ x^α e^{-x} dx = Γ(α + 1)`, exact for integer α."
laguerre_mass(α::Integer) = factorial(big(α))
laguerre_mass(α) = SpecialFunctions.gamma(big(α) + 1)
