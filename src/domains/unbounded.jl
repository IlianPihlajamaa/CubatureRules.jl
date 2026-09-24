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
    RealSpace{D}()
    RealLine()

The whole of `R^D`, and `RealLine()` for the line `(-∞, ∞)` it specialises to. Useful as the
base of a [`WeightedDomain`](@ref); see [`HermiteLine`](@ref) and [`GaussianSpace`](@ref).
"""
struct RealSpace{D,T} <: Domain{D,T} end
RealSpace{D}() where {D} = RealSpace{D,Int}()

"The whole line: `RealSpace{1}`."
const RealLine = RealSpace{1}

isreference(::HalfLine) = true
endpoints(::HalfLine) = (0, Inf)
isreference(::RealSpace) = true
endpoints(::RealSpace{1}) = (-Inf, Inf)
reference(::HalfLine) = HalfLine()
reference(::RealSpace{D}) where {D} = RealSpace{D}()
convert_domain(::Type{S}, ::HalfLine) where {S} = HalfLine{S}()
convert_domain(::Type{S}, ::RealSpace{D}) where {S,D} = RealSpace{D,S}()
measure(::HalfLine) = Inf
measure(::RealSpace) = Inf
indomain(x, ::HalfLine; tol = 0) = only(x) >= -tol
indomain(x, ::RealSpace; tol = 0) = all(isfinite, x)
isinterior(x, ::HalfLine; tol = 0) = only(x) > tol
isinterior(x, ::RealSpace; tol = 0) = all(isfinite, x)
Base.show(io::IO, ::HalfLine) = print(io, "HalfLine()")
Base.show(io::IO, ::RealSpace{D}) where {D} = D == 1 ? print(io, "RealLine()") :
                                              print(io, "RealSpace{", D, "}()")

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

The Hermite weight `e^{-|x|²}`, with mass `π^{D/2}` over `R^D` — `√π` on the line. Radial,
so it reads the same in any dimension.
"""
struct GaussianWeight end
(::GaussianWeight)(x) = exp(-sum(abs2, x))

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
HermiteLine() = GaussianSpace(1)

"""
    GaussianSpace(D)

`R^D` with the weight `e^{-|x|²}`: Stroud's `E_n^{r²}`, and the home of the Gaussian product
rules. `GaussianSpace(1)` is [`HermiteLine`](@ref).
"""
GaussianSpace(D::Integer) = WeightedDomain(RealSpace{D}(), GaussianWeight())

const LaguerreDomain = WeightedDomain{1,<:Any,<:HalfLine,<:ExponentialWeight}
const GaussianDomain{D} = WeightedDomain{D,<:Any,<:RealSpace{D},<:GaussianWeight}
"The one-dimensional Gaussian space, where the Gauss–Hermite rules live."
const HermiteDomain = GaussianDomain{1}

# The measure of a weighted unbounded domain is the weight's mass — finite, and what a
# rule's weights sum to.
measure(d::LaguerreDomain) = laguerre_mass(d.weight.α)
measure(::GaussianDomain{D}) where {D} = big(π)^(D // 2)

"""
    gaussian_moment(D, α)

`∫_{R^D} x^α e^{-|x|²} dx`: zero unless every exponent is even, and otherwise the product of
`Γ((αᵢ+1)/2)`, the weight being a product of one-dimensional ones.
"""
function gaussian_moment(D::Integer, α)
    length(α) == D || throw(ArgumentError("expected $D exponents, got $(length(α))"))
    any(isodd, α) && return zero(BigFloat)
    return prod(SpecialFunctions.gamma((big(a) + 1) / 2) for a in α)
end

monomial_moment(d::GaussianDomain{D}, α::NTuple{D,<:Integer}) where {D} =
    isreference(d) ? gaussian_moment(D, α) :
    throw(ArgumentError("moments are defined on the reference space"))

"`∫₀^∞ x^α e^{-x} dx = Γ(α + 1)`, exact for integer α."
laguerre_mass(α::Integer) = factorial(big(α))
laguerre_mass(α) = SpecialFunctions.gamma(big(α) + 1)
