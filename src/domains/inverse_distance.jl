# The weakly singular kernel of boundary-element methods on a triangle (PLAN §0.3, §6 Tier 5;
# notes/singular-bem.md):
#
#     ∫_T f(y) / |y − x₀| dy,    x₀ in the closed triangle T.
#
# The kernel is the weight of a `WeightedDomain` on the triangle, so a rule for it is an
# ordinary rule whose weights carry 1/|y − x₀|, exact on polynomials f up to its degree.

"""
    InverseDistance(x₀)

The weight `1 / |y − x₀|` on a triangle, the weakly singular kernel of boundary-element
methods, as the weight of a [`WeightedDomain`](@ref) on a `Simplex` in the plane. `x₀` must
lie in the closed triangle: at a vertex, on an edge or inside. It is taken exactly as given:
`0.3` is the binary number nearest 0.3, and `3//10` or `big"0.3"` the decimal.

```julia
T = Simplex((0, 0), (1, 0), (0, 1))
dom = WeightedDomain(T, InverseDistance((1//3, 1//3)))
r = rule(dom; degree = 10)
integrate(y -> 1 + y[1]^2, r)            # ∫_T (1 + y₁²) / |y − x₀| dy
```

Rules come from [`DuffyGauss`](@ref). Their weights carry the kernel, so `integrate(f, r)`
needs only the smooth factor `f`, and the rule is exact when `f` is a polynomial of degree up
to the rule's degree.
"""
struct InverseDistance{D,P<:Real}
    x0::SVector{D,P}
    function InverseDistance{D,P}(x0::SVector{D,P}) where {D,P<:Real}
        all(isfinite, x0) || throw(ArgumentError("the singular point must be finite, got $(Tuple(x0))"))
        return new{D,P}(x0)
    end
end
InverseDistance(x0::SVector{D,P}) where {D,P<:Real} = InverseDistance{D,P}(x0)
InverseDistance(x0::Union{Tuple,AbstractVector}) = InverseDistance(SVector{length(x0)}(promote(x0...)))

Base.show(io::IO, w::InverseDistance) = print(io, "InverseDistance(", Tuple(w.x0), ")")

"""
    InverseDistanceDomain

A triangle in the plane carrying an [`InverseDistance`](@ref) kernel.
"""
const InverseDistanceDomain = WeightedDomain{2,<:Any,<:Simplex{2},<:InverseDistance{2}}

# The kernel is not invariant under affine maps (only under similarities), so the rule is
# built where it is asked for: such a domain is its own reference.
isreference(::InverseDistanceDomain) = true
reference(d::InverseDistanceDomain) = d

_exact(x::Integer) = Rational{BigInt}(x)
_exact(x::Rational) = Rational{BigInt}(x)
_exact(x::AbstractFloat) = Rational{BigInt}(x)      # exact: a float is a binary fraction
_exact(x::Real) = Rational{BigInt}(BigFloat(x))
_cross(a, b) = a[1] * b[2] - a[2] * b[1]

"""
    singular_patches(dom) -> (patches, place)

The sub-triangles `(x₀, p, q)` of positive area into which `x₀` cuts the triangle, as pairs
`(p, q)` of exact rational vertices in the triangle's orientation, and where `x₀` lies
(`:vertex`, `:edge` or `:interior`). Decided in exact arithmetic, so a point on an edge is
on it and not near it. Throws if the triangle is degenerate or `x₀` is outside it.
"""
function singular_patches(dom::InverseDistanceDomain)
    v = [SVector{2,Rational{BigInt}}(_exact(p[1]), _exact(p[2])) for p in vertices(dom.base)]
    x0 = SVector{2,Rational{BigInt}}(_exact(dom.weight.x0[1]), _exact(dom.weight.x0[2]))
    A = _cross(v[2] - v[1], v[3] - v[1])
    iszero(A) && throw(ArgumentError("the triangle $(dom.base) is degenerate"))
    # λᵢ: the share of the area opposite vertex i, i.e. of the sub-triangle on edge (i+1, i+2)
    edges = ((2, 3), (3, 1), (1, 2))
    λ = [_cross(v[j] - x0, v[k] - x0) / A for (j, k) in edges]
    any(<(0), λ) && throw(ArgumentError(
        "the singular point $(Tuple(dom.weight.x0)) lies outside the triangle $(dom.base): the integral is not " *
        "singular there and an ordinary rule applies; rules for near-singular points are not available yet"))
    patches = [(v[j], v[k]) for ((j, k), l) in zip(edges, λ) if l > 0]
    place = length(patches) == 1 ? :vertex : length(patches) == 2 ? :edge : :interior
    return patches, place
end

"The number of sub-triangles: 1, 2 or 3 as `x₀` is a vertex, on an edge or inside."
function npatches(dom::InverseDistanceDomain)
    v = [SVector{2,Rational{BigInt}}(_exact(p[1]), _exact(p[2])) for p in vertices(dom.base)]
    x0 = SVector{2,Rational{BigInt}}(_exact(dom.weight.x0[1]), _exact(dom.weight.x0[2]))
    return count(((j, k),) -> !iszero(_cross(v[j] - x0, v[k] - x0)), ((2, 3), (3, 1), (1, 2)))
end

# On the sub-triangle (x₀, p, q) in polar coordinates about x₀ the kernel cancels the
# Jacobian, ∫∫ f (1/r) r dr dψ, and along the edge pq, at the foot f of the perpendicular from
# x₀ at distance h, the angle is ψ = atan(τ/h), dψ = h dτ/(h² + τ²) and the ray has length
# R = √(h² + τ²). So the patch contributes ∫ dτ h/√(h² + τ²) ∫₀¹ f(x₀ + ρ (e(τ) − x₀)) dρ;
# with f = 1 that is h [asinh(τ/h)] between the ends of the edge.
"The foot, unit direction, distance `h` and edge parameters `τ_p, τ_q` of a patch, in `S`."
function patch_frame(x0, p, q, ::Type{S}) where {S}
    x0, p, q = SVector{2,S}(x0), SVector{2,S}(p), SVector{2,S}(q)
    u = (q - p) / norm(q - p)
    f = p + dot(x0 - p, u) * u
    return f, u, norm(x0 - f), dot(p - f, u), dot(q - f, u)
end

"`∫_T 1/|y − x₀| dy`, in closed form, in the ambient `BigFloat` precision."
function measure(d::InverseDistanceDomain)
    patches, _ = singular_patches(d)
    x0 = SVector{2,BigFloat}(d.weight.x0)
    return sum(patches) do (p, q)
        _, _, h, τp, τq = patch_frame(x0, p, q, BigFloat)
        h * (asinh(τq / h) - asinh(τp / h))
    end
end
