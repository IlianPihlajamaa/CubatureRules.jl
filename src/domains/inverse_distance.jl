# The weakly singular kernel of boundary-element methods on a triangle (PLAN §0.3, §6 Tier 5;
# notes/singular-bem.md):
#
#     ∫_T f(y) / |y − x₀| dy,    x₀ in the closed triangle T,
#
# on a triangle in the plane (`Simplex{2}`) or in space (`SurfaceTriangle`). The kernel is the
# weight of a `WeightedDomain`, so a rule for it is an ordinary rule whose weights carry
# 1/|y − x₀|, exact on polynomials f up to its degree.

"""
    InverseDistance(x₀)

The weight `1 / |y − x₀|` on a triangle, the weakly singular kernel of boundary-element
methods, as the weight of a [`WeightedDomain`](@ref) on a triangle in the plane (a `Simplex`
with two-dimensional vertices) or in space (a [`SurfaceTriangle`](@ref)). `x₀` must lie in
the closed triangle: at a vertex, on an edge or inside.

```julia
T = Simplex((0, 0), (1, 0), (0, 1))
dom = WeightedDomain(T, InverseDistance((1//3, 1//3)))
r = rule(dom; degree = 10)
integrate(y -> 1 + y[1]^2, r)            # ∫_T (1 + y₁²) / |y − x₀| dy
```

Rules come from [`DuffyGauss`](@ref). Their weights carry the kernel, so `integrate(f, r)`
needs only the smooth factor `f`, and the rule is exact when `f` is a polynomial of degree up
to the rule's degree.

Exact coordinates (integers and rationals) are taken as given. Floating-point ones carry
rounding, and a point meant to lie on the triangle rarely does exactly: the centroid of a
triangle in space, computed in `Float64`, is off its plane by about an ulp. So `x₀` within 64
units in the last place (of the largest coordinate) of the triangle's plane, or of its
boundary, is moved onto it, and the rule is built for that point; the provenance says how
far it moved. A point further off is refused: integrals near, but not at, the singularity
are not covered yet.
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

A triangle in the plane (`Simplex{2}`) or in space ([`SurfaceTriangle`](@ref)) carrying an
[`InverseDistance`](@ref) kernel.
"""
const InverseDistanceDomain = Union{WeightedDomain{2,<:Any,<:Simplex{2},<:InverseDistance{2}},
                                    WeightedDomain{3,<:Any,<:SurfaceTriangle,<:InverseDistance{3}}}

# The kernel is not invariant under affine maps (only under similarities), so the rule is
# built where it is asked for: such a domain is its own reference.
isreference(::InverseDistanceDomain) = true
reference(d::InverseDistanceDomain) = d

_exact(x::Integer) = Rational{BigInt}(x)
_exact(x::Rational) = Rational{BigInt}(x)
_exact(x::AbstractFloat) = Rational{BigInt}(x)      # exact: a float is a binary fraction
_exact(x::Real) = Rational{BigInt}(BigFloat(x))
_ulp(x::Real) = Rational{BigInt}(value_ulp(x))

# twice the signed area of (o, a, b) along the triangle's orientation `n` (in space), or in
# the plane; and twice the unsigned area of the parallelogram on a and b
_signed_area2(a::SVector{2}, b::SVector{2}, _) = a[1] * b[2] - a[2] * b[1]
_signed_area2(a::SVector{3}, b::SVector{3}, n) = dot(n, cross(a, b))
_area2(a::SVector{2}, b::SVector{2}) = abs(a[1] * b[2] - a[2] * b[1])
_area2(a::SVector{3}, b::SVector{3}) = norm(cross(a, b))

"""
    singular_point(dom) -> (x₀, patches, place, moved)

Where the rule is built: the singular point `x₀` as exact rationals, the sub-triangles
`(x₀, p, q)` of positive area into which it cuts the triangle (as pairs `(p, q)` in the
triangle's orientation), where it lies (`:vertex`, `:edge` or `:interior`), and how far it
was moved onto the triangle's plane or boundary to absorb rounding (0 when it was not; see
[`InverseDistance`](@ref)). Decided in exact arithmetic. Throws if the triangle is
degenerate or `x₀` is not on it.
"""
function singular_point(dom::InverseDistanceDomain)
    D = length(dom.weight.x0)
    V = [SVector{D,Rational{BigInt}}(map(_exact, p)) for p in vertices(dom.base)]
    x = SVector{D,Rational{BigInt}}(map(_exact, dom.weight.x0))
    given = x
    tol = 64 * maximum(_ulp, vcat(collect(dom.weight.x0), (collect(p) for p in vertices(dom.base))...))
    e1, e2 = V[2] - V[1], V[3] - V[1]
    n = D == 3 ? cross(e1, e2) : nothing
    total = _signed_area2(e1, e2, n)
    iszero(total) && throw(ArgumentError("the triangle $(dom.base) is degenerate"))
    if D == 3                                       # onto the plane, if within rounding of it
        s = dot(n, x - V[1])
        nn = dot(n, n)
        s^2 <= tol^2 * nn || throw(ArgumentError(
            "the singular point $(Tuple(dom.weight.x0)) is not in the plane of $(dom.base) (it is " *
            "$(Float64(abs(s)) / sqrt(Float64(nn))) away): integrals near, but not at, the singularity are " *
            "not available yet"))
        x -= (s / nn) * n
    end
    edges = ((2, 3), (3, 1), (1, 2))
    λ(x) = [_signed_area2(V[j] - x, V[k] - x, n) / total for (j, k) in edges]
    if tol > 0 || any(<(0), λ(x))
        # the nearest point of the boundary
        c, δ2 = nothing, nothing
        for (j, k) in edges
            d = V[k] - V[j]
            τ = clamp(dot(x - V[j], d) / dot(d, d), 0, 1)
            p = V[j] + τ * d
            q = dot(x - p, x - p)
            (δ2 === nothing || q < δ2) && ((c, δ2) = (p, q))
        end
        if δ2 <= tol^2
            x = c
            for v in V                               # and onto a vertex, if as close to one
                dot(x - v, x - v) <= tol^2 && (x = v)
            end
        elseif any(<(0), λ(x))
            throw(ArgumentError(
                "the singular point $(Tuple(dom.weight.x0)) lies outside the triangle $(dom.base): the integral " *
                "is not singular there and an ordinary rule applies; integrals near, but not at, the " *
                "singularity are not available yet"))
        end
    end
    patches = [(V[j], V[k]) for ((j, k), l) in zip(edges, λ(x)) if l > 0]
    place = length(patches) == 1 ? :vertex : length(patches) == 2 ? :edge : :interior
    return x, patches, place, sqrt(Float64(dot(x - given, x - given)))
end

"The number of sub-triangles: 1, 2 or 3 as `x₀` is a vertex, on an edge or inside."
function npatches(dom::InverseDistanceDomain)
    try
        return length(singular_point(dom)[2])
    catch err
        err isa ArgumentError || rethrow()
        return 3                                     # refused when the rule is built
    end
end

# On the sub-triangle (x₀, p, q) in polar coordinates about x₀ the kernel cancels the
# Jacobian, ∫∫ f (1/r) r dr dψ, and along the edge pq, at the foot f of the perpendicular from
# x₀ at distance h, the angle is ψ = atan(τ/h), dψ = h dτ/(h² + τ²) and the ray has length
# R = √(h² + τ²). So the patch contributes ∫ dτ h/√(h² + τ²) ∫₀¹ f(x₀ + ρ (e(τ) − x₀)) dρ;
# with f = 1 that is h [asinh(τ/h)] between the ends of the edge.
"The foot, unit direction, distance `h` and edge parameters `τ_p, τ_q` of a patch, in `S`."
function patch_frame(x0, p, q, ::Type{S}) where {S}
    D = length(x0)
    x0, p, q = SVector{D,S}(x0), SVector{D,S}(p), SVector{D,S}(q)
    u = (q - p) / norm(q - p)
    f = p + dot(x0 - p, u) * u
    return f, u, norm(x0 - f), dot(p - f, u), dot(q - f, u)
end

"`∫_T 1/|y − x₀| dy`, in closed form, in the ambient `BigFloat` precision."
function measure(d::InverseDistanceDomain)
    x0, patches, _, _ = singular_point(d)
    return sum(patches) do (p, q)
        _, _, h, τp, τq = patch_frame(x0, p, q, BigFloat)
        h * (asinh(τq / h) - asinh(τp / h))
    end
end
