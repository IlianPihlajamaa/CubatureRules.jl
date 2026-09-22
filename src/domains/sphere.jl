# Spheres (PLAN §2.1, v0.4).
#
# `Sphere{D}` is the (D-1)-sphere embedded in R^D: its points carry D coordinates, so
# `Sphere{3}()` is the familiar surface in three-space and `Sphere{2}()` is the circle. The
# type parameter counts coordinates, as it does for every other domain here, rather than
# the dimension of the surface.
#
# A sphere has no boundary, which the rest of the package has to be told: `isinterior` is
# `indomain`, and the registry's `interior = true` filter must not throw away every rule on
# it. Nodes are meant to satisfy |x| = 1 exactly, so both predicates carry a tolerance that
# defaults to the resolution of the coordinate type rather than to zero.

"""
    Sphere{D}()
    Sphere(centre, radius)

The sphere of dimension `D - 1` embedded in `R^D`: `Sphere{3}()` is the unit sphere in
three-space, `Sphere{2}()` the unit circle. `Sphere{D}()` is the reference sphere, centred
at the origin with radius 1.

Its measure is the surface area `2 π^{D/2} / Γ(D/2) · r^{D-1}` — `2π` for the circle, `4π`
for the sphere — and a rule's weights sum to it.
"""
struct Sphere{D,T} <: Domain{D,T}
    centre::SVector{D,T}
    radius::T
    function Sphere{D,T}(centre, radius) where {D,T}
        radius > 0 || throw(ArgumentError("a sphere needs a positive radius, got $radius"))
        return new{D,T}(SVector{D,T}(centre), T(radius))
    end
end

function Sphere(centre, radius)
    D = length(centre)
    T = promote_type(eltype(centre), typeof(radius))
    return Sphere{D,T}(SVector{D,T}(centre), T(radius))
end
Sphere{D}() where {D} = Sphere{D,Int}(SVector{D,Int}(ntuple(_ -> 0, D)), 1)

isreference(d::Sphere) = all(iszero, d.centre) && isone(d.radius)
reference(::Sphere{D}) where {D} = Sphere{D}()
convert_domain(::Type{S}, d::Sphere{D}) where {S,D} =
    Sphere{D,S}(SVector{D,S}(map(S, d.centre)), S(d.radius))

"Surface area of the unit `(D-1)`-sphere, `2 π^{D/2} / Γ(D/2)`."
unit_sphere_area(D::Integer) = 2 * big(π)^(D / 2) / SpecialFunctions.gamma(big(D) / 2)
measure(d::Sphere{D}) where {D} = unit_sphere_area(D) * big(d.radius)^(D - 1)

# A point is on the sphere, not inside it: the default tolerance is the resolution of the
# coordinates, since |x| = 1 is only ever attained to rounding.
#
# The resolution is the coarser of two things, and both matter. A `BigFloat` node carries
# its own precision, which bounds how well |x| = 1 can hold in the stored numbers; and the
# radius is recomputed here in whatever precision is ambient, which bounds how well it can
# be seen to hold. Verification runs at a precision set by the rule's *weights*, which need
# not match its nodes', so taking either alone rejects the nodes of otherwise perfect rules.
_coord_eps(xi::BigFloat) = ldexp(BigFloat(1), -precision(xi))
_coord_eps(xi::T) where {T<:AbstractFloat} = eps(T)
_coord_eps(xi) = 0                              # exact coordinates resolve everything
_sphere_tol(x, r) = 64 * max(maximum(_coord_eps, x), _coord_eps(r))

function indomain(x, d::Sphere; tol = nothing)
    r = sqrt(sum(abs2, x .- d.centre))
    t = tol === nothing ? _sphere_tol(x, r) : tol
    return abs(r - d.radius) <= t * max(d.radius, one(r))
end
# The surface has no boundary, so every point of it counts as interior. Without this the
# registry's `interior = true` filter would reject every rule on a sphere.
isinterior(x, d::Sphere; tol = nothing) = indomain(x, d; tol)

function Base.show(io::IO, d::Sphere{D}) where {D}
    isreference(d) ? print(io, "Sphere{", D, "}()") :
    print(io, "Sphere(", Tuple(d.centre), ", ", d.radius, ")")
end

"""
    sphere_moment(D, α)

`∫_{S^{D-1}} x^α dσ` over the unit sphere: zero unless every exponent is even, and
otherwise `2 ∏ᵢ Γ((αᵢ+1)/2) / Γ((Σα + D)/2)`. Returned as a `BigFloat` at the working
precision — these moments are rational multiples of `π^{⌈D/2⌉}` and so are never exact in
rational arithmetic.
"""
function sphere_moment(D::Integer, α)
    length(α) == D || throw(ArgumentError("expected $D exponents, got $(length(α))"))
    any(isodd, α) && return zero(BigFloat)
    num = sum(SpecialFunctions.loggamma((big(a) + 1) / 2) for a in α)
    den = SpecialFunctions.loggamma((big(sum(α)) + D) / 2)
    return 2 * exp(num - den)
end

monomial_moment(d::Sphere{D}, α::NTuple{D,<:Integer}) where {D} =
    isreference(d) ? sphere_moment(D, α) :
    throw(ArgumentError("moments are defined on the reference sphere; map the rule instead"))
