# Balls (PLAN §2.1, v0.4): the solid counterpart of `Sphere{D}`, and the disk when D = 2.
#
# Everything about a ball follows from the sphere by one radial integral. Its volume is the
# sphere's area over D, and its monomial moments are the sphere's over |α| + D:
#
#     ∫_{B^D} x^α dx = ∫₀¹ r^{|α|+D-1} dr ∫_{S^{D-1}} ω^α dσ(ω) = (1/(|α|+D)) ∫_{S^{D-1}} ω^α dσ
#
# so the moments the sphere already computes serve here unchanged, and the same separation
# gives the rules: a radial rule against a rule on the sphere.

"""
    Ball{D}()
    Ball(centre, radius)

The solid ball of dimension `D`: `Ball{3}()` is the unit ball in three-space, `Ball{2}()`
the unit disk. `Ball{D}()` is the reference ball, centred at the origin with radius 1, and
its measure is the volume `π^{D/2} / Γ(D/2 + 1) · r^D`.
"""
struct Ball{D,T} <: Domain{D,T}
    centre::SVector{D,T}
    radius::T
    function Ball{D,T}(centre, radius) where {D,T}
        radius > 0 || throw(ArgumentError("a ball needs a positive radius, got $radius"))
        return new{D,T}(SVector{D,T}(centre), T(radius))
    end
end

function Ball(centre, radius)
    D = length(centre)
    T = promote_type(eltype(centre), typeof(radius))
    return Ball{D,T}(SVector{D,T}(centre), T(radius))
end
Ball{D}() where {D} = Ball{D,Int}(SVector{D,Int}(ntuple(_ -> 0, D)), 1)

"`Disk()` — the unit disk, `Ball{2}()`."
Disk() = Ball{2}()

isreference(d::Ball) = all(iszero, d.centre) && isone(d.radius)
reference(::Ball{D}) where {D} = Ball{D}()
convert_domain(::Type{S}, d::Ball{D}) where {S,D} =
    Ball{D,S}(SVector{D,S}(map(S, d.centre)), S(d.radius))

"Volume of the unit `D`-ball, `π^{D/2} / Γ(D/2 + 1)` — the sphere's area over `D`."
unit_ball_volume(D::Integer) = unit_sphere_area(D) / D
measure(d::Ball{D}) where {D} = unit_ball_volume(D) * big(d.radius)^D

indomain(x, d::Ball; tol = nothing) =
    sum(abs2, x .- d.centre) <= (d.radius * (1 + (tol === nothing ? _ball_tol(x) : tol)))^2
isinterior(x, d::Ball; tol = nothing) =
    sum(abs2, x .- d.centre) < (d.radius * (1 - (tol === nothing ? _ball_tol(x) : tol)))^2
_ball_tol(x) = 64 * maximum(_coord_eps, x)

function Base.show(io::IO, d::Ball{D}) where {D}
    isreference(d) ? print(io, "Ball{", D, "}()") : print(io, "Ball(", Tuple(d.centre), ", ", d.radius, ")")
end

"""
    ball_moment(D, α)

`∫_{B^D} x^α dx` over the unit ball: the sphere's moment divided by `|α| + D`, and zero
unless every exponent is even.
"""
ball_moment(D::Integer, α) = sphere_moment(D, α) / (sum(α) + D)

monomial_moment(d::Ball{D}, α::NTuple{D,<:Integer}) where {D} =
    isreference(d) ? ball_moment(D, α) :
    throw(ArgumentError("moments are defined on the reference ball; map the rule instead"))
