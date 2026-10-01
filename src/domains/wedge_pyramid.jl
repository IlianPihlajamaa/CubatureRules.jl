# Wedges (triangular prisms) and pyramids (PLAN §6, v0.6): the remaining finite-element
# shapes in three dimensions.
#
# Both are given by their vertices and must be affine images of their reference shapes. A
# rule's polynomial degree survives an affine map and no other map, so that is what keeps
# `map_to`, mesh integration and verification exact: a wedge whose top is not a translate
# of its bottom, or a pyramid whose base is not a parallelogram, is refused rather than
# integrated with a claim that no longer holds. `transform` takes those, with `NoClaim`.

const REF_WEDGE_VERTICES = SVector{6,SVector{3,Int}}((0, 0, -1), (1, 0, -1), (0, 1, -1), (0, 0, 1), (1, 0, 1), (0, 1, 1))
const REF_PYRAMID_VERTICES = SVector{5,SVector{3,Int}}((-1, -1, 0), (1, -1, 0), (1, 1, 0), (-1, 1, 0), (0, 0, 1))

# Two edge vectors that must agree for the shape to be affine: exactly in exact arithmetic,
# to rounding in floating point, where the vertices of a mesh carry it.
function _same_edge(a::SVector{3,T}, b::SVector{3,T}, scale) where {T}
    T <: AbstractFloat || return a == b
    return maximum(abs, a - b) <= 64 * eps(T) * max(scale, one(T))
end
_vertex_scale(vs) = maximum(v -> maximum(abs, v), vs)

# the type affine maps of a shape with coordinates of type T are computed in
_frametype(::Type{T}) where {T} = T <: Integer ? Rational{BigInt} : T <: Rational ? Rational{BigInt} : T

"""
    Wedge()
    Wedge(v₁, v₂, v₃, v₄, v₅, v₆)

A wedge, or triangular prism: the triangle `v₁ v₂ v₃` swept to `v₄ v₅ v₆`. `Wedge()` is the
reference wedge, the reference triangle `Simplex{2}()` times `[-1, 1]`, with vertices
`(0,0,-1)`, `(1,0,-1)`, `(0,1,-1)`, `(0,0,1)`, `(1,0,1)`, `(0,1,1)` and volume 1.

The top must be a translate of the bottom, `v₄ - v₁ = v₅ - v₂ = v₆ - v₃`, as for the elements
of an extruded mesh, so that the wedge is an affine image of the reference one. A rule keeps
its polynomial degree under an affine map and under no other, so a wedge with a twisted or
tapered top is refused; [`transform`](@ref) integrates over it with `NoClaim`.
"""
struct Wedge{T} <: Domain{3,T}
    vertices::SVector{6,SVector{3,T}}
    function Wedge{T}(v::SVector{6,SVector{3,T}}) where {T}
        e, s = v[4] - v[1], _vertex_scale(v)
        (_same_edge(v[5] - v[2], e, s) && _same_edge(v[6] - v[3], e, s)) ||
            throw(ArgumentError("the top of a wedge must be a translate of its bottom (v₄ - v₁ = v₅ - v₂ = v₆ - v₃) " *
                                "for it to be an affine image of the reference wedge; use `transform` for other shapes"))
        w = new{T}(v)
        iszero(det(affine_frame(w, _frametype(T))[1])) && throw(ArgumentError("degenerate wedge: zero volume"))
        return w
    end
end
Wedge() = Wedge{Int}(REF_WEDGE_VERTICES)
function Wedge(vs::Union{Tuple,AbstractVector}...)
    length(vs) == 6 || throw(ArgumentError("a wedge has 6 vertices, got $(length(vs))"))
    all(v -> length(v) == 3, vs) || throw(ArgumentError("wedge vertices need 3 coordinates each"))
    T = promote_type((eltype(v) for v in vs)...)
    return Wedge{T}(SVector{6,SVector{3,T}}(ntuple(i -> SVector{3,T}(Tuple(vs[i])), 6)))
end

"""
    Pyramid()
    Pyramid(v₁, v₂, v₃, v₄, v₅)

A pyramid with the quadrilateral base `v₁ v₂ v₃ v₄`, in order around it, and apex `v₅`.
`Pyramid()` is the reference pyramid, with base `[-1, 1]²` at `z = 0` and apex `(0, 0, 1)`,
volume 4/3.

The base must be a parallelogram, `v₁ + v₃ = v₂ + v₄`, for the pyramid to be an affine image
of the reference one; the apex can be anywhere off the plane of the base. As for
[`Wedge`](@ref), other shapes are refused and left to [`transform`](@ref).
"""
struct Pyramid{T} <: Domain{3,T}
    vertices::SVector{5,SVector{3,T}}
    function Pyramid{T}(v::SVector{5,SVector{3,T}}) where {T}
        _same_edge(v[2] - v[1], v[3] - v[4], _vertex_scale(v)) ||
            throw(ArgumentError("the base of a pyramid must be a parallelogram (v₁ + v₃ = v₂ + v₄) for it to be " *
                                "an affine image of the reference pyramid; use `transform` for other shapes"))
        p = new{T}(v)
        iszero(det(affine_frame(p, _frametype(T))[1])) && throw(ArgumentError("degenerate pyramid: zero volume"))
        return p
    end
end
Pyramid() = Pyramid{Int}(REF_PYRAMID_VERTICES)
function Pyramid(vs::Union{Tuple,AbstractVector}...)
    length(vs) == 5 || throw(ArgumentError("a pyramid has 5 vertices, got $(length(vs))"))
    all(v -> length(v) == 3, vs) || throw(ArgumentError("pyramid vertices need 3 coordinates each"))
    T = promote_type((eltype(v) for v in vs)...)
    return Pyramid{T}(SVector{5,SVector{3,T}}(ntuple(i -> SVector{3,T}(Tuple(vs[i])), 5)))
end

const WedgeOrPyramid = Union{Wedge,Pyramid}

vertices(d::WedgeOrPyramid) = d.vertices
isreference(d::Wedge) = d.vertices == REF_WEDGE_VERTICES
isreference(d::Pyramid) = d.vertices == REF_PYRAMID_VERTICES
reference(::Wedge) = Wedge()
reference(::Pyramid) = Pyramid()
convert_domain(::Type{S}, d::Wedge) where {S} = Wedge{S}(SVector{6,SVector{3,S}}(map(v -> SVector{3,S}(map(S, v)), d.vertices)))
convert_domain(::Type{S}, d::Pyramid) where {S} = Pyramid{S}(SVector{5,SVector{3,S}}(map(v -> SVector{3,S}(map(S, v)), d.vertices)))

"""
    affine_frame(dom, S) -> (A, b)

The affine map `x ↦ A x + b` from the reference shape onto a [`Wedge`](@ref) or
[`Pyramid`](@ref), in type `S`.
"""
function affine_frame(d::Wedge, ::Type{S}) where {S}
    v = map(p -> SVector{3,S}(map(S, p)), d.vertices)
    return hcat(v[2] - v[1], v[3] - v[1], (v[4] - v[1]) / 2), (v[1] + v[4]) / 2
end
function affine_frame(d::Pyramid, ::Type{S}) where {S}
    v = map(p -> SVector{3,S}(map(S, p)), d.vertices)
    c = (v[1] + v[3]) / 2
    return hcat((v[2] - v[1]) / 2, (v[4] - v[1]) / 2, v[5] - c), c
end

_reference_measure(::Wedge) = 1 // 1
_reference_measure(::Pyramid) = 4 // 3
function measure(d::WedgeOrPyramid)
    isreference(d) && return _reference_measure(d)
    S = _frametype(eltype(eltype(d.vertices)))
    return _reference_measure(d) * abs(det(affine_frame(d, S)[1]))
end

# reference coordinates of a point, in the arithmetic of the point and the vertices
function _reference_point(d::WedgeOrPyramid, x)
    S = _frametype(promote_type(eltype(x), eltype(eltype(d.vertices))))
    ξ = SVector{3,S}(map(S, Tuple(x)))
    isreference(d) && return ξ
    A, b = affine_frame(d, S)
    return A \ (ξ - b)
end

function indomain(x, d::Wedge; tol = 0)
    ξ = _reference_point(d, x)
    return ξ[1] >= -tol && ξ[2] >= -tol && 1 - ξ[1] - ξ[2] >= -tol && -1 - tol <= ξ[3] <= 1 + tol
end
function isinterior(x, d::Wedge; tol = 0)
    ξ = _reference_point(d, x)
    return ξ[1] > tol && ξ[2] > tol && 1 - ξ[1] - ξ[2] > tol && -1 + tol < ξ[3] < 1 - tol
end
function indomain(x, d::Pyramid; tol = 0)
    ξ = _reference_point(d, x)
    return -tol <= ξ[3] <= 1 + tol && abs(ξ[1]) <= 1 - ξ[3] + tol && abs(ξ[2]) <= 1 - ξ[3] + tol
end
function isinterior(x, d::Pyramid; tol = 0)
    ξ = _reference_point(d, x)
    return ξ[3] > tol && abs(ξ[1]) < 1 - ξ[3] - tol && abs(ξ[2]) < 1 - ξ[3] - tol
end

function Base.show(io::IO, d::WedgeOrPyramid)
    name = d isa Wedge ? "Wedge" : "Pyramid"
    isreference(d) ? print(io, name, "()") : print(io, name, "(", join((Tuple(v) for v in d.vertices), ", "), ")")
end
