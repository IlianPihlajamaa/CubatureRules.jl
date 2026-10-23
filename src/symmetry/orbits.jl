# Orbit algebra for the full symmetry group S_{D+1} of a D-simplex.
#
# The group permutes barycentric coordinates, so an orbit is described by the *pattern* of
# equal coordinates: a multiset of multiplicities summing to D+1. On the triangle
#
#     [3]       centroid         (1/3, 1/3, 1/3)             1 point,  weight only
#     [2, 1]    vertex type      (a, a, 1-2a)                3 points, (w, a)
#     [1, 1, 1] general          (a, b, 1-a-b)               6 points, (w, a, b)
#
# and on the tetrahedron [4], [3,1], [2,2], [2,1,1], [1,1,1,1]. Nothing below is specific to
# D = 2: the tetrahedron needs data, not code.
#
# A pattern with r distinct values has r-1 free coordinates (the last is fixed by Σλ = 1)
# and one weight, which is the weight *per point*.

"""
    OrbitPattern(mult, N)

An `S_N`-orbit type on the `(N-1)`-simplex, given by the multiplicities of the distinct
barycentric coordinates, e.g. `OrbitPattern([2, 1], 3)` for the 3-point vertex-type orbit
on the triangle.
"""
struct OrbitPattern
    mult::Vector{Int}
    N::Int
    labels::Vector{Vector{Int}}   # distinct arrangements of value labels over the N slots
end

function OrbitPattern(mult::AbstractVector{<:Integer}, N::Integer)
    all(>=(1), mult) || throw(ArgumentError("multiplicities must be positive"))
    sum(mult) == N || throw(ArgumentError("multiplicities $(mult) do not sum to $N"))
    base = reduce(vcat, [fill(i, m) for (i, m) in enumerate(mult)])
    labels = sort!(unique!([base[p] for p in permutations_of(N)]))
    return OrbitPattern(collect(Int, mult), Int(N), labels)
end

Base.:(==)(a::OrbitPattern, b::OrbitPattern) = a.mult == b.mult && a.N == b.N
Base.hash(p::OrbitPattern, h::UInt) = hash(p.mult, hash(p.N, h))
Base.show(io::IO, p::OrbitPattern) = print(io, "OrbitPattern(", p.mult, ")")

"Number of points in the orbit: N! / Π mᵢ!."
orbit_size(p::OrbitPattern) = length(p.labels)
"Free coordinate parameters (the weight excluded)."
ncoords(p::OrbitPattern) = length(p.mult) - 1
"Unknowns contributed to a moment system: one weight plus the free coordinates."
nunknowns(p::OrbitPattern) = 1 + ncoords(p)

"All permutations of 1:N in lexicographic order (N is at most 5 here)."
function permutations_of(N::Integer)
    N == 1 && return [[1]]
    out = Vector{Vector{Int}}()
    for first in 1:N, rest in permutations_of(N - 1)
        push!(out, vcat(first, [r >= first ? r + 1 : r for r in rest]))
    end
    return out
end

"""
    pattern_values(p, v)

The distinct barycentric values of pattern `p` given its free coordinates `v`; the last is
`(1 - Σ mᵢ vᵢ) / m_r`.
"""
function pattern_values(p::OrbitPattern, v::AbstractVector{S}) where {S}
    r = length(p.mult)
    length(v) == r - 1 || throw(ArgumentError("pattern $(p.mult) has $(r-1) free coordinates"))
    last = (one(S) - sum((p.mult[i] * v[i] for i in 1:(r - 1)); init = zero(S))) / p.mult[r]
    return vcat(v, last)
end

"""
    OrbitStructure

The shape of a symmetric rule before any numbers: a list of orbits in its field `orbits`,
whose parameter blocks `[w, p₁, …]` — the weight of each point, then the orbit's own
parameters — make up the rule's parameter vector `θ` in order. An orbit type defines
`orbit_size` and `nunknowns`; the structure's counts and block offsets follow, here, once for
every domain's structure (the simplex's, the box's and round domains', the octahedral one on
the sphere, the pyramid's and the wedge's).
"""
abstract type OrbitStructure end

npoints(s::OrbitStructure) = sum(orbit_size, s.orbits; init = 0)
nunknowns(s::OrbitStructure) = sum(nunknowns, s.orbits; init = 0)

"""
    margins(s::OrbitStructure, θ) -> NamedTuple

How far the parameters `θ` of the structure `s` are from failing to describe a positive
interior rule of that shape: `weight`, the smallest weight; `interior`, the smallest distance
of a node's coordinate from the edge of its range (a barycentric value from 0, a coordinate
from 0 and 1, a node from the boundary of the disk or ball); and `gap`, the smallest difference
between two of an orbit's values, where an orbit has several (on the sphere, the distance to
an orbit degeneracy, and there is no `interior`). One function for every domain's structure,
read by the searches and the eliminations.
"""
function margins end

"""
    moment_system(s::OrbitStructure, n, S, basis = moment_basis(s, n))

The defining equations of a fully symmetric rule with orbit structure `s` and exactness
degree `n`, in arithmetic `S`. Callable: `sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)`.
The seed tables record its residual, and the certificates quote it.
"""
function moment_system end

"What the degree-`n` moment systems of `s` in every arithmetic share, computed once: nothing, except on the simplex."
moment_basis(::OrbitStructure, n) = nothing

"""
    refinement_system(s::OrbitStructure, n, S, basis = moment_basis(s, n))

The system Gauss–Newton solves for the rule: the [`moment_system`](@ref) itself, except where
another set of equations with the same roots is better conditioned.
"""
refinement_system(s::OrbitStructure, n, ::Type{S}, basis = moment_basis(s, n)) where {S} = moment_system(s, n, S, basis)

"Offsets of each orbit's block in the parameter vector."
function param_offsets(s::OrbitStructure)
    offs = Vector{Int}(undef, length(s.orbits))
    k = 0
    for (i, o) in enumerate(s.orbits)
        offs[i] = k
        k += nunknowns(o)
    end
    return offs
end

"""
    SymmetricStructure(patterns)

A fully symmetric rule's orbit structure: an ordered list of [`OrbitPattern`](@ref)s. The
parameter vector of the rule is the concatenation, per orbit, of `[w, v₁, …, v_{r-1}]`.
"""
struct SymmetricStructure <: OrbitStructure
    N::Int
    orbits::Vector{OrbitPattern}
end
function SymmetricStructure(patterns::AbstractVector, N::Integer)
    return SymmetricStructure(Int(N), [p isa OrbitPattern ? p : OrbitPattern(p, N) for p in patterns])
end

dimension(s::SymmetricStructure) = s.N - 1
Base.show(io::IO, s::SymmetricStructure) =
    print(io, "SymmetricStructure(", join((string(o.mult) for o in s.orbits), ", "), ")")

"""
    expand(s::SymmetricStructure, θ) -> (λs, ws)

Barycentric coordinates (vectors of length N) and weights of every node of the symmetric
rule with parameters `θ`. Exact under the group action: the points of an orbit are
permutations of the *same* computed values.
"""
function expand(s::SymmetricStructure, θ::AbstractVector{S}) where {S}
    λs = Vector{Vector{S}}()
    ws = Vector{S}()
    for (o, off) in zip(s.orbits, param_offsets(s))
        w = θ[off + 1]
        vals = pattern_values(o, θ[(off + 2):(off + nunknowns(o))])
        for lab in o.labels
            push!(λs, [vals[l] for l in lab])
            push!(ws, w)
        end
    end
    return λs, ws
end

"""
    orbit_pattern_types(N)

Every orbit pattern of `S_N` on the `(N-1)`-simplex (the integer partitions of `N`),
in a fixed order: larger parts first.
"""
function orbit_pattern_types(N::Integer)
    # integer partitions of N, largest parts first, as multiplicity vectors
    parts = Vector{Vector{Int}}()
    function rec(rem, maxpart, acc)
        rem == 0 && (push!(parts, copy(acc)); return)
        for k in min(rem, maxpart):-1:1
            push!(acc, k); rec(rem - k, k, acc); pop!(acc)
        end
    end
    rec(N, N, Int[])
    return [OrbitPattern(p, N) for p in parts]
end
