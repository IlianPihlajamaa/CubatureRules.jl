# Domains (PLAN §2.1). Separate from rule families and from rules.

"""
    Domain{D,T}

An integration domain of dimension `D` whose geometry is described in number type `T`.
"""
abstract type Domain{D,T} end

dimension(::Domain{D}) where {D} = D
dimension(::Type{<:Domain{D}}) where {D} = D

"""
    measure(dom)

The total mass of the domain: its length, area or volume, or on a weighted domain the
integral of the weight. The weights of a rule on `dom` add up to this. Unbounded unweighted
domains return `Inf`.
"""
function measure end

"""
    isreference(dom) -> Bool

Whether `dom` is the reference domain of its kind — `[-1, 1]`, the unit simplex, the unit
sphere and so on — on which families construct their rules.
"""
function isreference end

"""
    reference(dom)

The reference domain of the same kind as `dom`. A rule requested on `dom` is built on
`reference(dom)` and mapped. A domain carrying a [`MomentWeight`](@ref) is its own reference,
because its moments are stated for that interval.
"""
function reference end

"""
    vertices(s::Simplex)

The vertices of the simplex, as a tuple of `SVector`s.
"""
function vertices end

# ---------------------------------------------------------------------------------------
# Interval

"""
    Interval(a, b)
    Interval()

The closed interval `[a, b]`. `Interval()` is the reference interval `[-1, 1]`.
"""
struct Interval{T} <: Domain{1,T}
    a::T
    b::T
    function Interval{T}(a, b) where {T}
        a < b || throw(ArgumentError("an interval needs a < b, got [$a, $b]"))
        return new{T}(a, b)
    end
end
Interval(a::T, b::T) where {T} = Interval{T}(a, b)
Interval(a, b) = Interval(promote(a, b)...)
Interval() = Interval(-1, 1)

isreference(d::Interval) = d.a == -1 && d.b == 1
reference(::Interval) = Interval()
measure(d::Interval) = d.b - d.a
convert_domain(::Type{T}, d::Interval) where {T} = Interval{T}(T(d.a), T(d.b))

Base.show(io::IO, d::Interval) = isreference(d) ? print(io, "Interval()") :
                                 print(io, "Interval(", d.a, ", ", d.b, ")")

# ---------------------------------------------------------------------------------------
# Simplex

"""
    Simplex{D}()
    Simplex(v₀, v₁, …, v_D)

A `D`-dimensional simplex. `Simplex{D}()` is the reference simplex with vertices
`0, e₁, …, e_D` and measure `1/D!`. Vertices may be given as tuples or `SVector`s.

Barycentric coordinates `λ` are ordered to match the vertices, so on the reference
simplex `λ = (1 - Σx, x₁, …, x_D)`.
"""
struct Simplex{D,T,N} <: Domain{D,T}
    vertices::SVector{N,SVector{D,T}}
    function Simplex{D,T,N}(v::SVector{N,SVector{D,T}}) where {D,T,N}
        N == D + 1 || throw(ArgumentError("a $D-simplex has $(D+1) vertices, got $N"))
        return new{D,T,N}(v)
    end
end

function Simplex{D}() where {D}
    D isa Int && D >= 1 || throw(ArgumentError("simplex dimension must be a positive Int"))
    v = SVector{D + 1,SVector{D,Int}}(ntuple(i -> SVector{D,Int}(ntuple(j -> Int(j == i - 1), D)), D + 1))
    return Simplex{D,Int,D + 1}(v)
end

function Simplex(vs::Union{Tuple,AbstractVector}...)
    N = length(vs)
    D = N - 1
    all(v -> length(v) == D, vs) ||
        throw(ArgumentError("a simplex with $N vertices needs $D coordinates per vertex"))
    T = promote_type((eltype(v) for v in vs)...)
    sv = SVector{N,SVector{D,T}}(ntuple(i -> SVector{D,T}(Tuple(vs[i])), N))
    return Simplex{D,T,N}(sv)
end

vertices(s::Simplex) = s.vertices
@inline _ref_vertex(i::Int, ::Val{D}) where {D} = SVector{D,Int}(ntuple(j -> Int(j == i - 1), Val(D)))
function isreference(s::Simplex{D}) where {D}
    for i in 1:(D + 1)
        s.vertices[i] == _ref_vertex(i, Val(D)) || return false
    end
    return true
end
reference(::Simplex{D}) where {D} = Simplex{D}()

function convert_domain(::Type{S}, s::Simplex{D,T,N}) where {S,D,T,N}
    return Simplex{D,S,N}(SVector{N,SVector{D,S}}(map(v -> SVector{D,S}(map(S, v)), s.vertices)))
end

"The edge matrix `[v₁-v₀ … v_D-v₀]` of a simplex."
@inline edge_matrix(s::Simplex{D,T}) where {D,T} =
    SMatrix{D,D}(ntuple(k -> s.vertices[(k - 1) ÷ D + 2][(k - 1) % D + 1] - s.vertices[1][(k - 1) % D + 1], Val(D * D)))

function measure(s::Simplex{D}) where {D}
    isreference(s) && return 1 // factorial(D)
    A = edge_matrix(s)
    return abs(det(A)) / factorial(D)
end

function Base.show(io::IO, s::Simplex{D}) where {D}
    if isreference(s)
        print(io, "Simplex{", D, "}()")
    else
        print(io, "Simplex(", join((Tuple(v) for v in s.vertices), ", "), ")")
    end
end

"""
    barycentric(s::Simplex, x)

Barycentric coordinates of the point `x` with respect to `s`, in generic arithmetic.
"""
function barycentric(s::Simplex{D}, x::AbstractVector) where {D}
    if isreference(s)
        λ0 = one(eltype(x)) - sum(x)
        return SVector{D + 1}(λ0, x...)
    end
    μ = edge_matrix(s) \ (SVector{D}(x) - s.vertices[1])
    return SVector{D + 1}(one(eltype(μ)) - sum(μ), μ...)
end

"""
    cartesian(s::Simplex, λ)

Cartesian point with barycentric coordinates `λ` with respect to `s`.
"""
cartesian(s::Simplex{D}, λ::AbstractVector) where {D} =
    sum(λ[i] * s.vertices[i] for i in 1:(D + 1))

# ---------------------------------------------------------------------------------------
# Point-in-domain predicates.
#
# Tolerance policy: membership is decided in the arithmetic of the point, with an explicit
# absolute tolerance on the barycentric coordinates (default zero). `isinterior` is strict:
# every barycentric coordinate must exceed `tol`. Verification (§8) uses tol = 0, so a node
# is "interior" only if it is interior in its own number type.

"""
    indomain(x, dom; tol = 0)

Whether `x` lies in the closed domain, allowing barycentric coordinates down to `-tol`.
"""
indomain(x, s::Simplex; tol = 0) = all(>=(-tol), barycentric(s, x))
indomain(x, d::Interval; tol = 0) = d.a - tol <= only(x) <= d.b + tol

"""
    isinterior(x, dom; tol = 0)

Whether `x` lies strictly inside the domain: every barycentric coordinate `> tol`.
"""
isinterior(x, s::Simplex; tol = 0) = all(>(tol), barycentric(s, x))
isinterior(x, d::Interval; tol = 0) = d.a + tol < only(x) < d.b - tol

# ---------------------------------------------------------------------------------------
# Orthotope (box).

"""
    Orthotope(lo, hi)
    Orthotope{D}()

The box `[lo₁, hi₁] × ⋯ × [lo_D, hi_D]`. `Orthotope{D}()` is the reference box `[-1, 1]^D`,
the natural home of tensor-product rules.
"""
struct Orthotope{D,T} <: Domain{D,T}
    lo::SVector{D,T}
    hi::SVector{D,T}
    function Orthotope{D,T}(lo, hi) where {D,T}
        all(lo .< hi) || throw(ArgumentError("an orthotope needs lo < hi in every direction"))
        return new{D,T}(lo, hi)
    end
end

function Orthotope(lo, hi)
    length(lo) == length(hi) || throw(ArgumentError("lo and hi need the same length"))
    D = length(lo)
    T = promote_type(eltype(lo), eltype(hi))
    return Orthotope{D,T}(SVector{D,T}(lo), SVector{D,T}(hi))
end
Orthotope{D}() where {D} = Orthotope{D,Int}(SVector{D,Int}(ntuple(_ -> -1, D)), SVector{D,Int}(ntuple(_ -> 1, D)))

isreference(d::Orthotope) = all(d.lo .== -1) && all(d.hi .== 1)
reference(::Orthotope{D}) where {D} = Orthotope{D}()
measure(d::Orthotope) = prod(d.hi - d.lo)
convert_domain(::Type{S}, d::Orthotope{D}) where {S,D} =
    Orthotope{D,S}(SVector{D,S}(map(S, d.lo)), SVector{D,S}(map(S, d.hi)))
indomain(x, d::Orthotope; tol = 0) = all(d.lo .- tol .<= x .<= d.hi .+ tol)
isinterior(x, d::Orthotope; tol = 0) = all(d.lo .+ tol .< x .< d.hi .- tol)

function Base.show(io::IO, d::Orthotope{D}) where {D}
    isreference(d) ? print(io, "Orthotope{", D, "}()") : print(io, "Orthotope(", Tuple(d.lo), ", ", Tuple(d.hi), ")")
end

# ---------------------------------------------------------------------------------------
# Weighted domains.

"""
    WeightedDomain(base, weight)

The domain `base` equipped with the measure `weight(x) dx`. The base is stored with its
concrete type so that access costs nothing.
"""
struct WeightedDomain{D,T,B<:Domain{D,T},W} <: Domain{D,T}
    base::B
    weight::W
end

isreference(d::WeightedDomain) = isreference(d.base)
reference(d::WeightedDomain) = WeightedDomain(reference(d.base), d.weight)
convert_domain(::Type{S}, d::WeightedDomain) where {S} = WeightedDomain(convert_domain(S, d.base), d.weight)
Base.show(io::IO, d::WeightedDomain) = print(io, "WeightedDomain(", d.base, ", ", d.weight, ")")

"""
    JacobiWeight(α, β)

The weight `(1 - x)^α (1 + x)^β` on the reference interval `[-1, 1]`.
"""
struct JacobiWeight{T}
    α::T
    β::T
end
JacobiWeight(α, β) = JacobiWeight(promote(α, β)...)
(w::JacobiWeight)(x) = (1 - x)^w.α * (1 + x)^w.β

measure(d::WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight}) = jacobi_mass(d.weight.α, d.weight.β)
indomain(x, d::WeightedDomain; kw...) = indomain(x, d.base; kw...)
isinterior(x, d::WeightedDomain; kw...) = isinterior(x, d.base; kw...)

# ---------------------------------------------------------------------------------------
# Domains named by the plan but not yet implemented. They exist so that the selector can
# say *when* they arrive rather than failing with a MethodError.

struct NotYetImplemented <: Exception
    what::String
    release::String
end
Base.showerror(io::IO, e::NotYetImplemented) =
    print(io, e.what, " is not implemented yet; it is scheduled for CubatureRules ", e.release, ".")

for (name, rel) in (
                    (:Polytope, "v0.6"), (:Wedge, "v0.6"), (:Pyramid, "v0.6"))
    @eval begin
        @doc "`$($(string(name)))` — not yet implemented (scheduled for $($rel))." struct $name{D,T} <: Domain{D,T}
            $name{D,T}(args...) where {D,T} = throw(NotYetImplemented($(string(name)) * " domains", $rel))
        end
        $name(args...) = throw(NotYetImplemented($(string(name)) * " domains", $rel))
        $name{D}(args...) where {D} = throw(NotYetImplemented($(string(name)) * " domains", $rel))
    end
end
