# Domains. Separate from rule families and from rules.

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

The value is exact where it can be: a `Rational` (or an integer) for domains with integer or
rational coordinates and weights with rational moments. Otherwise it is a `BigFloat` at the
ambient precision, which callers set with `setprecision` — a complex one for an oscillatory
weight — or, for coordinates in a floating-point type, a number of that type. Compare with a
rule's `sum(weights(r))` after converting both to one type.
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

"""
The type of a domain's coordinates given as tuples or vectors: each entry promoted, so that
`(0, 0.5)` gives `Float64` and not the `Real` that `eltype` of the tuple would.
"""
_coordinate_type(vs...) = promote_type((typeof(x) for v in vs for x in v)...)

# the type affine maps of a shape with coordinates of type T are computed in
_frametype(::Type{T}) where {T} = T <: Integer ? Rational{BigInt} : T <: Rational ? Rational{BigInt} : T

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
        s = new{D,T,N}(v)
        iszero(det(SMatrix{D,D,_frametype(T)}(edge_matrix(s)))) &&
            throw(ArgumentError("degenerate simplex: its vertices $(Tuple(Tuple.(v))) span no volume"))
        return s
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
    T = _coordinate_type(vs...)
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

"""
    affine_frame(dom, S) -> (A, b, J)

The affine map `ξ ↦ A ξ + b` from the reference domain of `dom` onto `dom`, in type `S`, and
the factor `J` by which it scales the measure: `|det A|`; for a triangle in space its area
element, for a sphere `ρ^(D−1)`. In one dimension `A`, `b` and `J` are numbers. The one
description of where a domain lies, from which [`map_to`](@ref), `integrate` over a domain,
[`affine_map`](@ref) and verification's map back to the reference are all made.
"""
function affine_frame end

@inline function affine_frame(s::Simplex{D}, ::Type{S}) where {D,S}
    A = SMatrix{D,D,S}(edge_matrix(s))
    return A, SVector{D,S}(s.vertices[1]), abs(det(A))
end
@inline function affine_frame(d::Interval, ::Type{S}) where {S}
    h = (S(d.b) - S(d.a)) / 2
    return h, S(d.a) + h, abs(h)
end

function measure(s::Simplex{D,T}) where {D,T}
    isreference(s) && return 1 // factorial(D)
    # in exact arithmetic for integer and rational vertices, so the measure keeps their exactness
    A = SMatrix{D,D,_frametype(T)}(edge_matrix(s))
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
# every barycentric coordinate must exceed `tol`. Verification uses tol = 0, so a node
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
    T = _coordinate_type(lo, hi)
    return Orthotope{D,T}(SVector{D,T}(lo), SVector{D,T}(hi))
end
Orthotope{D}() where {D} = Orthotope{D,Int}(SVector{D,Int}(ntuple(_ -> -1, D)), SVector{D,Int}(ntuple(_ -> 1, D)))

isreference(d::Orthotope) = all(d.lo .== -1) && all(d.hi .== 1)
@inline function affine_frame(d::Orthotope{D}, ::Type{S}) where {D,S}
    h = SVector{D,S}((d.hi - d.lo) ./ 2)
    return SMatrix{D,D,S}(Diagonal(h)), SVector{D,S}(d.lo) + h, abs(prod(h))
end
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
concrete type so that access costs nothing. A weight that only makes sense on one kind of
domain — a [`JacobiWeight`](@ref) on an interval, say — is refused on any other here, by
[`check_weight`](@ref), rather than further down.
"""
struct WeightedDomain{D,T,B<:Domain{D,T},W} <: Domain{D,T}
    base::B
    weight::W
    function WeightedDomain(base::B, weight::W) where {D,T,B<:Domain{D,T},W}
        check_weight(base, weight)
        return new{D,T,B,W}(base, weight)
    end
end

"""
    check_weight(base, weight)

Throw an `ArgumentError` if `weight` is not a weight on domains like `base`. Every weight of
this package declares where it lives; a weight defined elsewhere is accepted anywhere unless
it adds a method.
"""
check_weight(base::Domain, weight) = nothing
_wrong_base(weight, what, base) =
    throw(ArgumentError("$(nameof(typeof(weight))) is a weight on $(what), not on $(base)"))

"""
    same_domain(a, b) -> Bool

Whether two domains are the same up to the number types of their parameters:
`JacobiWeight(1, 1)` and `JacobiWeight(1.0, 1.0)` are, `Interval()` and a weighted interval
are not. Fields are compared recursively, numbers by value; anything else, a function for
instance, by `==`.
"""
function same_domain(a, b)
    a === b && return true
    (a isa Number && b isa Number) && return a == b
    if (a isa Union{Tuple,AbstractArray}) && (b isa Union{Tuple,AbstractArray})
        return length(a) == length(b) && all(same_domain(x, y) for (x, y) in zip(a, b))
    end
    Ta, Tb = typeof(a), typeof(b)
    (Base.typename(Ta) === Base.typename(Tb) && isstructtype(Ta) && fieldcount(Ta) > 0) || return a == b
    return all(same_domain(getfield(a, i), getfield(b, i)) for i in 1:fieldcount(Ta))
end

isreference(d::WeightedDomain) = isreference(d.base)
reference(d::WeightedDomain) = WeightedDomain(reference(d.base), d.weight)
convert_domain(::Type{S}, d::WeightedDomain) where {S} = WeightedDomain(convert_domain(S, d.base), d.weight)
Base.show(io::IO, d::WeightedDomain) = print(io, "WeightedDomain(", d.base, ", ", d.weight, ")")

"""
    JacobiWeight(α, β)

The weight `(1 - x)^α (1 + x)^β` on the reference interval `[-1, 1]`, and `(b - x)^α (x - a)^β`
on `Interval(a, b)`, as for [`PrincipalValue`](@ref) and [`FunctionWeight`](@ref). A rule on
`[-1, 1]` maps onto `[a, b]` with its weights scaled by `h^(α+β+1)`, `h = (b - a)/2`.
Calling the weight, `w(x)`, evaluates it on `[-1, 1]`.
"""
struct JacobiWeight{T}
    α::T
    β::T
end
JacobiWeight(α, β) = JacobiWeight(promote(α, β)...)
(w::JacobiWeight)(x) = (1 - x)^w.α * (1 + x)^w.β

"An interval with a [`JacobiWeight`](@ref)."
const JacobiInterval = WeightedDomain{1,<:Any,<:Interval,<:JacobiWeight}

check_weight(base::Interval, ::JacobiWeight) = nothing
check_weight(base::Domain, w::JacobiWeight) = _wrong_base(w, "an Interval", base)

function measure(d::JacobiInterval)
    m = jacobi_mass(d.weight.α, d.weight.β)
    isreference(d.base) && return m
    h = d.base.b isa Integer && d.base.a isa Integer ? (d.base.b - d.base.a) // 2 : (d.base.b - d.base.a) / 2
    return m * h^(d.weight.α + d.weight.β + 1)
end
# the frame of the base, with the weight's own scaling: (b - x)^α (x - a)^β dx = h^(α+β+1) × reference
@inline function affine_frame(d::JacobiInterval, ::Type{S}) where {S}
    A, b, _ = affine_frame(d.base, S)
    return A, b, abs(A)^(d.weight.α + d.weight.β + 1)
end
indomain(x, d::WeightedDomain; kw...) = indomain(x, d.base; kw...)
isinterior(x, d::WeightedDomain; kw...) = isinterior(x, d.base; kw...)
