# Rules: the construction-time record and the runtime (static) form (PLAN §2.4).

"""
    QuadratureRule{D,T,Dom,C,S,W}

A quadrature rule of dimension `D` in number type `T` on domain `Dom`, with exactness
claim `C`. Storage is a type parameter: `S` is the node container (`Vector{T}` in 1D,
`Vector{SVector{D,T}}` otherwise) and `W` the weight container.

This is the *construction-time record*: it carries everything needed to explain, verify and
cite itself. For hot loops convert it with [`static`](@ref).

Fields are accessed through [`nodes`](@ref), [`weights`](@ref), [`domain`](@ref),
[`exactness`](@ref), [`provenance`](@ref) and [`certificate`](@ref).
"""
struct QuadratureRule{D,T,Dom<:Domain{D},C<:ExactnessClaim,S,W}
    nodes::S
    weights::W
    domain::Dom
    exactness::C
    provenance::Provenance
    certificate::Union{Nothing,Certificate}
    function QuadratureRule{D,T,Dom,C,S,W}(nodes, weights, domain, claim, prov, cert) where {D,T,Dom,C,S,W}
        length(nodes) == length(weights) ||
            throw(DimensionMismatch("$(length(nodes)) nodes but $(length(weights)) weights"))
        return new{D,T,Dom,C,S,W}(nodes, weights, domain, claim, prov, cert)
    end
end

function QuadratureRule(nodes::AbstractVector{<:SVector{D,T}}, weights::AbstractVector{T}, domain::Domain{D},
                        claim::ExactnessClaim, prov::Provenance, cert = nothing) where {D,T}
    D == 1 && return QuadratureRule(map(only, nodes), weights, domain, claim, prov, cert)
    return QuadratureRule{D,T,typeof(domain),typeof(claim),typeof(nodes),typeof(weights)}(
        nodes, weights, domain, claim, prov, cert)
end

function QuadratureRule(nodes::AbstractVector{T}, weights::AbstractVector{T}, domain::Domain{1},
                        claim::ExactnessClaim, prov::Provenance, cert = nothing) where {T<:Number}
    return QuadratureRule{1,T,typeof(domain),typeof(claim),typeof(nodes),typeof(weights)}(
        nodes, weights, domain, claim, prov, cert)
end

nodes(r::QuadratureRule) = r.nodes
weights(r::QuadratureRule) = r.weights
domain(r::QuadratureRule) = r.domain
exactness(r::QuadratureRule) = r.exactness
provenance(r::QuadratureRule) = r.provenance
certificate(r::QuadratureRule) = r.certificate
npoints(r::QuadratureRule) = length(r.weights)
degree(r::QuadratureRule) = degree(r.exactness)
family(r::QuadratureRule) = r.provenance.family
derivation(r::QuadratureRule) = r.provenance.derivation
Base.eltype(::Type{<:QuadratureRule{D,T}}) where {D,T} = T
Base.length(r::QuadratureRule) = npoints(r)

"Node `i` as an `SVector{D}` (also in 1D, where [`nodes`](@ref) stores scalars)."
node_vector(r::QuadratureRule{1}, i) = SVector(r.nodes[i])
node_vector(r::QuadratureRule, i) = r.nodes[i]

# ---------------------------------------------------------------------------------------
# Equality and hashing are by content — nodes, weights, domain, claim, family — so user
# code can memoise rules in a Dict. Provenance paths and certificates are not part of
# identity: two constructions of the same rule are the same rule.

Base.:(==)(a::QuadratureRule, b::QuadratureRule) =
    a.nodes == b.nodes && a.weights == b.weights && a.domain == b.domain &&
    a.exactness == b.exactness && family(a) == family(b)

Base.hash(r::QuadratureRule, h::UInt) =
    hash(r.nodes, hash(r.weights, hash(r.domain, hash(r.exactness, hash(family(r), hash(:QuadratureRule, h))))))

function Base.isapprox(a::QuadratureRule, b::QuadratureRule; kw...)
    npoints(a) == npoints(b) || return false
    a.domain == b.domain || return false
    return all(i -> isapprox(a.nodes[i], b.nodes[i]; kw...), eachindex(a.nodes)) &&
           isapprox(a.weights, b.weights; kw...)
end

"""
    rule_hash(r)

A SHA-256 content hash of the rule's nodes, weights, precision, domain and claim, stable
across Julia versions and platforms (unlike `Base.hash`). Used by the reproducibility CI
job (PLAN §5 item 9).
"""
function rule_hash(r::QuadratureRule)
    io = IOBuffer()
    print(io, family(r), '|', r.domain, '|', describe(r.exactness), '|', eltype(r))
    eltype(r) === BigFloat && print(io, '|', precision(first(r.weights)))
    for i in eachindex(r.weights)
        print(io, '|')
        for c in node_vector(r, i)
            print(io, _exact_string(c), ',')
        end
        print(io, _exact_string(r.weights[i]))
    end
    return bytes2hex(sha256(take!(io)))
end

_exact_string(x::BigFloat) = string(x)                      # MPFR: round-trips at its precision
_exact_string(x::AbstractFloat) = repr(x)
_exact_string(x::Rational) = string(numerator(x), "//", denominator(x))
_exact_string(x) = string(x)

# ---------------------------------------------------------------------------------------
# Runtime form.

"""
    StaticQuadratureRule{D,T,N,Dom,C,F}

The runtime form of a rule: nodes and weights in `SVector`s, no provenance or
certificate (the family name is kept as the type parameter `F`). `isbits` whenever `T`
and the domain are, so it can live in registers and on a GPU. Produced by
[`static`](@ref).
"""
struct StaticQuadratureRule{D,T,N,Dom<:Domain{D},C<:ExactnessClaim,F,S}
    nodes::S
    weights::SVector{N,T}
    domain::Dom
    exactness::C
end

"""
    static(rule)

Convert a construction-time [`QuadratureRule`](@ref) to its runtime form
[`StaticQuadratureRule`](@ref). Intended for small rules in hot loops.
"""
function static(r::QuadratureRule{D,T}) where {D,T}
    N = npoints(r)
    F = Symbol(family(r))
    w = SVector{N,T}(r.weights)
    x = D == 1 ? SVector{N,T}(r.nodes) : SVector{N,SVector{D,T}}(r.nodes)
    return StaticQuadratureRule{D,T,N,typeof(r.domain),typeof(r.exactness),F,typeof(x)}(x, w, r.domain, r.exactness)
end

nodes(r::StaticQuadratureRule) = r.nodes
weights(r::StaticQuadratureRule) = r.weights
domain(r::StaticQuadratureRule) = r.domain
exactness(r::StaticQuadratureRule) = r.exactness
npoints(::StaticQuadratureRule{D,T,N}) where {D,T,N} = N
degree(r::StaticQuadratureRule) = degree(r.exactness)
family(::StaticQuadratureRule{D,T,N,Dom,C,F}) where {D,T,N,Dom,C,F} = String(F)
Base.eltype(::Type{<:StaticQuadratureRule{D,T}}) where {D,T} = T

function Base.show(io::IO, r::StaticQuadratureRule{D,T,N}) where {D,T,N}
    print(io, "StaticQuadratureRule{", D, ",", T, "}(", family(r), ", ", N, " points, ",
          describe(r.exactness), ", on ", r.domain, ")")
end

const AnyRule = Union{QuadratureRule,StaticQuadratureRule}
