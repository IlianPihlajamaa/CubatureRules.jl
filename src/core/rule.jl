# Rules.
#
# There is one rule type. A separate SVector-backed runtime form existed during development
# and was removed after measurement: 1.0–1.4× faster at best even on a near-free integrand, slower
# for large rules, compile time growing to seconds per rule size, and no allocation saving —
# integrating with this type already allocates nothing.

"""
    QuadratureRule{D,T,Dom,C,S,W}

A quadrature rule of dimension `D` in number type `T` on domain `Dom`, with exactness
claim `C`. Storage is a type parameter: `S` is the node container (`Vector{T}` in 1D,
`Vector{SVector{D,T}}` otherwise) and `W` the weight container.

A rule carries everything needed to explain, verify and cite itself, and integrating with
it allocates nothing when the integrand does not. Code that needs an `isbits` rule, such as
a GPU kernel, can take `SVector{N}(nodes(r))` and `SVector{N}(weights(r))` directly.

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

# Real nodes with complex weights: the rules for an oscillatory weight e^{iωx}. The type
# parameter T stays the real node type, which is what a caller asks for with `T = Float64`.
function QuadratureRule(nodes::AbstractVector{T}, weights::AbstractVector{Complex{T}}, domain::Domain{1},
                        claim::ExactnessClaim, prov::Provenance, cert = nothing) where {T<:Real}
    return QuadratureRule{1,T,typeof(domain),typeof(claim),typeof(nodes),typeof(weights)}(
        nodes, weights, domain, claim, prov, cert)
end

"The rule `r` with the provenance `p`: nodes, weights, domain, claim and certificate unchanged."
with_provenance(r::QuadratureRule, p::Provenance) =
    QuadratureRule(r.nodes, r.weights, r.domain, r.exactness, p, r.certificate)

"""
    nodes(r)

The nodes of the rule: a vector of numbers in one dimension, otherwise a vector of
`SVector{D}`. For a rule built with `rule` on a mapped domain, the nodes are on that domain.
"""
nodes(r::QuadratureRule) = r.nodes

"""
    weights(r)

The weights of the rule, in the same order as [`nodes`](@ref). On a weighted domain they
include the weight, so they add up to its total mass.
"""
weights(r::QuadratureRule) = r.weights

"""
    domain(r)

The domain the rule integrates over.
"""
domain(r::QuadratureRule) = r.domain

"""
    exactness(r) -> ExactnessClaim

The rule's exactness claim: a [`PolynomialDegree`](@ref), a [`SpanOf`](@ref) or
[`NoClaim`](@ref). Use [`degree`](@ref) for the degree itself.
"""
exactness(r::QuadratureRule) = r.exactness

"""
    provenance(r) -> Provenance

How the rule was made: the family, the seed, the steps of the construction, the citations,
the licence of its data, and why the selector chose it. See
[`Provenance`](@ref CubatureRules.Provenance).
"""
provenance(r::QuadratureRule) = r.provenance

"""
    certificate(r) -> Union{Certificate,Nothing}

The record of how the rule's defining equations were solved: their residual, the condition
number, the guard digits and the number of iterations. This describes the solver, not the
correctness of the rule; for that, use [`check`](@ref). See
[`Certificate`](@ref CubatureRules.Certificate).
"""
certificate(r::QuadratureRule) = r.certificate
npoints(r::QuadratureRule) = length(r.weights)
degree(r::QuadratureRule) = degree(r.exactness)

"""
    family(r) -> String

The name of the family that built the rule.
"""
family(r::QuadratureRule) = r.provenance.family
derivation(r::QuadratureRule) = r.provenance.derivation

"""
    numbertype(r) -> Type

The number type of the rule's nodes and weights: `Float64`, `BigFloat`, `Rational{BigInt}`, …
(the weights of an oscillatory rule may be complex over it).
"""
numbertype(::QuadratureRule{D,T}) where {D,T} = T

# A rule is a collection of (node, weight) pairs: `for (x, w) in r`.
Base.length(r::QuadratureRule) = npoints(r)
Base.eltype(::Type{<:QuadratureRule{D,T,Dom,C,S,W}}) where {D,T,Dom,C,S,W} = Tuple{eltype(S),eltype(W)}
Base.iterate(r::QuadratureRule, i::Int = 1) = i > npoints(r) ? nothing : ((r.nodes[i], r.weights[i]), i + 1)

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
job.
"""
function rule_hash(r::QuadratureRule)
    io = IOBuffer()
    print(io, family(r), '|')
    _canonical(io, r.domain)
    print(io, '|', describe(r.exactness), '|', numbertype(r))
    numbertype(r) === BigFloat && print(io, '|', precision(real(first(r.weights))))
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

# A domain as `rule_hash` sees it: the name of its type and its fields, recursively, numbers
# written exactly. `show` is for people and may change; this may not, without changing every
# hash. A weight given by a function, which has no stable name, is written as its label.
_canonical(io::IO, x::Number) = print(io, _exact_string(x))     # the value, not its type: Int is Int32 on some platforms
_canonical(io::IO, x::Union{Symbol,AbstractString}) = print(io, repr(x))
_canonical(io::IO, ::Nothing) = print(io, "nothing")
function _canonical(io::IO, x::Union{Tuple,AbstractArray})
    print(io, '[')
    foreach(y -> (_canonical(io, y); print(io, ',')), x)
    print(io, ']')
end
_canonical(io::IO, w::Union{FunctionWeight,MomentWeight}) = print(io, nameof(typeof(w)), '(', repr(w.label), ')')
function _canonical(io::IO, x)
    T = typeof(x)
    (isstructtype(T) && !(x isa Function)) || throw(ArgumentError("rule_hash cannot write a $T canonically"))
    print(io, nameof(T), '(')
    for i in 1:fieldcount(T)
        _canonical(io, getfield(x, i))
        print(io, ',')
    end
    print(io, ')')
end


_is_exact_type(::Type{T}) where {T} = T <: Rational || T <: Integer
"The precision of a rule's numbers in bits: 0 for exact rationals, a BigFloat's own precision."
_rule_precision(r::QuadratureRule) = _is_exact_type(numbertype(r)) ? 0 :
                                     numbertype(r) === BigFloat ? precision(real(first(r.weights))) : bits_of(numbertype(r))
