# Tensor-product rules on boxes (PLAN §2.2, §6 Tier 1).
#
# The claim is the conservative one: a tensor rule is exact on a tensor-product polynomial
# space, which is larger than the total-degree space, but the selector ranks candidates by
# node count against a *total* degree, so that is what is reported. The decision is recorded
# here so `show` never prints something arbitrary.

"""
    TensorProduct(families::Tuple)
    TensorProduct(family, D)

Tensor product of 1D rules on an [`Orthotope`](@ref), one family per axis. With `m_i` points
on axis `i` it has `∏ m_i` points and claims the **total** degree `min_i (2m_i - 1)`, not the
larger tensor-product space it is in fact exact on.
"""
struct TensorProduct{Fs<:Tuple} <: CombinatorFamily
    families::Fs
end
TensorProduct(f::RuleFamily, D::Integer) = TensorProduct(ntuple(_ -> f, Int(D)))

derivation(::Type{<:TensorProduct}) = Derived()
derivation(f::TensorProduct) = any(fi -> derivation(fi) isa Seeded, f.families) ? Seeded() : Derived()
family_name(::TensorProduct) = "TensorProduct"
function describe_family(f::TensorProduct)
    names = unique(describe_family.(f.families))
    return "TensorProduct(" * (length(names) == 1 ? only(names) : join(names, " ⊗ ")) * ")"
end

function candidates(::Type{<:TensorProduct}, dom::Orthotope{D}, c::PolynomialDegree) where {D}
    isreference(dom) || return TensorProduct[]
    out = TensorProduct[]
    for F in leaf_families(), inner in candidates(F, Interval(), c)
        push!(out, TensorProduct(inner, D))
    end
    return out
end

npoints(f::TensorProduct, dom::Orthotope, degree::Integer) =
    prod(npoints(fi, Interval(), degree) for fi in f.families)
claimed_degree(f::TensorProduct, dom, degree) = minimum(claimed_degree(fi, Interval(), degree) for fi in f.families)
function degree_range(f::TensorProduct, dom::Orthotope)
    lo = maximum(first(degree_range(fi, Interval())) for fi in f.families)
    hi = minimum(last(degree_range(fi, Interval())) for fi in f.families)
    return lo:hi
end
degree_range(::TensorProduct, dom) = 1:0
supports_type(f::TensorProduct, ::Type{T}) where {T} = all(fi -> supports_type(fi, T), f.families)
function properties(f::TensorProduct, dom, degree)
    ps = [properties(fi, Interval(), degree) for fi in f.families]
    return (positive = all(p -> p.positive, ps), interior = all(p -> p.interior, ps),
            symmetry = :none, nested = all(p -> p.nested, ps))
end
cost_estimate(f::TensorProduct, dom, degree, T) =
    sum(cost_estimate(fi, Interval(), degree, T) for fi in f.families) * _precision_factor(T)

function build(f::TensorProduct, dom::Orthotope{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    length(f.families) == D ||
        throw(ArgumentError("TensorProduct with $(length(f.families)) factors cannot build on $(dom)"))
    parts = [build(fi, Interval(), degree, ctx) for fi in f.families]
    ms = [npoints(p) for p in parts]
    xs = Vector{SVector{D,T}}(undef, prod(ms))
    ws = Vector{T}(undef, prod(ms))
    # the product of the factors' weights is formed with guard digits and rounded once, so
    # the result carries the precision of the factors rather than the ambient one
    with_bits(ctx.bits + 32) do
        for (k, I) in enumerate(CartesianIndices(Tuple(ms)))
            xs[k] = SVector{D,T}(ntuple(j -> nodes(parts[j])[I[j]], D))
            ws[k] = finalize_number(ctx, prod(ntuple(j -> big(weights(parts[j])[I[j]]), D)))
        end
    end
    claim = PolynomialDegree(minimum(exactness(p).d for p in parts))   # `degree` is the argument here
    certs = [certificate(p) for p in parts if certificate(p) !== nothing]
    cert = isempty(certs) ? nothing :
           Certificate(equations = "per axis: " * first(certs).equations,
                       residual = maximum(c -> c.residual, certs), residual_bits = first(certs).residual_bits,
                       digits = target_digits(ctx), guard_digits = maximum(c -> c.guard_digits, certs),
                       cond = maximum(c -> c.cond, certs), iterations = maximum(c -> c.iterations, certs))
    prov = Provenance(family = "TensorProduct", derivation = derivation(f),
                      path = vcat(["tensor product of $D one-dimensional rules"],
                                  ["axis $j: " * describe_family(f.families[j]) for j in 1:D],
                                  ["claim: conservative total degree $(claim.d), not the tensor-product space"]),
                      seed_source = join(unique(provenance(p).seed_source for p in parts), "; "),
                      citations = unique(reduce(vcat, [provenance(p).citations for p in parts])),
                      license = "MIT")
    return QuadratureRule(xs, ws, Orthotope{D}(), claim, prov, cert)
end

"""
    r₁ ⊗ r₂

Tensor product of two rules on intervals or boxes, giving a rule on the product box. The
claim is the conservative total degree `min(d₁, d₂)` (PLAN §2.3).
"""
function ⊗(a::QuadratureRule{Da,Ta}, b::QuadratureRule{Db,Tb}) where {Da,Ta,Db,Tb}
    (_is_box(a.domain) && _is_box(b.domain)) ||
        throw(ArgumentError("⊗ combines rules on intervals or boxes; got $(a.domain) and $(b.domain)"))
    D = Da + Db
    T = promote_type(Ta, Tb)
    xs = Vector{SVector{D,T}}(undef, npoints(a) * npoints(b))
    ws = Vector{T}(undef, npoints(a) * npoints(b))
    k = 0
    for i in 1:npoints(a), j in 1:npoints(b)
        k += 1
        xs[k] = SVector{D,T}(ntuple(l -> l <= Da ? T(node_vector(a, i)[l]) : T(node_vector(b, j)[l - Da]), D))
        ws[k] = T(weights(a)[i]) * T(weights(b)[j])
    end
    lo = SVector{D,T}(ntuple(l -> l <= Da ? T(_box_lo(a.domain)[l]) : T(_box_lo(b.domain)[l - Da]), D))
    hi = SVector{D,T}(ntuple(l -> l <= Da ? T(_box_hi(a.domain)[l]) : T(_box_hi(b.domain)[l - Da]), D))
    claim = _tensor_claim(a.exactness, b.exactness)
    prov = Provenance(family = "TensorProduct", derivation = derivation(a) isa Seeded || derivation(b) isa Seeded ?
                                                             Seeded() : Derived(),
                      path = ["⊗ of $(family(a)) and $(family(b))",
                              "claim: conservative total degree, not the tensor-product space"],
                      seed_source = provenance(a).seed_source * "; " * provenance(b).seed_source,
                      citations = unique(vcat(provenance(a).citations, provenance(b).citations)))
    return QuadratureRule(xs, ws, Orthotope{D,T}(lo, hi), claim, prov, nothing)
end

_is_box(::Union{Interval,Orthotope}) = true
_is_box(::Domain) = false
_box_lo(d::Interval) = (d.a,)
_box_hi(d::Interval) = (d.b,)
_box_lo(d::Orthotope) = d.lo
_box_hi(d::Orthotope) = d.hi
_tensor_claim(a::PolynomialDegree, b::PolynomialDegree) = PolynomialDegree(min(a.d, b.d))
_tensor_claim(a::ExactnessClaim, b::ExactnessClaim) = NoClaim()
