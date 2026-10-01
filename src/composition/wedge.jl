# Product rules on the wedge (PLAN §6, v0.6): a triangle rule times a rule on [-1, 1].
#
# A polynomial of total degree d in (x, y, z) is a sum of p(x, y) zʳ with deg p + r ≤ d, so a
# triangle rule of degree d times a line rule of degree d integrates it exactly. The claim
# is that total degree, the smaller of the two, as for `TensorProduct`.

"""
    WedgeProduct(triangle, line)
    WedgeProduct()

Product rules on a [`Wedge`](@ref): a rule on the triangle times a rule on `[-1, 1]`, exact to
the smaller of their two degrees. `rule(Wedge(); degree)` takes the triangle family with the
fewest points at that degree; `WedgeProduct()` is Xiao–Gimbutas times Gauss–Legendre.
"""
struct WedgeProduct{A<:RuleFamily,B<:RuleFamily} <: CombinatorFamily
    tri::A
    line::B
end
WedgeProduct() = WedgeProduct(XiaoGimbutas(), GaussLegendre())

derivation(::Type{<:WedgeProduct}) = Derived()
derivation(f::WedgeProduct) = derivation(f.tri) isa Seeded || derivation(f.line) isa Seeded ? Seeded() : Derived()
family_name(::WedgeProduct) = "WedgeProduct"
describe_family(f::WedgeProduct) = "WedgeProduct(" * describe_family(f.tri) * " × " * describe_family(f.line) * ")"

# every triangle family, and the conical one, which is itself a combinator and so is not a
# leaf; on the line Gauss–Legendre, which no other line family beats on points
function candidates(::Type{<:WedgeProduct}, dom::Wedge, c::PolynomialDegree)
    isreference(dom) || return WedgeProduct[]
    tris = RuleFamily[]
    for F in leaf_families(), f in candidates(F, Simplex{2}(), c)
        push!(tris, f)
    end
    append!(tris, candidates(ConicalProduct, Simplex{2}(), c))
    return [WedgeProduct(t, GaussLegendre()) for t in tris]
end

npoints(f::WedgeProduct, dom, degree::Integer) =
    npoints(f.tri, Simplex{2}(), degree) * npoints(f.line, Interval(), degree)
claimed_degree(f::WedgeProduct, dom, degree) =
    min(claimed_degree(f.tri, Simplex{2}(), degree), claimed_degree(f.line, Interval(), degree))
degree_range(f::WedgeProduct, dom::Wedge) =
    intersect(degree_range(f.tri, Simplex{2}()), degree_range(f.line, Interval()))
degree_range(::WedgeProduct, dom) = 1:0
supports_type(f::WedgeProduct, ::Type{T}) where {T} = supports_type(f.tri, T) && supports_type(f.line, T)
function properties(f::WedgeProduct, dom, degree)
    pt, pl = properties(f.tri, Simplex{2}(), degree), properties(f.line, Interval(), degree)
    return (positive = pt.positive && pl.positive, interior = pt.interior && pl.interior, symmetry = :none, nested = false)
end
cost_estimate(f::WedgeProduct, dom, degree, T) =
    cost_estimate(f.tri, Simplex{2}(), degree, T) + cost_estimate(f.line, Interval(), degree, T)

function build(f::WedgeProduct, dom::Wedge, degree::Int, ctx::BuildContext{T}) where {T}
    t = build(f.tri, Simplex{2}(), degree, ctx)
    l = build(f.line, Interval(), degree, ctx)
    nt, nl = npoints(t), npoints(l)
    xs = Vector{SVector{3,T}}(undef, nt * nl)
    ws = Vector{T}(undef, nt * nl)
    # the products are formed with guard digits and rounded once, as in TensorProduct
    with_bits(ctx.bits + 32) do
        k = 0
        for i in 1:nt, j in 1:nl
            k += 1
            p = node_vector(t, i)
            xs[k] = SVector{3,T}(p[1], p[2], nodes(l)[j])
            ws[k] = finalize_number(ctx, big(weights(t)[i]) * big(weights(l)[j]))
        end
    end
    claim = PolynomialDegree(min(exactness(t).d, exactness(l).d))
    certs = [c for c in (certificate(t), certificate(l)) if c !== nothing]
    cert = isempty(certs) ? nothing :
           Certificate(equations = join(("$(n): " * c.equations for (n, c) in zip(("triangle", "line"), certs)), "; "),
                       residual = maximum(c -> c.residual, certs), residual_bits = maximum(c -> c.residual_bits, certs),
                       digits = target_digits(ctx), guard_digits = minimum(c -> c.guard_digits, certs),
                       cond = maximum(c -> c.cond, certs), iterations = maximum(c -> c.iterations, certs))
    prov = Provenance(family = "WedgeProduct", derivation = derivation(f),
                      path = vcat(["product of a triangle rule and a rule on [-1, 1]",
                                   "triangle: $(describe_family(f.tri)), $(nt) points, degree $(exactness(t).d)",
                                   "line: $(describe_family(f.line)), $(nl) points, degree $(exactness(l).d)"],
                                  ["triangle " * s for s in provenance(t).path]),
                      seed_source = provenance(t).seed_source,
                      citations = unique(vcat(provenance(t).citations, provenance(l).citations)),
                      license = provenance(t).license)
    return QuadratureRule(xs, ws, Wedge(), claim, prov, cert)
end
