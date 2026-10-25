# Product rules on the wedge: a triangle rule times a rule on [-1, 1].
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
inner_families(f::WedgeProduct) = (f.tri, f.line)

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
    return (positive = pt.positive && pl.positive, interior = pt.interior && pl.interior, symmetry = NoSymmetry(), nested = false)
end
cost_estimate(f::WedgeProduct, dom, degree, T) =
    cost_estimate(f.tri, Simplex{2}(), degree, T) + cost_estimate(f.line, Interval(), degree, T)

function build(f::WedgeProduct, dom::Wedge, degree::Int, ctx::BuildContext{T}) where {T}
    # Both factors with guard digits, as in TensorProduct, so each product of weights is rounded
    # once. A seeded triangle family at Float64 precision still delivers its stored table: its
    # weights are the table's correctly rounded values, and refining them would turn an
    # instant rule into a Newton solve for half an ulp.
    fctx = _factor_context(ctx)
    t = _stamp_license(f.tri, build(f.tri, Simplex{2}(), degree, fctx(f.tri)))
    l = _stamp_license(f.line, build(f.line, Interval(), degree, fctx(f.line)))
    nt, nl = npoints(t), npoints(l)
    xs = Vector{SVector{3,T}}(undef, nt * nl)
    ws = Vector{T}(undef, nt * nl)
    with_bits(ctx.bits + 32) do
        k = 0
        for i in 1:nt, j in 1:nl
            k += 1
            p = node_vector(t, i)
            xs[k] = SVector{3,T}(finalize_number(ctx, p[1]), finalize_number(ctx, p[2]), finalize_number(ctx, nodes(l)[j]))
            ws[k] = finalize_number(ctx, big(weights(t)[i]) * big(weights(l)[j]))
        end
    end
    claim = PolynomialDegree(min(exactness(t).d, exactness(l).d))
    named = [(n, c) for (n, c) in (("triangle", certificate(t)), ("line", certificate(l))) if c !== nothing]
    cert = combine_certificates(last.(named), join(("$(n): " * c.equations for (n, c) in named), "; "), target_digits(ctx))
    prov = Provenance(family = "WedgeProduct", derivation = derivation(f),
                      path = vcat(["product of a triangle rule and a rule on [-1, 1]",
                                   "triangle: $(describe_family(f.tri)), $(nt) points, degree $(exactness(t).d)",
                                   "line: $(describe_family(f.line)), $(nl) points, degree $(exactness(l).d)"],
                                  ["triangle " * s for s in provenance(t).path]),
                      seed_source = join(unique((provenance(t).seed_source, provenance(l).seed_source)), "; "),
                      citations = unique(vcat(provenance(t).citations, provenance(l).citations)),
                      comparisons = unique(vcat(provenance(t).comparisons, provenance(l).comparisons)),
                      license = join_licenses((provenance(t).license, provenance(l).license)))
    return QuadratureRule(xs, ws, Wedge(), claim, prov, cert)
end
