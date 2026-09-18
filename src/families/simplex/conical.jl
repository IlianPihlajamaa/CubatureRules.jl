# Conical (Stroud / collapsed-coordinate) product rules on the D-simplex: Gauss–Jacobi in
# collapsed coordinates. The always-available fallback at any degree (PLAN §6 Tier 1).
#
#     x₁ = u₁,  x_j = u_j Π_{i<j} (1 - u_i),    dx = Π_i (1 - u_i)^{D-i} du
#
# so direction i uses the Gauss–Jacobi weight (1-u)^{D-i} on [0, 1].

"""
    ConicalProduct(inner)

Conical product rule on a simplex built from the 1D family `inner` (currently
[`GaussJacobi`](@ref)). With `m` points per direction it has `m^D` points, all interior,
all positive, and polynomial degree `2m - 1`.
"""
struct ConicalProduct{F<:RuleFamily} <: CombinatorFamily
    inner::F
end
ConicalProduct() = ConicalProduct(GaussJacobi())

derivation(::Type{<:ConicalProduct}) = Derived()
family_name(::ConicalProduct) = "ConicalProduct"
describe_family(f::ConicalProduct) = "ConicalProduct(" * family_name(f.inner) * ")"

# Combinators recurse into the leaf families (PLAN §2.5), to nesting depth 1.
function candidates(::Type{<:ConicalProduct}, dom::Simplex, c::PolynomialDegree)
    isreference(dom) || return ConicalProduct[]
    out = ConicalProduct[]
    for F in leaf_families()
        conical_compatible(F) || continue
        for inner in candidates(F, Interval(), c)
            push!(out, ConicalProduct(inner))
        end
    end
    return out
end

conical_points(degree) = max(1, cld(degree + 1, 2))
npoints(::ConicalProduct, dom::Simplex{D}, degree::Integer) where {D} = conical_points(degree)^D
claimed_degree(::ConicalProduct, dom, degree) = 2conical_points(degree) - 1
degree_range(::ConicalProduct, dom) = 0:typemax(Int)
properties(::ConicalProduct, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

const STROUD_1971 = Citation(key = "Stroud1971", authors = ["Arthur H. Stroud"],
                             title = "Approximate Calculation of Multiple Integrals",
                             journal = "Prentice-Hall", year = 1971)

"""
    conical_work(D, m, bits) -> (nodes, weights, iterations, onedim_nodes)

Conical product nodes (as `Vector{Vector{BigFloat}}`) and weights on the reference
`D`-simplex with `m` Gauss–Jacobi points per direction, at precision `bits`.
"""
function conical_work(D::Int, m::Int, bits::Int)
    oned = [gauss_jacobi_work(m, D - i, 0, bits) for i in 1:D]
    iters = maximum(t -> t[3], oned)
    with_bits(bits) do
        us = [(1 .+ t[1]) ./ 2 for t in oned]
        ws = [t[2] ./ big(2)^(D - i + 1) for (i, t) in enumerate(oned)]
        X = Vector{Vector{BigFloat}}()
        W = BigFloat[]
        for I in CartesianIndices(ntuple(_ -> m, D))
            x = Vector{BigFloat}(undef, D)
            scale = one(BigFloat)
            w = one(BigFloat)
            for j in 1:D
                u = us[j][I[j]]
                x[j] = u * scale
                scale *= 1 - u
                w *= ws[j][I[j]]
            end
            push!(X, x)
            push!(W, w)
        end
        X, W, iters, [t[1] for t in oned]
    end
end

function build(f::ConicalProduct, dom::Simplex{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    f.inner isa GaussJacobi ||
        throw(ArgumentError("ConicalProduct is implemented over GaussJacobi, got $(describe_family(f.inner))"))
    isexact(ctx) && throw(ArgumentError("conical product nodes are irrational; $(T) is not supported"))
    m = conical_points(degree)
    guard = gj_guard_bits(m) + 8
    X, W, iters, u1d = conical_work(D, m, ctx.bits + guard)
    xs = [SVector{D,T}(ntuple(j -> finalize_number(ctx, x[j]), D)) for x in X]
    ws = [finalize_number(ctx, w) for w in W]
    res = maximum(i -> gj_residual(finalize_number.(Ref(ctx), u1d[i]),
                                   D - i, 0, 2ctx.bits + guard), 1:D)
    cert = Certificate(equations = "P_$m^(D-i,0)(u) = 0 in each collapsed direction (Newton correction)",
                       residual = BigFloat(res; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = 1.0, iterations = iters)
    prov = Provenance(family = "ConicalProduct", derivation = Derived(),
                      path = ["Gauss–Jacobi ($m points) in each of $D collapsed coordinates",
                              "tensor product mapped by the Duffy/Stroud collapse"],
                      seed_source = "Golub–Welsch (Float64) per direction",
                      citations = [STROUD_1971, GOLUB_WELSCH_1969])
    return QuadratureRule(xs, ws, Simplex{D}(), PolynomialDegree(2m - 1), prov, cert)
end
