# Conical (Stroud / collapsed-coordinate) product rules on the D-simplex: Gauss–Jacobi in
# collapsed coordinates. The always-available fallback at any degree (PLAN §6 Tier 1).
#
#     x₁ = u₁,  x_j = u_j Π_{i<j} (1 - u_i),    dx = Π_i (1 - u_i)^{D-i} du
#
# so direction i uses the Gauss–Jacobi weight (1-u)^{D-i} on [0, 1].

"""
    ConicalProduct(inner)

Conical product rule on a simplex or a [`Pyramid`](@ref), built from the 1D family `inner`
(currently [`GaussJacobi`](@ref)): a product rule on a cube carried over by collapsed
coordinates, with Gauss–Jacobi weights that absorb the Jacobian. With `m` points per
direction it has `m^D` points, all interior, all positive, and polynomial degree `2m - 1`.
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

# --- the pyramid ---------------------------------------------------------------------------
#
#     x = ξ (1 − ζ),  y = η (1 − ζ),  z = ζ,     dx dy dz = (1 − ζ)² dξ dη dζ
#
# carries the cube [-1, 1]² × [0, 1] onto the reference pyramid. A monomial xᵃ yᵇ zᶜ becomes
# ξᵃ ηᵇ (1 − ζ)ᵃ⁺ᵇ ζᶜ: degree a and b in ξ and η, and at most a + b + c in ζ with the weight
# (1 − ζ)². So Gauss–Legendre in ξ and η and Gauss–Jacobi(2, 0) in ζ, m points each, are
# exact to total degree 2m − 1.

function candidates(::Type{<:ConicalProduct}, dom::Pyramid, c::PolynomialDegree)
    isreference(dom) || return ConicalProduct[]
    return [ConicalProduct(inner) for F in leaf_families() if conical_compatible(F)
            for inner in candidates(F, Interval(), c)]
end
npoints(::ConicalProduct, dom::Pyramid, degree::Integer) = conical_points(degree)^3

function build(f::ConicalProduct, dom::Pyramid, degree::Int, ctx::BuildContext{T}) where {T}
    f.inner isa GaussJacobi ||
        throw(ArgumentError("ConicalProduct is implemented over GaussJacobi, got $(describe_family(f.inner))"))
    isexact(ctx) && throw(ArgumentError("conical product nodes are irrational; $(T) is not supported"))
    m = conical_points(degree)
    guard = gj_guard_bits(m) + 8
    bits = ctx.bits + guard
    xg, wg, it1 = gauss_jacobi_work(m, 0, 0, bits)
    xz, wz, it2 = gauss_jacobi_work(m, 2, 0, bits)
    xs = SVector{3,T}[]
    ws = T[]
    with_bits(bits) do
        ζ = (1 .+ xz) ./ 2
        ωz = wz ./ 8                                   # ∫₀¹ g (1 − ζ)² dζ = ⅛ ∫ g (1 − u)² du
        for i in 1:m, j in 1:m, k in 1:m
            s = 1 - ζ[k]
            push!(xs, SVector{3,T}(finalize_number(ctx, xg[i] * s), finalize_number(ctx, xg[j] * s),
                                   finalize_number(ctx, ζ[k])))
            push!(ws, finalize_number(ctx, wg[i] * wg[j] * ωz[k]))
        end
    end
    res = max(gj_residual(finalize_number.(Ref(ctx), xg), 0, 0, 2ctx.bits + guard),
              gj_residual(finalize_number.(Ref(ctx), xz), 2, 0, 2ctx.bits + guard))
    cert = Certificate(equations = "P_$m(ξ) = 0 and P_$m^(2,0)(u) = 0 in the collapsed coordinates (Newton correction)",
                       residual = BigFloat(res; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = 1.0, iterations = max(it1, it2))
    prov = Provenance(family = "ConicalProduct", derivation = Derived(),
                      path = ["Gauss–Legendre ($m points) in ξ and η, Gauss–Jacobi(2, 0) ($m points) in ζ",
                              "collapsed onto the pyramid: x = ξ(1 − ζ), y = η(1 − ζ), z = ζ"],
                      seed_source = "Golub–Welsch (Float64) per direction",
                      citations = [STROUD_1971, GOLUB_WELSCH_1969])
    return QuadratureRule(xs, ws, Pyramid(), PolynomialDegree(2m - 1), prov, cert)
end
