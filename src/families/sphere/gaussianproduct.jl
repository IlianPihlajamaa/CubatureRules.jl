# Rules on Gaussian-weighted R^D — Stroud's E_n^{r²} (PLAN §6 Tier 1, v0.4).
#
# Stroud tabulates these at fixed degrees, and several of his entries go negative above some
# dimension: the classical 2n²+1 degree-5 rule has a negative axis weight for n ≥ 5, and the
# 2^n-corner variant a negative weight at the origin for n ≥ 4. Rather than ship a patchwork
# of dimension-restricted tables, this family derives the rule, at any degree and dimension,
# with positive weights throughout — from the separation that the weight is radial:
#
#     ∫_{R^D} f e^{-|x|²} dx = ∫₀^∞ r^{D-1} e^{-r²} ∫_{S^{D-1}} f(rω) dσ(ω) dr.
#
# The substitution u = r² turns the radial factor into a generalised Laguerre integral,
#
#     ∫₀^∞ r^{D-1} e^{-r²} g(r) dr = ½ ∫₀^∞ u^{(D-2)/2} e^{-u} g(√u) du,
#
# which is exactly Gauss–Laguerre with α = (D-2)/2, halved.
#
# Half-integer powers of u would spoil that — g(r) = r^d is u^{d/2} — but they never arise.
# A monomial of odd total degree has some odd exponent, and a sphere moment with an odd
# exponent is zero, so the *angular* factor annihilates every odd degree whatever the radial
# rule does. Only even d reaches the radial factor, where u^{d/2} is a polynomial, so ⌈(d+2)/4⌉
# Laguerre points suffice: a degree-11 rule in three dimensions needs three of them.

"""
    GaussianProduct(angular = SphereProduct())

Rules on [`GaussianSpace`](@ref)`(D)`, the whole of `R^D` weighted by `e^{-|x|²}`: a
generalised Gauss–Laguerre rule in the radius against the sphere family `angular` in the
angle. Stroud's `E_n^{r²}`, derived rather than tabulated — available at every degree and
dimension, with positive weights throughout.

The radial factor needs only `⌈(d+2)/4⌉` points, because odd degrees are annihilated by the
angular factor and never reach it.
"""
struct GaussianProduct{F<:RuleFamily} <: CombinatorFamily
    angular::F
end
GaussianProduct() = GaussianProduct(SphereProduct())

derivation(::Type{<:GaussianProduct}) = Derived()
family_name(::GaussianProduct) = "GaussianProduct"
describe_family(f::GaussianProduct) = "GaussianProduct(" * family_name(f.angular) * ")"

function candidates(::Type{<:GaussianProduct}, dom::GaussianDomain{D}, c::PolynomialDegree) where {D}
    (isreference(dom) && D >= 2) || return GaussianProduct[]
    out = GaussianProduct[]
    for F in leaf_families(), angular in candidates(F, Sphere{D}(), c)
        push!(out, GaussianProduct(angular))
    end
    return out
end

"Laguerre points in u = r²: only even degrees reach the radial factor."
gaussian_radial_points(degree) = max(1, cld(degree + 2, 4))

npoints(f::GaussianProduct, dom::GaussianDomain{D}, degree::Integer) where {D} =
    gaussian_radial_points(degree) * npoints(f.angular, Sphere{D}(), degree)
claimed_degree(f::GaussianProduct, dom::GaussianDomain{D}, degree) where {D} =
    min(2 * (2gaussian_radial_points(degree) - 1) + 1, claimed_degree(f.angular, Sphere{D}(), degree))
degree_range(f::GaussianProduct, dom::GaussianDomain{D}) where {D} = degree_range(f.angular, Sphere{D}())
degree_range(::GaussianProduct, dom) = 1:0
properties(f::GaussianProduct, dom::GaussianDomain, degree) =
    (positive = true, interior = true, symmetry = :none, nested = false)
cost_estimate(f::GaussianProduct, dom::GaussianDomain{D}, degree, T) where {D} =
    float(gaussian_radial_points(degree)) * cost_estimate(f.angular, Sphere{D}(), degree, T)

"""
    gaussian_radial_work(D, n, bits) -> (r, w)

The `n`-point rule for the radial weight `r^{D-1} e^{-r²}` on `[0, ∞)`, at precision `bits`:
Gauss–Laguerre with `α = (D-2)/2` in `u = r²`, with nodes `√u` and weights halved.
"""
function gaussian_radial_work(D::Int, n::Int, bits::Int)
    rec = laguerre_recurrence((D - 2) // 2)
    u, w, _ = gauss_from_recurrence(n, rec, bits; name = "GaussianProduct")
    return with_bits(bits) do
        [sqrt(ui) for ui in u], [wi / 2 for wi in w]
    end
end

function build(f::GaussianProduct, dom::GaussianDomain{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    isexact(ctx) && throw(ArgumentError("GaussianProduct nodes are irrational; $(T) is not supported"))
    n = gaussian_radial_points(degree)
    guard = gj_guard_bits(n) + 16
    wbits = ctx.bits + guard
    ang = rule(f.angular, Sphere{D}(); digits = floor(Int, wbits * log10(2)), degree, cancel = ctx.cancel)
    r, wr = gaussian_radial_work(D, n, wbits)
    xs, ws = with_bits(wbits) do
        X = Vector{SVector{D,BigFloat}}()
        W = BigFloat[]
        for i in 1:n, k in 1:npoints(ang)
            push!(X, SVector{D,BigFloat}(r[i] .* BigFloat.(nodes(ang)[k])))
            push!(W, wr[i] * BigFloat(weights(ang)[k]))
        end
        X, W
    end
    xout = [SVector{D,T}(ntuple(j -> finalize_number(ctx, x[j]), D)) for x in xs]
    wout = [finalize_number(ctx, w) for w in ws]
    resid = with_bits(2wbits) do
        abs(sum(BigFloat.(wout)) - measure(GaussianSpace(D)))
    end
    cert = Certificate(equations = "radial Gauss–Laguerre (α = (D-2)/2 in u = r²) × " *
                                   "$(family_name(f.angular)) on the sphere " *
                                   "(deviation of Σw from π^(D/2))",
                       residual = BigFloat(resid; precision = 64), residual_bits = 2wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = certificate(ang).cond, iterations = certificate(ang).iterations)
    prov = Provenance(family = "GaussianProduct", derivation = Derived(),
                      path = ["radial: $(n)-point Gauss–Laguerre, α = $((D - 2)//2), in u = r²",
                              "angular: " * describe_family(f.angular) * ", $(npoints(ang)) points",
                              "odd degrees vanish through the angular factor, so only even ones " *
                              "reach the radial rule"],
                      seed_source = provenance(ang).seed_source,
                      citations = vcat(provenance(ang).citations, [STROUD_1971]), symmetry = :none)
    return QuadratureRule(xout, wout, GaussianSpace(D), PolynomialDegree(claimed_degree(f, dom, degree)),
                          prov, cert)
end
