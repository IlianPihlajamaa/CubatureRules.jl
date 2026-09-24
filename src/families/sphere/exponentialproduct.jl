# Rules on R^D with the weight e^{-|x|} — Stroud's E_n^r (PLAN §6 Tier 1, v0.5).
#
# The same separation as for the Gaussian weight (see gaussianproduct.jl):
#
#     ∫_{R^D} f e^{-|x|} dx = ∫₀^∞ r^{D-1} e^{-r} ∫_{S^{D-1}} f(rω) dσ(ω) dr,
#
# and the same observation that odd degrees are annihilated by the angular factor, so only
# even powers of r reach the radial one. With u = r² the radial factor becomes
#
#     ∫₀^∞ r^{D-1} e^{-r} g(r) dr = ½ ∫₀^∞ u^{(D-2)/2} e^{-√u} g(√u) du,
#
# and a Gauss rule in u needs half as many points as one in r — ⌈(d+2)/4⌉ for degree d, as for
# the Gaussian weight. The difference is that u^{(D-2)/2} e^{-√u} has no classical family. It
# does have exact moments,
#
#     ∫₀^∞ u^{j + (D-2)/2} e^{-√u} du = 2 Γ(2j + D) = 2 (2j + D - 1)!,
#
# so its modified moments against the Laguerre polynomials of the same α are computed exactly,
# in rational arithmetic, and the radial rule comes from ModifiedChebyshev. The alternative —
# Gauss–Laguerre in r with α = D - 1 — needs no moments but twice the radial points.

"""
    ExponentialProduct(angular = SphereProduct())

Rules on [`ExponentialSpace`](@ref)`(D)`, the whole of `R^D` weighted by `e^{-|x|}`: a radial
Gauss rule against the sphere family `angular`. Stroud's `E_n^r`, derived rather than
tabulated — available at every degree and dimension `D ≥ 2`, with positive weights.

The radial factor needs only `⌈(d+2)/4⌉` points: odd degrees vanish through the angular
factor, and in `u = r²` the even ones are polynomials, integrated by a Gauss rule for
`u^{(D-2)/2} e^{-√u}`. That weight has no classical family; its rule is built from exact
modified moments by [`ModifiedChebyshev`](@ref).
"""
struct ExponentialProduct{F<:RuleFamily} <: CombinatorFamily
    angular::F
end
ExponentialProduct() = ExponentialProduct(SphereProduct())

derivation(::Type{<:ExponentialProduct}) = Derived()
family_name(::ExponentialProduct) = "ExponentialProduct"
describe_family(f::ExponentialProduct) = "ExponentialProduct(" * family_name(f.angular) * ")"

function candidates(::Type{<:ExponentialProduct}, dom::ExponentialDomain{D}, c::PolynomialDegree) where {D}
    (isreference(dom) && D >= 2) || return ExponentialProduct[]
    out = ExponentialProduct[]
    for F in leaf_families(), angular in candidates(F, Sphere{D}(), c)
        push!(out, ExponentialProduct(angular))
    end
    return out
end

npoints(f::ExponentialProduct, dom::ExponentialDomain{D}, degree::Integer) where {D} =
    gaussian_radial_points(degree) * npoints(f.angular, Sphere{D}(), degree)
claimed_degree(f::ExponentialProduct, dom::ExponentialDomain{D}, degree) where {D} =
    min(2 * (2gaussian_radial_points(degree) - 1) + 1, claimed_degree(f.angular, Sphere{D}(), degree))
degree_range(f::ExponentialProduct, dom::ExponentialDomain{D}) where {D} = degree_range(f.angular, Sphere{D}())
degree_range(::ExponentialProduct, dom) = 1:0
properties(f::ExponentialProduct, dom::ExponentialDomain, degree) =
    (positive = true, interior = true, symmetry = :none, nested = false)
cost_estimate(f::ExponentialProduct, dom::ExponentialDomain{D}, degree, T) where {D} =
    float(gaussian_radial_points(degree)) * cost_estimate(f.angular, Sphere{D}(), degree, T)

"Generalised binomial coefficient `C(x, m)` for rational `x`, as an exact product."
gbinomial(x, m::Integer) = m == 0 ? one(x) : prod((x - m + i) // i for i in 1:m)

"""
    exponential_radial_weight(D) -> MomentWeight

The radial weight of `e^{-|x|}` on `R^D` after `u = r²`: `u^β e^{-√u}` on `[0, ∞)` with
`β = (D-2)/2`, described by its exact modified moments against the monic Laguerre polynomials
`πₖ(u) = Σⱼ (-1)^(k+j) k!/j! C(k+β, k-j) uʲ` of the same `β`.
"""
function exponential_radial_weight(D::Integer)
    β = Rational{BigInt}(D - 2, 2)
    aux = MonicRecurrence((k, T) -> T(2k + 1 + β), (k, T) -> k == 0 ? zero(T) : T(k * (k + β)))
    moment(k) = sum((-1)^(k + j) * factorial(big(k)) // factorial(big(j)) * gbinomial(k + β, k - j) *
                    2 * factorial(big(2j + D - 1)) for j in 0:k)
    return MomentWeight(aux, (k, T) -> T(moment(k)); label = "u^$β e^-√u on [0, ∞)", support = (0, Inf))
end

"""
    exponential_radial_work(D, n, bits) -> (r, w, cond)

The `n`-point rule for the radial weight `r^{D-1} e^{-r}` on `[0, ∞)`, exact for even powers
up to `r^{4n-2}`: the Gauss rule for `u^{(D-2)/2} e^{-√u}` with nodes `√u` and weights halved.
"""
function exponential_radial_work(D::Int, n::Int, bits::Int)
    u, w, _, _, _, cond = modified_chebyshev_work(exponential_radial_weight(D), n, bits)
    r, wr = with_bits(bits) do
        [sqrt(ui) for ui in u], [wi / 2 for wi in w]
    end
    return r, wr, cond
end

function build(f::ExponentialProduct, dom::ExponentialDomain{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    isexact(ctx) && throw(ArgumentError("ExponentialProduct nodes are irrational; $(T) is not supported"))
    n = gaussian_radial_points(degree)
    guard = gj_guard_bits(n) + 16
    wbits = ctx.bits + guard
    ang = rule(f.angular, Sphere{D}(); digits = floor(Int, wbits * log10(2)), degree, cancel = ctx.cancel)
    r, wr, rcond = exponential_radial_work(D, n, wbits)
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
        abs(sum(BigFloat.(wout)) - measure(ExponentialSpace(D)))
    end
    cert = Certificate(equations = "radial Gauss rule for u^((D-2)/2) e^(-√u) in u = r², from exact " *
                                   "modified moments × $(family_name(f.angular)) on the sphere " *
                                   "(deviation of Σw from |S^(D-1)| (D-1)!)",
                       residual = BigFloat(resid; precision = 64), residual_bits = 2wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = max(rcond, certificate(ang).cond), iterations = certificate(ang).iterations)
    prov = Provenance(family = "ExponentialProduct", derivation = Derived(),
                      path = ["radial: $(n)-point Gauss rule for u^$((D - 2)//2) e^-√u in u = r², by " *
                              "Wheeler's algorithm from exact modified moments against Laguerre",
                              "angular: " * describe_family(f.angular) * ", $(npoints(ang)) points",
                              "odd degrees vanish through the angular factor, so only even ones " *
                              "reach the radial rule"],
                      seed_source = provenance(ang).seed_source,
                      citations = vcat(provenance(ang).citations, [STROUD_1971, WHEELER_1974]), symmetry = :none)
    return QuadratureRule(xout, wout, ExponentialSpace(D), PolynomialDegree(claimed_degree(f, dom, degree)),
                          prov, cert)
end
