# Product rules on a ball: a radial rule against a rule on the
# sphere, which is the whole content of the separation
#
#     ∫_{B^D} f dx = ∫₀¹ r^{D-1} ∫_{S^{D-1}} f(rω) dσ(ω) dr.
#
# A polynomial of degree d is Σ_k r^k p_k(ω) with p_k homogeneous of degree k, and a sphere
# moment of odd degree vanishes, so an angular rule exact to degree d annihilates every odd k
# and only the even powers r^k = u^{k/2}, k ≤ d, reach the radial factor. With u = r²,
#
#     ∫₀¹ r^{D-1} g(r) dr = ½ ∫₀¹ u^{(D-2)/2} g(√u) du,
#
# a Jacobi weight in u, whose n-point Gauss rule is exact to degree 2n − 1 in u, so
# ⌈(d+1)/4⌉ radial points suffice, as for the Gaussian weight (gaussianproduct.jl), where a
# rule in r would need ⌈(d+1)/2⌉.
#
# The angular factor is any sphere family, which is why this is a combinator: with Lebedev
# inside it inherits Lebedev's node counts, and the selector ranks the combinations for you.

"""
    BallProduct(angular = SphereProduct())

Product rules on a [`Ball`](@ref): a Gauss–Jacobi rule in the radius against the sphere
family `angular` on the surface. The radial rule is in `u = r²`: odd degrees vanish through
the angular factor, so `n = ⌈(d+1)/4⌉` radial points suffice, and the rule has `n` times the
angular rule's node count, all interior, all positive.

`BallProduct(LebedevRule())` is the cheap one wherever a Lebedev seed exists; the selector
offers both and ranks them by node count.
"""
struct BallProduct{F<:RuleFamily} <: CombinatorFamily
    angular::F
end
BallProduct() = BallProduct(SphereProduct())

derivation(::Type{<:BallProduct}) = Derived()
family_name(::BallProduct) = "BallProduct"
describe_family(f::BallProduct) = "BallProduct(" * describe_family(f.angular) * ")"
inner_families(f::BallProduct) = (f.angular,)

# Combinators recurse into the families that live on the matching sphere.
function candidates(::Type{<:BallProduct}, dom::Ball{D}, c::PolynomialDegree) where {D}
    isreference(dom) || return BallProduct[]
    out = BallProduct[]
    for F in leaf_families(), angular in candidates(F, Sphere{D}(), c)
        push!(out, BallProduct(angular))
    end
    return out
end

"Gauss–Jacobi points in u = r²: only even degrees reach the radial factor."
ball_radial_points(degree) = max(1, cld(degree + 1, 4))
npoints(f::BallProduct, dom::Ball{D}, degree::Integer) where {D} =
    ball_radial_points(degree) * npoints(f.angular, Sphere{D}(), degree)
claimed_degree(f::BallProduct, dom::Ball{D}, degree) where {D} =
    min(2 * (2ball_radial_points(degree) - 1) + 1, claimed_degree(f.angular, Sphere{D}(), degree))
degree_range(f::BallProduct, dom::Ball{D}) where {D} = degree_range(f.angular, Sphere{D}())
degree_range(::BallProduct, dom) = 1:0
properties(f::BallProduct, dom::Ball{D}, degree) where {D} =
    (positive = true, interior = true, symmetry = NoSymmetry(), nested = false)
cost_estimate(f::BallProduct, dom::Ball{D}, degree, T) where {D} =
    float(ball_radial_points(degree)) * cost_estimate(f.angular, Sphere{D}(), degree, T)

"""
    ball_radial_work(D, n, bits) -> (r, w)

The `n`-point rule for the radial weight `r^{D-1}` on `[0, 1]`, at precision `bits`, in
`u = r²`: `∫₀¹ r^{D-1} g(r) dr = Σ wᵢ g(rᵢ)` for every even polynomial `g` of degree
`≤ 4n - 2`. It is the Gauss–Jacobi rule with `α = 0`, `β = (D-2)/2` mapped by
`u = (1+t)/2`, with nodes `√u` and weights scaled by `2^{-β-2}`.
"""
function ball_radial_work(D::Int, n::Int, bits::Int)
    # an integer exponent for even D; for odd D a half-integer, which as a BigFloat is exact,
    # where a Rational would send the Gamma functions of the weight's mass through Float64
    β = iseven(D) ? (D - 2) ÷ 2 : with_bits(() -> BigFloat(D - 2) / 2, bits)
    t, v, _ = gauss_jacobi_work(n, 0, β, bits)
    return with_bits(bits) do
        r = [sqrt((1 + ti) / 2) for ti in t]
        scale = BigFloat(2)^(-(β + 2))
        w = [vi * scale for vi in v]
        r, w
    end
end

function build(f::BallProduct, dom::Ball{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    isexact(ctx) && throw(ArgumentError("BallProduct nodes are irrational; $(T) is not supported"))
    n = ball_radial_points(degree)
    guard = gj_guard_bits(n) + 16
    wbits = ctx.bits + guard
    # the angular rule is built by the ordinary machinery, at the working precision
    ang = rule(f.angular, Sphere{D}(); digits = floor(Int, wbits * log10(2)), degree,
               cancel = ctx.cancel)
    r, wr = ball_radial_work(D, n, wbits)
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
    # the defining equation of the construction: the weights integrate the volume, and the
    # radial nodes satisfy the Jacobi polynomial they came from
    resid = with_bits(2wbits) do
        abs(sum(BigFloat.(wout)) - measure(Ball{D}()))
    end
    cert = Certificate(equations = "radial Gauss–Jacobi (β = (D-2)/2 in u = r²) × $(family_name(f.angular)) " *
                                   "on the sphere (deviation of Σw from the volume)",
                       residual = BigFloat(resid; precision = 64), residual_bits = 2wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = certificate(ang).cond, iterations = certificate(ang).iterations)
    prov = Provenance(family = "BallProduct", derivation = Derived(),
                      path = ["radial: $(n)-point Gauss–Jacobi in u = r² for r^$(D-1) on [0,1]",
                              "angular: " * describe_family(f.angular) * ", $(npoints(ang)) points",
                              "claim: min(4n-1, angular degree) = $(claimed_degree(f, dom, degree))"],
                      seed_source = provenance(ang).seed_source,
                      citations = vcat(provenance(ang).citations, [STROUD_1971]),
                      comparisons = provenance(ang).comparisons, symmetry = NoSymmetry(),
                      license = provenance(ang).license)
    return QuadratureRule(xout, wout, Ball{D}(), PolynomialDegree(claimed_degree(f, dom, degree)), prov, cert)
end
