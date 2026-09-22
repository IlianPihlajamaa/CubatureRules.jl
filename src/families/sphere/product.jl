# Product rules on spheres (PLAN §6 Tier 1, v0.4): the fallback that exists at every degree
# and every precision, as the conical product does on simplices.
#
# On S², with x = (sinθ cosφ, sinθ sinφ, cosθ) and dσ = sinθ dθ dφ, the substitution
# t = cosθ turns the integral into
#
#     ∫_{S²} f dσ = ∫_{-1}^{1} ∫_0^{2π} f(t, φ) dφ dt
#
# which is a Gauss–Legendre rule in t against the trapezoid rule in φ. The trapezoid rule is
# spectrally exact on the circle: m equally spaced points integrate e^{ikφ} exactly for
# |k| < m, so m = d+1 covers every harmonic of degree ≤ d, and n = ⌈(d+1)/2⌉ Gauss points
# cover the t-dependence. That is (d+1)(d+2)/2 points, about twice what a Lebedev rule of
# the same degree needs — the price of a rule that exists at every degree with no table.
#
# On the circle the trapezoid rule alone is the whole story, and there it is optimal: d+1
# points for degree d, which no rule can beat.

"""
    SphereProduct()

Product rules on a sphere: Gauss–Legendre in `cos θ` against the trapezoid rule in `φ` on
[`Sphere`](@ref)`{3}()`, and the trapezoid rule alone on `Sphere{2}()`, the circle.

Available at every degree and every precision, with positive weights throughout. On the
circle it is optimal (`d+1` points for degree `d`); on the sphere it costs
`⌈(d+1)/2⌉ (d+1)` points, roughly twice a minimal rule, which is what it buys for needing
no table.
"""
struct SphereProduct <: RuleFamily end

derivation(::Type{SphereProduct}) = Derived()
describe_family(::SphereProduct) = "SphereProduct"
family_name(::SphereProduct) = "SphereProduct"

candidates(::Type{SphereProduct}, dom::Sphere{D}, ::PolynomialDegree) where {D} =
    (D == 2 || D == 3) && isreference(dom) ? [SphereProduct()] : SphereProduct[]

"Gauss points in cos θ, and equally spaced points in φ, for degree `d` on `S^(D-1)`."
sphere_product_shape(::Sphere{3}, d::Integer) = (cld(d + 1, 2), d + 1)
sphere_product_shape(::Sphere{2}, d::Integer) = (1, d + 1)

npoints(::SphereProduct, dom::Sphere, degree::Integer) = prod(sphere_product_shape(dom, degree))
claimed_degree(::SphereProduct, dom::Sphere, degree) = degree
degree_range(::SphereProduct, dom::Sphere{D}) where {D} = (D == 2 || D == 3) ? (0:typemax(Int)) : (1:0)
properties(::SphereProduct, dom::Sphere, degree) =
    (positive = true, interior = true, symmetry = :none, nested = false)

const ATKINSON_1982 = Citation(key = "Atkinson1982", authors = ["Kendall Atkinson"],
                               title = "Numerical integration on the sphere",
                               journal = "Journal of the Australian Mathematical Society, Series B", year = 1982,
                               volume = "23", pages = "332--347", doi = "10.1017/S0334270000000278")

function build(f::SphereProduct, dom::Sphere{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    isexact(ctx) && throw(ArgumentError("SphereProduct nodes are irrational; $(T) is not supported"))
    (D == 2 || D == 3) || throw(ArgumentError("SphereProduct is defined on the circle and the 2-sphere"))
    n, m = sphere_product_shape(reference(dom), degree)
    guard = gj_guard_bits(max(n, m)) + 8
    xs, ws, cert, path = with_bits(ctx.bits + guard) do
        two_π = 2 * BigFloat(π)
        if D == 2
            φ = [two_π * (k - 1) / m for k in 1:m]
            x = [[cos(φk), sin(φk)] for φk in φ]
            w = fill(two_π / m, m)
            res = maximum(xi -> abs(xi[1]^2 + xi[2]^2 - 1), x)
            x, w, res, ["trapezoid rule, $m equally spaced points"]
        else
            t, wt, _ = gauss_from_recurrence(n, jacobi_recurrence(0, 0), ctx.bits + guard; name = "SphereProduct")
            φ = [two_π * (k - 1) / m for k in 1:m]
            x = Vector{Vector{BigFloat}}(undef, n * m)
            w = Vector{BigFloat}(undef, n * m)
            i = 0
            for j in 1:n
                ρ = sqrt(max(zero(BigFloat), (1 - t[j]) * (1 + t[j])))    # sin θ, without cancellation
                for k in 1:m
                    i += 1
                    x[i] = [ρ * cos(φ[k]), ρ * sin(φ[k]), t[j]]
                    w[i] = wt[j] * two_π / m
                end
            end
            res = maximum(xi -> abs(xi[1]^2 + xi[2]^2 + xi[3]^2 - 1), x)
            x, w, res, ["Gauss–Legendre in cos θ, $n points", "trapezoid rule in φ, $m points"]
        end
    end
    xout = [SVector{D,T}(map(xi -> finalize_number(ctx, xi), x)) for x in xs]
    wout = [finalize_number(ctx, w) for w in ws]
    certificate = Certificate(equations = "|x| = 1 at every node (deviation of the delivered nodes)",
                              residual = BigFloat(cert; precision = 64), residual_bits = ctx.bits + guard,
                              digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = "SphereProduct", derivation = Derived(), path = path,
                      seed_source = "none (closed form × Gauss–Legendre)",
                      citations = [ATKINSON_1982, GOLUB_WELSCH_1969], symmetry = :none)
    return QuadratureRule(xout, wout, Sphere{D}(), PolynomialDegree(degree), prov, certificate)
end
