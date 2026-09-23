# Product rules on spheres (PLAN §6 Tier 1, v0.4): the fallback that exists at every degree,
# every dimension and every precision, as the conical product does on simplices.
#
# In hyperspherical coordinates on S^{D-1} ⊂ R^D,
#
#     x₁ = cos θ₁,  x₂ = sin θ₁ cos θ₂,  …,  x_{D-1} = sin θ₁ ⋯ sin θ_{D-2} cos φ,
#     x_D = sin θ₁ ⋯ sin θ_{D-2} sin φ,      dσ = ∏_{k=1}^{D-2} sin^{D-1-k}(θ_k) dθ_k dφ,
#
# and the substitution t_k = cos θ_k turns each polar factor into a Jacobi weight:
#
#     sin^{D-1-k}(θ_k) dθ_k  =  (1 - t_k²)^{(D-2-k)/2} dt_k,
#
# which is Gauss–Jacobi with α = β = (D-2-k)/2. The azimuth keeps the trapezoid rule, which
# is spectrally exact on a circle: m equally spaced points integrate e^{ikφ} exactly for
# |k| < m. So degree d costs ⌈(d+1)/2⌉^{D-2} (d+1) points.
#
# Two dimensions fall out of the same formula rather than being special cases: on the circle
# there are no polar factors and the trapezoid rule alone is optimal, d+1 points for degree
# d; on S² there is one polar factor with α = β = 0, which is Gauss–Legendre in cos θ.

"""
    SphereProduct()

Product rules on a [`Sphere`](@ref): Gauss–Jacobi in each polar angle against the trapezoid
rule in the azimuth. Available at every degree, in every dimension, at every precision, with
positive weights throughout.

Costs `⌈(d+1)/2⌉^{D-2} (d+1)` points on `Sphere{D}`, so it is the fallback rather than the
economical choice — on `Sphere{3}` roughly twice a Lebedev rule of the same degree, and it
grows with dimension. On the circle it is optimal: `d+1` points for degree `d`.
"""
struct SphereProduct <: RuleFamily end

derivation(::Type{SphereProduct}) = Derived()
describe_family(::SphereProduct) = "SphereProduct"
family_name(::SphereProduct) = "SphereProduct"

candidates(::Type{SphereProduct}, dom::Sphere{D}, ::PolynomialDegree) where {D} =
    (D >= 2 && isreference(dom)) ? [SphereProduct()] : SphereProduct[]

"""
    sphere_product_shape(dom, d) -> (n, m)

Gauss–Jacobi points per polar angle and equally spaced points in the azimuth, for degree `d`
on `S^{D-1}`. There are `D-2` polar angles, so the node count is `n^{D-2} m`.
"""
sphere_product_shape(::Sphere{D}, d::Integer) where {D} = (max(1, cld(d + 1, 2)), max(1, d + 1))

function npoints(::SphereProduct, dom::Sphere{D}, degree::Integer) where {D}
    n, m = sphere_product_shape(dom, degree)
    return n^(D - 2) * m
end
claimed_degree(::SphereProduct, dom::Sphere, degree) = degree
degree_range(::SphereProduct, dom::Sphere{D}) where {D} = D >= 2 ? (0:typemax(Int)) : (1:0)
degree_range(::SphereProduct, dom) = 1:0
properties(::SphereProduct, dom::Sphere, degree) =
    (positive = true, interior = true, symmetry = :none, nested = false)

const ATKINSON_1982 = Citation(key = "Atkinson1982", authors = ["Kendall Atkinson"],
                               title = "Numerical integration on the sphere",
                               journal = "Journal of the Australian Mathematical Society, Series B", year = 1982,
                               volume = "23", pages = "332--347", doi = "10.1017/S0334270000000278")

"""
    sphere_product_work(D, n, m, bits) -> (nodes, weights)

The product rule on `S^{D-1}` with `n` Gauss–Jacobi points in each of the `D-2` polar angles
and `m` equally spaced azimuths, at precision `bits`. Nodes are `Vector{BigFloat}` of length
`D`; the weights sum to the sphere's area.
"""
function sphere_product_work(D::Int, n::Int, m::Int, bits::Int)
    # one Gauss–Jacobi rule per polar angle: α = β = (D-2-k)/2 for k = 1 … D-2
    polar = [gauss_jacobi_work(n, (D - 2 - k) // 2, (D - 2 - k) // 2, bits) for k in 1:(D - 2)]
    return with_bits(bits) do
        two_π = 2 * BigFloat(π)
        φ = [two_π * (j - 1) / m for j in 1:m]
        X = Vector{Vector{BigFloat}}()
        W = BigFloat[]
        # every combination of one node per polar angle, times every azimuth
        for idx in Iterators.product(ntuple(_ -> 1:n, D - 2)...), j in 1:m
            x = Vector{BigFloat}(undef, D)
            w = two_π / m
            radial = one(BigFloat)          # ∏ sin θ accumulated so far
            for k in 1:(D - 2)
                t = polar[k][1][idx[k]]
                w *= polar[k][2][idx[k]]
                x[k] = radial * t
                radial *= sqrt(max(zero(BigFloat), (1 - t) * (1 + t)))   # sin θ, no cancellation
            end
            x[D - 1] = radial * cos(φ[j])
            x[D] = radial * sin(φ[j])
            push!(X, x)
            push!(W, w)
        end
        X, W
    end
end

function build(f::SphereProduct, dom::Sphere{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    isexact(ctx) && throw(ArgumentError("SphereProduct nodes are irrational; $(T) is not supported"))
    D >= 2 || throw(ArgumentError("SphereProduct needs at least the circle, Sphere{2}"))
    n, m = sphere_product_shape(reference(dom), degree)
    guard = gj_guard_bits(max(n, m)) + 8
    xs, ws = sphere_product_work(D, n, m, ctx.bits + guard)
    res = with_bits(ctx.bits + guard) do
        maximum(x -> abs(sum(abs2, x) - 1), xs)
    end
    xout = [SVector{D,T}(ntuple(j -> finalize_number(ctx, x[j]), D)) for x in xs]
    wout = [finalize_number(ctx, w) for w in ws]
    path = D == 2 ? ["trapezoid rule, $m equally spaced points"] :
           ["Gauss–Jacobi in each of the $(D - 2) polar angles, $n points each",
            "trapezoid rule in the azimuth, $m points"]
    cert = Certificate(equations = "|x| = 1 at every node (deviation of the delivered nodes)",
                       residual = BigFloat(res; precision = 64), residual_bits = ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = "SphereProduct", derivation = Derived(), path = path,
                      seed_source = "none (closed form × Gauss–Jacobi)",
                      citations = [ATKINSON_1982, GOLUB_WELSCH_1969], symmetry = :none)
    return QuadratureRule(xout, wout, Sphere{D}(), PolynomialDegree(degree), prov, cert)
end
