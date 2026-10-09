# Rules for the kernel 1/|y − x₀| on a curved (quadratic) triangle, as a weight on the
# reference triangle (notes/singular-bem.md, stage 6).
#
# With x₀ = χ(ξ₀) and the Duffy map ξ = ξ₀ + s v(t), v = a + t (b − a), on each piece of the
# reference triangle cut at ξ₀, χ(ξ) − χ(ξ₀) = s W(s, t) with W = Dχ(ξ₀) v + s Q(v) exactly
# (Q(v) = ½ vᵀ D²χ v, constant for a quadratic map), so
#
#     J(ξ)/|χ(ξ) − χ(ξ₀)| dξ = |a × b| J(ξ₀ + s v)/|W(s, t)| ds dt :
#
# the singularity is gone, and along each ray the weight s ↦ J/|W| is smooth and positive (for
# an element without folds). Its Gauss rule, from the Stieltjes procedure on a Gauss–Legendre
# discretisation, makes the radial direction exact for polynomials in ξ. Across the rays the
# integrand is analytic, behaving like 1/|Dχ(ξ₀) v(t)| as the curvature terms vanish: the sinh
# substitution of q₀(t) = |Dχ(ξ₀) v(t)|² (the flat case's q) takes that out, and Gauss–Legendre
# after it is refined until the Dubiner moments agree between two resolutions, as in DuffySinh.

"""
    DuffyCurved()

Rules for the weakly singular kernel on a curved triangle, the
[`CurvedInverseDistance`](@ref) weight `J(ξ)/|χ(ξ) − χ(ξ₀)|` on the reference triangle,
with `x₀ = χ(ξ₀)` at a vertex, on an edge or inside. The reference triangle is cut at `ξ₀`
and mapped by the Duffy map; each ray gets the Gauss rule of its own weight, so the radial
direction is exact for polynomials in `ξ`, and across the rays Gauss–Legendre after a sinh
substitution is refined until every polynomial of the degree is integrated to the requested
precision. Positive weights, nodes inside the reference triangle; `npoints` gives the
smallest size, which grows with the precision.
"""
struct DuffyCurved <: RuleFamily end

derivation(::Type{DuffyCurved}) = Derived()
family_name(::DuffyCurved) = "DuffyCurved"

candidates(::Type{DuffyCurved}, dom::CurvedKernelDomain, ::PolynomialDegree) =
    isreference(dom.base) && curved_geometry(dom).place !== :near ? [DuffyCurved()] : DuffyCurved[]
npoints(::DuffyCurved, dom, degree::Integer) =
    length(curved_geometry(dom).patches) * duffy_m(degree) * (duffy_m(degree) + 2)
claimed_degree(::DuffyCurved, dom, degree) = 2duffy_m(degree) - 1
degree_range(::DuffyCurved, dom) = 0:typemax(Int)
properties(::DuffyCurved, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

"The `m`-point Gauss rule on `[0, 1]` of `s ↦ J(ξ₀ + s v)/|Dχ(ξ₀) v + s Q(v)|`, to `bits`."
function curved_ray_rule(el, X0, v, Av, m, bits)
    Qv = second_term(el, v)
    discretize(M, wb) = with_bits(wb) do
        x, w = gauss_jacobi_work(M, 0, 0, wb)[1:2]
        s = (1 .+ x) ./ 2
        s, [wi / 2 * area_element(el, X0 + si * v) / norm(Av + si * Qv) for (si, wi) in zip(s, w)]
    end
    return converged_stieltjes(discretize, m, bits; name = "DuffyCurved", failure = (M, gap) ->
        "the Gauss rule of a ray's weight on $(el) did not converge to $bits bits with $M points " *
        "(the recurrence still changed by $gap): is the element folded?")
end

function build(f::DuffyCurved, dom::CurvedKernelDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("nodes for the curved kernel are irrational; $(T) is not supported"))
    isreference(dom.base) || throw(ArgumentError("CurvedInverseDistance is a weight on the reference triangle Simplex{2}()"))
    g = curved_geometry(dom)
    g.place === :near && throw(ArgumentError("DuffyCurved is for ξ₀ on the reference triangle; $(Tuple(dom.weight.ξ0)) is off it"))
    el = dom.weight.element
    m = duffy_m(degree)
    K = 2m - 1
    wbits = ctx.bits + 32
    tol = ldexp(BigFloat(1), -(ctx.bits + 4))
    function assemble(nt)
        xs, ws = SVector{2,BigFloat}[], BigFloat[]
        with_bits(wbits) do
            u, wu = gauss_jacobi_work(nt, 0, 0, wbits)[1:2]
            X0 = SVector{2,BigFloat}(g.x0)
            A = jacobian_matrix(el, X0)
            for (p, q) in g.patches
                a, b = SVector{2,BigFloat}(p) - X0, SVector{2,BigFloat}(q) - X0
                d = b - a
                # q₀(t) = |A v(t)|² = |A d|² ((t − t*)² + η²)
                Aa, Ad = A * a, A * d
                ts = -dot(Aa, Ad) / dot(Ad, Ad)
                η = sqrt(max(dot(Aa, Aa) * dot(Ad, Ad) - dot(Aa, Ad)^2, zero(BigFloat))) / dot(Ad, Ad)
                v0, v1 = asinh(-ts / η), asinh((1 - ts) / η)
                vm, hv = (v0 + v1) / 2, (v1 - v0) / 2
                J0 = abs(a[1] * b[2] - a[2] * b[1])
                for (uk, wk) in zip(u, wu)
                    t = ts + η * sinh(vm + hv * uk)
                    dtdu = hv * η * cosh(vm + hv * uk)
                    v = a + t * d
                    s, ω, _ = curved_ray_rule(el, X0, v, A * v, m, wbits)
                    for (si, ωi) in zip(s, ω)
                        push!(xs, X0 + si * v)
                        push!(ws, wk * dtdu * ωi * J0)
                    end
                end
            end
        end
        return xs, ws
    end
    moments(xs, ws) = with_bits(wbits) do
        dub = DubinerBasis{BigFloat}(K; normalize = true)
        total = zeros(BigFloat, dubiner_length(K))
        for (x, w) in zip(xs, ws)
            dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, x[1], x[2])
            total .+= w .* dub.φ
        end
        total
    end
    # refined until two resolutions agree; the coarser is delivered (see DuffySinh)
    nt = m + 2
    xs, ws = assemble(nt)
    prev = moments(xs, ws)
    gap = BigFloat(Inf)
    for _ in 1:12
        checkcancel(ctx.cancel)
        nt2 = ceil(Int, 1.5nt)
        xs2, ws2 = assemble(nt2)
        cur = moments(xs2, ws2)
        gap = with_bits(wbits) do
            maximum(abs, cur - prev) / max(one(BigFloat), maximum(abs, cur))
        end
        ctx.verbose >= 1 && @info @sprintf("DuffyCurved: %d against %d rays per piece, moments differ by %.1e", nt, nt2, Float64(gap))
        gap <= tol && break
        nt, xs, ws, prev = nt2, xs2, ws2, cur
    end
    gap <= tol || throw(RefinementError("DuffyCurved", "the angular resolution did not converge to $(ctx.bits) bits"))
    nodes = [SVector{2,T}(finalize_number(ctx, x[1]), finalize_number(ctx, x[2])) for x in xs]
    weights = [finalize_number(ctx, w) for w in ws]
    cert = Certificate(equations = "integrals of the orthonormal Dubiner polynomials to degree $K against the weight, " *
                                   "converged between two angular resolutions",
                       residual = BigFloat(gap; precision = 64), residual_bits = wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
    prov = Provenance(family = "DuffyCurved", derivation = Derived(),
                      path = ["kernel: 1/|y − x₀| on $(el), x₀ = χ$(Tuple(dom.weight.ξ0)) " *
                              (g.place === :vertex ? "a vertex" : g.place === :edge ? "on an edge" : "inside") *
                              ", as the weight J(ξ)/|χ(ξ) − χ(ξ₀)| on the reference triangle: " *
                              "$(length(g.patches)) piece$(length(g.patches) == 1 ? "" : "s") with a vertex at ξ₀",
                              "Duffy map ξ = ξ₀ + s (a + t (b − a)) on each: χ(ξ) − χ(ξ₀) = s (Dχ(ξ₀) v + s ½ vᵀD²χ v)",
                              "s: $m-point Gauss rule of each ray's weight J/|W|, by the Stieltjes procedure",
                              "t: $nt-point Gauss–Legendre after the sinh substitution of |Dχ(ξ₀) v(t)|², refined until " *
                              "the Dubiner moments to degree $K agreed to $(ctx.bits) bits"],
                      seed_source = "none (derived)", citations = [DUFFY_1982, GAUTSCHI_1982, JOHNSTON_ELLIOTT_2005])
    return QuadratureRule(nodes, weights, dom, PolynomialDegree(K), prov, cert)
end
