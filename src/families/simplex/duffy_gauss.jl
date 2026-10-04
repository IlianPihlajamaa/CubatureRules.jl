# Rules for the weakly singular kernel 1/|y − x₀| on a triangle (notes/singular-bem.md).
#
# On a sub-triangle (x₀, p, q) with a = p − x₀ and b = q − x₀, the Duffy map
#
#     y = x₀ + s (a + t (b − a)),   (s, t) ∈ [0, 1]²,   dy = s |a × b| ds dt
#
# gives |y − x₀| = s √q(t), q(t) = |a + t (b − a)|², so f(y)/|y − x₀| dy = f · |a × b| / √q(t)
# ds dt: the s cancels, and a polynomial f of degree d is a polynomial of degree ≤ d in s and
# in t. Gauss–Legendre in s times the Gauss rule of the weight 1/√q(t) is therefore exact.
#
# That weight comes from the Stieltjes procedure on a discretisation of it. q has the complex
# roots t* ± iη, and with t = t* + η sinh u it disappears altogether, dt/√q = du/|b − a|, so
# Gauss–Legendre in u is a discretisation that converges geometrically however thin the
# sub-triangle, on an interval in u that grows only like log(1/η).

"""
    DuffyGauss()

Rules for the kernel `1/|y − x₀|` on a triangle, an [`InverseDistance`](@ref) weight. The
triangle is cut at `x₀` into one, two or three sub-triangles with a vertex at `x₀` (as `x₀` is
a vertex, on an edge or inside). On each, the Duffy map `y = x₀ + s (a + t (b − a))` cancels
the singularity and leaves the weight `1/√q(t)`, `q(t) = |a + t (b − a)|²`, in the collapsed
direction; the rule is `m`-point Gauss–Legendre in `s` times the `m`-point Gauss rule of that
weight in `t`, with `m = ⌈(d + 1)/2⌉`. It integrates `f(y)/|y − x₀|` exactly for every
polynomial `f` of degree `2m − 1`, with `m²` points per sub-triangle, positive weights and
interior nodes. The weight depends on the shape of each sub-triangle, so every rule is built
for its triangle and point; nothing is tabulated.

```julia
dom = WeightedDomain(Simplex((0, 0), (1, 0), (0, 1)), InverseDistance((0, 0)))
r = rule(dom; degree = 9)                # 25 points
integrate(y -> 1.0, r)                   # √2 log(1 + √2)
```
"""
struct DuffyGauss <: RuleFamily end

derivation(::Type{DuffyGauss}) = Derived()
family_name(::DuffyGauss) = "DuffyGauss"

const DUFFY_1982 = Citation(key = "Duffy1982", authors = ["M. G. Duffy"],
                            title = "Quadrature over a pyramid or cube of integrands with a singularity at a vertex",
                            journal = "SIAM Journal on Numerical Analysis", year = 1982, volume = "19",
                            pages = "1260–1262", doi = "10.1137/0719090")

duffy_m(degree) = max(1, cld(degree + 1, 2))
candidates(::Type{DuffyGauss}, dom::InverseDistanceDomain, ::PolynomialDegree) = [DuffyGauss()]
npoints(::DuffyGauss, dom, degree::Integer) = npatches(dom) * duffy_m(degree)^2
claimed_degree(::DuffyGauss, dom, degree) = 2duffy_m(degree) - 1
degree_range(::DuffyGauss, dom) = 0:typemax(Int)
degree_for_npoints(::DuffyGauss, dom, n::Integer) = 2isqrt(n ÷ npatches(dom)) - 1
properties(::DuffyGauss, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

"""
    inverse_sqrt_quadratic_rule(a, b, m, bits, gl) -> (t, λ, M, gap, usedbits)

The `m`-point Gauss rule on `[0, 1]` of the weight `1/|a + t (b − a)|`, accurate to `bits`,
from the Stieltjes procedure on the sinh-substituted Gauss–Legendre discretisation. `gl(M,
bits)` returns the `M`-point Gauss–Legendre rule on `[−1, 1]`, so that callers can share it.
"""
function inverse_sqrt_quadratic_rule(a, b, m::Integer, bits::Integer, gl)
    function discretize(M, wbits)
        with_bits(wbits) do
            A, B = SVector{2,BigFloat}(a), SVector{2,BigFloat}(b)
            d = B - A
            dd = dot(d, d)
            ts = -dot(A, d) / dd                       # q(t) = |d|² ((t − t*)² + η²)
            η = abs(_cross(A, B)) / dd
            u0, u1 = asinh(-ts / η), asinh((1 - ts) / η)
            x, w = gl(M, wbits)
            c, h = (u0 + u1) / 2, (u1 - u0) / 2
            [ts + η * sinh(c + h * xi) for xi in x], [h * wi / sqrt(dd) for wi in w]
        end
    end
    return converged_stieltjes(discretize, m, bits; name = "DuffyGauss", failure = (M, gap) ->
        "the Gauss rule of 1/√q(t) on a sub-triangle with edge vectors $(Float64.(a)), $(Float64.(b)) did not " *
        "converge to $bits bits with $M points (the recurrence still changed by $gap)")
end

function build(f::DuffyGauss, dom::InverseDistanceDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("nodes for the kernel 1/|y − x₀| are irrational; $(T) is not supported"))
    patches, place = singular_patches(dom)
    m = duffy_m(degree)
    x0 = SVector{2,Rational{BigInt}}(_exact(dom.weight.x0[1]), _exact(dom.weight.x0[2]))
    ctx.verbose >= 1 && @info @sprintf("DuffyGauss: %d sub-triangle%s, %d × %d points each, target %d bits",
                                       length(patches), length(patches) == 1 ? "" : "s", m, m, ctx.bits)
    # Gauss–Legendre rules on [−1, 1], shared by the s-direction and every discretisation
    cache = Dict{Tuple{Int,Int},Tuple{Vector{BigFloat},Vector{BigFloat}}}()
    gl(M, bits) = get!(() -> gauss_jacobi_work(M, 0, 0, bits)[1:2], cache, (M, bits))
    trules = [inverse_sqrt_quadratic_rule(p - x0, q - x0, m, ctx.bits + 8, gl) for (p, q) in patches]
    wbits = maximum(r[5] for r in trules)
    xs, ws = SVector{2,T}[], T[]
    with_bits(wbits) do
        sx, sw = gl(m, wbits)
        X0 = SVector{2,BigFloat}(x0)
        for ((p, q), (t, λ, _, _, _)) in zip(patches, trules)
            a, b = SVector{2,BigFloat}(p - x0), SVector{2,BigFloat}(q - x0)
            J = abs(_cross(a, b))
            for (si, swi) in zip(sx, sw), (tj, λj) in zip(t, λ)
                s, w = (1 + si) / 2, swi / 2
                y = X0 + s * (a + tj * (b - a))
                push!(xs, SVector{2,T}(finalize_number(ctx, y[1]), finalize_number(ctx, y[2])))
                push!(ws, finalize_number(ctx, w * λj * J))
            end
        end
    end
    Ms = [r[3] for r in trules]
    gap = maximum(r[4] for r in trules)
    cert = Certificate(equations = "recurrence of the discretised weight 1/√q(t) on each sub-triangle: " *
                                   "$(join(Ms, ", ")) against half as many points",
                       residual = BigFloat(gap; precision = 64), residual_bits = wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
    prov = Provenance(family = "DuffyGauss", derivation = Derived(),
                      path = ["kernel: 1/|y − x₀| on $(dom.base), x₀ = $(Tuple(dom.weight.x0)) " *
                              (place === :vertex ? "a vertex" : place === :edge ? "on an edge" : "inside") *
                              ": $(length(patches)) sub-triangle$(length(patches) == 1 ? "" : "s") with a vertex at x₀",
                              "Duffy map y = x₀ + s (a + t (b − a)) on each: f(y)/|y − x₀| dy = f |a × b| / √q(t) ds dt",
                              "s: $m-point Gauss–Legendre on [0, 1]",
                              "t: $m-point Gauss rule of 1/√q(t), from the Stieltjes procedure on a sinh-substituted " *
                              "Gauss–Legendre discretisation ($(join(Ms, ", ")) points), doubled until the recurrence " *
                              "agreed to $(ctx.bits + 8) bits; nodes by Newton on the recurrence"],
                      seed_source = "none (derived)", citations = [DUFFY_1982, GAUTSCHI_1982])
    return QuadratureRule(xs, ws, dom, PolynomialDegree(2m - 1), prov, cert)
end
