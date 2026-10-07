# Rules for the kernel 1/|y − x₀| on a triangle with x₀ off it: near-singular integrals
# (notes/singular-bem.md, stage 3).
#
# The triangle is cut at the point c nearest to x₀, at the distance H = |x₀ − c| > 0. On a
# sub-triangle (c, p, q), with a = p − c, b = q − c and the Duffy map y = c + s e(t),
# e(t) = a + t (b − a),
#
#     f(y) / |y − x₀| dy = f |a × b| s / √Q(s) ds dt,   Q(s) = q(t) s² + 2ℓ(t) s + H²,
#
# with q(t) = |e(t)|² and ℓ(t) = e(t)·(c − x₀). For each t, s/√Q is a weight on [0, 1] whose
# moments ∫ sᵏ/√Q ds follow a three-term recurrence from a closed form; its Gauss rule makes
# the s-direction exact for polynomials, however close x₀ is. What is left in t is analytic
# but not polynomial, and behaves like 1/√q(t) as x₀ approaches the triangle: t = t* + η sinh v,
# with q's roots t* ± iη, takes that out, and Gauss–Legendre in v is refined until the rule's
# integrals of every polynomial of its degree agree between two resolutions to the target
# precision. So the rule is exact to working precision, like every other rule here, but its
# size grows with that precision.

"""
    DuffySinh()

Rules for the kernel `1/|y − x₀|` on a triangle, an [`InverseDistance`](@ref) weight, when
`x₀` is off the triangle: near-singular integrals, where an ordinary rule needs many points
once `x₀` is close. The triangle is cut at its point `c` nearest to `x₀` and mapped by the
Duffy map from `c`. In the radial direction each line gets the Gauss rule of its own weight
`s / |y(s) − x₀|`, from the moments of that weight, so it is exact for polynomials however
close `x₀` is; in the angular direction, Gauss–Legendre after a sinh substitution, with as
many points as it takes for the rule to integrate `f(y)/|y − x₀|` for every polynomial `f` of
its degree to the requested precision. Positive weights, interior nodes; the number of points
grows with the precision asked for, and `npoints` gives the smallest.

```julia
T = SurfaceTriangle((0, 0, 0), (1, 0, 0), (0, 1, 0))
r = rule(WeightedDomain(T, InverseDistance((0.3, 0.3, 1e-3))); degree = 9)
integrate(y -> 1.0, r)
```
"""
struct DuffySinh <: RuleFamily end

derivation(::Type{DuffySinh}) = Derived()
family_name(::DuffySinh) = "DuffySinh"

const JOHNSTON_ELLIOTT_2005 = Citation(key = "JohnstonElliott2005", authors = ["P. R. Johnston", "D. Elliott"],
                                       title = "A sinh transformation for evaluating nearly singular boundary element integrals",
                                       journal = "International Journal for Numerical Methods in Engineering", year = 2005,
                                       volume = "62", pages = "564–578", doi = "10.1002/nme.1208")

candidates(::Type{DuffySinh}, dom::InverseDistanceDomain, ::PolynomialDegree) =
    kernel_geometry(dom).place === :near ? [DuffySinh()] : DuffySinh[]
npoints(::DuffySinh, dom, degree::Integer) = npatches(dom) * duffy_m(degree) * (duffy_m(degree) + 2)
claimed_degree(::DuffySinh, dom, degree) = 2duffy_m(degree) - 1
degree_range(::DuffySinh, dom) = 0:typemax(Int)
properties(::DuffySinh, dom, degree) = (positive = true, interior = true, symmetry = :none, nested = false)

"""
    radial_kernel_moments(q, ℓ, H², K) -> Vector

`Jₖ = ∫₀¹ sᵏ / √(q s² + 2ℓ s + H²) ds` for `k = 0 … K`, in the arithmetic of the arguments:
`J₀` in closed form, then `q k Jₖ = √Q(1) − [k = 1] √Q(0) − (2k − 1) ℓ Jₖ₋₁ − (k − 1) H² Jₖ₋₂`,
from `d/ds (sᵏ⁻¹ √Q)`. The forward recurrence amplifies errors by up to `2|ℓ|/q + H²/q` a step:
callers carry guard bits for it.
"""
function radial_kernel_moments(q, ℓ, H2, K::Integer)
    J = Vector{typeof(q)}(undef, K + 1)
    Q1, Q0 = sqrt(q + 2ℓ + H2), sqrt(H2)
    s0 = -ℓ / q                                        # Q = q ((s − s0)² + κ²)
    κ2 = H2 / q - s0^2
    J[1] = κ2 > 0 ? (asinh((1 - s0) / sqrt(κ2)) - asinh(-s0 / sqrt(κ2))) / sqrt(q) :
                    log((1 - s0) / -s0) / sqrt(q)      # κ = 0: x₀ on the line of the ray, behind c
    K >= 1 && (J[2] = (Q1 - Q0 - ℓ * J[1]) / q)
    for k in 2:K
        J[k + 1] = (Q1 - (2k - 1) * ℓ * J[k] - (k - 1) * H2 * J[k - 1]) / (q * k)
    end
    return J
end

"""
    radial_kernel_rule(q, ℓ, H², m, bits) -> (s, ω)

The `m`-point Gauss rule on `[0, 1]` of the weight `s / √(q s² + 2ℓ s + H²)`, as BigFloats of
precision `bits`: its moments `J₁ … J₂ₘ` from [`radial_kernel_moments`](@ref), at a precision
raised for the recurrence's amplification and for the Hankel conditioning of ordinary moments,
and checked at a higher one; then Wheeler's algorithm and Newton on the recurrence. Accepted
when the rule reproduces all `2m` moments to `bits`, which is exactness in `s` to degree
`2m − 1`; otherwise the guard is doubled.
"""
function radial_kernel_rule(q, ℓ, H2, m::Integer, bits::Integer)
    growth = with_bits(64) do
        log2(2 + 2abs(BigFloat(ℓ)) / BigFloat(q) + BigFloat(H2) / BigFloat(q))
    end
    guard = 64 + 10m + ceil(Int, 2m * Float64(growth))
    for _ in 1:6
        pbits = bits + guard
        μ = with_bits(pbits) do
            radial_kernel_moments(BigFloat(q), BigFloat(ℓ), BigFloat(H2), 2m + 2)[2:end]      # Jₖ₊₁: the weight's moments
        end
        check = with_bits(pbits + 64) do
            radial_kernel_moments(BigFloat(q), BigFloat(ℓ), BigFloat(H2), 2m + 2)[2:end]
        end
        tol = ldexp(BigFloat(1), -(bits + 8))
        good = with_bits(pbits) do
            all(k -> abs(μ[k] - check[k]) <= tol * abs(check[k]), eachindex(μ))
        end
        if good
            local s, ω
            ok = try
                α, β = with_bits(pbits) do
                    wheeler(MomentWeight(monomial_recurrence(), (k, T) -> T(μ[k + 1])), m + 1, BigFloat)   # the Gauss driver reads a₀ … aₘ
                end
                s, ω, _ = gauss_from_recurrence(m, moment_recurrence(α, β), pbits; name = "DuffySinh")
                with_bits(pbits) do
                    all(k -> abs(sum(ω .* s .^ (k - 1)) - μ[k]) <= tol * μ[k], 1:2m)
                end
            catch err
                err isa MomentBreakdownError || err isa RefinementError || rethrow()
                false
            end
            ok && return s, ω
        end
        guard *= 2
    end
    throw(RefinementError("DuffySinh", "the radial Gauss rule for q = $(Float64(q)), ℓ = $(Float64(ℓ)), " *
                                       "H² = $(Float64(H2)) did not reproduce its moments to $bits bits"))
end

function build(f::DuffySinh, dom::InverseDistanceDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("nodes for the kernel 1/|y − x₀| are irrational; $(T) is not supported"))
    g = kernel_geometry(dom)
    g.place === :near || throw(ArgumentError("DuffySinh is for x₀ off the triangle; $(Tuple(dom.weight.x0)) is on it: use DuffyGauss"))
    D = length(g.x0)
    m = duffy_m(degree)
    K = 2m - 1                                         # the degree the rule is exact to
    wbits = ctx.bits + 32
    tol = ldexp(BigFloat(1), -(ctx.bits + 4))
    # Duffy geometry of each sub-triangle, exact where it can be
    geo = map(g.patches) do (p, q)
        a, b = p - g.c, q - g.c
        (a, b, b - a)
    end
    # the rule with nt Gauss–Legendre points in v on each sub-triangle
    function assemble(nt)
        xs, ws = SVector{D,BigFloat}[], BigFloat[]
        with_bits(wbits) do
            v, wv = gauss_jacobi_work(nt, 0, 0, wbits)[1:2]
            c, x0 = SVector{D,BigFloat}(g.c), SVector{D,BigFloat}(g.x0)
            H2 = BigFloat(g.H2)
            for (ar, br, dr) in geo
                a, b, d = SVector{D,BigFloat}(ar), SVector{D,BigFloat}(br), SVector{D,BigFloat}(dr)
                dd = dot(d, d)
                ts, η = -dot(a, d) / dd, _area2(a, b) / dd    # q(t) = |d|² ((t − t*)² + η²)
                v0, v1 = asinh(-ts / η), asinh((1 - ts) / η)
                vm, hv = (v0 + v1) / 2, (v1 - v0) / 2
                J = _area2(a, b)
                for (vk, wk) in zip(v, wv)
                    t = ts + η * sinh(vm + hv * vk)
                    e = a + t * d
                    qt, ℓt = dot(e, e), dot(e, c - x0)
                    s, ω = radial_kernel_rule(qt, ℓt, H2, m, wbits)
                    dtdv = hv * sqrt(qt / dd)                 # dt = η cosh v dv = √(q/|d|²) dv
                    for (si, ωi) in zip(s, ω)
                        push!(xs, c + si * e)
                        push!(ws, wk * dtdv * ωi * J)
                    end
                end
            end
        end
        return xs, ws
    end
    # the integrals of the triangle's orthonormal Dubiner polynomials to degree K
    V = [SVector{D,BigFloat}(map(BigFloat, p)) for p in vertices(dom.base)]
    E = hcat(V[2] - V[1], V[3] - V[1])
    A = inv(E' * E) * E'
    moments(xs, ws) = with_bits(wbits) do
        dub = DubinerBasis{BigFloat}(K; normalize = true)
        total = zeros(BigFloat, dubiner_length(K))
        for (x, w) in zip(xs, ws)
            ξ = A * (x - V[1])
            dubiner!(dub.φ, dub.gx, dub.gy, dub.ws, ξ[1], ξ[2])
            total .+= w .* dub.φ
        end
        total
    end
    # Refine until two resolutions agree, and deliver the coarser: the angular integrand is
    # analytic, so the error falls geometrically and the finer rule is far more accurate than
    # the difference, which then bounds the coarser one's error.
    nt = m + 2
    xs, ws = assemble(nt)
    prev = moments(xs, ws)
    gap = BigFloat(Inf)
    for _ in 1:12
        nt2 = ceil(Int, 1.5nt)
        xs2, ws2 = assemble(nt2)
        cur = moments(xs2, ws2)
        gap = with_bits(wbits) do
            maximum(abs, cur - prev) / max(one(BigFloat), maximum(abs, cur))
        end
        ctx.verbose >= 1 && @info @sprintf("DuffySinh: %d against %d points per sub-triangle in t, moments differ by %.1e",
                                           nt, nt2, Float64(gap))
        gap <= tol && break
        nt, xs, ws, prev = nt2, xs2, ws2, cur
        nt > 2000 && throw(RefinementError("DuffySinh", "the angular resolution did not converge"))
    end
    gap <= tol || throw(RefinementError("DuffySinh", "the angular resolution did not converge to $(ctx.bits) bits"))
    nodes = [SVector{D,T}(ntuple(i -> finalize_number(ctx, x[i]), D)) for x in xs]
    weights = [finalize_number(ctx, w) for w in ws]
    H = sqrt(Float64(g.H2))
    cert = Certificate(equations = "integrals of the orthonormal Dubiner polynomials to degree $K, " *
                                   "converged between two angular resolutions",
                       residual = BigFloat(gap; precision = 64), residual_bits = wbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
    prov = Provenance(family = "DuffySinh", derivation = Derived(),
                      path = [@sprintf("kernel: 1/|y − x₀| on %s, x₀ = %s off the triangle, at %.3e from its nearest point; %d sub-triangle%s with a vertex there",
                                       string(dom.base), string(Tuple(dom.weight.x0)), H, length(g.patches),
                                       length(g.patches) == 1 ? "" : "s"),
                              "Duffy map y = c + s (a + t (b − a)) on each: f(y)/|y − x₀| dy = f |a × b| s / √Q(s) ds dt",
                              "s: $m-point Gauss rule of s/√Q(s) on each line, from its moments (closed form and a " *
                              "three-term recurrence) by Wheeler's algorithm, checked against them",
                              "t: $nt-point Gauss–Legendre after t = t* + η sinh v, refined until the integrals of the " *
                              "Dubiner polynomials to degree $K agreed to $(ctx.bits) bits"],
                      seed_source = "none (derived)", citations = [DUFFY_1982, JOHNSTON_ELLIOTT_2005])
    return QuadratureRule(nodes, weights, dom, PolynomialDegree(K), prov, cert)
end
