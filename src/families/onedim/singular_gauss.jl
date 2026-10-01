# Gauss-type rules for Cauchy principal-value and Hadamard finite-part integrals (PLAN §6
# Tier 2, v0.5). The kernels and the Hilbert transform are in domains/kernel.jl.
#
# With q_k the orthogonal polynomials of the Jacobi weight w (orthonormal up to the common
# factor that makes q₀ = 1, as in the Gauss driver), the integrals
#
#     ρ_k(t) = ⨍ w q_k / (x − t),       σ_k(t) = ⨎ w q_k / (x − t)² = ρ_k'(t)
#
# satisfy b_{k+1} ρ_{k+1} = (t − a_k) ρ_k − b_k ρ_{k−1} + μ₀ δ_{k0} and its t-derivative, from
# ρ₀ and σ₀. Run forward on (−1, 1), where no solution of the recurrence dominates another,
# this is stable.
#
# Principal value (Hunter 1972). For f of degree ≤ 2n, (f(x) − f(t))/(x − t) has degree
# ≤ 2n − 1 and the n-point Gauss rule integrates it exactly, so
#
#     ⨍ w f / (x − t) = Σ λ_k f(x_k)/(x_k − t) + f(t) [ρ₀ − Σ λ_k/(x_k − t)],
#
# n + 1 nodes and degree 2n. The bracket is ρ_n(t)/q_n(t) (apply the rule to q_n), which has
# no cancellation in it.
#
# Finite part. Every f of degree ≤ 2n is a + b(x − t) + (x − t)² p with p of degree ≤ 2n − 2.
# A rule with nodes y_1 … y_n and t integrates the last term exactly when λ_k = ω_k (y_k − t)²
# is an n-point rule for w exact to degree 2n − 2: the y_k are the zeros of q_n − c q_{n−1}
# for some c, and λ_k their Christoffel weights. The term b(x − t) has finite part ρ₀, and
# Σ λ_k/(y_k − t) = ρ₀ holds exactly when ρ_n − c ρ_{n−1} = 0 (apply the principal-value
# rule above, on these nodes, to the node polynomial). That fixes c, and the weight at t
# takes care of a. So the nodes are the zeros of
#
#     Q(x) = ρ_{n−1}(t) q_n(x) − ρ_n(t) q_{n−1}(x),
#
# and t is never one of them: Q(t) is the Casoratian of the solutions q and ρ, equal to
# −μ₀/b_n ≠ 0. All n zeros are real, but one may lie outside the interval; on (−1, 1) that
# happens for about half of all n, in runs of length π/(2θ) with t = cos θ. Then a larger n
# is taken, up to 2n₀, and past that the interpolatory rule on 2n₀ Gauss nodes and t, of
# degree 2n₀, which is never larger.

"""
    SingularGauss()

Rules for the Cauchy principal value and the Hadamard finite part,
[`PrincipalValue`](@ref) and [`FinitePart`](@ref), with a Jacobi weight on an interval.
Each rule has `n + 1` nodes, `t` and `n` others, and is exact to degree `2n`. The
principal-value rule is Hunter's: the Gauss–Jacobi nodes and `t`. The finite-part rule
places its `n` nodes at the zeros of `ρ_{n-1}(t) q_n(x) − ρ_n(t) q_{n-1}(x)`, which makes it
exact to degree `2n` without the derivative `f'(t)` that Gauss-type finite-part rules
usually need; near the ends of the interval it sometimes needs a larger `n` to keep its
nodes inside.

The weights are signed and large near `t`, and the rule evaluates the integrand at `t`. Both
are inherent to the integrals: `Σ |w_i|` grows as `t` approaches a node, so the
construction moves to `n + 1` points when `t` falls within an eighth of a gap of a node.
"""
struct SingularGauss <: RuleFamily end

derivation(::Type{SingularGauss}) = Derived()
family_name(::SingularGauss) = "SingularGauss"

const HUNTER_1972 = Citation(key = "Hunter1972", authors = ["D. B. Hunter"],
                             title = "Some Gauss-type formulae for the evaluation of Cauchy principal values of integrals",
                             journal = "Numerische Mathematik", year = 1972, volume = "19",
                             pages = "419--424", doi = "10.1007/BF01404924")
const GAUTSCHI_WIMP_1987 = Citation(key = "GautschiWimp1987", authors = ["Walter Gautschi", "Jet Wimp"],
                                    title = "Computing the Hilbert transform of a Jacobi weight function",
                                    journal = "BIT", year = 1987, volume = "27", pages = "203--215")

candidates(::Type{SingularGauss}, dom::KernelDomain, ::PolynomialDegree) = [SingularGauss()]
npoints(::SingularGauss, dom, degree::Integer) = cld(degree, 2) + 1
claimed_degree(::SingularGauss, dom, degree) = 2cld(degree, 2)
degree_range(::SingularGauss, dom) = 0:typemax(Int)
degree_for_npoints(::SingularGauss, dom, n::Integer) = 2(n - 1)
properties(::SingularGauss, dom, degree) = (positive = false, interior = true, symmetry = :none, nested = false)

# ρ_k and σ_k for k = 0 … K (index k + 1)
function kernel_recurrence(K, t, a, b, μ, ρ0, σ0)
    ρ, σ = Vector{BigFloat}(undef, K + 1), Vector{BigFloat}(undef, K + 1)
    ρ[1], σ[1] = ρ0, σ0
    for k in 0:(K - 1)
        r = (t - a[k + 1]) * ρ[k + 1] + (k == 0 ? μ : -b[k] * ρ[k])
        s = ρ[k + 1] + (t - a[k + 1]) * σ[k + 1] - (k == 0 ? zero(t) : b[k] * σ[k])
        ρ[k + 2], σ[k + 2] = r / b[k + 1], s / b[k + 1]
    end
    return ρ, σ
end

# q_{n−1}(x), q_n(x), their derivatives, and Σ_{k<n} q_k(x)²
function q_pair(n, x, a, b)
    q0, q1, d0, d1, s = zero(x), one(x), zero(x), zero(x), one(x)
    for k in 0:(n - 1)
        bk = k == 0 ? zero(x) : b[k]
        q2 = ((x - a[k + 1]) * q1 - bk * q0) / b[k + 1]
        d2 = (q1 + (x - a[k + 1]) * d1 - bk * d0) / b[k + 1]
        q0, q1, d0, d1 = q1, q2, d1, d2
        k < n - 1 && (s += q1^2)
    end
    return q0, q1, d0, d1, s
end

# How close t comes to a node, as a fraction of the gap around it (the interval ends bound
# the outer gaps); ≤ 1/2 inside the node range
function gap_ratio(x, t)
    j = searchsortedlast(x, t)
    lo = j == 0 ? -one(t) : x[j]
    hi = j == length(x) ? one(t) : x[j + 1]
    δ = min(j == 0 ? hi - lo : t - x[j], j == length(x) ? hi - lo : x[j + 1] - t)
    return Float64(δ / (hi - lo))
end

const NEAR_NODE = 1 / 8

# Hunter's principal-value rule on n Gauss nodes: nodes, weights (t last)
function hunter_rule(n, t, ρ, a, b, rec, bits)
    x, λ, _ = gauss_from_recurrence(n, rec, bits)
    ωt = ρ[n + 1] / q_pair(n, t, a, b)[2]
    return x, λ ./ (x .- t), ωt
end

# the interpolatory finite-part rule on n Gauss nodes and t
function interpolatory_fp_rule(n, t, ρ, σ, a, b, μ, rec, bits)
    x, λ, _ = gauss_from_recurrence(n, rec, bits)
    ω = map(x, λ) do xj, λj
        # f(t) + (x − t) g(x), with g interpolated on the Gauss nodes in the q_k
        q0, q1, s = zero(xj), one(xj), ρ[1]
        for k in 0:(n - 2)
            q2 = ((xj - a[k + 1]) * q1 - (k == 0 ? zero(xj) : b[k] * q0)) / b[k + 1]
            q0, q1 = q1, q2
            s += ρ[k + 2] * q1
        end
        λj / (xj - t) * s / μ
    end
    ωt = σ[n + 1] / q_pair(n, t, a, b)[2]           # the Lagrange polynomial of t is q_n/q_n(t)
    return x, ω, ωt
end

# the zeros of ρ_{n−1} q_n − ρ_n q_{n−1}, or nothing when one lies outside (−1, 1)
function finite_part_nodes(n, ρ, a, b, μ)
    iszero(ρ[n]) && return nothing
    c = b[n] * ρ[n + 1] / ρ[n]                       # the node polynomial is π_n − c π_{n−1}
    d = Float64.(a[1:n])
    d[n] += Float64(c)
    isfinite(d[n]) || return nothing
    seed = n == 1 ? d : eigvals(SymTridiagonal(d, Float64.(b[1:(n - 1)])))
    all(s -> -1 < s < 1, seed) || return nothing
    sort!(seed)
    y = BigFloat.(seed)
    tol = ldexp(BigFloat(1), -(precision(BigFloat) - 6))
    for i in eachindex(y)
        yi, last_small = y[i], false
        for it in 1:100
            p0, p1, d0, d1, _ = q_pair(n, yi, a, b)
            δ = (ρ[n] * p1 - ρ[n + 1] * p0) / (ρ[n] * d1 - ρ[n + 1] * d0)
            yi -= δ
            small = abs(δ) <= tol
            small && last_small && break
            last_small = small
            it == 100 && throw(RefinementError("SingularGauss", "Newton did not converge at finite-part node $i of $n"))
        end
        y[i] = yi
    end
    (all(yi -> -1 < yi < 1, y) && allunique(y)) || return nothing
    λ = [μ / q_pair(n, yi, a, b)[5] for yi in y]
    return y, λ
end

function finite_part_rule(n, t, ρ, σ, a, b, μ)
    r = finite_part_nodes(n, ρ, a, b, μ)
    r === nothing && return nothing
    y, λ = r
    p0, p1, _, _, _ = q_pair(n, t, a, b)
    ωt = (ρ[n] * σ[n + 1] - ρ[n + 1] * σ[n]) / (ρ[n] * p1 - ρ[n + 1] * p0)
    return y, λ ./ (y .- t) .^ 2, ωt
end

# The highest degree to which the rule is exact, from `d` up: a rule centred in a symmetric
# weight integrates every odd polynomial, and t at a zero of ρ_n gives the interpolatory rule
# degree 2n + 1, so the construction's own count can be low.
function probe_degree(nodes, ω, m, a, b, d, bits)
    K = length(m) - 1
    tol = ldexp(BigFloat(1), -(bits - 24))
    Σ, scale = zeros(BigFloat, K + 1), zeros(BigFloat, K + 1)
    for (x, w) in zip(nodes, ω)
        q0, q1 = zero(x), one(x)
        for k in 0:K
            Σ[k + 1] += w * q1
            scale[k + 1] += abs(w * q1)
            k == K && break
            q2 = ((x - a[k + 1]) * q1 - (k == 0 ? zero(x) : b[k] * q0)) / b[k + 1]
            q0, q1 = q1, q2
        end
    end
    while d < K && abs(Σ[d + 2] - m[d + 2]) <= tol * scale[d + 2]
        d += 1
    end
    return d
end

function build(f::SingularGauss, dom::KernelDomain, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("principal-value and finite-part nodes are irrational; $(T) is not supported"))
    check_kernel(dom)
    k = dom.weight
    n0 = cld(degree, 2)
    nmax = 2n0 + 1
    K = 2nmax + 3                                     # moments to probe the degree with
    wbits = ctx.bits + gj_guard_bits(nmax) + 16
    ρ0, σ0 = kernel_hilbert(k, dom.base.a, dom.base.b, wbits)
    rec = jacobi_recurrence(k.α, k.β)
    with_bits(wbits) do
        A, B = BigFloat(dom.base.a), BigFloat(dom.base.b)
        h, c = (B - A) / 2, (B + A) / 2
        t = (BigFloat(k.t) - c) / h
        a = [BigFloat(rec.a(j, BigFloat)) for j in 0:K]
        b = [BigFloat(rec.b(j, BigFloat)) for j in 1:(K + 1)]
        μ = BigFloat(rec.mass(BigFloat))
        ρ, σ = kernel_recurrence(K, t, a, b, μ, BigFloat(ρ0), BigFloat(σ0))
        quasi = false                                # nodes from the finite-part polynomial, not Gauss
        if n0 == 0
            x, ω, ωt = BigFloat[], BigFloat[], k.order == 1 ? ρ[1] : σ[1]
            how = "the single node t, weighted by the integral of 1"
        elseif k.order == 1
            choice = hunter_rule(n0, t, ρ, a, b, rec, wbits)
            if gap_ratio(choice[1], t) < NEAR_NODE
                alt = hunter_rule(n0 + 1, t, ρ, a, b, rec, wbits)
                gap_ratio(alt[1], t) > gap_ratio(choice[1], t) && (choice = alt)
            end
            x, ω, ωt = choice
            how = "the $(length(x)) Gauss–Jacobi nodes and t; weights λ_k/(x_k − t), and ρ_n(t)/q_n(t) at t"
        else
            choice = nothing
            for n in n0:(2n0)
                r = finite_part_rule(n, t, ρ, σ, a, b, μ)
                r === nothing && continue
                gap_ratio(r[1], t) >= NEAR_NODE && (choice = r; break)
            end
            if choice !== nothing
                x, ω, ωt = choice
                quasi = true
                how = "the $(length(x)) zeros of ρ_{n−1}(t) q_n − ρ_n(t) q_{n−1} and t; weights λ_k/(y_k − t)², " *
                      "with λ_k the Christoffel weights of the zeros"
            else
                choice = interpolatory_fp_rule(2n0, t, ρ, σ, a, b, μ, rec, wbits)
                if gap_ratio(choice[1], t) < NEAR_NODE
                    alt = interpolatory_fp_rule(2n0 + 1, t, ρ, σ, a, b, μ, rec, wbits)
                    gap_ratio(alt[1], t) > gap_ratio(choice[1], t) && (choice = alt)
                end
                x, ω, ωt = choice
                how = "the $(length(x)) Gauss–Jacobi nodes and t, interpolatory (no n ≤ $(2n0) placed every " *
                      "node of the degree-2n rule inside the interval)"
            end
        end
        m = k.order == 1 ? ρ : σ
        nodes = vcat(x, t)
        weights = vcat(ω, ωt)
        d = probe_degree(nodes, weights, m, a, b, 0, wbits)
        d >= degree || throw(RefinementError("SingularGauss", "the rule for $(k) came out exact to degree $d only, below $degree"))
        # to the interval asked for; t itself is taken as given, not mapped there and back
        s = kernel_scale(k, h)
        p = sortperm(nodes)
        xs = [i == length(nodes) ? finalize_number(ctx, BigFloat(k.t)) : finalize_number(ctx, c + h * nodes[i]) for i in p]
        ws = [finalize_number(ctx, s * weights[i]) for i in p]
        # the defining equations at the delivered nodes, back on [−1, 1]
        n = length(x)
        res = maximum(x; init = zero(t)) do xi
            xh = (BigFloat(finalize_number(ctx, c + h * xi)) - c) / h
            q = q_pair(n, xh, a, b)
            quasi ? abs((ρ[n] * q[2] - ρ[n + 1] * q[1]) / (ρ[n] * q[4] - ρ[n + 1] * q[3])) : abs(q[2] / q[4])
        end
        cert = Certificate(equations = n == 0 ? "none (one node)" :
                                       quasi ? "ρ_$(n - 1)(t) q_$n(y_i) − ρ_$n(t) q_$(n - 1)(y_i) = 0 (Newton correction)" :
                                       "q_$n(x_i) = 0 (Newton correction |q/q'|)",
                           residual = BigFloat(res; precision = 64), residual_bits = wbits,
                           digits = target_digits(ctx), guard_digits = floor(Int, (wbits - ctx.bits) * log10(2)))
        prov = Provenance(family = "SingularGauss", derivation = Derived(),
                          path = ["kernel: $(k) on $(dom.base)",
                                  "Hilbert transform of the weight: closed form (₂F₁ and logarithmic series), " *
                                  "checked at two precisions",
                                  "ρ_k(t)$(k.order == 2 ? " and σ_k(t)" : ""): forward recurrence at $wbits bits",
                                  "nodes: " * how,
                                  "degree: $d, the largest k with the rule exact on q_0 … q_k at $wbits bits"],
                          seed_source = "none (derived)",
                          citations = k.order == 1 ? [HUNTER_1972, GAUTSCHI_WIMP_1987] : [GAUTSCHI_WIMP_1987])
        QuadratureRule(xs, ws, dom, PolynomialDegree(d), prov, cert)
    end
end
