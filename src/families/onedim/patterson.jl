# Gauss–Patterson rules (PLAN §6, v0.5): the nested sequence of 1, 3, 7, 15, … points on
# [-1, 1] in which each rule keeps every node of the one before and adds the nodes that raise
# the degree most — Kronrod's extension, iterated (Patterson 1968).
#
# If the previous rule has n nodes with node polynomial Π, the new nodes are the roots of the
# polynomial q of degree n + 1 for which Π q is orthogonal to every polynomial of degree ≤ n.
# The 2n + 1 nodes together are then exact to degree 3n + 1 (interpolation gives 2n, the
# orthogonality n + 1 more) and, being symmetric, to 3n + 2.
#
# Π q is odd, of degree 2n + 1 and orthogonal to degree ≤ n, so it is a combination of the odd
# Legendre polynomials from n + 2 to 2n + 1. With the top coefficient fixed, the other
# (n - 1)/2 follow from requiring Π q to vanish at the positive old nodes (zero comes free
# with oddness): one linear system per level. The new nodes are the remaining roots. They
# interlace with the old ones, one between each adjacent pair and one beyond the last, so each
# is bracketed, and the sign of q on either side is known in advance. The root search never
# evaluates q at an old node, where it is 0/0.
#
# The construction loses digits as the levels go up: about 6 at 63 points, 17 at 127, 43 at
# 255 and 95 at 511. Newton on the new nodes directly, with q kept in product form, loses the
# same at 255, so the loss comes from the problem rather than this formulation. It tracks how
# nearly each rule reaches the next degree: the 127-point rule misses degree 192 by only
# 1e-20, the 255-point rule degree 384 by 6e-41, the 511-point rule degree 768 by 3e-81, and
# nodes that almost satisfy one more equation are poorly pinned down by the ones they must.
# The working precision is therefore raised until two precisions agree, as for
# ModifiedChebyshev. The starting guard comes from those measurements, so one confirming run
# is normally enough.
#
# The weights are interpolatory: ∫ ℓᵢ over [-1, 1] for the Lagrange basis polynomial ℓᵢ,
# computed exactly by a Gauss–Legendre rule with more than half as many points.

"""
    GaussPatterson()

Gauss–Patterson rules on [`Interval`](@ref): the nested sequence of 1, 3, 7, 15, 31, 63,
127, 255 and 511 points in which each rule keeps every node of the one before and adds the
nodes that raise the degree most. The rule with `2^(k+1) - 1` points has degree `3·2^k - 1`.

Because the rules are nested, a caller that moves to the next level reuses every function
value it already has, and an [`embedded`](@ref) pair estimates the error from the same values.
The cost is extra points: the 255-point rule has degree 383, which Gauss–Legendre reaches
with 192. All weights are positive and all nodes interior. The rules are derived at any
precision.

From 127 points on, a rule misses the next degree by less than `Float64` rounding (by `1e-20`
at 127 points), so in `Float64` it is exact one degree higher to within rounding.
[`verify`](@ref) reports sharpness there as not resolvable, with the miss measured during
construction, rather than as an understated claim.
"""
struct GaussPatterson <: RuleFamily end

derivation(::Type{GaussPatterson}) = Derived()
describe_family(::GaussPatterson) = "GaussPatterson"

const PATTERSON_1968 = Citation(key = "Patterson1968", authors = ["T. N. L. Patterson"],
                                title = "The optimum addition of points to quadrature formulae",
                                journal = "Mathematics of Computation", year = 1968, volume = "22",
                                pages = "847--856", doi = "10.1090/S0025-5718-68-99866-9")

"The last level offered. Level 9 (1023 points) also computes, with interlacing nodes and positive
weights, but takes 90 s in Float64 and is beyond the published tables, so it is not claimed."
const PATTERSON_MAXLEVEL = 8

patterson_npoints(k::Integer) = 2^(k + 1) - 1
patterson_degree(k::Integer) = k == 0 ? 1 : 3 * 2^k - 1

"The lowest level whose degree reaches `d`."
function patterson_level(d::Integer)
    k = 0
    while patterson_degree(k) < d
        k += 1
    end
    return k
end

candidates(::Type{GaussPatterson}, dom::Interval, c::PolynomialDegree) =
    isreference(dom) && patterson_level(c.d) <= PATTERSON_MAXLEVEL ? [GaussPatterson()] : GaussPatterson[]

npoints(::GaussPatterson, dom, degree::Integer) = patterson_npoints(patterson_level(degree))
claimed_degree(::GaussPatterson, dom, degree) = patterson_degree(patterson_level(degree))
degree_range(::GaussPatterson, dom::Interval) = 0:patterson_degree(PATTERSON_MAXLEVEL)
degree_range(::GaussPatterson, dom) = 1:0
properties(::GaussPatterson, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = true)

function degree_for_npoints(::GaussPatterson, dom, n::Integer)
    for k in 0:PATTERSON_MAXLEVEL
        patterson_npoints(k) == n && return patterson_degree(k)
    end
    throw(ArgumentError("Gauss–Patterson rules have 1, 3, 7, 15, …, " *
                        "$(patterson_npoints(PATTERSON_MAXLEVEL)) points, not $n"))
end

"""
    patterson_combination(e, S, x) -> (F, F′)

`F(x) = Σᵢ eᵢ P_{Sᵢ}(x)` and its derivative, for ascending Legendre degrees `S`, by the
three-term recurrence.
"""
function patterson_combination(e, S, x)
    p0, p1, d0, d1 = one(x), x, zero(x), one(x)
    F, dF = zero(x), zero(x)
    k = 1
    if S[1] == 1
        F, dF, k = e[1] * p1, e[1] * d1, 2
    end
    for j in 1:S[end]-1
        p2 = ((2j + 1) * x * p1 - j * p0) / (j + 1)
        d2 = d0 + (2j + 1) * p1
        if k <= length(S) && S[k] == j + 1
            F += e[k] * p2
            dF += e[k] * d2
            k += 1
        end
        p0, p1, d0, d1 = p1, p2, d1, d2
    end
    return F, dF
end

"""
    patterson_root(e, S, old, lo, hi, negleft, tol) -> x or nothing

The new node in `(lo, hi)`: the root of `q = F/Π` there, where `F` is the Legendre combination
`(e, S)` and `Π` the node polynomial of the positive nodes `old` and zero. `negleft` says
whether `q < 0` to the left of the root, which interlacing fixes in advance. Safeguarded
Newton; it stops once the steps stop shrinking, which at a high level happens well above the
rounding unit because of the digits the construction loses. `nothing` if it does not converge.
"""
function patterson_root(e, S, old, lo, hi, negleft, tol)
    x = (lo + hi) / 2
    prevstep = hi - lo
    # Well into convergence a Newton step at least halves the last one, until rounding noise
    # takes over. That noise sits near cond × eps, which at 511 points is 1e-88 against a
    # rounding unit of 1e-185, so "small" has to be relative to the bracket, not to `tol`.
    small = (hi - lo) / 1_000_000
    for _ in 1:4096
        F, dF = patterson_combination(e, S, x)
        s = inv(x) + sum((2x / (x^2 - t^2) for t in old); init = zero(x))        # Π′/Π
        Πneg = isodd(count(t -> x < t, old))                                     # sign of Π for x > 0
        if (signbit(F) != Πneg) == negleft
            lo = x
        else
            hi = x
        end
        xn = x - F / (dF - F * s)
        step = abs(xn - x)
        # before the bracket test: a converged step can be zero, or land exactly on the end of
        # the bracket just moved to x, and bisecting from there would throw the root away
        step <= tol && return x
        if lo < xn < hi
            step < small && step >= prevstep / 2 && return xn
            prevstep = step
        elseif prevstep < small
            return x              # converged; the noise now points out of the tightened bracket
        else
            xn = (lo + hi) / 2
        end
        hi - lo <= tol && return xn
        x = xn
    end
    return nothing
end

"""
    patterson_nodes(K, bits) -> Vector{BigFloat} or nothing

The positive nodes of the level-`K` Gauss–Patterson rule, ascending, computed at `bits`
through every level below it. `nothing` if a root search fails, which more precision cures.
"""
function patterson_nodes(K::Integer, bits::Integer; cancel = nothing)
    return with_bits(bits) do
        tol = ldexp(BigFloat(1), -bits + 4)
        pos = BigFloat[]
        for _ in 1:K
            checkcancel(cancel)
            n = 2length(pos) + 1
            m = (n + 1) ÷ 2
            S = collect(n+2:2:2n+1)
            A = zeros(BigFloat, m - 1, m - 1)
            b = zeros(BigFloat, m - 1)
            for (r, t) in enumerate(pos)
                p0, p1, c = one(t), t, 1
                for j in 1:2n
                    p2 = ((2j + 1) * t * p1 - j * p0) / (j + 1)
                    if c <= m && S[c] == j + 1
                        c < m ? (A[r, c] = p2) : (b[r] = -p2)
                        c += 1
                    end
                    p0, p1 = p1, p2
                end
            end
            e = vcat(m > 1 ? A \ b : BigFloat[], one(BigFloat))
            edges = vcat([zero(BigFloat)], pos, [one(BigFloat)])
            new = BigFloat[]
            for i in 1:m
                y = patterson_root(e, S, pos, edges[i], edges[i+1], isodd(m - i + 1), tol)
                y === nothing && return nothing
                push!(new, y)
            end
            pos = sort!(vcat(pos, new))
        end
        pos
    end
end

"""
    patterson_rule(pos, bits) -> (x, w)

The symmetric rule on the nodes `±pos` and zero with interpolatory weights: `wᵢ = ∫ ℓᵢ`,
integrated exactly by Gauss–Legendre as `Σₗ gₗ Ω(xₗ) / ((xₗ - xᵢ) Ω′(xᵢ))`, with `Ω` the node
polynomial.
"""
function patterson_rule(pos::Vector{BigFloat}, bits::Integer)
    return with_bits(bits) do
        x = vcat(-reverse(pos), [zero(BigFloat)], pos)
        N = length(x)
        N == 1 && return x, [BigFloat(2)]
        # 2M - 1 ≥ N - 1; M = (N + 1)/2 is even from three points on, so zero is not a
        # Gauss–Legendre node and the division below never meets the centre node
        M = (N + 1) ÷ 2
        xg, wg, _ = gauss_jacobi_work(M, 0, 0, bits)
        Ω = [prod(xl - xi for xi in x) for xl in xg]
        w = similar(x)
        c = (N + 1) ÷ 2
        for i in c:N
            dΩ = prod((x[i] - x[j] for j in 1:N if j != i); init = one(BigFloat))
            s = sum(xg[l] == x[i] ? wg[l] * dΩ : wg[l] * Ω[l] / (xg[l] - x[i]) for l in 1:M)
            w[i] = s / dΩ
            w[N+1-i] = w[i]
        end
        x, w
    end
end

"""
    patterson_guard_bits(K)

Starting guard for level `K`: the digits the construction was measured to lose there, plus
the usual margin. The agreement test decides; this only saves re-runs.
"""
patterson_guard_bits(K::Integer) = 32 + (0, 0, 0, 0, 4, 24, 64, 150, 330, 720)[K+1]

"""
    patterson_work(K, bits; verbose, cancel) -> (x, w, usedbits, cond)

The level-`K` Gauss–Patterson rule accurate to `bits` bits, as BigFloats. The working
precision rises until two runs 64 or more bits apart agree to `bits`; `cond` is the
amplification that measured (the discrepancy at the lower precision over its rounding unit).
"""
function patterson_work(K::Integer, bits::Integer; verbose::Integer = 0, cancel = nothing,
                        maxbits::Integer = 16 * bits + 8192)
    K == 0 && return (with_bits(() -> [zero(BigFloat)], bits), with_bits(() -> [BigFloat(2)], bits), bits, 1.0)
    base = bits + patterson_guard_bits(K)
    tol = ldexp(BigFloat(1), -(bits + 4))
    function compute(b)
        pos = patterson_nodes(K, b; cancel)
        pos === nothing && return nothing
        return patterson_rule(pos, b)
    end
    prev, prev_bits = compute(base), base
    extra = 64
    gap = BigFloat(NaN)
    while true
        trial = prev_bits + extra
        trial > maxbits && throw(RefinementError("GaussPatterson", """
            the $(patterson_npoints(K))-point Gauss–Patterson rule did not reach $(bits) bits, \
            even at $(prev_bits) bits of working precision\
            $(isnan(gap) ? "" : " (successive precisions still disagree by $(Float64(gap)))")."""))
        cur = compute(trial)
        if cur === nothing
            verbose >= 1 && @info @sprintf("  %d bits: a root search failed", trial)
            extra *= 2
            continue
        end
        if prev !== nothing
            gap = with_bits(trial) do
                rule_discrepancy(prev[1], prev[2], cur[1], cur[2])
            end
            verbose >= 1 && @info @sprintf("  %d bits: disagrees with %d bits by %.3e (want %.3e)",
                                           trial, prev_bits, Float64(gap), Float64(tol))
            cond = gap > 0 ? exp2(prev_bits + Float64(log2(gap))) : 1.0
            gap <= tol && return (cur[1], cur[2], trial, max(cond, 1.0))
            # the gap measured the amplification; go straight to the precision it asks for
            need = bits + 36 + ceil(Int, log2(max(cond, 1.0)))
            if need > trial
                prev, prev_bits, extra = compute(need), need, 64
                continue
            end
        end
        prev, prev_bits = cur, trial
        extra = 64
    end
end

function build(f::GaussPatterson, dom::Interval, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("Gauss–Patterson nodes are irrational; $(T) is not supported"))
    K = patterson_level(degree)
    N = patterson_npoints(K)
    d = patterson_degree(K)
    ctx.verbose >= 1 && @info @sprintf("Gauss–Patterson: level %d, %d points, degree %d, target %d bits",
                                       K, N, d, ctx.bits)
    x, w, used, cond = patterson_work(K, ctx.bits; verbose = ctx.verbose, cancel = ctx.cancel)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    # the claim itself, on the delivered numbers: Σ wᵢ Pⱼ(xᵢ) = 2δⱼ₀ for j ≤ d
    res = with_bits(used) do
        s = zeros(BigFloat, d + 1)
        for i in 1:N
            xi, wi = BigFloat(xs[i]), BigFloat(ws[i])
            p0, p1 = one(xi), xi
            s[1] += wi
            d >= 1 && (s[2] += wi * p1)
            for j in 1:d-1
                p0, p1 = p1, ((2j + 1) * xi * p1 - j * p0) / (j + 1)
                s[j+2] += wi * p1
            end
        end
        s[1] -= 2
        maximum(abs, s) / 2
    end
    # The rule misses degree d + 1 by little: 1e-20 at 127 points, 6e-41 at 255, and the
    # construction's amplification is about the inverse of that. Measure it on the unrounded
    # rule, with the test function verification uses (orthonormal p_a², a = (d + 1)/2), so a
    # rule delivered at lower precision can be reported as unresolvable there, not understated.
    miss = with_bits(used) do
        a = (d + 1) ÷ 2
        abs(sum(w[i] * (sqrt(BigFloat(2a + 1) / 2) * legendre_value(a, x[i]))^2 for i in 1:N) - 1)
    end
    cert = Certificate(equations = "Σᵢ wᵢ Pⱼ(xᵢ) = 2δⱼ₀ for j ≤ $d (relative to the mass)",
                       residual = BigFloat(res; precision = 64), residual_bits = used,
                       digits = target_digits(ctx), guard_digits = floor(Int, (used - ctx.bits) * log10(2)),
                       cond = cond, next_error = BigFloat(miss; precision = 64))
    prov = Provenance(family = "GaussPatterson", derivation = Derived(),
                      path = ["$K Kronrod extensions from the midpoint rule, at $used bits",
                              "each: Legendre coefficients of Π q from one linear system, new nodes " *
                              "by bracketed Newton between the old ones",
                              "weights: interpolatory, integrated by Gauss–Legendre",
                              "precision: raised until two runs agreed to $(ctx.bits) bits"],
                      seed_source = "none (derived)", citations = [PATTERSON_1968], symmetry = :reflection)
    return QuadratureRule(xs, ws, Interval(), PolynomialDegree(d), prov, cert)
end

"Legendre `Pⱼ(x)` by the three-term recurrence."
function legendre_value(j::Integer, x)
    j == 0 && return one(x)
    p0, p1 = one(x), x
    for k in 1:j-1
        p0, p1 = p1, ((2k + 1) * x * p1 - k * p0) / (k + 1)
    end
    return p1
end
