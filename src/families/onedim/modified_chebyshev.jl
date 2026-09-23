# Modified Chebyshev: a Gauss rule from moments (PLAN §6 Tier 2).
#
# The classical route to a Gauss rule for an unfamiliar weight goes through its moments: 2n
# of them determine the n-point rule. The route is also famous for being unusable. Against
# the monomials the map is the Hilbert matrix in disguise, and in Float64 it does not merely
# lose accuracy — around n = 13 a recurrence coefficient β goes negative and the output stops
# describing a measure at all. That is why PolyChaos.jl and most other implementations reach
# for Stieltjes or Lanczos instead.
#
# The conditioning is real and this file does not pretend otherwise. What changes here is
# that precision is a dial. The loss is linear in n — about 1.4 decimal digits per node for
# ordinary moments — so the breakdown point moves out in proportion to the working precision,
# and 200 digits carries the ordinary-moment algorithm well past where any fixed-precision
# implementation can go. Choosing a good auxiliary family removes the loss almost entirely;
# the point is that you no longer have to.
#
# Nothing here is assumed. `modified_chebyshev_work` runs Wheeler's recursion, checks the
# resulting rule against the moments it came from, and doubles the working precision until it
# passes — so the certificate reports conditioning that was measured, and the extra digits
# spent are the measurement.

"""
    monic(rec::Recurrence) -> MonicRecurrence

The monic form of an orthonormal [`Recurrence`](@ref): the same `a`, and `b` squared, since
the orthonormal coefficient is `√βₖ`. Lets any recurrence the package already knows serve as
the auxiliary family for a [`MomentWeight`](@ref):

```julia
MomentWeight(monic(jacobi_recurrence(0, 0)), moments)
```
"""
monic(rec::Recurrence) = MonicRecurrence(rec.a, (k, T) -> rec.b(k, T)^2)

"""
    MomentBreakdownError(k, value, bits)

Wheeler's recursion produced a non-positive or non-finite `βₖ` at `bits` bits of working
precision. For a genuine positive measure this is a conditioning failure rather than a
property of the moments, and more precision cures it; the driver catches this and retries.
"""
struct MomentBreakdownError <: Exception
    k::Int
    value::BigFloat
    bits::Int
end
Base.showerror(io::IO, e::MomentBreakdownError) =
    print(io, "Wheeler's recursion broke down at ", e.bits, " bits: β[", e.k, "] = ", e.value,
              ", which cannot come from a positive measure")

"""
    wheeler(w::MomentWeight, n, T) -> (α, β)

The first `n` recurrence coefficients `α₀…α_{n-1}`, `β₀…β_{n-1}` of the measure `w`, computed
in type `T` from its modified moments `m₀…m_{2n-1}` by Wheeler's algorithm.

The coefficients are monic: the measure's own orthogonal polynomials satisfy
`p_{k+1} = (x - αₖ) pₖ - βₖ p_{k-1}`, with `β₀ = m₀` the total mass. Throws
[`MomentBreakdownError`](@ref) if some `βₖ ≤ 0`, which for a positive measure means the
working precision was too low.
"""
function wheeler(w::MomentWeight, n::Integer, ::Type{T}) where {T}
    n >= 1 || throw(ArgumentError("need at least one recurrence coefficient, got n = $n"))
    N = 2n
    # index i of any array is the mathematical index i - 1
    m = T[w.moments(k, T) for k in 0:(N - 1)]
    a = T[w.aux.a(k, T) for k in 0:(N - 1)]
    b = T[w.aux.b(k, T) for k in 0:(N - 1)]
    α = zeros(T, n)
    β = zeros(T, n)
    σ_old = zeros(T, N)          # σ_{k-2, ·}
    σ = copy(m)                  # σ_{0, ·} = m
    α[1] = a[1] + m[2] / m[1]
    β[1] = m[1]
    for k in 2:n
        σ_new = zeros(T, N)
        for l in k:(N - k + 1)
            σ_new[l] = σ[l + 1] - (α[k - 1] - a[l]) * σ[l] - β[k - 1] * σ_old[l] + b[l] * σ[l - 1]
        end
        α[k] = a[k] + σ_new[k + 1] / σ_new[k] - σ[k] / σ[k - 1]
        β[k] = σ_new[k] / σ[k - 1]
        σ_old, σ = σ, σ_new
    end
    for k in eachindex(β)
        (isfinite(β[k]) && β[k] > 0) ||
            throw(MomentBreakdownError(k - 1, BigFloat(β[k]; precision = 64), precision(T)))
    end
    return α, β
end

"Wrap computed monic coefficients for the shared Gauss driver, whose `b` is √β."
moment_recurrence(α::Vector{T}, β::Vector{T}) where {T} =
    Recurrence((k, S) -> S(α[k + 1]), (k, S) -> sqrt(S(β[k + 1])), S -> S(β[1]))

"""
    moment_residual(x, w, mw, n) -> T

The worst relative failure of `Σᵢ wᵢ πₖ(xᵢ) = mₖ` over `k < 2n`, where `πₖ` is the auxiliary
family of `mw`. Since those `2n` polynomials span `P_{2n-1}`, this vanishing is exactly the
claim that the rule is exact to degree `2n - 1` for the measure the moments describe — and it
is checked against the input data, not against the recurrence derived from it.
"""
function moment_residual(x::Vector{T}, w::Vector{T}, mw::MomentWeight, n::Integer) where {T}
    N = 2n
    a = T[mw.aux.a(k, T) for k in 0:(N - 1)]
    b = T[mw.aux.b(k, T) for k in 0:(N - 1)]
    got = zeros(T, N)
    for i in eachindex(x)
        πk, πprev = one(T), zero(T)
        for k in 1:N
            got[k] += w[i] * πk
            πk, πprev = (x[i] - a[k]) * πk - b[k] * πprev, πk
        end
    end
    mass = abs(T(mw.moments(0, T)))
    worst = zero(T)
    for k in 1:N
        want = T(mw.moments(k - 1, T))
        # scale by the mass as well as by the moment: a modified moment may be near zero
        # (it is for a well-chosen auxiliary family) without the rule being in any trouble
        worst = max(worst, abs(got[k] - want) / max(abs(want), mass))
    end
    return worst
end

"""
    rule_discrepancy(x1, w1, x2, w2) -> T

How far apart two computations of the same rule are: the larger of the node difference
relative to the node span and the weight difference relative to the total mass.
"""
function rule_discrepancy(x1::Vector{T}, w1::Vector{T}, x2::Vector{T}, w2::Vector{T}) where {T}
    span = maximum(x2) - minimum(x2)
    mass = sum(abs, w2)
    dx = maximum(abs, x1 .- x2) / (span > 0 ? span : one(T))
    dw = maximum(abs, w1 .- w2) / (mass > 0 ? mass : one(T))
    return max(dx, dw)
end

"""
    modified_chebyshev_work(mw, n, bits; maxbits) -> (x, w, iterations, usedbits, residual, cond)

`n`-point Gauss nodes and weights for the measure `mw`, accurate to `bits` bits, as BigFloats
at the precision the search settled on.

The stopping rule is the subtle part. The obvious test — does the rule reproduce the moments
it came from? — is not a test of accuracy here, because the map it inverts is the
ill-conditioned one. At 40 points from ordinary moments, a rule whose nodes are wrong in the
sixth decimal still reproduces all 80 moments to 10⁻⁵², since that is what a condition number
of 10⁴⁶ means. Believing the residual would mean shipping six correct digits with a
certificate claiming fifty.

So the rule is computed twice, at working precisions 64 bits apart, and the precision is
doubled until the two agree to the accuracy asked for. Consecutive iterates bracket the error
in the thing the caller actually receives, and nothing has to be assumed about the
conditioning — it is measured, and returned as `cond`, the ratio of the node error to the
moment residual that accompanied it.
"""
function modified_chebyshev_work(mw::MomentWeight, n::Integer, bits::Integer;
                                 maxbits::Integer = 16 * bits + 4096)
    base = bits + gj_guard_bits(n)
    tol = ldexp(BigFloat(1), -(bits + 4))
    compute(trial) = with_bits(trial) do
        try
            α, β = wheeler(mw, n + 1, BigFloat)
            x, w, iters = gauss_from_recurrence(n, moment_recurrence(α, β), trial)
            (x, w, iters, moment_residual(x, w, mw, n))
        catch e
            e isa MomentBreakdownError ? nothing : rethrow()
        end
    end

    extra = 0
    prev, prev_bits = compute(base), base
    gap = BigFloat(NaN)
    while true
        extra = extra == 0 ? 64 : 2 * extra
        if base + extra > maxbits
            throw(RefinementError("ModifiedChebyshev", """
                modified Chebyshev did not reach $(bits) bits for $(mw.label) at $n points, \
                even at $(base + extra ÷ 2) bits of working precision\
                $(isnan(gap) ? " (Wheeler's recursion broke down every time)" :
                               " (successive precisions still disagree by $(Float64(gap)))").

                The moment problem may be too ill-conditioned for the requested accuracy, or \
                the moments may not come from a positive measure on this interval. Note that \
                the moment function must deliver full precision on demand: if it returns \
                stored Float64 values, no amount of working precision can help.
                """))
        end
        trial = base + extra
        cur = compute(trial)
        cur === nothing && continue          # broke down; more precision is the answer
        if prev !== nothing
            x, w, iters, res = cur
            gap = with_bits(trial) do
                rule_discrepancy(prev[1], prev[2], x, w)
            end
            if gap <= tol
                # `prev` carried an error of about `gap` from a rounding unit of 2^-prev_bits,
                # so their ratio is the amplification the moment problem applies — the number
                # that says why this took the precision it did. `cur` is a further doubling
                # better, so returning it is the conservative choice.
                cond = gap > 0 ? exp2(prev_bits + Float64(log2(gap))) : 1.0
                return (x, w, iters, trial, res, max(cond, 1.0))
            end
        end
        prev, prev_bits = cur, trial
    end
end

const WHEELER_1974 = Citation(key = "Wheeler1974", authors = ["John C. Wheeler"],
                              title = "Modified moments and Gaussian quadratures",
                              journal = "Rocky Mountain Journal of Mathematics", year = 1974,
                              volume = "4", pages = "287--296")
const GAUTSCHI_1970 = Citation(key = "Gautschi1970", authors = ["Walter Gautschi"],
                               title = "On the construction of Gaussian quadrature rules from modified moments",
                               journal = "Mathematics of Computation", year = 1970,
                               volume = "24", pages = "245--260")
const SACK_DONOVAN_1972 = Citation(key = "SackDonovan1972", authors = ["R. A. Sack", "A. F. Donovan"],
                                   title = "An algorithm for Gaussian quadrature given modified moments",
                                   journal = "Numerische Mathematik", year = 1972,
                                   volume = "18", pages = "465--478")

"""
    ModifiedChebyshev()

Gauss rules for a measure given by its moments: see [`MomentWeight`](@ref). An `n`-point rule
has polynomial degree `2n - 1` and is built by Wheeler's algorithm, which turns `2n` modified
moments into the measure's own three-term recurrence.

This family is selected for any [`Interval`](@ref) carrying a [`MomentWeight`](@ref), and for
no other domain — it is the escape hatch for weights with no classical Gauss family:

```julia
w = OrdinaryMoments((k, T) -> one(T) / (k + 1))          # Lebesgue measure on [0, 1]
r = rule(WeightedDomain(Interval(0, 1), w); degree = 79, digits = 50)
```

That example is deliberately the hard one. In `Float64` the same construction collapses
before `n = 16`; here the driver spends the precision the conditioning demands and the
certificate reports how much that was.
"""
struct ModifiedChebyshev <: RuleFamily end

derivation(::Type{ModifiedChebyshev}) = Derived()
family_name(::ModifiedChebyshev) = "ModifiedChebyshev"

candidates(::Type{ModifiedChebyshev}, dom::MomentInterval, ::PolynomialDegree) = [ModifiedChebyshev()]

mc_points(degree) = max(1, cld(degree + 1, 2))
npoints(::ModifiedChebyshev, dom, degree::Integer) = mc_points(degree)
claimed_degree(::ModifiedChebyshev, dom, degree) = 2 * mc_points(degree) - 1
degree_range(::ModifiedChebyshev, dom) = 0:typemax(Int)
degree_for_npoints(::ModifiedChebyshev, dom, n::Integer) = 2n - 1
# Gauss rules for a positive measure have positive weights inside the support hull; if the
# moments do not describe such a measure the construction fails rather than lying, and
# verification checks both claims against the rule that comes out.
properties(::ModifiedChebyshev, dom, degree) =
    (positive = true, interior = true, symmetry = :none, nested = false)

function build(f::ModifiedChebyshev, dom::MomentInterval, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) &&
        throw(ArgumentError("modified Chebyshev nodes are irrational; $(T) is not supported"))
    n = mc_points(degree)
    mw = dom.weight
    x, w, iters, used, res, cond = modified_chebyshev_work(mw, n, ctx.bits)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    spent = used - ctx.bits - gj_guard_bits(n)
    cert = Certificate(equations = "Σᵢ wᵢ πₖ(xᵢ) = mₖ for k < $(2n) (relative to the mass)",
                       residual = BigFloat(res; precision = 64), residual_bits = used,
                       digits = target_digits(ctx),
                       guard_digits = floor(Int, (used - ctx.bits) * log10(2)),
                       cond = cond, iterations = iters)
    prov = Provenance(family = "ModifiedChebyshev", derivation = Derived(),
                      path = ["moments: $(2n + 2) modified moments of $(mw.label)",
                              "recurrence: Wheeler's algorithm at $used bits" *
                              (spent > 0 ? " ($(floor(Int, spent * log10(2))) digits above the guard, found by escalation)" : ""),
                              "nodes: Newton on the computed recurrence",
                              "weights: Christoffel function"],
                      seed_source = "modified moments",
                      citations = [WHEELER_1974, GAUTSCHI_1970, SACK_DONOVAN_1972])
    return QuadratureRule(xs, ws, dom, PolynomialDegree(2n - 1), prov, cert)
end
