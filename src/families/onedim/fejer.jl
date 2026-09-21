# Fejér rules (PLAN §6 Tier 1): interpolatory rules on Chebyshev nodes, with closed-form
# weights, in generic arithmetic at any precision.
#
#   type 1   x_k = cos((2k-1)π/2n),  k = 1..n      (the Chebyshev points of the first kind)
#   type 2   x_k = cos(kπ/(n+1)),    k = 1..n      (interior Chebyshev points of the second
#                                                   kind — Clenshaw–Curtis without endpoints)
#
# Both have positive weights at every order and all nodes interior, which is what makes them
# useful where Newton–Cotes fails. Both node sets are symmetric about the origin with
# symmetric weights, so an odd number of points integrates one degree more than
# interpolation alone guarantees.

"""
    Fejer(kind = 1)

Fejér quadrature on Chebyshev nodes: `kind = 1` uses the Chebyshev points of the first kind,
`kind = 2` the interior points of the second kind. `n` points integrate degree `n - 1`, or
`n` when `n` is odd. Positive weights and interior nodes at every order.
"""
struct Fejer <: RuleFamily
    kind::Int
    function Fejer(kind::Integer = 1)
        kind in (1, 2) || throw(ArgumentError("Fejér quadrature has kind 1 or 2, got $kind"))
        return new(Int(kind))
    end
end

derivation(::Type{Fejer}) = Derived()
describe_family(f::Fejer) = "Fejer($(f.kind))"

const FEJER_1933 = Citation(key = "Fejer1933", authors = ["Leopold Fejér"],
                            title = "Mechanische Quadraturen mit positiven Cotesschen Zahlen",
                            journal = "Mathematische Zeitschrift", year = 1933, volume = "37",
                            pages = "287--309", doi = "10.1007/BF01474575")

candidates(::Type{Fejer}, dom::Interval, ::PolynomialDegree) =
    isreference(dom) ? [Fejer(1), Fejer(2)] : Fejer[]

"An odd number of points, which buys one degree for free."
fejer_points(d::Integer) = max(1, isodd(d) ? Int(d) : Int(d) + 1)
fejer_degree(n::Integer) = isodd(n) ? n : n - 1

npoints(::Fejer, dom, degree::Integer) = fejer_points(degree)
claimed_degree(::Fejer, dom, degree) = fejer_degree(fejer_points(degree))
degree_range(::Fejer, dom::Interval) = 0:typemax(Int)
degree_range(::Fejer, dom) = 1:0
properties(::Fejer, dom, degree) = (positive = true, interior = true, symmetry = :reflection, nested = false)

"""
    fejer_work(n, kind, bits) -> (x, w)

Fejér nodes (ascending) and weights on `[-1, 1]` as BigFloats of precision `bits`, from the
closed forms

    type 1:  w_k = (2/n) (1 − 2 Σ_{j=1}^{⌊n/2⌋} cos(2jθ_k) / (4j² − 1))
    type 2:  w_k = (4 sin θ_k / (n+1)) Σ_{j=1}^{⌈(n+1)/2⌉} sin((2j−1)θ_k) / (2j − 1)
"""
function fejer_work(n::Integer, kind::Int, bits::Integer)
    n >= 1 || throw(ArgumentError("need at least one node"))
    return with_bits(bits) do
        π_ = BigFloat(π)
        θ = kind == 1 ? [(2k - 1) * π_ / (2n) for k in 1:n] : [k * π_ / (n + 1) for k in 1:n]
        x = cos.(θ)
        w = Vector{BigFloat}(undef, n)
        for k in 1:n
            if kind == 1
                s = sum((cos(2j * θ[k]) / (4j^2 - 1) for j in 1:(n ÷ 2)); init = zero(π_))
                w[k] = 2 * (1 - 2s) / n
            else
                s = sum((sin((2j - 1) * θ[k]) / (2j - 1) for j in 1:cld(n + 1, 2)); init = zero(π_))
                w[k] = 4 * sin(θ[k]) * s / (n + 1)
            end
        end
        p = sortperm(x)
        x[p], w[p]
    end
end

function build(f::Fejer, dom::Interval, degree::Int, ctx::BuildContext{T}) where {T}
    isexact(ctx) && throw(ArgumentError("Fejér nodes are irrational; $(T) is not supported"))
    n = fejer_points(degree)
    guard = 32 + 4 * ceil(Int, log2(n + 1))
    x, w = fejer_work(n, f.kind, ctx.bits + guard)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    # the closed form is exact, so the only error is the rounding of the delivered numbers
    xhi, whi = fejer_work(n, f.kind, 2ctx.bits + guard)
    err = with_bits(2ctx.bits + guard) do
        maximum(max(abs(BigFloat(xs[i]) - xhi[i]), abs(BigFloat(ws[i]) - whi[i]) / whi[i]) for i in 1:n)
    end
    cert = Certificate(equations = "closed-form Fejér $(f.kind) nodes and weights (rounding error of the delivered rule)",
                       residual = BigFloat(err; precision = 64), residual_bits = 2ctx.bits + guard,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)))
    prov = Provenance(family = "Fejer", derivation = Derived(),
                      path = ["$(n) Chebyshev points of the $(f.kind == 1 ? "first" : "second") kind",
                              "weights: closed form at $(ctx.bits + guard) bits, rounded once to $T"],
                      seed_source = "none (closed form)", citations = [FEJER_1933], symmetry = :reflection)
    return QuadratureRule(xs, ws, Interval(), PolynomialDegree(fejer_degree(n)), prov, cert)
end
