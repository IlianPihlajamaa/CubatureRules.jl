# Newton–Cotes rules (PLAN §6 Tier 1): equally spaced nodes, weights from exact integration
# of the Lagrange basis. Everything is done in Rational{BigInt} and rounded once, so
# `T = Rational{BigInt}` gives the exact rule — the form people want when writing papers.
#
# An n-point rule integrates degree n-1 by construction, and one more when n is odd, because
# the node set is symmetric and the extra odd moment vanishes.

"""
    NewtonCotes(; open = false)

Newton–Cotes rules on an interval: `n` equally spaced nodes, closed (endpoints included) or
open (endpoints excluded). Weights are exact rationals. Positive only for small `n` —
closed up to 8 points, open up to 2 — after which the weights alternate in sign and the
rules stop being usable in practice.
"""
Base.@kwdef struct NewtonCotes <: RuleFamily
    open::Bool = false
end

derivation(::Type{NewtonCotes}) = Derived()
describe_family(f::NewtonCotes) = f.open ? "NewtonCotes(open)" : "NewtonCotes(closed)"
supports_type(::NewtonCotes, ::Type{T}) where {T} = T <: AbstractFloat || T <: Rational

const NEWTON_COTES = Citation(key = "DavisRabinowitz1984", authors = ["Philip J. Davis", "Philip Rabinowitz"],
                              title = "Methods of Numerical Integration", journal = "Academic Press",
                              year = 1984)

candidates(::Type{NewtonCotes}, dom::Interval, ::PolynomialDegree) =
    isreference(dom) ? [NewtonCotes(false), NewtonCotes(true)] : NewtonCotes[]

"Points needed for degree `d`: an odd count, which buys one degree for free."
nc_points(f::NewtonCotes, d::Integer) = max(f.open ? 1 : 2, isodd(d) ? Int(d) : Int(d) + 1)
nc_degree(f::NewtonCotes, n::Integer) = isodd(n) ? n : n - 1

npoints(f::NewtonCotes, dom, degree::Integer) = nc_points(f, degree)
claimed_degree(f::NewtonCotes, dom, degree) = nc_degree(f, nc_points(f, degree))
degree_range(::NewtonCotes, dom::Interval) = 0:typemax(Int)
degree_range(::NewtonCotes, dom) = 1:0
properties(f::NewtonCotes, dom, degree) =
    (positive = nc_points(f, degree) <= (f.open ? 2 : 8), interior = f.open, symmetry = :reflection,
     nested = false)

"""
    newton_cotes_exact(n, open) -> (x, w)

Exact nodes and weights on `[-1, 1]`, as `Rational{BigInt}`. The weights are the integrals
of the Lagrange basis, obtained by expanding each basis polynomial by convolution.
"""
function newton_cotes_exact(n::Int, open::Bool)
    n >= 1 || throw(ArgumentError("need at least one node"))
    (open || n >= 2) || throw(ArgumentError("a closed Newton–Cotes rule needs at least two nodes"))
    x = open ? [-1 + 2 * big(i) // (n + 1) for i in 1:n] : (n == 1 ? [big(0) // 1] :
        [-1 + 2 * big(i - 1) // (n - 1) for i in 1:n])
    w = Vector{Rational{BigInt}}(undef, n)
    for i in 1:n
        # coefficients of Π_{j≠i} (t - x_j), lowest order first
        c = [one(Rational{BigInt})]
        for j in 1:n
            j == i && continue
            c = vcat(zero(Rational{BigInt}), c) .- vcat(x[j] .* c, zero(Rational{BigInt}))
        end
        denom = prod(x[i] - x[j] for j in 1:n if j != i; init = one(Rational{BigInt}))
        w[i] = sum(isodd(k - 1) ? zero(Rational{BigInt}) : c[k] * 2 // k for k in eachindex(c)) / denom
    end
    return x, w
end

function build(f::NewtonCotes, dom::Interval, degree::Int, ctx::BuildContext{T}) where {T}
    n = nc_points(f, degree)
    x, w = newton_cotes_exact(n, f.open)
    xs = [finalize_number(ctx, xi) for xi in x]
    ws = [finalize_number(ctx, wi) for wi in w]
    err = isexact(ctx) ? big(0.0) :
          with_bits(4ctx.bits) do
              maximum(abs(BigFloat(ws[i]) - BigFloat(w[i])) / max(abs(BigFloat(w[i])), big(1e-40)) for i in 1:n)
          end
    cert = Certificate(equations = "∫ ℓᵢ over [-1, 1] in exact rational arithmetic (relative rounding error)",
                       residual = BigFloat(err; precision = 64), residual_bits = isexact(ctx) ? 0 : 4ctx.bits,
                       digits = target_digits(ctx), guard_digits = 0)
    prov = Provenance(family = "NewtonCotes", derivation = Derived(),
                      path = ["$(n) equally spaced nodes, $(f.open ? "open" : "closed")",
                              "weights: exact integration of the Lagrange basis in Rational{BigInt}",
                              isexact(ctx) ? "exact output" : "rounded once to $T"],
                      seed_source = "none (closed form)", citations = [NEWTON_COTES], symmetry = :reflection)
    return QuadratureRule(xs, ws, Interval(), PolynomialDegree(nc_degree(f, n)), prov, cert)
end
