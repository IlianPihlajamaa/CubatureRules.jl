# Grundmann–Möller rules on the D-simplex (PLAN §6 Tier 1): explicit formula, exact
# rationals, any dimension, odd degree d = 2s + 1.
#
#     ∫ f ≈ Σ_{i=0}^{s} (-1)^i 2^{-2s} (d+D-2i)^d / (i! (d+D-i)!)
#               Σ_{|β| = s-i} f(λ = (2β + 1) / (d + D - 2i))
#
# over the reference simplex (measure 1/D!), β ranging over (D+1)-part compositions.
# Weights alternate in sign for s ≥ 1.

"""
    GrundmannMöller()

Grundmann–Möller rules on a simplex of any dimension, at odd degree `2s + 1`, with exact
rational nodes and weights (`T = Rational{BigInt}` is supported). Fully symmetric, all
nodes interior; weights are negative for degree ≥ 3. Also available as `GrundmannMoeller`.
"""
struct GrundmannMöller <: RuleFamily end
const GrundmannMoeller = GrundmannMöller

derivation(::Type{GrundmannMöller}) = Derived()
supports_type(::GrundmannMöller, ::Type{T}) where {T} = T <: AbstractFloat || T <: Rational

const GRUNDMANN_MOLLER_1978 = Citation(
    key = "GrundmannMoller1978", authors = ["Axel Grundmann", "H. M. Möller"],
    title = "Invariant integration formulas for the n-simplex by combinatorial methods",
    journal = "SIAM Journal on Numerical Analysis", year = 1978, volume = "15", pages = "282--290",
    doi = "10.1137/0715019")

candidates(::Type{GrundmannMöller}, dom::Simplex, ::PolynomialDegree) =
    isreference(dom) ? [GrundmannMöller()] : GrundmannMöller[]

gm_s(degree) = max(0, cld(degree - 1, 2))
npoints(::GrundmannMöller, dom::Simplex{D}, degree::Integer) where {D} =
    sum(i -> binomial(gm_s(degree) - i + D, D), 0:gm_s(degree))
claimed_degree(::GrundmannMöller, dom, degree) = 2gm_s(degree) + 1
degree_range(::GrundmannMöller, dom::Simplex) = 0:typemax(Int)
properties(::GrundmannMöller, dom, degree) =
    (positive = gm_s(degree) == 0, interior = true, symmetry = Symbol("S", dimension(dom) + 1), nested = false)

"All compositions of `total` into `parts` non-negative parts, in lexicographic order."
function compositions(total::Int, parts::Int)
    parts == 1 && return [[total]]
    out = Vector{Vector{Int}}()
    for first in total:-1:0, rest in compositions(total - first, parts - 1)
        push!(out, vcat(first, rest))
    end
    return out
end

"""
    grundmann_moller_exact(D, s) -> (λs, ws)

Exact barycentric nodes and weights (`Rational{BigInt}`) of the degree-`2s+1` rule on the
reference `D`-simplex.
"""
function grundmann_moller_exact(D::Int, s::Int)
    d = 2s + 1
    λs = Vector{Vector{Rational{BigInt}}}()
    ws = Rational{BigInt}[]
    for i in 0:s
        den = big(d + D - 2i)
        w = (isodd(i) ? -1 : 1) * den^d // (big(2)^(2s) * factorial(big(i)) * factorial(big(d + D - i)))
        for β in compositions(s - i, D + 1)
            push!(λs, [(2b + 1) // den for b in β])
            push!(ws, w)
        end
    end
    return λs, ws
end

function build(::GrundmannMöller, dom::Simplex{D}, degree::Int, ctx::BuildContext{T}) where {D,T}
    s = gm_s(degree)
    λs, ws = grundmann_moller_exact(D, s)
    xs = [SVector{D,T}(ntuple(j -> finalize_number(ctx, λ[j + 1]), D)) for λ in λs]
    wt = [finalize_number(ctx, w) for w in ws]
    err = isexact(ctx) ? big(0.0) :
          with_bits(4ctx.bits) do
              maximum(abs(BigFloat(wt[i]) - BigFloat(ws[i])) / abs(BigFloat(ws[i])) for i in eachindex(ws))
          end
    cert = Certificate(equations = "closed form (relative rounding error of the weights)",
                       residual = BigFloat(err; precision = 64), residual_bits = isexact(ctx) ? 0 : 4ctx.bits,
                       digits = target_digits(ctx), guard_digits = 0)
    prov = Provenance(family = "GrundmannMöller", derivation = Derived(),
                      path = ["closed form in exact Rational{BigInt}", isexact(ctx) ? "exact output" : "rounded once to $T"],
                      seed_source = "none (explicit formula)", citations = [GRUNDMANN_MOLLER_1978], symmetry = Symbol("S", D + 1))
    return QuadratureRule(xs, wt, Simplex{D}(), PolynomialDegree(2s + 1), prov, cert)
end
