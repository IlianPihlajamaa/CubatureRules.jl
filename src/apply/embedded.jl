# Embedded rules: an error estimate at no extra integrand evaluations (PLAN §3.3).
#
# A pair of rules that share nodes gives two estimates from one set of function values, and
# their difference estimates the error of the finer one. That is a property of the rules
# themselves, so it belongs here.
#
# What does not belong here is a driver that picks rules until an answer converges. One
# existed until v0.5 (`integrate(f, domain; rtol)`, walking `RuleSequence`s) and was
# removed: it assumed that one family forms a clean sequence of rules, which seeded families
# with gaps and upper limits do not, it built and refined a new rule at every step, and it
# was the one place where construction happened out of sight. Adaptive integration is
# served elsewhere — QuadGK.jl in one dimension, HCubature.jl on boxes and
# HAdaptiveIntegration.jl on simplices.

"""
    EmbeddedRule(fine, coarse)

Two rules that share nodes: a fine rule and a coarser one whose nodes are a subset. The
difference of the two estimates the error of the fine one at **no extra integrand
evaluations**, which is how Gauss–Kronrod works and what nested families give for free.

`integrate(f, e)` returns the fine value; `integrate(f, e; error = true)` returns an
[`IntegrationResult`](@ref).
"""
struct EmbeddedRule{D,T,R<:QuadratureRule{D,T},C<:ExactnessClaim}
    fine::R
    coarse_weights::Vector{T}          # aligned with the fine rule's nodes; zero off the coarse rule
    coarse_claim::C                    # a claim, not a degree: the coarse rule may have none
end

function EmbeddedRule(fine::QuadratureRule{D,T}, coarse::QuadratureRule{D}) where {D,T}
    w = zeros(T, npoints(fine))
    tol = eltype(fine) <: AbstractFloat ? 1024 * eps(real(float(T))) : zero(T)
    for j in 1:npoints(coarse)
        xj = node_vector(coarse, j)
        i = findfirst(k -> _same_node(node_vector(fine, k), xj, tol), 1:npoints(fine))
        i === nothing && throw(ArgumentError("the coarse rule's node $(xj) is not a node of the fine rule, " *
                                             "so the two are not embedded"))
        w[i] = T(weights(coarse)[j])
    end
    return EmbeddedRule{D,T,typeof(fine),typeof(exactness(coarse))}(fine, w, exactness(coarse))
end

# Relative, not absolute: a double-exponential rule on an unbounded domain spreads its nodes
# over sixty orders of magnitude, and an absolute tolerance makes every node near the origin
# look like every other one — which would silently attach a coarse weight to the wrong node.
_same_node(a, b, tol) = all(i -> abs(a[i] - b[i]) <= tol * max(abs(a[i]), abs(b[i])), eachindex(a))

nodes(e::EmbeddedRule) = nodes(e.fine)
weights(e::EmbeddedRule) = weights(e.fine)
domain(e::EmbeddedRule) = domain(e.fine)
exactness(e::EmbeddedRule) = exactness(e.fine)
npoints(e::EmbeddedRule) = npoints(e.fine)
degree(e::EmbeddedRule) = degree(e.fine)
"The fine rule's degree, or -1 when it claims no polynomial degree."
_result_degree(r) = exactness(r) isa PolynomialDegree ? degree(r) : -1
family(e::EmbeddedRule) = family(e.fine)
provenance(e::EmbeddedRule) = provenance(e.fine)
Base.show(io::IO, e::EmbeddedRule) =
    print(io, "EmbeddedRule(", family(e), ", ", npoints(e), " points, ", describe(exactness(e)),
          ", embedding a rule with ", describe(e.coarse_claim), ")")

"""
    IntegrationResult

What an error-estimating call returns: the value, an estimate of its error, the number of
integrand evaluations, and the family and degree of the rule that produced it.
"""
struct IntegrationResult{V,E}
    value::V
    error_estimate::E
    neval::Int
    degree::Int
    family::String
end

function Base.show(io::IO, ::MIME"text/plain", r::IntegrationResult)
    println(io, "IntegrationResult")
    println(io, "  value      : ", r.value)
    println(io, "  error est. : ", r.error_estimate)
    println(io, "  evaluations: ", r.neval)
    print(io, "  rule       : ", r.family, ", degree ", r.degree)
end
Base.show(io::IO, r::IntegrationResult) = print(io, "IntegrationResult(", r.value, " ± ", r.error_estimate, ")")

function integrate(f::F, e::EmbeddedRule; error::Bool = false) where {F}
    x, w, wc = nodes(e.fine), weights(e.fine), e.coarse_weights
    @inbounds v = f(x[1])
    fine = w[1] * v
    coarse = wc[1] * v
    @inbounds for i in 2:length(w)
        v = f(x[i])
        fine += w[i] * v
        coarse += wc[i] * v
    end
    error || return fine
    return IntegrationResult(fine, abs(fine - coarse), npoints(e), _result_degree(e.fine), family(e))
end

"""
    embedded(family, domain; degree, T, digits)

The rule of `family` at `degree` together with an embedded coarser rule of the same family,
sharing nodes. Defined where a family is nested; Gauss–Kronrod carries its Gauss subset by
construction.
"""
function embedded(f::GaussKronrod, dom::Interval; degree::Integer, T = nothing, digits = nothing)
    Tout, bits = resolve_precision(T, digits)
    n = kronrod_halves(degree)
    guard = 32 + 4 * ceil(Int, log2(n + 1))
    fine = rule(f, dom; degree, T, digits)
    gauss = with_bits(bits + guard) do
        xg, wg = QuadGK.gauss(BigFloat, n)
        xg, wg
    end
    xs = [finalize_number(BuildContext{Tout}(bits), xi) for xi in gauss[1]]
    ws = [finalize_number(BuildContext{Tout}(bits), wi) for wi in gauss[2]]
    prov = Provenance(family = "GaussJacobi", derivation = Derived(), path = ["Gauss subset of the Kronrod rule"],
                      seed_source = "QuadGK.jl", citations = [QUADGK_JL], symmetry = :reflection)
    coarse = QuadratureRule(xs, ws, Interval(), PolynomialDegree(2n - 1), prov)
    return EmbeddedRule(fine, coarse)
end

"""
    integrate(f, domain; kwargs...)

Not supported: `integrate` takes a rule, built with [`rule`](@ref). This method exists only
to say so, and to point to packages that integrate adaptively to a tolerance.
"""
function integrate(f, dom::Domain; kwargs...)
    throw(ArgumentError(
        "integrate takes a rule, not a domain: build one first, e.g. " *
        "integrate(f, rule($(dom); degree = 20)). CubatureRules does not integrate adaptively " *
        "to a tolerance; for that, use QuadGK.jl (intervals), HCubature.jl (boxes) or " *
        "HAdaptiveIntegration.jl (simplices), or compare rules of two degrees, or use an " *
        "EmbeddedRule for an error estimate."))
end
