# Rule sequences, embedded error estimates and accuracy-targeted integration
# (PLAN §3.3, §3.4, §5 item 7).
#
# One abstraction buys three things: convergence studies, error estimates at no extra
# integrand evaluations, and the `rtol` interface below. Note what this is *not*: it is
# order-adaptive, walking a sequence of rules of increasing degree over the whole domain,
# not space-adaptive. For integrands with localised features, QuadGK.jl and HCubature.jl
# subdivide space and will do better; this fills the gap in dimensions ≥ 2, where nothing
# offers order-adaptive high-precision integration.

"""
    RuleSequence(family, domain; degrees, T, digits)

A lazy sequence of rules of increasing degree. Iterating constructs each rule in turn:

```julia
for r in RuleSequence(ConicalProduct(), Simplex{2}())
    @show degree(r), integrate(f, r)
end
```

The default schedule roughly doubles the degree each step, clipped to what the family
offers. Requested degrees that the family would answer with the same rule are dropped, so
the sequence never yields one rule twice.
"""
struct RuleSequence{F<:RuleFamily,D<:Domain}
    family::F
    domain::D
    degrees::Vector{Int}
    T::Type
    bits::Int
end

function RuleSequence(f::RuleFamily, dom::Domain; degrees = nothing, T = nothing, digits = nothing,
                      maxdegree::Integer = 64)
    Tout, bits = resolve_precision(T, digits)
    ds = degrees === nothing ? default_schedule(f, dom, maxdegree) : collect(Int, degrees)
    # A family rounds a requested degree up to one it has, so two requests can land on the
    # same rule — GaussKronrod answers both 1 and 3 with its 3-point rule. Walking that rule
    # twice makes the successive difference exactly zero, and a tolerance-driven integrate
    # then reports convergence on the first rule it built: exp on [-1,1] came back 6.5e-5
    # wrong with an error estimate of 0.0. Keep only requests that reach a new rule.
    ref = reference_domain(dom)
    if ref !== nothing
        ds = unique(d -> claimed_degree(f, ref, d), ds)
    end
    isempty(ds) && throw(ArgumentError("$(describe_family(f)) offers no degrees on $(dom)"))
    return RuleSequence{typeof(f),typeof(dom)}(f, dom, ds, Tout, bits)
end

"Degrees to walk: 1, 3, 7, 15, … within the family's range."
function default_schedule(f::RuleFamily, dom::Domain, maxdegree::Integer)
    rg = degree_range(f, reference_domain(dom))
    ds = Int[]
    d = max(1, first(rg))
    while d <= min(maxdegree, last(rg))
        push!(ds, d)
        d = 2d + 1
    end
    return ds
end

"""
    LevelSequence(family, domain; levels, T, digits)

The counterpart of [`RuleSequence`](@ref) for a family parameterised by a level rather than
a degree — the double-exponential rules. The sequence walks `TanhSinh(m)`, `ExpSinh(m)` or
`SinhSinh(m)` for increasing `m`, which roughly doubles the point count each step and keeps
every node of the previous step.
"""
struct LevelSequence{F<:RuleFamily,D<:Domain}
    families::Vector{F}
    domain::D
    T::Type
    bits::Int
end

function LevelSequence(f::RuleFamily, dom::Domain; levels = 2:8, T = nothing, digits = nothing)
    needs_degree(f) && throw(ArgumentError("$(describe_family(f)) is parameterised by a degree; " *
                                           "use RuleSequence"))
    Tout, bits = resolve_precision(T, digits)
    F = typeof(f)
    return LevelSequence{F,typeof(dom)}([F(m) for m in levels], dom, Tout, bits)
end

Base.length(s::LevelSequence) = length(s.families)
Base.eltype(::Type{<:LevelSequence}) = QuadratureRule
function Base.iterate(s::LevelSequence, i::Int = 1)
    i > length(s.families) && return nothing
    r = rule(s.families[i], s.domain; T = s.T === BigFloat ? nothing : s.T,
             digits = s.T === BigFloat ? floor(Int, s.bits * log10(2)) : nothing)
    return r, i + 1
end
Base.show(io::IO, s::LevelSequence) =
    print(io, "LevelSequence(", describe_family(first(s.families)), " on ", s.domain,
          ", levels ", [f.level for f in s.families], ")")

Base.length(s::RuleSequence) = length(s.degrees)
Base.eltype(::Type{<:RuleSequence}) = QuadratureRule
function Base.iterate(s::RuleSequence, i::Int = 1)
    i > length(s.degrees) && return nothing
    r = rule(s.family, s.domain; degree = s.degrees[i], T = s.T === BigFloat ? nothing : s.T,
             digits = s.T === BigFloat ? floor(Int, s.bits * log10(2)) : nothing)
    return r, i + 1
end
Base.show(io::IO, s::RuleSequence) =
    print(io, "RuleSequence(", describe_family(s.family), " on ", s.domain, ", degrees ", s.degrees, ")")

# ---------------------------------------------------------------------------------------

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

What an accuracy-targeted or error-estimating call returns: never a bare number.
"""
struct IntegrationResult{V,E}
    value::V
    error_estimate::E
    neval::Int
    degree::Int
    family::String
    converged::Bool
end

function Base.show(io::IO, ::MIME"text/plain", r::IntegrationResult)
    println(io, "IntegrationResult", r.converged ? "" : " (tolerance not met)")
    println(io, "  value     : ", r.value)
    println(io, "  error est.: ", r.error_estimate)
    println(io, "  evaluations: ", r.neval)
    print(io, "  rule      : ", r.family, ", degree ", r.degree)
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
    return IntegrationResult(fine, abs(fine - coarse), npoints(e), _result_degree(e.fine), family(e), true)
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

# ---------------------------------------------------------------------------------------

"""
    integrate(f, domain; rtol, atol = 0, family, T, digits, maxdegree = 64)

Accuracy-targeted integration: walk a [`RuleSequence`](@ref) until successive estimates
agree to the tolerance, and return an [`IntegrationResult`](@ref).

This is the one signature that takes a domain rather than a rule, and the exception is
principled: a tolerance inherently requires a *sequence* of rules, so construction is the
algorithm rather than an accidental cost (PLAN §3.1).

It is **order**-adaptive, not space-adaptive: for integrands with localised features use
`QuadGK.jl` or `HCubature.jl` instead.
"""
function integrate(f::F, dom::Domain; rtol = nothing, atol = 0, family = nothing, T = nothing,
                   digits = nothing, maxdegree::Integer = 64) where {F}
    rtol === nothing && throw(ArgumentError(
        "integrate(f, domain) needs an accuracy target: integrate(f, domain; rtol = 1e-12). " *
        "For a fixed rule, build it first — integrate(f, rule(domain; degree = 20)) — so that the " *
        "construction is visible and can be hoisted out of a loop."))
    fam = family === nothing ? adaptive_family(dom, maxdegree) : family
    # a family with no degree to walk is walked by level instead — the unbounded domains
    # have nothing else, and on them this is the whole of tolerance-driven integration
    seq = needs_degree(fam) ? RuleSequence(fam, dom; T, digits, maxdegree) :
          LevelSequence(fam, dom; T, digits)
    prev = nothing
    neval = 0
    best = nothing
    lasterr = nothing
    for r in seq
        v = integrate(f, r)
        neval += npoints(r)
        if prev !== nothing
            err = abs(v - prev)
            lasterr = err
            if err <= max(atol, rtol * abs(v))
                # `family` is the keyword here, so the rule reports its own name
                return IntegrationResult(v, err, neval, _result_degree(r), provenance(r).family, true)
            end
        end
        prev = v
        best = (v, _result_degree(r), provenance(r).family)
    end
    # not converged: report the last difference achieved rather than nothing useful
    v, d, fname = best
    return IntegrationResult(v, lasterr === nothing ? oftype(abs(v), Inf) : lasterr, neval, d, fname, false)
end

"The family an `rtol` call walks: the cheapest candidate at a mid-range degree, nested first."
function adaptive_family(dom::Domain, maxdegree::Integer)
    cands = gather(dom, min(15, maxdegree), Float64)
    if isempty(cands)
        # no family answers a degree here; a double-exponential one may still converge
        cl = claimless_instances(dom)
        isempty(cl) && throw(NoRuleError(no_degree_message(dom)))
        return first(cl)
    end
    ref = reference_domain(dom)
    nested = filter(c -> properties(c.family, ref, 15).nested, cands)
    return first(isempty(nested) ? cands : nested).family
end
