# The registry and selector (PLAN §2.5). Discovery is reflective; instantiation is explicit.
#
# `subtypes` is called inside functions, at run time, on every query. It must never be
# hoisted into a `const`: that would be evaluated at precompile time and silently omit every
# family a downstream package defines. The test suite adds a family from a separately
# precompiled package and asserts that it is found.

"""
    families() -> Vector

Every concrete rule-family type currently loaded (from this package or any other), found
by walking `subtypes(RuleFamily)` through abstract layers. Sorted by name, because
`subtypes`' order is unspecified.
"""
function families()
    out = Any[]
    stack = Any[RuleFamily]
    while !isempty(stack)
        T = pop!(stack)
        for S in InteractiveUtils.subtypes(T)
            isabstracttype(S) ? push!(stack, S) : push!(out, S)
        end
    end
    return sort!(out; by = F -> string(F))
end

"Leaf families: every loaded family that is not a combinator (combinator depth is 1)."
leaf_families() = filter(F -> !(F <: CombinatorFamily), families())

"""
    Candidate

One constructible option for a request, with the metadata used to rank and explain it.
"""
struct Candidate
    family::RuleFamily
    name::String
    npoints::Int
    degree::Int
    derivation::DerivationClass
    positive::Bool
    interior::Bool
    symmetry::Symbol
    cost::Float64
end

function Candidate(f::RuleFamily, dom, degree::Integer, T)
    p = properties(f, dom, degree)
    return Candidate(f, describe_family(f), npoints(f, dom, degree), claimed_degree(f, dom, degree),
                     derivation(f), p.positive, p.interior, p.symmetry, cost_estimate(f, dom, degree, T))
end

"Documented total order: node count, then derived before seeded, then family name."
rank_key(c::Candidate) = (c.npoints, c.derivation isa Derived ? 0 : 1, c.name)

function reference_domain(dom::Domain)
    applicable(reference, dom) || return nothing
    return reference(dom)
end

"All candidates for `claim` on the reference of `dom`, ranked, before filtering."
function gather(dom::Domain, degree::Integer, T)
    ref = reference_domain(dom)
    ref === nothing && return Candidate[]
    claim = PolynomialDegree(degree)
    out = Candidate[]
    for F in families()
        for f in candidates(F, ref, claim)
            degree in degree_range(f, ref) || continue
            push!(out, Candidate(f, ref, degree, T))
        end
    end
    return sort!(out; by = rank_key)
end

passes(c::Candidate; positive, interior) = (!positive || c.positive) && (!interior || c.interior)

"""
    NoRuleError

Thrown when no loaded family can satisfy a request. The message says why, and what is
nearby.
"""
struct NoRuleError <: Exception
    msg::String
end
Base.showerror(io::IO, e::NoRuleError) = print(io, e.msg)

# ---------------------------------------------------------------------------------------
# rule(...)

"""
    rule(domain; degree, T = Float64, digits, positive = false, interior = false, family, cancel, seed)
    rule(family, domain; degree, npoints, T, digits, cancel, seed)

Construct a quadrature rule of at least polynomial degree `degree` on `domain`.

Without `family`, every loaded family is asked for candidates, which are ranked by node
count, then derived-before-seeded, then name; the first one satisfying the filters is
built, and the choice is recorded in `provenance(rule).selection`.

Precision: `digits = 200` gives `BigFloat` at 200 decimal digits; `T = BigFloat` uses the
ambient BigFloat precision at the time of the call; any other `T` its own precision.
`T = Rational{BigInt}` selects families with exact rational rules.

`cancel` takes a [`CancellationToken`](@ref); `seed` a seed source for seeded families.

There is no default degree: `rule(domain)` errors, listing what is available.
"""
function rule(dom::Domain; degree = nothing, npoints = nothing, T = nothing, digits = nothing,
              positive::Bool = false, interior::Bool = false, family = nothing, cancel = nothing, seed = nothing)
    if family !== nothing
        return rule(family, dom; degree, npoints, T, digits, cancel, seed)
    end
    npoints === nothing || throw(ArgumentError("`npoints` requires an explicit `family`"))
    degree === nothing && throw(NoRuleError(no_degree_message(dom)))
    degree >= 0 || throw(ArgumentError("degree must be non-negative"))
    Tout, bits = resolve_precision(T, digits)
    cands = gather(dom, degree, Tout)
    ok = filter(c -> passes(c; positive, interior) && supports_type(c.family, Tout), cands)
    isempty(ok) && throw(NoRuleError(unsatisfiable_message(dom, degree, Tout, positive, interior, cands)))
    chosen = first(ok)
    selection = "selected by rule(): " * join(("$(c.name) ($(c.npoints) points)" for c in ok), " < ") *
                "; ranked by (npoints, derived before seeded, family name)" *
                (positive ? "; positive = true" : "") * (interior ? "; interior = true" : "")
    r = _build(chosen.family, dom, degree, Tout, bits, cancel, seed)
    r = QuadratureRule(r.nodes, r.weights, r.domain, r.exactness, with_selection(r.provenance, selection), r.certificate)
    return isreference(dom) ? r : map_to(r, dom)
end

function rule(f::RuleFamily, dom::Domain; degree = nothing, npoints = nothing, T = nothing, digits = nothing,
              cancel = nothing, seed = nothing, positive::Bool = false, interior::Bool = false)
    ref = reference_domain(dom)
    ref === nothing && throw(NoRuleError("$(describe_family(f)) does not support $(dom)"))
    if npoints !== nothing
        degree === nothing || throw(ArgumentError("give either `degree` or `npoints`, not both"))
        applicable(degree_for_npoints, f, ref, npoints) ||
            throw(ArgumentError("$(describe_family(f)) does not take an `npoints` request"))
        degree = degree_for_npoints(f, ref, npoints)
    end
    if degree === nothing
        needs_degree(f) && throw(NoRuleError(no_degree_message(dom; only = f)))
        degree = 0          # the family is parameterised by something other than degree
    end
    Tout, bits = resolve_precision(T, digits)
    dep = missing_dependency(f)
    dep === nothing ||
        throw(NoRuleError("$(describe_family(f)) builds its rules with $(dep), which is not loaded; " *
                          "run `using $(chopsuffix(dep, ".jl"))` first"))
    if needs_degree(f)
        isempty(candidates(typeof(f), ref, PolynomialDegree(degree))) &&
            throw(NoRuleError(unsatisfiable_message(dom, degree, Tout, false, false, Candidate[]; only = f)))
        degree in degree_range(f, ref) ||
            throw(NoRuleError(unsatisfiable_message(dom, degree, Tout, false, false, Candidate[]; only = f)))
    end
    supports_type(f, Tout) ||
        throw(NoRuleError("$(describe_family(f)) cannot deliver nodes and weights in $Tout"))
    c = Candidate(f, ref, degree, Tout)
    (positive && !c.positive) && throw(NoRuleError("$(c.name) at degree $degree does not have positive weights"))
    (interior && !c.interior) && throw(NoRuleError("$(c.name) at degree $degree has boundary nodes"))
    r = _build(f, dom, degree, Tout, bits, cancel, seed)
    r = QuadratureRule(r.nodes, r.weights, r.domain, r.exactness,
                       with_selection(r.provenance, "family chosen explicitly by the caller"), r.certificate)
    return isreference(dom) ? r : map_to(r, dom)
end

function _build(f, dom, degree, T, bits, cancel, seed)
    ctx = BuildContext{T}(bits; cancel)
    ref = reference(dom)
    return seed === nothing ? build(f, ref, Int(degree), ctx) : build(f, ref, Int(degree), ctx; seed)
end

# ---------------------------------------------------------------------------------------
# Diagnostics: say precisely why, and what is nearby.

_range_string(r) = isempty(r) ? "none" : last(r) == typemax(Int) ? "$(first(r))–∞" : "$(first(r))–$(last(r))"

"Every family instance applicable to the domain at some degree, for diagnostics."
function applicable_instances(dom::Domain)
    ref = reference_domain(dom)
    ref === nothing && return Any[]
    out = Any[]
    for F in families(), f in candidates(F, ref, PolynomialDegree(0))
        push!(out, f)
    end
    # seeded families may not cover degree 0 but still exist on the domain
    for F in families(), f in candidates(F, ref, PolynomialDegree(1))
        any(g -> g == f, out) || push!(out, f)
    end
    return sort!(out; by = describe_family)
end

"""
    claimless_instances(dom) -> families

Families that live on `dom` but are parameterised by something other than a degree. They
are exact on no polynomial space, so they never appear among the candidates for a degree —
which is exactly why an explanation of "no rule found" has to go and look for them.
"""
function claimless_instances(dom::Domain)
    ref = reference_domain(dom)
    ref === nothing && return Any[]
    out = Any[]
    for F in families()
        hasmethod(F, Tuple{}) || continue
        # this runs while building an error message, so a downstream family that cannot be
        # default-constructed must not replace the error the caller is about to see
        f = try
            F()
        catch
            continue
        end
        (!needs_degree(f) && home_domain(f) == ref) && push!(out, f)
    end
    return sort!(out; by = describe_family)
end

# Named, not offered: a NoClaim family can still be what the caller wants.
function _claimless_hint(io::IO, dom::Domain)
    fs = claimless_instances(dom)
    isempty(fs) && return
    names = [string(nameof(typeof(f))) for f in fs]
    one = length(fs) == 1
    print(io, "\n", join(names, ", "), one ? " lives on " : " live on ", dom,
          one ? " but claims" : " but claim", " no polynomial degree, so ",
          one ? "it is" : "they are", " never offered for one. Ask by name: ",
          "rule(", first(names), "(4), ", dom, ").")
end

function no_degree_message(dom::Domain; only = nothing)
    io = IOBuffer()
    print(io, "rule(", dom, ") needs a `degree`: there is no defensible default, since the right degree ",
          "depends on the integrand.")
    insts = only === nothing ? applicable_instances(dom) : [only]
    ref = reference_domain(dom)
    if isempty(insts) || ref === nothing
        print(io, " No loaded family answers a degree request on ", dom, ".", _planned(dom))
        _claimless_hint(io, dom)
    else
        print(io, " Available on ", dom, ":")
        for f in insts
            print(io, "\n  ", rpad(describe_family(f), 28), " degrees ", _range_string(degree_range(f, ref)),
                  "  (", derivation(f) isa Derived ? "derived" : "seeded", ")")
        end
        print(io, "\nExample: rule(", dom, "; degree = 10)")
    end
    return String(take!(io))
end

_planned(dom) = ""
_planned(::Union{Ball,Polytope,Wedge,Pyramid}) = " That domain is scheduled for a later release."

function unsatisfiable_message(dom, degree, T, positive, interior, cands; only = nothing)
    io = IOBuffer()
    reqs = String[]
    positive && push!(reqs, "positive weights")
    interior && push!(reqs, "interior nodes")
    T <: Rational && push!(reqs, "exact $(T) nodes and weights")
    what = isempty(reqs) ? "rule" : "rule with " * join(reqs, " and ")
    who = only === nothing ? "" : " from $(describe_family(only))"
    print(io, "No ", what, who, " exists on ", dom, " at degree ", degree, " among the loaded families.")
    ref = reference_domain(dom)
    if ref === nothing
        print(io, _planned(dom))
        return String(take!(io))
    end
    # reasons for candidates at this degree that were filtered out
    for c in cands
        why = String[]
        (positive && !c.positive) && push!(why, "has negative weights")
        (interior && !c.interior) && push!(why, "has boundary nodes")
        supports_type(c.family, T) || push!(why, "cannot produce $(T)")
        isempty(why) || print(io, "\n  ", c.name, " (", c.npoints, " points) ", join(why, " and "), ".")
    end
    # nearest alternatives at other degrees
    near = String[]
    insts = only === nothing ? applicable_instances(dom) : [only]
    for f in insts
        rg = degree_range(f, ref)
        (isempty(rg) || !supports_type(f, T)) && continue
        # the closest degree (lower first on ties) at which this family meets the filters
        c = nothing
        for δ in 1:64, d in (degree - δ, degree + δ)
            d in rg || continue
            cd = Candidate(f, ref, d, T)
            passes(cd; positive, interior) && (c = cd; break)
        end
        c === nothing && continue
        push!(near, "$(c.name) degree $(c.degree) ($(c.npoints) points" *
                    (last(rg) == typemax(Int) ? ", always available" : "") * ")")
    end
    isempty(near) || print(io, "\nNearest: ", join(near, ", or "), ".")
    only === nothing && _claimless_hint(io, dom)
    return String(take!(io))
end

# ---------------------------------------------------------------------------------------
# available / compare

"""
    CandidateTable

The result of [`available`](@ref) and [`compare`](@ref): a vector of rows (NamedTuples)
with a readable `show`.
"""
struct CandidateTable
    title::String
    rows::Vector{NamedTuple}
end
Base.length(t::CandidateTable) = length(t.rows)
Base.getindex(t::CandidateTable, i) = t.rows[i]
Base.iterate(t::CandidateTable, s...) = iterate(t.rows, s...)
Base.isempty(t::CandidateTable) = isempty(t.rows)

function Base.show(io::IO, ::MIME"text/plain", t::CandidateTable)
    println(io, t.title)
    isempty(t.rows) && return print(io, "  (none)")
    ks = keys(first(t.rows))
    cols = vcat([[string(k) for k in ks]], [[_cell(r[k]) for k in ks] for r in t.rows])
    widths = [maximum(length(c[j]) for c in cols) for j in eachindex(ks)]
    for (i, c) in enumerate(cols)
        print(io, "  ", join((rpad(c[j], widths[j]) for j in eachindex(ks)), "  "))
        i < length(cols) && println(io)
    end
end
_cell(x::Bool) = x ? "yes" : "no"
_cell(x::DerivationClass) = x isa Derived ? "derived" : "seeded"
_cell(x::UnitRange) = _range_string(x)
_cell(x) = string(x)

"""
    available(domain; degree, positive = false, interior = false, T = Float64)

Without `degree`: every loaded family applicable to `domain`, with its degree range.
With `degree`: the candidates that satisfy the request, ranked exactly as [`rule`](@ref)
would rank them — nothing is constructed.
"""
function available(dom::Domain; degree = nothing, positive::Bool = false, interior::Bool = false, T = Float64)
    ref = reference_domain(dom)
    if degree === nothing
        rows = NamedTuple[]
        ref === nothing && return CandidateTable("No families for $(dom).", rows)
        for f in applicable_instances(dom)
            p = properties(f, ref, max(first(degree_range(f, ref)), 1))
            push!(rows, (family = describe_family(f), degrees = degree_range(f, ref), derivation = derivation(f),
                         symmetry = p.symmetry))
        end
        return CandidateTable("Families available on $(dom):", rows)
    end
    cands = filter(c -> passes(c; positive, interior) && supports_type(c.family, T), gather(dom, degree, T))
    rows = NamedTuple[(family = c.name, npoints = c.npoints, degree = c.degree, derivation = c.derivation,
                       positive = c.positive, interior = c.interior, symmetry = c.symmetry) for c in cands]
    return CandidateTable("Candidates on $(dom) for degree ≥ $(degree), ranked:", rows)
end

"""
    compare(domain, degree)

Side-by-side properties of every candidate at `degree`, including those `available` would
filter out, with a rough construction-cost estimate.
"""
function compare(dom::Domain, degree::Integer; T = Float64)
    cands = gather(dom, degree, T)
    rows = NamedTuple[(family = c.name, npoints = c.npoints, degree = c.degree, derivation = c.derivation,
                       positive = c.positive, interior = c.interior, symmetry = c.symmetry,
                       exact_rational = supports_type(c.family, Rational{BigInt}),
                       cost = @sprintf("%.3g", c.cost)) for c in cands]
    return CandidateTable("Comparison on $(dom) at degree $(degree):", rows)
end
