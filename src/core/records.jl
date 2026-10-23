# Metadata carried by a constructed rule: where it came from (Provenance), whether the
# construction converged (Certificate) and whether it is what it claims
# (Verification). Certification and verification are different checks with different types.

"""
    Citation

A bibliographic reference. `cite(rule)` renders these. `kind` is `:article` (a journal
paper), `:book`, or `:software`, a package or data set, whose `journal` names where it is
published (GitHub, Zenodo) and `url` where; left at `:auto` it is an article when it has a
volume or pages, and a book otherwise.
"""
Base.@kwdef struct Citation
    key::String
    authors::Vector{String}
    title::String
    journal::String = ""
    year::Int
    volume::String = ""
    pages::String = ""
    doi::String = ""
    kind::Symbol = :auto
    url::String = ""
end

# The package itself: the source of every rule in its own seed tables, which its searches
# produced. CITATION.cff at the repository root says the same for GitHub.
const CUBATURERULES_JL = Citation(key = "Pihlajamaa_CubatureRules", authors = ["Ilian Pihlajamaa"],
                                  title = "CubatureRules.jl: quadrature and cubature rules computed on demand, " *
                                          "at any degree and precision",
                                  journal = "GitHub", year = 2026, kind = :software,
                                  url = "https://github.com/IlianPihlajamaa/CubatureRules.jl")

"""
    Derived()
    Seeded()

The two derivation classes. A *derived* rule is constructed from first
principles at any order; a *seeded* rule is refined from a starting point that exists only
at tabulated orders.
"""
abstract type DerivationClass end

"A family whose rules are computed from a formula or a well-conditioned algorithm, at any degree."
struct Derived <: DerivationClass end

"A family whose rules are refined from stored starting values, and so exist only at tabulated degrees."
struct Seeded <: DerivationClass end
Base.show(io::IO, ::Derived) = print(io, "Derived()")
Base.show(io::IO, ::Seeded) = print(io, "Seeded()")

"""
    Provenance

Where a rule came from: family, derivation class and path, seed source, citations, the
licence of any data used, and — when chosen by the selector — the selection record.
`symmetry` is the [`SymmetryGroup`](@ref) the rule claims (`NoSymmetry()`, `SymmetricGroup{3}()`, …), which
verification checks.

`citations` are the rule's sources: the publications or packages its numbers or its method
come from. `comparisons` are published rules of the same class that it was compared with,
for example by point count, but not made from: further reading, not sources.
"""
Base.@kwdef struct Provenance
    family::String
    derivation::DerivationClass
    path::Vector{String} = String[]
    seed_source::String = "none"
    citations::Vector{Citation} = Citation[]
    comparisons::Vector{Citation} = Citation[]
    license::String = "MIT"
    selection::String = ""
    symmetry::SymmetryGroup = NoSymmetry()
end

"""
    Provenance(p::Provenance; fields...)

A copy of `p` with the fields given as keywords replaced, e.g. `Provenance(p; license = l)`.
"""
function Provenance(p::Provenance; kw...)
    fields = (family = p.family, derivation = p.derivation, path = p.path, seed_source = p.seed_source,
              citations = p.citations, comparisons = p.comparisons, license = p.license,
              selection = p.selection, symmetry = p.symmetry)
    return Provenance(; merge(fields, values(kw))...)
end
with_selection(p::Provenance, s::AbstractString) = Provenance(p; selection = String(s))
with_step(p::Provenance, s::AbstractString; symmetry::SymmetryGroup = p.symmetry) =
    Provenance(p; path = vcat(p.path, String(s)), symmetry)

"""
    Certificate

The measured outcome of refinement: *did Newton converge on the system it
was given?* Not to be confused with [`Verification`](@ref).

- `equations` — which defining equations the residual refers to
- `residual` — the defining-equation residual at the delivered (rounded) rule
- `residual_bits` — the precision, in bits, in which that residual was evaluated
- `digits` — decimal digits delivered
- `guard_digits` — extra digits carried during construction
- `cond` — estimated condition number of the Newton Jacobian (1 if not applicable)
- `iterations` — Newton iterations used (0 for closed-form rules)
- `next_error` — the rule's error at one degree above its claim, measured at working
  precision; recorded (otherwise `nothing`) by families whose rules come so close to the next
  degree that a rule delivered at ordinary precision cannot show the difference
"""
Base.@kwdef struct Certificate
    equations::String
    residual::BigFloat
    residual_bits::Int
    digits::Int
    guard_digits::Int
    cond::Float64 = 1.0
    iterations::Int = 0
    next_error::Union{Nothing,BigFloat} = nothing
end

"""
    combine_certificates(certs, equations, digits) -> Certificate or nothing

The certificate of a rule assembled from parts that each have one (a product rule): the
largest residual and condition number, the fewest guard digits and residual bits, the most
iterations. `nothing` when no part has a certificate.
"""
function combine_certificates(certs, equations::AbstractString, digits::Integer)
    isempty(certs) && return nothing
    return Certificate(; equations = String(equations), residual = maximum(c -> c.residual, certs),
                       residual_bits = minimum(c -> c.residual_bits, certs), digits = Int(digits),
                       guard_digits = minimum(c -> c.guard_digits, certs), cond = maximum(c -> c.cond, certs),
                       iterations = maximum(c -> c.iterations, certs))
end

"""
    Verification

The outcome of checking a rule against a basis of its claimed space: *is this
rule what it claims to be?*

- `basis` — the basis used, e.g. "orthonormal Dubiner"
- `degree` — degree up to which exactness was tested
- `max_residual`, `tolerance` — largest basis residual and the tolerance it was held to
- `exact` — exactness passed
- `sharp` — not exact at `degree + 1` (`nothing` if not tested, or not resolvable; see `sharp_note`)
- `sharp_residual` — the degree-`d+1` residual
- `weights_sum_ok` — the weights sum to the measure of the domain
- `interior` — every node lies inside the domain
- `positive` — no weight is negative (a weight that underflowed to zero is not)
- `symmetric` — the rule is invariant under the symmetry its provenance records (`nothing`
  when it records none)
- `method` — `:exact_integration` or `:convergence_sweep`
- `empirical` — whether the result is empirical rather than certified
- `sharp_note` — why sharpness was left undecided, when it was tested but could not be resolved
- `precision_bits` — arithmetic precision of the check
- `errors`, `decreasing` — for a convergence sweep: the error of each rule, and whether they
  never grew except at the floor (`exact` then also requires the last to meet the tolerance)
"""
Base.@kwdef struct Verification
    basis::String
    degree::Int
    max_residual::BigFloat
    tolerance::BigFloat
    exact::Bool
    sharp::Union{Bool,Nothing}
    sharp_residual::BigFloat
    weights_sum_ok::Bool
    interior::Bool
    positive::Bool
    symmetric::Union{Bool,Nothing}
    method::Symbol = :exact_integration
    empirical::Bool = false
    precision_bits::Int
    sharp_note::String = ""
    errors::Vector{BigFloat} = BigFloat[]
    decreasing::Union{Nothing,Bool} = nothing
end

"""
    passed(v::Verification)

Whether exactness, sharpness (when tested), the weight sum and symmetry (when recorded) all
hold. Positivity and interiority are reported, not required: a rule with negative weights or
nodes on the boundary can be exactly what it claims to be.
"""
passed(v::Verification) = v.exact && v.sharp !== false && v.weights_sum_ok && v.symmetric !== false
