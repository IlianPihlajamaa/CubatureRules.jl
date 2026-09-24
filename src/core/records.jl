# Metadata carried by a constructed rule: where it came from (Provenance), whether the
# construction converged (Certificate, §1) and whether it is what it claims (Verification,
# §8). Certification and verification are different checks with different types.

"""
    Citation

A bibliographic reference. `cite(rule)` renders these.
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
end

"""
    Derived()
    Seeded()

The two derivation classes (PLAN §1). A *derived* rule is constructed from first
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
`symmetry` names the symmetry group the rule claims (`:none`, `:S3`, …), which
verification checks.
"""
Base.@kwdef struct Provenance
    family::String
    derivation::DerivationClass
    path::Vector{String} = String[]
    seed_source::String = "none"
    citations::Vector{Citation} = Citation[]
    license::String = "MIT"
    selection::String = ""
    symmetry::Symbol = :none
end

with_selection(p::Provenance, s::AbstractString) =
    Provenance(p.family, p.derivation, p.path, p.seed_source, p.citations, p.license, String(s), p.symmetry)
with_step(p::Provenance, s::AbstractString; symmetry::Symbol = p.symmetry) =
    Provenance(p.family, p.derivation, vcat(p.path, String(s)), p.seed_source, p.citations, p.license, p.selection, symmetry)

"""
    Certificate

The measured outcome of refinement (PLAN §1, §2.7): *did Newton converge on the system it
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
    Verification

The outcome of checking a rule against a basis of its claimed space (PLAN §8): *is this
rule what it claims to be?*

- `basis` — the basis used, e.g. "orthonormal Dubiner"
- `degree` — degree up to which exactness was tested
- `max_residual`, `tolerance` — largest basis residual and the tolerance it was held to
- `exact` — exactness passed
- `sharp` — not exact at `degree + 1` (`nothing` if not tested, or not resolvable; see `sharp_note`)
- `sharp_residual` — the degree-`d+1` residual
- `weights_sum_ok`, `interior`, `positive`, `symmetric` — structural invariants
  (`nothing` where the invariant is not claimed / not applicable)
- `method` — `:exact_integration` or `:convergence_sweep`
- `empirical` — whether the result is empirical rather than certified
- `sharp_note` — why sharpness was left undecided, when it was tested but could not be resolved
- `precision_bits` — arithmetic precision of the check
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
end

"""
    passed(v::Verification)

Whether exactness, sharpness (when tested), the weight sum and symmetry (when claimed) all
hold. Positivity and interiority are reported but only fail a check when claimed.
"""
passed(v::Verification) = v.exact && v.sharp !== false && v.weights_sum_ok && v.symmetric !== false
