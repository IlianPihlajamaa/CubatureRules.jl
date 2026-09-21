# The rule-family interface (PLAN §2.5).
#
# A family is one type `F <: RuleFamily` plus methods of the functions below. Nothing is
# registered anywhere: families are discovered with `subtypes` at call time.

"""
    RuleFamily

Supertype of every rule family. To add a family, define `struct MyFamily <: RuleFamily`
and methods of:

- [`candidates`](@ref)`(::Type{MyFamily}, domain, claim)` — instances applicable on the domain
- [`build`](@ref)`(f::MyFamily, domain, degree, ctx)` — construct the rule on a reference domain
- [`npoints`](@ref)`(f, domain, degree)`, [`properties`](@ref), [`degree_range`](@ref),
  [`derivation`](@ref), and optionally [`cost_estimate`](@ref), [`supports_type`](@ref)

The family is then found by [`rule`](@ref) and [`available`](@ref) with no registration step.
"""
abstract type RuleFamily end

"""
    CombinatorFamily <: RuleFamily

A family built from other families (conical, tensor, Smolyak products). The registry
expands combinators over the leaf families to nesting depth 1.
"""
abstract type CombinatorFamily <: RuleFamily end

"""
    candidates(::Type{F}, domain, claim) -> Vector

The instances of family `F` that can satisfy `claim` on `domain`. An empty vector means
"not applicable" — there is no separate `supports` predicate to disagree with it.
"""
candidates(::Type{<:RuleFamily}, ::Domain, ::ExactnessClaim) = RuleFamily[]

"""
    build(f::RuleFamily, domain, degree::Int, ctx::BuildContext) -> QuadratureRule

Construct the rule of family instance `f` with at least polynomial degree `degree` on the
*reference* domain `domain`, in the number type and precision of `ctx`.
"""
function build end

"""
    npoints(f, domain, degree)

Number of nodes the rule of `f` at `degree` would have, without constructing it.
"""
function npoints end

"""
    properties(f, domain, degree) -> NamedTuple

`(positive, interior, symmetry, nested)` for the rule of `f` at `degree`, without
constructing it.
"""
function properties end

"""
    degree_range(f, domain) -> UnitRange{Int}

Degrees for which `f` can construct a rule on `domain`. Unbounded ranges end at
`typemax(Int)`.
"""
function degree_range end

"""
    derivation(f) -> Derived() | Seeded()

Whether the family is constructed from first principles at any order (`Derived`) or
refined from starting points available only at tabulated orders (`Seeded`).
"""
function derivation end
derivation(f::RuleFamily) = derivation(typeof(f))

"""
    claimed_degree(f, domain, degree)

The degree actually claimed by the rule built for a request of `degree` (may exceed it,
e.g. Grundmann–Möller only exists at odd degree).
"""
claimed_degree(f::RuleFamily, dom, degree) = degree

"""
    cost_estimate(f, domain, degree, T)

A rough relative cost of constructing the rule, in arbitrary units comparable across
families. Used for planning, never for correctness.
"""
cost_estimate(f::RuleFamily, dom, degree, T) = float(npoints(f, dom, degree)) * _precision_factor(T)
_precision_factor(::Type{T}) where {T} = T === BigFloat ? precision(BigFloat) / 53 * 20 : 1.0

"""
    supports_type(f, T)

Whether `f` can deliver nodes and weights in type `T`. Families with irrational nodes do
not support `Rational`.
"""
supports_type(f::RuleFamily, ::Type{T}) where {T} = T <: AbstractFloat

"Human-readable family name used in provenance, ranking and display."
family_name(f::RuleFamily) = string(nameof(typeof(f)))

"Instance description used in tables, e.g. `ConicalProduct(GaussJacobi)`."
describe_family(f::RuleFamily) = family_name(f)

"`true` for families that can act as the 1D factor of a conical product."
conical_compatible(::Type{<:RuleFamily}) = false

"""
    missing_dependency(f) -> nothing or a package name

The package a family needs before it can build anything, when that package is not loaded.
The selector uses it to say so instead of reporting the family as nonexistent. Families
whose rules come from an upstream package define this.
"""
missing_dependency(::RuleFamily) = nothing
