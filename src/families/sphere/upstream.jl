# Interop with an installed Lebedev.jl (PLAN §0.1, §0.2), which is GPL-3 licensed.
#
# This family ships no data. It hands back, in this package's types, the rule that the
# caller's own copy of Lebedev.jl produced, and it says so in the provenance: the family
# name, the licence and the citation all point at Lebedev.jl rather than at us.
#
# Three deliberate restrictions keep that honest, and they are enforced in code rather than
# described in a comment:
#
#   * `selectable` is false, so the selector never returns one of these from a request that
#     said only "degree 19". A caller receives a GPL-licensed rule only by naming the family.
#   * only `Float64` is supported. Refining these values to 200 digits would produce numbers
#     derived from a GPL table while looking like ours, which is exactly the laundering this
#     package's seed policy exists to prevent. The refinement machinery is deliberately not
#     wired up here.
#   * nothing from here is ever written to `src/data`. The shipped tables stay in-house.
#
# What this *is* for: using rules you already have the right to use, inside the same
# interface as everything else — `integrate`, `check`, `map_to`, `static` — and verifying
# them against spherical harmonics with machinery that owes nothing to their source.

"""
    UpstreamLebedev()

Lebedev rules taken from an installed [Lebedev.jl](https://github.com/stefabat/Lebedev.jl),
wrapped in this package's types.

Requires `using Lebedev`, and is **never chosen by the selector**: ask for it by name,

```julia
rule(UpstreamLebedev(), Sphere{3}(); degree = 29)
```

The rule that comes back is Lebedev.jl's, not ours. Its provenance records that, including
the GPL-3 licence, which this package cannot relicense — a rule from this family is not
covered by CubatureRules.jl's MIT licence and may not be redistributed under it.

Only `Float64` is offered. These are `Float64` tables, and refining them here would make the
result a derived work of a copyleft-licensed table wearing this package's provenance. For
arbitrary precision either use [`Lebedev`](@ref), whose seeds are generated in-house, or
supply numbers you may license yourself through [`ExplicitSeed`](@ref).
"""
struct UpstreamLebedev <: RuleFamily end

derivation(::Type{UpstreamLebedev}) = Seeded()
family_name(::UpstreamLebedev) = "Lebedev.jl"
describe_family(::UpstreamLebedev) = "UpstreamLebedev (Lebedev.jl, GPL-3)"
# nothing once Lebedev.jl is loaded and the extension has defined `build` for this family,
# the same test the delegated 1D families use
missing_dependency(f::UpstreamLebedev) =
    applicable(build, f, Sphere{3}(), 1, BuildContext{Float64}(53)) ? nothing : "Lebedev.jl"

# Never offered by `rule(domain; degree)` or listed by `available`: a caller must name it,
# because the rule it returns carries terms this package cannot pass on.
selectable(::UpstreamLebedev) = false

properties(::UpstreamLebedev, dom, degree) =
    (positive = true, interior = true, symmetry = :Oh, nested = false)

const LEBEDEV_LAIKOV_1999 = Citation(
    key = "LebedevLaikov1999", authors = ["V. I. Lebedev", "D. N. Laikov"],
    title = "A quadrature formula for the sphere of the 131st algebraic order of accuracy",
    journal = "Doklady Mathematics", year = 1999, volume = "59", pages = "477--481")

"The licence a rule from this family carries, recorded in its provenance."
const UPSTREAM_LEBEDEV_LICENSE =
    "GPL-3.0, via Lebedev.jl. This rule is not covered by CubatureRules.jl's MIT licence " *
    "and may not be redistributed under it. Lebedev & Laikov (1999) ask to be cited."

# Without the extension loaded there is nothing to offer; `rule` reports the missing
# dependency before it gets this far.
candidates(::Type{UpstreamLebedev}, dom::Domain, ::PolynomialDegree) = UpstreamLebedev[]
degree_range(::UpstreamLebedev, dom) = 1:0
