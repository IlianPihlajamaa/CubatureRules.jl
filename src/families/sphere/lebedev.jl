# Lebedev-style O_h-symmetric sphere rules, refined to arbitrary precision
# (PLAN §6 Tier 3, v0.4).
#
# Seeds. Nothing is copied from any published table (PLAN §0.2). The orbit structures and
# Float64 seeds in `src/data/lebedev_seeds.toml` were found in-house by multistart
# Levenberg–Marquardt and Gauss–Newton on the O_h-invariant moment system
# (`scripts/generate_lebedev_seeds.jl`), searching point counts upward until a positive,
# non-degenerate rule appears. The published counts are used only to annotate the log.
#
# The counts found agree with Lebedev's at every degree the search covers but one: at degree
# 13 his 74-point rule has a negative weight, so the search — which requires positive
# weights — reports the 78-point rule instead. That is not a disagreement about the
# mathematics, only about what is being asked for.
#
# Only odd degrees are tabulated. Every O_h orbit is centrally symmetric, so the odd
# harmonics integrate to zero whatever the parameters and a rule of degree 2k is
# automatically of degree 2k+1; an even-degree request is answered by the odd rule above it.

"Supertype of the sources a [`LebedevRule`](@ref) rule can be seeded from."
abstract type LebedevSeeds end

"""
    InHouseSeeds

The seeds shipped in `src/data/lebedev_seeds.toml`, generated in-house and MIT-licensed.
"""
struct InHouseSeeds <: LebedevSeeds end

"""
    LebedevJLSeeds

The caller's own installation of Lebedev.jl, which is GPL-3. Available only once that
package is loaded, never chosen for you, and the licence travels with every rule built from
it. Spelled `UpstreamLebedev()`; see [`LebedevRule`](@ref).
"""
struct LebedevJLSeeds <: LebedevSeeds end

"""
    LebedevRule()

Octahedrally symmetric, positive-weight rules on `Sphere{3}()`, refined by Gauss–Newton on
the `O_h`-invariant moment system to any requested precision. Seeded: available at the
degrees in the shipped seed table.

The point counts are the smallest the in-house search reached subject to positive weights,
which is why degree 13 is answered by a 78-point rule rather than the published 74-point one
— that rule has a negative weight. Minimal rules are not unique, so a shipped rule need not
coincide node for node with any published table.

The type parameter says where the seeds came from. `LebedevRule()` is
`LebedevRule{InHouseSeeds}()`, the table shipped here; `UpstreamLebedev()` is
`LebedevRule{LebedevJLSeeds}()`, the caller's own Lebedev.jl, which reaches degree 125 and
carries that package's GPL-3 licence into every rule built from it.
"""
struct LebedevRule{S<:LebedevSeeds} <: RuleFamily end
LebedevRule() = LebedevRule{InHouseSeeds}()

# `UpstreamLebedev()` is `LebedevRule{LebedevJLSeeds}()`: Lebedev rules from the caller's own
# Lebedev.jl, which is GPL-3. Documented on `Lebedev` above rather than here, because a
# docstring on an alias of a parametric type documents the same binding twice, which
# Documenter reports as a duplicate.
#
# This package ships none of those numbers. A rule built from them records in its provenance
# that it came from Lebedev.jl under GPL-3 and that this package's MIT licence does not cover
# it — including after refinement, which derives from their table and says so. It is never
# selected on its own: `rule(Sphere{3}(); degree = 29)` returns an in-house rule unless
# `copyleft = true` is passed, while `available` lists the variant with its licence in the
# family name.
const UpstreamLebedev = LebedevRule{LebedevJLSeeds}

derivation(::Type{<:LebedevRule}) = Seeded()
family_name(::LebedevRule) = "Lebedev"
describe_family(::LebedevRule{InHouseSeeds}) = "Lebedev"
describe_family(::LebedevRule{LebedevJLSeeds}) = "Lebedev (Lebedev.jl, GPL-3)"
family_license(::LebedevRule{InHouseSeeds}) = ""
family_license(::LebedevRule{LebedevJLSeeds}) = UPSTREAM_LEBEDEV_LICENSE

# Never chosen on its own: a rule carrying terms this package cannot pass on reaches a
# caller only through `copyleft = true` or by naming the family.
selectable(::LebedevRule{LebedevJLSeeds}) = false

# nothing once Lebedev.jl is loaded and its extension has supplied the table hook
missing_dependency(::LebedevRule{LebedevJLSeeds}) =
    applicable(upstream_lebedev_table, 1) ? nothing : "Lebedev.jl"

"""
    upstream_lebedev_table(degree) -> (order, nodes, weights)

The rule of at least `degree` from the caller's Lebedev.jl, with weights already scaled to
sum to the sphere's area. Defined by the package extension; calling it without Lebedev.jl
loaded is an error that `missing_dependency` reports first.
"""
function upstream_lebedev_table end

"The terms a rule seeded from Lebedev.jl carries, recorded in its provenance."
const UPSTREAM_LEBEDEV_LICENSE =
    "GPL-3.0, via Lebedev.jl. This rule is seeded from that package's table and is not " *
    "covered by CubatureRules.jl's MIT licence; it may not be redistributed under it. " *
    "Lebedev & Laikov (1999) ask to be cited."

const LEBEDEV_LAIKOV_1999 = Citation(
    key = "LebedevLaikov1999", authors = ["V. I. Lebedev", "D. N. Laikov"],
    title = "A quadrature formula for the sphere of the 131st algebraic order of accuracy",
    journal = "Doklady Mathematics", year = 1999, volume = "59", pages = "477--481")

const LEBEDEV_1976 = Citation(
    key = "Lebedev1976", authors = ["V. I. Lebedev"],
    title = "Quadratures on a sphere",
    journal = "USSR Computational Mathematics and Mathematical Physics", year = 1976, volume = "16",
    pages = "10--24", doi = "10.1016/0041-5553(76)90100-2")

# --- the seed table ---------------------------------------------------------------------

"""
    OctahedralSeedEntry

One row of the sphere seed table: the degree it serves, its point count, the orbit structure
and the `Float64` parameters found for it.
"""
struct OctahedralSeedEntry
    degree::Int
    npoints::Int
    structure::OctahedralStructure
    seed::Vector{Float64}
    status::String
    note::String
end

const LEBEDEV_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "lebedev_seeds.toml")
isfile(LEBEDEV_SEED_FILE) && include_dependency(LEBEDEV_SEED_FILE)

function load_octahedral_seeds(path)
    isfile(path) || return OctahedralSeedEntry[]
    t = TOML.parsefile(path)
    out = OctahedralSeedEntry[]
    for e in get(t, "rule", Any[])
        s = OctahedralStructure(Symbol.(e["structure"]))
        push!(out, OctahedralSeedEntry(e["degree"], e["npoints"], s, Float64.(e["seed"]),
                                       get(e, "status", "ok"), get(e, "note", "")))
    end
    sort!(out; by = e -> e.degree)
    return out
end

# Static package data, parsed once at precompile time (not a cache of rules: see §7).
const LEBEDEV_SEEDS = load_octahedral_seeds(LEBEDEV_SEED_FILE)

lebedev_entries() = filter(e -> e.status == "ok", LEBEDEV_SEEDS)

"The entry serving a request for `degree`: fewest points, then highest degree."
function lebedev_entry_for(degree::Integer)
    es = filter(e -> e.degree >= degree, lebedev_entries())
    isempty(es) && return nothing
    return first(sort(es; by = e -> (e.npoints, -e.degree)))
end

function _lebedev_entry(degree)
    e = lebedev_entry_for(degree)
    e === nothing && throw(ArgumentError("no Lebedev seed covers degree $degree " *
                                         "(shipped range $(degree_range(LebedevRule(), Sphere{3}())))"))
    return e
end

function _inhouse_candidates(dom::Sphere{3}, c::PolynomialDegree)
    (isreference(dom) && lebedev_entry_for(c.d) !== nothing) || return LebedevRule[]
    return [LebedevRule()]
end
_inhouse_candidates(dom::Domain, ::PolynomialDegree) = LebedevRule[]

# `families()` yields the type without its parameter, so one method answers for both
# variants: the in-house seeds always, the upstream ones only when the extension is loaded.
function candidates(::Type{<:LebedevRule}, dom::Domain, c::PolynomialDegree)
    out = LebedevRule[]
    append!(out, _inhouse_candidates(dom, c))
    append!(out, upstream_lebedev_candidates(dom, c))
    return out
end

"Candidates backed by the caller's Lebedev.jl; empty until its extension is loaded."
upstream_lebedev_candidates(dom, c) = LebedevRule[]

npoints(::LebedevRule{InHouseSeeds}, dom, degree::Integer) = _lebedev_entry(degree).npoints
claimed_degree(::LebedevRule{InHouseSeeds}, dom, degree) = _lebedev_entry(degree).degree
function degree_range(::LebedevRule{InHouseSeeds}, dom::Sphere{3})
    es = lebedev_entries()
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::LebedevRule{InHouseSeeds}, dom) = 1:0
properties(::LebedevRule, dom, degree) = (positive = true, interior = true, symmetry = :Oh, nested = false)
cost_estimate(f::LebedevRule{InHouseSeeds}, dom, degree, T) =
    float(npoints(f, dom, degree)) * length(invariant_exponents(claimed_degree(f, dom, degree))) *
    _precision_factor(T)

# --- construction -------------------------------------------------------------------------

function build(f::LebedevRule{InHouseSeeds}, dom::Sphere{3}, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {T}
    isexact(ctx) && throw(ArgumentError("Lebedev nodes are irrational; $(T) is not supported"))
    e = _lebedev_entry(degree)
    n = e.degree
    structure, θ64 = e.structure, e.seed
    seed_desc = "stored Float64 table src/data/lebedev_seeds.toml " *
                "(generated in-house by multistart from the orbit structure)"
    if seed isa MultistartSeed
        found = multistart_octahedral(structure, n; nstarts = seed.nstarts, rng_seed = seed.rng_seed,
                                      first_only = false, cancel = ctx.cancel)
        isempty(found) && throw(RefinementError("Lebedev",
            "multistart found no valid degree-$n rule with structure $structure in $(seed.nstarts) starts"))
        best = argmax(t -> (t[2], t[3]), found)
        θ64 = best[1]
        seed_desc = describe(seed) * "; $(length(found)) valid rule(s) found, largest minimum weight kept"
    elseif seed isa ExplicitSeed
        length(seed.θ) == nunknowns(structure) || throw(ArgumentError(
            "an explicit seed for degree $n needs $(nunknowns(structure)) parameters (structure $structure)"))
        θ64, seed_desc = seed.θ, describe(seed)
    end
    θ, res, guard = refine_octahedral(structure, n, θ64, ctx.bits; cancel = ctx.cancel, verbose = ctx.verbose)
    res.converged || throw(RefinementError("Lebedev",
        "Gauss–Newton did not converge at degree $n (residual $(Float64(res.residual)) after " *
        "$(res.iterations) iterations); the seed is not returned unrefined"))
    wbits = ctx.bits + guard
    xs, ws = with_bits(() -> expand(structure, θ), wbits)
    wmin, dmin = octahedral_margins(structure, θ)
    (wmin > 0 && dmin > 0) || throw(RefinementError("Lebedev",
        "refined degree-$n rule left the admissible set (min weight $(Float64(wmin)), margin $(Float64(dmin)))"))
    nodes = [SVector{3,T}(ntuple(j -> finalize_number(ctx, x[j]), 3)) for x in xs]
    wt = [finalize_number(ctx, w) for w in ws]
    θhat = [finalize_number(ctx, t) for t in θ]
    rbits = 2ctx.bits + guard
    resid = with_bits(rbits) do
        sys = OctahedralMomentSystem(structure, n, BigFloat)
        r, _ = sys(BigFloat.(θhat))
        maximum(abs, r)
    end
    cert = Certificate(equations = "O_h-invariant moment system in orbit parameters " *
                                   "(invariants p₄^a p₆^b, 4a + 6b ≤ $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = res.cond, iterations = res.iterations)
    prov = Provenance(family = "Lebedev", derivation = Seeded(),
                      path = ["seed: " * seed_desc,
                              "structure: " * string(structure),
                              "refine: Gauss–Newton with rank-revealing solve at $(wbits) bits",
                              @sprintf("guard: %d bits from measured cond(J) = %.2e", guard, res.cond_max)],
                      seed_source = seed_desc, citations = [LEBEDEV_1976], symmetry = :Oh,
                      license = "MIT (seeds generated in-house; no published table was used)")
    return QuadratureRule(nodes, wt, Sphere{3}(), PolynomialDegree(n), prov, cert)
end
