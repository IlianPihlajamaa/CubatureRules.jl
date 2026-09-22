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

"""
    Lebedev()

Octahedrally symmetric, positive-weight rules on `Sphere{3}()`, refined by Gauss–Newton on
the `O_h`-invariant moment system to any requested precision. Seeded: available at the
degrees in the shipped seed table.

The point counts are the smallest the in-house search reached subject to positive weights,
which is why degree 13 is answered by a 78-point rule rather than the published 74-point one
— that rule has a negative weight. Minimal rules are not unique, so a shipped rule need not
coincide node for node with any published table.
"""
struct Lebedev <: RuleFamily end

derivation(::Type{Lebedev}) = Seeded()

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
                                         "(shipped range $(degree_range(Lebedev(), Sphere{3}())))"))
    return e
end

function candidates(::Type{Lebedev}, dom::Sphere{3}, c::PolynomialDegree)
    (isreference(dom) && lebedev_entry_for(c.d) !== nothing) || return Lebedev[]
    return [Lebedev()]
end
candidates(::Type{Lebedev}, dom::Domain, ::PolynomialDegree) = Lebedev[]

npoints(::Lebedev, dom, degree::Integer) = _lebedev_entry(degree).npoints
claimed_degree(::Lebedev, dom, degree) = _lebedev_entry(degree).degree
function degree_range(::Lebedev, dom::Sphere{3})
    es = lebedev_entries()
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::Lebedev, dom) = 1:0
properties(::Lebedev, dom, degree) = (positive = true, interior = true, symmetry = :Oh, nested = false)
cost_estimate(f::Lebedev, dom, degree, T) =
    float(npoints(f, dom, degree)) * length(invariant_exponents(claimed_degree(f, dom, degree))) *
    _precision_factor(T)

# --- construction -------------------------------------------------------------------------

function build(f::Lebedev, dom::Sphere{3}, degree::Int, ctx::BuildContext{T};
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
    θ, res, guard = refine_octahedral(structure, n, θ64, ctx.bits; cancel = ctx.cancel)
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
