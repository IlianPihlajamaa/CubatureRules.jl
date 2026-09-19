# Fully symmetric, positive, interior tetrahedron rules refined to arbitrary precision
# (PLAN §6 Tier 3, v0.2).
#
# Unlike the triangle table, nothing here comes from a paper: the point counts are the
# smallest the in-house search reached (`scripts/generate_tetrahedron_seeds.jl`), walking
# upwards from a count below which no fully symmetric rule can exist. They are not proven
# minimal and are not attributed to any published family. The seeds are MIT-licensed.

"""
    FullySymmetric()

Fully symmetric (`S₄`), positive-weight, interior-node rules on the tetrahedron, refined by
Gauss–Newton on the `S₄`-invariant moment system to any requested precision. Seeded:
available at the degrees in the shipped seed table. Point counts are the smallest found
by the in-house orbit-structure search, not proven minimal.

Keyword `seed` of [`rule`](@ref) selects the seed source, as for [`XiaoGimbutas`](@ref).
"""
struct FullySymmetric <: RuleFamily end

derivation(::Type{FullySymmetric}) = Seeded()

const TETRAHEDRON_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "tetrahedron_s4_seeds.toml")
isfile(TETRAHEDRON_SEED_FILE) && include_dependency(TETRAHEDRON_SEED_FILE)
const TETRAHEDRON_SEEDS = load_symmetric_seeds(TETRAHEDRON_SEED_FILE, 4)

tet_entries() = filter(e -> e.status == "ok", TETRAHEDRON_SEEDS)
tet_entry_for(degree::Integer) = seed_entry_for(tet_entries(), degree)

function candidates(::Type{FullySymmetric}, dom::Simplex{3}, c::PolynomialDegree)
    (isreference(dom) && tet_entry_for(c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end

function _tet_entry(degree)
    e = tet_entry_for(degree)
    e === nothing && throw(ArgumentError("no FullySymmetric tetrahedron seed covers degree $degree " *
                                         "(shipped range $(degree_range(FullySymmetric(), Simplex{3}())))"))
    return e
end

npoints(::FullySymmetric, dom, degree::Integer) = _tet_entry(degree).npoints
claimed_degree(::FullySymmetric, dom, degree) = _tet_entry(degree).degree
function degree_range(::FullySymmetric, dom::Simplex{3})
    es = tet_entries()
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::FullySymmetric, dom) = 1:0
properties(::FullySymmetric, dom, degree) = (positive = true, interior = true, symmetry = :S4, nested = false)
cost_estimate(f::FullySymmetric, dom, degree, T) =
    float(npoints(f, dom, degree)) * tet_length(claimed_degree(f, dom, degree)) / 10 * _precision_factor(T)

function build(f::FullySymmetric, dom::Simplex{3}, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed())
    return build_symmetric("FullySymmetric", _tet_entry(degree), ctx; seed,
                           table = "src/data/tetrahedron_s4_seeds.toml", citations = Citation[],
                           license = "MIT (seeds and point counts generated in-house)")
end
