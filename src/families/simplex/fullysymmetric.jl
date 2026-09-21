# Fully symmetric, positive, interior tetrahedron rules refined to arbitrary precision
# (PLAN §6 Tier 3, v0.2).
#
# Nothing here is taken from a paper: the point counts are the smallest the in-house search
# reached (`scripts/generate_tetrahedron_seeds.jl`), walking upwards from a count below
# which no fully symmetric rule can exist, then node elimination. The seeds are
# MIT-licensed.
#
# They have since been compared with the published counts for the same class (fully
# symmetric, positive weights, interior nodes). Degrees 1–10 agree exactly with Witherden &
# Vincent (2015), Table 1; degrees 7, 9 and 11–14 are smaller than Zhang, Cui & Liu (2009),
# Table 4.2 (35 vs 36, 59 vs 61, 102 vs 109, 124 vs 140, 145 vs 171, 179 vs 236). Xiao &
# Gimbutas and Jaśkowiec & Sukumar report fewer points still, but their tetrahedral rules
# are not fully symmetric, so they are not the same class.

"""
    FullySymmetric()

Fully symmetric (`S₄`), positive-weight, interior-node rules on the tetrahedron, refined by
Gauss–Newton on the `S₄`-invariant moment system to any requested precision. Seeded:
available at the degrees in the shipped seed table.

Point counts are the smallest found by the in-house search, not proven minima. They agree
with Witherden & Vincent (2015) for degrees 1–10 and are smaller than Zhang, Cui & Liu
(2009) at degrees 7, 9 and 11–14.

Keyword `seed` of [`rule`](@ref) selects the seed source, as for [`XiaoGimbutas`](@ref).
"""
struct FullySymmetric <: RuleFamily end

derivation(::Type{FullySymmetric}) = Seeded()

const WITHERDEN_VINCENT_2015 = Citation(
    key = "WitherdenVincent2015", authors = ["Freddie D. Witherden", "Peter E. Vincent"],
    title = "On the identification of symmetric quadrature rules for finite element methods",
    journal = "Computers & Mathematics with Applications", year = 2015, volume = "69", pages = "1232--1241",
    doi = "10.1016/j.camwa.2015.03.017")
const ZHANG_CUI_LIU_2009 = Citation(
    key = "ZhangCuiLiu2009", authors = ["Linbo Zhang", "Tao Cui", "Hui Liu"],
    title = "A set of symmetric quadrature rules on triangles and tetrahedra",
    journal = "Journal of Computational Mathematics", year = 2009, volume = "27", pages = "89--96")

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
    e = _tet_entry(degree)
    return build_symmetric("FullySymmetric", e, ctx; seed, lower = seed_entry_below(tet_entries(), e.degree),
                           table = "src/data/tetrahedron_s4_seeds.toml",
                           citations = [WITHERDEN_VINCENT_2015, ZHANG_CUI_LIU_2009],
                           license = "MIT (seeds and point counts generated in-house; the citations are the " *
                                     "published rules of the same class, for comparison, not a source)")
end
