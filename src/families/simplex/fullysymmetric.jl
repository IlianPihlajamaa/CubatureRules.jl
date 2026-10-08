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

Positive-weight, interior-node rules invariant under the whole symmetry group of their
domain, refined by Gauss–Newton on the invariant moment system to any requested precision.
Seeded: available at the degrees in the shipped seed tables.

- On the tetrahedron (`Simplex{3}()`), the 24 permutations of the barycentric coordinates
  (`S₄`). Point counts agree with Witherden & Vincent (2015) for degrees 1–10 and are smaller
  than Zhang, Cui & Liu (2009) at degrees 7, 9 and 11–14.
- On the 4-simplex (`Simplex{4}()`), the 120 permutations of the barycentric coordinates
  (`S₅`).
- On the square, the cube and the 4-cube (`Orthotope{D}()`, `D = 2, 3, 4`), the 8, 48 or 384
  signed permutations of the coordinates. At the same total degree they need far fewer
  points than the tensor Gauss rule, the more so the higher the dimension; the tensor rule
  remains the one for tensor-product (`Q_k`) integrands, which it integrates exactly with
  fewer points.
- On the disk and the ball (`Disk()`, `Ball{3}()`), the same signed permutations, 8 or 48.
- On the pyramid (`Pyramid()`), the 8 symmetries of its square base (`C₄ᵥ`).
- On the wedge (`Wedge()`), the permutations of the triangle's barycentric coordinates and
  the reflection `z ↦ −z`, 12 in all (`D₃ₕ`).

Point counts are the smallest found by the in-house search, not proven minima. Keyword
`seed` of [`rule`](@ref) selects the seed source on the tetrahedron and the 4-simplex, as
for [`XiaoGimbutas`](@ref).
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
const SIMPLEX4_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "simplex4_s5_seeds.toml")
isfile(TETRAHEDRON_SEED_FILE) && include_dependency(TETRAHEDRON_SEED_FILE)
isfile(SIMPLEX4_SEED_FILE) && include_dependency(SIMPLEX4_SEED_FILE)
const TETRAHEDRON_SEEDS = load_symmetric_seeds(TETRAHEDRON_SEED_FILE, 4)
const SIMPLEX4_SEEDS = load_symmetric_seeds(SIMPLEX4_SEED_FILE, 5)

simplex_entries(D::Integer) = filter(e -> e.status == "ok", D == 3 ? TETRAHEDRON_SEEDS : D == 4 ? SIMPLEX4_SEEDS :
                                                          empty(TETRAHEDRON_SEEDS))
tet_entries() = simplex_entries(3)
tet_entry_for(degree::Integer) = seed_entry_for(tet_entries(), degree)
simplex_table(D) = D == 3 ? "src/data/tetrahedron_s4_seeds.toml" : "src/data/simplex4_s5_seeds.toml"

function candidates(::Type{FullySymmetric}, dom::Simplex{D}, c::PolynomialDegree) where {D}
    (D in (3, 4) && isreference(dom) && seed_entry_for(simplex_entries(D), c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end

function _simplex_entry(D, degree)
    e = seed_entry_for(simplex_entries(D), degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the $(D)-simplex covers degree $degree " *
                                         "(shipped range $(degree_range(FullySymmetric(), Simplex{D}())))"))
    return e
end

npoints(::FullySymmetric, ::Simplex{D}, degree::Integer) where {D} = _simplex_entry(D, degree).npoints
claimed_degree(::FullySymmetric, ::Simplex{D}, degree) where {D} = _simplex_entry(D, degree).degree
function degree_range(::FullySymmetric, ::Simplex{D}) where {D}
    es = simplex_entries(D)
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::FullySymmetric, dom) = 1:0
properties(::FullySymmetric, ::Simplex{D}, degree) where {D} =
    (positive = true, interior = true, symmetry = Symbol("S", D + 1), nested = false)
cost_estimate(f::FullySymmetric, dom::Simplex{D}, degree, T) where {D} =
    float(npoints(f, dom, degree)) * simplex_basis_length(D, claimed_degree(f, dom, degree)) / 10 * _precision_factor(T)

function build(f::FullySymmetric, dom::Simplex{D}, degree::Int, ctx::BuildContext;
               seed::SeedSource = TableSeed()) where {D}
    e = _simplex_entry(D, degree)
    return build_symmetric("FullySymmetric", e, ctx; seed, lower = seed_entry_below(simplex_entries(D), e.degree),
                           table = simplex_table(D),
                           citations = D == 3 ? [WITHERDEN_VINCENT_2015, ZHANG_CUI_LIU_2009] : Citation[],
                           license = D == 3 ? "MIT (seeds and point counts generated in-house; the citations are " *
                                              "the published rules of the same class, for comparison, not a source)" :
                                              "MIT (seeds and point counts generated in-house)")
end
