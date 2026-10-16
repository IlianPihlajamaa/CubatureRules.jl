# Fully symmetric, positive, interior rules on the pyramid (PLAN §6 Tier 3), refined to
# arbitrary precision: `FullySymmetric` on `Pyramid()`.
#
# The rules are invariant under the 8 symmetries of the square base (C₄ᵥ); their orbits are
# square orbits in the collapsed coordinates at a free height (symmetry/pyramid_orbits.jl). The
# seeds in src/data/pyramid_c4v_seeds.toml were found in-house (scripts/generate_orbit_seeds.jl):
# multistart at low degree, grow → eliminate above. No numbers were taken from a publication;
# the point counts are compared with Witherden & Vincent (2015), who tabulate the same class.
# The conical product needs (⌈(n + 1)/2⌉)³ points at degree n.

const PYRAMID_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "pyramid_c4v_seeds.toml")
isfile(PYRAMID_SEED_FILE) && include_dependency(PYRAMID_SEED_FILE)

# static package data, parsed once at precompile time
const PYRAMID_SEEDS = load_seeds(PYRAMID_SEED_FILE, PyramidStructure, PyramidStructure ∘ int_lists)

pyramid_entries() = filter(e -> e.status == "ok", PYRAMID_SEEDS)
function _pyramid_entry(degree)
    e = seed_entry_for(pyramid_entries(), degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the pyramid covers degree $degree"))
    return e
end
const PYRAMID_TABLE = "src/data/pyramid_c4v_seeds.toml"

function candidates(::Type{FullySymmetric}, dom::Pyramid, c::PolynomialDegree)
    (isreference(dom) && seed_entry_for(pyramid_entries(), c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end
npoints(::FullySymmetric, ::Pyramid, degree::Integer) = _pyramid_entry(degree).npoints
claimed_degree(::FullySymmetric, ::Pyramid, degree) = _pyramid_entry(degree).degree
function degree_range(::FullySymmetric, ::Pyramid)
    es = pyramid_entries()
    return isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
properties(::FullySymmetric, ::Pyramid, degree) = (positive = true, interior = true, symmetry = C4v(), nested = false)
cost_estimate(f::FullySymmetric, dom::Pyramid, degree, T) = float(npoints(f, dom, degree)) * _precision_factor(T)

function build(f::FullySymmetric, dom::Pyramid, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed())
    e = _pyramid_entry(degree)
    return build_seeded(e, ctx; seed, family = "FullySymmetric", domain = Pyramid(), symmetry = C4v(),
                        table = PYRAMID_TABLE, origin = ORBIT_TABLE_ORIGIN,
                        equations = "C₄ᵥ-invariant moment system on the pyramid (orthonormal invariant " *
                                    "polynomials in the collapsed coordinates, degree ≤ $(e.degree))",
                        citations = [WITHERDEN_VINCENT_2015],
                        license = "MIT (seeds generated in-house; point counts compared with Witherden & Vincent 2015)")
end
