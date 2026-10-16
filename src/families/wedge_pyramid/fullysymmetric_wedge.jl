# Fully symmetric, positive, interior rules on the wedge (PLAN §6 Tier 3), refined to arbitrary
# precision: `FullySymmetric` on `Wedge()`.
#
# The rules are invariant under the 12 symmetries of the prism (D₃ₕ): the permutations of the
# triangle's barycentric coordinates and z ↦ −z. Their orbits are triangle orbits at z = 0 or
# at a pair ±z (symmetry/wedge_orbits.jl). The seeds in src/data/wedge_d3h_seeds.toml were
# found in-house (scripts/generate_orbit_seeds.jl): multistart at low degree, grow → eliminate
# above. No numbers were taken from a publication; the point counts are compared with
# Witherden & Vincent (2015), who tabulate the same class. The product of a triangle rule and
# Gauss–Legendre needs several times as many points at the same total degree.

const WEDGE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "wedge_d3h_seeds.toml")
isfile(WEDGE_SEED_FILE) && include_dependency(WEDGE_SEED_FILE)

# static package data, parsed once at precompile time
const WEDGE_SEEDS = load_seeds(WEDGE_SEED_FILE, WedgeStructure, WedgeStructure ∘ int_lists)

wedge_entries() = filter(e -> e.status == "ok", WEDGE_SEEDS)
function _wedge_entry(degree)
    e = seed_entry_for(wedge_entries(), degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the wedge covers degree $degree"))
    return e
end
const WEDGE_TABLE = "src/data/wedge_d3h_seeds.toml"

function candidates(::Type{FullySymmetric}, dom::Wedge, c::PolynomialDegree)
    (isreference(dom) && seed_entry_for(wedge_entries(), c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end
npoints(::FullySymmetric, ::Wedge, degree::Integer) = _wedge_entry(degree).npoints
claimed_degree(::FullySymmetric, ::Wedge, degree) = _wedge_entry(degree).degree
function degree_range(::FullySymmetric, ::Wedge)
    es = wedge_entries()
    return isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
properties(::FullySymmetric, ::Wedge, degree) = (positive = true, interior = true, symmetry = D3h(), nested = false)
cost_estimate(f::FullySymmetric, dom::Wedge, degree, T) = float(npoints(f, dom, degree)) * _precision_factor(T)

function build(f::FullySymmetric, dom::Wedge, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed())
    e = _wedge_entry(degree)
    return build_seeded(e, ctx; seed, family = "FullySymmetric", domain = Wedge(), symmetry = D3h(),
                        table = WEDGE_TABLE, origin = ORBIT_TABLE_ORIGIN,
                        equations = "D₃ₕ-invariant moment system on the wedge (triangle invariants times even " *
                                    "Legendre polynomials in z, degree ≤ $(e.degree))",
                        citations = [WITHERDEN_VINCENT_2015],
                        license = "MIT (seeds generated in-house; point counts compared with Witherden & Vincent 2015)")
end
