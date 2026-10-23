# Fully symmetric, positive, interior rules on the square and the cube,
# refined to arbitrary precision: `FullySymmetric` on `Orthotope{2}()` and `Orthotope{3}()`.
#
# The seeds in src/data/square_d4_seeds.toml and src/data/cube_oh_seeds.toml were found
# in-house (scripts/generate_box_seeds.jl): multistart over the orbit structures of each point
# count at low degree, grow → eliminate above. No numbers were taken from a publication; the
# point counts are compared with Witherden & Vincent (2015), who tabulate the same class.
#
# At the same total degree they need far fewer points than the tensor Gauss rule, most of all
# on the cube. A tensor rule integrates the larger space of tensor-product polynomials, which
# is what the mass and stiffness matrices of tensor-product (Q_k) elements need; these are for
# total-degree integrands.

const SQUARE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "square_d4_seeds.toml")
const CUBE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "cube_oh_seeds.toml")
const HYPERCUBE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "hypercube_b4_seeds.toml")
isfile(SQUARE_SEED_FILE) && include_dependency(SQUARE_SEED_FILE)
isfile(CUBE_SEED_FILE) && include_dependency(CUBE_SEED_FILE)
isfile(HYPERCUBE_SEED_FILE) && include_dependency(HYPERCUBE_SEED_FILE)

"The seed table at `path` of fully symmetric rules on the `D`-cube, or on the `D`-ball with `geometry = RoundGeometry()`."
load_box_seeds(path, D; geometry::OrbitGeometry = BoxGeometry()) =
    load_seeds(path, BoxStructure{typeof(geometry)}, m -> BoxStructure(int_lists(m), D; geometry))

# static package data, parsed once at precompile time
const SQUARE_SEEDS = load_box_seeds(SQUARE_SEED_FILE, 2)
const CUBE_SEEDS = load_box_seeds(CUBE_SEED_FILE, 3)
const HYPERCUBE_SEEDS = load_box_seeds(HYPERCUBE_SEED_FILE, 4)

box_entries(D::Integer) = filter(e -> e.status == "ok", D == 2 ? SQUARE_SEEDS : D == 3 ? CUBE_SEEDS :
                                                      D == 4 ? HYPERCUBE_SEEDS : empty(SQUARE_SEEDS))
box_entry_for(D::Integer, degree::Integer) = seed_entry_for(box_entries(D), degree)
function _box_entry(D, degree)
    e = box_entry_for(D, degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the $(D)-cube covers degree $degree"))
    return e
end
box_table(D) = D == 2 ? "src/data/square_d4_seeds.toml" : D == 3 ? "src/data/cube_oh_seeds.toml" :
               "src/data/hypercube_b4_seeds.toml"
box_symmetry(D) = D == 2 ? D4() : D == 3 ? Oh() : B4()

function candidates(::Type{FullySymmetric}, dom::Orthotope{D}, c::PolynomialDegree) where {D}
    (D in (2, 3, 4) && isreference(dom) && box_entry_for(D, c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end
npoints(::FullySymmetric, dom::Orthotope{D}, degree::Integer) where {D} = _box_entry(D, degree).npoints
claimed_degree(::FullySymmetric, dom::Orthotope{D}, degree) where {D} = _box_entry(D, degree).degree
function degree_range(::FullySymmetric, dom::Orthotope{D}) where {D}
    es = box_entries(D)
    return isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
properties(::FullySymmetric, dom::Orthotope{D}, degree) where {D} =
    (positive = true, interior = true, symmetry = box_symmetry(D), nested = false)
cost_estimate(f::FullySymmetric, dom::Orthotope{D}, degree, T) where {D} =
    float(npoints(f, dom, degree)) * _precision_factor(T)

"Where the seed tables of the box, the ball, the pyramid and the wedge come from."
const ORBIT_TABLE_ORIGIN = "generated in-house: multistart, then grow → eliminate"

function build(f::FullySymmetric, dom::Orthotope{D}, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed()) where {D}
    e = _box_entry(D, degree)
    return build_seeded(e, ctx; seed, family = "FullySymmetric", domain = Orthotope{D}(), symmetry = box_symmetry(D),
                        table = box_table(D), origin = ORBIT_TABLE_ORIGIN,
                        equations = "B_$D-invariant moment system (symmetrised even Legendre products, degree ≤ $(e.degree))",
                        citations = [WITHERDEN_VINCENT_2015],
                        license = "MIT (seeds and point counts generated in-house; the citation is the published " *
                                  "rules of the same class, for comparison, not a source)")
end
