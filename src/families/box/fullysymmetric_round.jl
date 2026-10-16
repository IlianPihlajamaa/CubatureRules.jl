# Fully symmetric, positive, interior rules on the disk and the ball (PLAN §6 Tier 3), refined
# to arbitrary precision: `FullySymmetric` on `Disk()` and `Ball{3}()`.
#
# The rules are invariant under the signed permutations of the coordinates, the 8 symmetries of
# the square on the disk and the 48 of the cube (O_h) on the ball, so their orbits are those of
# the box (symmetry/box.jl); the moment equations are in the orthonormal invariant polynomials of
# the disk and the ball (symmetry/round.jl). The seeds in src/data/disk_d4_seeds.toml and
# src/data/ball_oh_seeds.toml were found in-house (scripts/generate_box_seeds.jl): multistart
# at low degree, grow → eliminate above. Product rules, radial Gauss–Jacobi times a rule on the
# circle or the sphere, need several times as many points at the same total degree.

const DISK_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "disk_d4_seeds.toml")
const BALL_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "ball_oh_seeds.toml")
isfile(DISK_SEED_FILE) && include_dependency(DISK_SEED_FILE)
isfile(BALL_SEED_FILE) && include_dependency(BALL_SEED_FILE)

# static package data, parsed once at precompile time
const DISK_SEEDS = load_box_seeds(DISK_SEED_FILE, 2; geometry = RoundGeometry())
const BALL_SEEDS = load_box_seeds(BALL_SEED_FILE, 3; geometry = RoundGeometry())

round_entries(D::Integer) = filter(e -> e.status == "ok", D == 2 ? DISK_SEEDS : D == 3 ? BALL_SEEDS : empty(DISK_SEEDS))
round_entry_for(D::Integer, degree::Integer) = seed_entry_for(round_entries(D), degree)
function _round_entry(D, degree)
    e = round_entry_for(D, degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the $(D)-ball covers degree $degree"))
    return e
end
round_table(D) = D == 2 ? "src/data/disk_d4_seeds.toml" : "src/data/ball_oh_seeds.toml"
round_symmetry(D) = D == 2 ? D4() : Oh()

function candidates(::Type{FullySymmetric}, dom::Ball{D}, c::PolynomialDegree) where {D}
    (D in (2, 3) && isreference(dom) && round_entry_for(D, c.d) !== nothing) || return FullySymmetric[]
    return [FullySymmetric()]
end
npoints(::FullySymmetric, dom::Ball{D}, degree::Integer) where {D} = _round_entry(D, degree).npoints
claimed_degree(::FullySymmetric, dom::Ball{D}, degree) where {D} = _round_entry(D, degree).degree
function degree_range(::FullySymmetric, dom::Ball{D}) where {D}
    es = round_entries(D)
    return isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
properties(::FullySymmetric, dom::Ball{D}, degree) where {D} =
    (positive = true, interior = true, symmetry = round_symmetry(D), nested = false)
cost_estimate(f::FullySymmetric, dom::Ball{D}, degree, T) where {D} =
    float(npoints(f, dom, degree)) * _precision_factor(T)

function build(f::FullySymmetric, dom::Ball{D}, degree::Int, ctx::BuildContext; seed::SeedSource = TableSeed()) where {D}
    e = _round_entry(D, degree)
    return build_seeded(e, ctx; seed, family = "FullySymmetric", domain = Ball{D}(), symmetry = round_symmetry(D),
                        table = round_table(D), origin = ORBIT_TABLE_ORIGIN,
                        equations = "B_$D-invariant moment system on the $(D)-ball (orthonormal invariant " *
                                    "polynomials, degree ≤ $(e.degree))",
                        license = "MIT (seeds and point counts generated in-house)")
end
