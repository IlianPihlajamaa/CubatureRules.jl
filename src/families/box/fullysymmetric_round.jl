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
const DISK_SEEDS = load_box_seeds(DISK_SEED_FILE, 2)
const BALL_SEEDS = load_box_seeds(BALL_SEED_FILE, 3)

round_entries(D::Integer) = filter(e -> e.status == "ok", D == 2 ? DISK_SEEDS : D == 3 ? BALL_SEEDS : OrbitSeedEntry{BoxStructure}[])
round_entry_for(D::Integer, degree::Integer) = seed_entry_for(round_entries(D), degree)
function _round_entry(D, degree)
    e = round_entry_for(D, degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the $(D)-ball covers degree $degree"))
    return e
end
round_table(D) = D == 2 ? "src/data/disk_d4_seeds.toml" : "src/data/ball_oh_seeds.toml"
round_symmetry(D) = D == 2 ? :D4 : :Oh

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

"""
    refine_round(structure, n, θ64, bits; cancel, verbose) -> (θ, result, guard_bits)

Refine a `Float64` seed of a fully symmetric rule on the disk or the ball to `bits` bits, by
[`refine_system`](@ref) on the [`RoundMomentSystem`](@ref).
"""
refine_round(structure::BoxStructure, n::Integer, θ64::Vector{Float64}, bits::Integer;
             cancel = nothing, verbose::Integer = 0) =
    refine_system(S -> RoundMomentSystem(structure, n, S), θ64, bits; cancel, verbose,
                  label = "refinement to degree $n on the $(structure.D)-ball")

function build(f::FullySymmetric, dom::Ball{D}, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {D,T}
    isexact(ctx) && throw(ArgumentError("fully symmetric nodes on the ball are irrational; $(T) is not supported"))
    seed isa TableSeed || throw(ArgumentError("FullySymmetric on the $(D)-ball takes only the stored table seed"))
    e = _round_entry(D, degree)
    n, s = e.degree, e.structure
    seed_desc = "stored Float64 table $(round_table(D)) (generated in-house: multistart, then grow → eliminate)"
    if ctx.bits <= 53 && isfinite(e.residual)
        # a converged Float64 rule whose residual was checked at high precision when the table
        # was written: at 53 bits or fewer Newton has nothing to add
        θ = e.seed
        steps = ["refine: none — the stored Float64 seed is the rule (residual $(@sprintf("%.1e", e.residual)) " *
                 "at $(e.residual_bits) bits, recorded in the table)"]
        guard, cond, iters = 0, 1.0, 0
    else
        θ, res, guard = refine_round(s, n, e.seed, ctx.bits; cancel = ctx.cancel, verbose = ctx.verbose)
        steps = [@sprintf("refine: Gauss–Newton at %d bits (%d guard), %d iterations, cond %.1e",
                          ctx.bits + guard, guard, res.iterations, res.cond_max)]
        cond, iters = res.cond, res.iterations
    end
    nodes, ws = with_bits(() -> expand(s, θ), max(ctx.bits + guard, 64))
    xs = [SVector{D,T}(ntuple(j -> finalize_number(ctx, x[j]), D)) for x in nodes]
    wt = [finalize_number(ctx, w) for w in ws]
    # the moment equations at the delivered parameters
    rbits = 2ctx.bits + 32
    resid = with_bits(rbits) do
        θd = BigFloat[]
        for (o, off) in zip(s.orbits, param_offsets(s))
            push!(θd, BigFloat(finalize_number(ctx, θ[off + 1])))
            append!(θd, (BigFloat(finalize_number(ctx, θ[off + 1 + i])) for i in 1:nparams(o)))
        end
        maximum(abs, RoundMomentSystem(s, n, BigFloat)(θd; jacobian = false)[1])
    end
    cert = Certificate(equations = "B_$D-invariant moment system on the $(D)-ball (orthonormal invariant " *
                                   "polynomials, degree ≤ $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)), cond = cond,
                       iterations = iters)
    prov = Provenance(family = "FullySymmetric", derivation = Seeded(),
                      path = vcat(["seed: " * seed_desc, "structure: " * string([o.mult for o in s.orbits])], steps),
                      seed_source = seed_desc, symmetry = round_symmetry(D),
                      license = "MIT (seeds and point counts generated in-house)")
    return QuadratureRule(xs, wt, Ball{D}(), PolynomialDegree(n), prov, cert)
end
