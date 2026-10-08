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
const PYRAMID_SEEDS = load_orbit_seeds(PYRAMID_SEED_FILE, PyramidStructure, PyramidStructure)

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
properties(::FullySymmetric, ::Pyramid, degree) = (positive = true, interior = true, symmetry = :C4v, nested = false)
cost_estimate(f::FullySymmetric, dom::Pyramid, degree, T) = float(npoints(f, dom, degree)) * _precision_factor(T)

"""
    refine_pyramid(structure, n, θ64, bits; cancel, verbose) -> (θ, result, guard_bits)

Refine a `Float64` seed of a fully symmetric rule on the pyramid to `bits` bits, by
[`refine_system`](@ref) on the [`PyramidMomentSystem`](@ref).
"""
refine_pyramid(structure::PyramidStructure, n::Integer, θ64::Vector{Float64}, bits::Integer;
               cancel = nothing, verbose::Integer = 0) =
    refine_system(S -> PyramidMomentSystem(structure, n, S), θ64, bits; cancel, verbose,
                  label = "refinement to degree $n on the pyramid")

function build(f::FullySymmetric, dom::Pyramid, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {T}
    isexact(ctx) && throw(ArgumentError("fully symmetric nodes on the pyramid are irrational; $(T) is not supported"))
    seed isa TableSeed || throw(ArgumentError("FullySymmetric on the pyramid takes only the stored table seed"))
    e = _pyramid_entry(degree)
    n, s = e.degree, e.structure
    seed_desc = "stored Float64 table $(PYRAMID_TABLE) (generated in-house: multistart, then grow → eliminate)"
    if ctx.bits <= 53 && isfinite(e.residual)
        # a converged Float64 rule whose residual was checked at high precision when the table
        # was written: at 53 bits or fewer Newton has nothing to add
        θ = e.seed
        steps = ["refine: none — the stored Float64 seed is the rule (residual $(@sprintf("%.1e", e.residual)) " *
                 "at $(e.residual_bits) bits, recorded in the table)"]
        guard, cond, iters = 0, 1.0, 0
    else
        θ, res, guard = refine_pyramid(s, n, e.seed, ctx.bits; cancel = ctx.cancel, verbose = ctx.verbose)
        steps = [@sprintf("refine: Gauss–Newton at %d bits (%d guard), %d iterations, cond %.1e",
                          ctx.bits + guard, guard, res.iterations, res.cond_max)]
        cond, iters = res.cond, res.iterations
    end
    nodes, ws = with_bits(() -> expand(s, θ), max(ctx.bits + guard, 64))
    xs = [SVector{3,T}(ntuple(j -> finalize_number(ctx, x[j]), 3)) for x in nodes]
    wt = [finalize_number(ctx, w) for w in ws]
    # the moment equations at the delivered parameters
    rbits = 2ctx.bits + 32
    resid = with_bits(rbits) do
        θd = [BigFloat(finalize_number(ctx, t)) for t in θ]
        maximum(abs, PyramidMomentSystem(s, n, BigFloat)(θd; jacobian = false)[1])
    end
    cert = Certificate(equations = "C₄ᵥ-invariant moment system on the pyramid (orthonormal invariant " *
                                   "polynomials in the collapsed coordinates, degree ≤ $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)), cond = cond,
                       iterations = iters)
    prov = Provenance(family = "FullySymmetric", derivation = Seeded(),
                      path = vcat(["seed: " * seed_desc, "structure: " * string([o.base.mult for o in s.orbits])], steps),
                      seed_source = seed_desc, symmetry = :C4v,
                      citations = [WITHERDEN_VINCENT_2015],
                      license = "MIT (seeds generated in-house; point counts compared with Witherden & Vincent 2015)")
    return QuadratureRule(xs, wt, Pyramid(), PolynomialDegree(n), prov, cert)
end
