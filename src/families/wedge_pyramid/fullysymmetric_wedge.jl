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
const WEDGE_SEEDS = load_orbit_seeds(WEDGE_SEED_FILE, WedgeStructure, WedgeStructure)

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
properties(::FullySymmetric, ::Wedge, degree) = (positive = true, interior = true, symmetry = :D3h, nested = false)
cost_estimate(f::FullySymmetric, dom::Wedge, degree, T) = float(npoints(f, dom, degree)) * _precision_factor(T)

"""
    refine_wedge(structure, n, θ64, bits; cancel, verbose) -> (θ, result, guard_bits)

Refine a `Float64` seed of a fully symmetric rule on the wedge to `bits` bits, by
[`refine_system`](@ref) on the [`WedgeMomentSystem`](@ref).
"""
refine_wedge(structure::WedgeStructure, n::Integer, θ64::Vector{Float64}, bits::Integer;
             cancel = nothing, verbose::Integer = 0) =
    refine_system(S -> WedgeMomentSystem(structure, n, S), θ64, bits; cancel, verbose,
                  label = "refinement to degree $n on the wedge")

function build(f::FullySymmetric, dom::Wedge, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {T}
    isexact(ctx) && throw(ArgumentError("fully symmetric nodes on the wedge are irrational; $(T) is not supported"))
    seed isa TableSeed || throw(ArgumentError("FullySymmetric on the wedge takes only the stored table seed"))
    e = _wedge_entry(degree)
    n, s = e.degree, e.structure
    seed_desc = "stored Float64 table $(WEDGE_TABLE) (generated in-house: multistart, then grow → eliminate)"
    if ctx.bits <= 53 && isfinite(e.residual)
        # a converged Float64 rule whose residual was checked at high precision when the table
        # was written: at 53 bits or fewer Newton has nothing to add
        θ = e.seed
        steps = ["refine: none — the stored Float64 seed is the rule (residual $(@sprintf("%.1e", e.residual)) " *
                 "at $(e.residual_bits) bits, recorded in the table)"]
        guard, cond, iters = 0, 1.0, 0
    else
        θ, res, guard = refine_wedge(s, n, e.seed, ctx.bits; cancel = ctx.cancel, verbose = ctx.verbose)
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
        maximum(abs, WedgeMomentSystem(s, n, BigFloat)(θd; jacobian = false)[1])
    end
    cert = Certificate(equations = "D₃ₕ-invariant moment system on the wedge (triangle invariants times even " *
                                   "Legendre polynomials in z, degree ≤ $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)), cond = cond,
                       iterations = iters)
    prov = Provenance(family = "FullySymmetric", derivation = Seeded(),
                      path = vcat(["seed: " * seed_desc, "structure: " * string(wedge_structure_lists(s))], steps),
                      seed_source = seed_desc, symmetry = :D3h,
                      citations = [WITHERDEN_VINCENT_2015],
                      license = "MIT (seeds generated in-house; point counts compared with Witherden & Vincent 2015)")
    return QuadratureRule(xs, wt, Wedge(), PolynomialDegree(n), prov, cert)
end
