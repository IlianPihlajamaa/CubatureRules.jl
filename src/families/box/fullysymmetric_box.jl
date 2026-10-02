# Fully symmetric, positive, interior rules on the square and the cube (PLAN §6 Tier 3, v0.6),
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

struct BoxSeedEntry
    degree::Int
    npoints::Int
    structure::BoxStructure
    seed::Vector{Float64}
    status::String
    note::String
    residual::Float64            # recorded by scripts/certify_tables.jl; NaN until checked
    residual_bits::Int
end

const SQUARE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "square_d4_seeds.toml")
const CUBE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "cube_oh_seeds.toml")
isfile(SQUARE_SEED_FILE) && include_dependency(SQUARE_SEED_FILE)
isfile(CUBE_SEED_FILE) && include_dependency(CUBE_SEED_FILE)

function load_box_seeds(path, D)
    isfile(path) || return BoxSeedEntry[]
    out = BoxSeedEntry[]
    for e in get(TOML.parsefile(path), "rule", Any[])
        s = BoxStructure([Int[x for x in m] for m in e["structure"]], D)
        push!(out, BoxSeedEntry(e["degree"], e["npoints"], s, Float64.(e["seed"]), get(e, "status", "ok"),
                                get(e, "note", ""), Float64(get(e, "residual", NaN)), Int(get(e, "residual_bits", 0))))
    end
    return sort!(out; by = e -> e.degree)
end

# static package data, parsed once at precompile time
const SQUARE_SEEDS = load_box_seeds(SQUARE_SEED_FILE, 2)
const CUBE_SEEDS = load_box_seeds(CUBE_SEED_FILE, 3)

box_entries(D::Integer) = filter(e -> e.status == "ok", D == 2 ? SQUARE_SEEDS : D == 3 ? CUBE_SEEDS : BoxSeedEntry[])
box_entry_for(D::Integer, degree::Integer) = seed_entry_for(box_entries(D), degree)
function _box_entry(D, degree)
    e = box_entry_for(D, degree)
    e === nothing && throw(ArgumentError("no FullySymmetric seed on the $(D)-cube covers degree $degree"))
    return e
end
box_table(D) = D == 2 ? "src/data/square_d4_seeds.toml" : "src/data/cube_oh_seeds.toml"
box_symmetry(D) = D == 2 ? :D4 : :Oh

function candidates(::Type{FullySymmetric}, dom::Orthotope{D}, c::PolynomialDegree) where {D}
    (D in (2, 3) && isreference(dom) && box_entry_for(D, c.d) !== nothing) || return FullySymmetric[]
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

"""
    refine_box(structure, n, θ64, bits; cancel, verbose) -> (θ, result, guard_bits)

Refine a `Float64` seed of a fully symmetric box rule to `bits` bits, by
[`refine_system`](@ref) on the [`BoxMomentSystem`](@ref).
"""
refine_box(structure::BoxStructure, n::Integer, θ64::Vector{Float64}, bits::Integer;
           cancel = nothing, verbose::Integer = 0) =
    refine_system(S -> BoxMomentSystem(structure, n, S), θ64, bits; cancel, verbose,
                  label = "box refinement to degree $n")

function build(f::FullySymmetric, dom::Orthotope{D}, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {D,T}
    isexact(ctx) && throw(ArgumentError("fully symmetric box nodes are irrational; $(T) is not supported"))
    seed isa TableSeed || throw(ArgumentError("FullySymmetric on the $(D)-cube takes only the stored table seed"))
    e = _box_entry(D, degree)
    n, s = e.degree, e.structure
    seed_desc = "stored Float64 table $(box_table(D)) (generated in-house: multistart, then grow → eliminate)"
    shipped = ctx.bits <= 53 && isfinite(e.residual)
    if shipped
        # a converged Float64 rule whose residual was checked at high precision when the table
        # was written: at 53 bits or fewer Newton has nothing to add
        θ = e.seed
        steps = ["refine: none — the stored Float64 seed is the rule (residual $(@sprintf("%.1e", e.residual)) " *
                 "at $(e.residual_bits) bits, recorded in the table)"]
        guard, cond, iters = 0, 1.0, 0
    else
        θ, res, guard = refine_box(s, n, e.seed, ctx.bits; cancel = ctx.cancel, verbose = ctx.verbose)
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
        maximum(abs, BoxMomentSystem(s, n, BigFloat)(θd; jacobian = false)[1])
    end
    cert = Certificate(equations = "B_$D-invariant moment system (symmetrised even Legendre products, degree ≤ $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)), cond = cond,
                       iterations = iters)
    prov = Provenance(family = "FullySymmetric", derivation = Seeded(),
                      path = vcat(["seed: " * seed_desc, "structure: " * string([o.mult for o in s.orbits])], steps),
                      seed_source = seed_desc, citations = [WITHERDEN_VINCENT_2015], symmetry = box_symmetry(D),
                      license = "MIT (seeds and point counts generated in-house; the citation is the published " *
                                "rules of the same class, for comparison, not a source)")
    return QuadratureRule(xs, wt, Orthotope{D}(), PolynomialDegree(n), prov, cert)
end
