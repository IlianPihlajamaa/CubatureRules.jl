# Minimal-point fully symmetric triangle rules, refined to arbitrary precision
# (PLAN §6 Tier 3 — the v0.1 cornerstone).
#
# Seeds. No numbers are copied from any published table (PLAN §0.2). The point counts and
# the fact that the rules are fully S_3-symmetric come from Xiao & Gimbutas (2010); the orbit
# structures and Float64 seeds in `src/data/triangle_s3_seeds.toml` were found in-house by
# multistart Gauss–Newton from the orbit structure alone (`scripts/generate_triangle_seeds.jl`)
# and are MIT-licensed. Where several rules share the point count, the one with the largest
# minimum barycentric coordinate is shipped, so a shipped rule need not coincide node for
# node with the published table.

"""
    XiaoGimbutas()

Fully symmetric, positive-weight, interior-node triangle rules with the minimal point counts
of Xiao & Gimbutas (2010), refined by Gauss–Newton on the `S_3`-invariant moment system to
any requested precision. Seeded: available at the degrees in the shipped seed table.

Keyword `seed` of [`rule`](@ref) selects the seed source: `TableSeed()` (default) or
`MultistartSeed()` to recover the rule from its orbit structure alone.
"""
struct XiaoGimbutas <: RuleFamily end

derivation(::Type{XiaoGimbutas}) = Seeded()

const XIAO_GIMBUTAS_2010 = Citation(
    key = "XiaoGimbutas2010", authors = ["Hong Xiao", "Zydrunas Gimbutas"],
    title = "A numerical algorithm for the construction of efficient quadrature rules in two and higher dimensions",
    journal = "Computers & Mathematics with Applications", year = 2010, volume = "59", pages = "663--676",
    doi = "10.1016/j.camwa.2009.10.027")

# --- the seed table ---------------------------------------------------------------------

struct SymmetricSeedEntry
    degree::Int
    npoints::Int
    structure::SymmetricStructure
    seed::Vector{Float64}
    status::String
    note::String
end

const TRIANGLE_SEED_FILE = joinpath(@__DIR__, "..", "..", "data", "triangle_s3_seeds.toml")
isfile(TRIANGLE_SEED_FILE) && include_dependency(TRIANGLE_SEED_FILE)

function load_symmetric_seeds(path, N)
    isfile(path) || return SymmetricSeedEntry[]
    t = TOML.parsefile(path)
    out = SymmetricSeedEntry[]
    for e in get(t, "rule", Any[])
        s = SymmetricStructure([Int.(m) for m in e["structure"]], N)
        push!(out, SymmetricSeedEntry(e["degree"], e["npoints"], s, Float64.(e["seed"]),
                                      get(e, "status", "ok"), get(e, "note", "")))
    end
    sort!(out; by = e -> e.degree)
    return out
end

# Static package data, parsed once at precompile time (not a cache of rules: see §7).
const TRIANGLE_SEEDS = load_symmetric_seeds(TRIANGLE_SEED_FILE, 3)

xg_entries() = filter(e -> e.status == "ok", TRIANGLE_SEEDS)

"The shipped entry serving a request for `degree`: fewest points, then highest degree."
function xg_entry_for(degree::Integer)
    es = filter(e -> e.degree >= degree, xg_entries())
    isempty(es) && return nothing
    return first(sort(es; by = e -> (e.npoints, -e.degree)))
end

function candidates(::Type{XiaoGimbutas}, dom::Simplex{2}, c::PolynomialDegree)
    (isreference(dom) && xg_entry_for(c.d) !== nothing) || return XiaoGimbutas[]
    return [XiaoGimbutas()]
end

npoints(::XiaoGimbutas, dom, degree::Integer) = _xg_entry(degree).npoints
claimed_degree(::XiaoGimbutas, dom, degree) = _xg_entry(degree).degree
function degree_range(::XiaoGimbutas, dom::Simplex{2})
    es = xg_entries()
    isempty(es) ? (1:0) : (0:maximum(e -> e.degree, es))
end
degree_range(::XiaoGimbutas, dom) = 1:0
properties(::XiaoGimbutas, dom, degree) = (positive = true, interior = true, symmetry = :S3, nested = false)
cost_estimate(f::XiaoGimbutas, dom, degree, T) =
    float(npoints(f, dom, degree)) * dubiner_length(claimed_degree(f, dom, degree)) / 10 * _precision_factor(T)

function _xg_entry(degree)
    e = xg_entry_for(degree)
    e === nothing && throw(ArgumentError("no XiaoGimbutas seed covers degree $degree " *
                                         "(shipped range $(degree_range(XiaoGimbutas(), Simplex{2}())))"))
    return e
end

# --- refinement ------------------------------------------------------------------------

"""
    refine_symmetric_triangle(structure, n, θ64, bits; cancel) -> (θ, result, guard_bits)

Refine a `Float64` seed of a fully symmetric triangle rule to `bits` bits: guard bits from
the condition number measured at the seed, Gauss–Newton at `bits + guard`, and one
re-run with more guard if the condition number found along the way demands it.
"""
function refine_symmetric_triangle(structure::SymmetricStructure, n::Integer, θ64::Vector{Float64},
                                   bits::Integer; cancel = nothing, basis = invariant_basis(3, n))
    sys64 = TriangleMomentSystem(structure, n, Float64, basis)
    r64, J64 = sys64(θ64)
    κ0 = lsq_step(J64, r64; rank_rtol = 1e-14)[2]
    guard = guard_bits_from_cond(κ0)
    for attempt in 1:2
        wbits = bits + guard
        res = with_bits(wbits) do
            sys = TriangleMomentSystem(structure, n, BigFloat, basis)
            θ0 = BigFloat.(θ64)
            gauss_newton(sys, θ0; step_tol = ldexp(BigFloat(1), -(bits + 16)),
                         res_floor = ldexp(BigFloat(1), -(wbits - 12)),
                         rank_rtol = ldexp(BigFloat(1), -(wbits ÷ 2)), maxiter = 60, cancel)
        end
        needed = guard_bits_from_cond(res.cond_max)
        if needed <= guard || attempt == 2
            return res.θ, res, guard
        end
        guard = needed
    end
end

function build(f::XiaoGimbutas, dom::Simplex{2}, degree::Int, ctx::BuildContext{T};
               seed::SeedSource = TableSeed()) where {T}
    isexact(ctx) && throw(ArgumentError("XiaoGimbutas nodes are irrational; $(T) is not supported"))
    e = _xg_entry(degree)
    n = e.degree
    basis = invariant_basis(3, n)
    structure, θ64, seed_desc = e.structure, e.seed, "stored Float64 table src/data/triangle_s3_seeds.toml " *
                                                     "(generated in-house by multistart from the orbit structure)"
    if seed isa MultistartSeed
        found = multistart(structure, n; nstarts = seed.nstarts, rng_seed = seed.rng_seed, basis,
                           first_only = false, cancel = ctx.cancel)
        isempty(found) && throw(RefinementError("XiaoGimbutas",
            "multistart found no valid degree-$n rule with structure $structure in $(seed.nstarts) starts"))
        best = argmax(t -> (t[3], -t[4]), found)
        θ64 = best[1]
        seed_desc = describe(seed) * "; $(length(found)) valid rule(s) found, largest minimum barycentric kept"
    end
    if seed isa ExplicitSeed
        length(seed.θ) == nunknowns(structure) || throw(ArgumentError(
            "an explicit seed for degree $n needs $(nunknowns(structure)) parameters (structure $structure)"))
        θ64, seed_desc = seed.θ, describe(seed)
    end
    θ, res, guard = refine_symmetric_triangle(structure, n, θ64, ctx.bits; cancel = ctx.cancel, basis)
    res.converged || throw(RefinementError("XiaoGimbutas",
        "Gauss–Newton did not converge at degree $n (residual $(Float64(res.residual)) after $(res.iterations) " *
        "iterations); the seed is not returned unrefined"))
    wbits = ctx.bits + guard
    λs, ws = with_bits(() -> expand(structure, θ), wbits)
    wmin = minimum(ws)
    λmin = minimum(minimum, λs)
    (wmin > 0 && λmin > 0) || throw(RefinementError("XiaoGimbutas",
        "refined degree-$n rule left the admissible set (min weight $(Float64(wmin)), min barycentric $(Float64(λmin)))"))
    # round once: each distinct barycentric value and weight is rounded, then permuted
    xs = [SVector{2,T}(finalize_number(ctx, λ[2]), finalize_number(ctx, λ[3])) for λ in λs]
    wt = [finalize_number(ctx, w) for w in ws]
    # certificate: defining-equation residual at the delivered parameters, at 2× precision
    θhat = [finalize_number(ctx, t) for t in θ]
    rbits = 2ctx.bits + guard
    resid = with_bits(rbits) do
        sys = TriangleMomentSystem(structure, n, BigFloat, basis)
        r, _ = sys(BigFloat.(θhat))
        maximum(abs, r)
    end
    cert = Certificate(equations = "S₃-invariant moment system in orbit parameters (orthonormal Dubiner basis, degree $n)",
                       residual = BigFloat(resid; precision = 64), residual_bits = rbits,
                       digits = target_digits(ctx), guard_digits = floor(Int, guard * log10(2)),
                       cond = res.cond, iterations = res.iterations)
    prov = Provenance(family = "XiaoGimbutas", derivation = Seeded(),
                      path = ["seed: " * seed_desc,
                              "structure: " * string(structure),
                              "refine: Gauss–Newton with rank-revealing SVD solve at $(wbits) bits",
                              @sprintf("guard: %d bits from measured cond(J) = %.2e", guard, res.cond_max)],
                      seed_source = seed_desc, citations = [XIAO_GIMBUTAS_2010], symmetry = :S3,
                      license = "MIT (seeds generated in-house; point counts from Xiao & Gimbutas 2010)")
    return QuadratureRule(xs, wt, Simplex{2}(), PolynomialDegree(n), prov, cert)
end
