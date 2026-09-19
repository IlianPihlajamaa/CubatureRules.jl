# Seed sources behind one interface (PLAN §1): a stored table, or an orbit structure plus
# multistart search. Both yield a Float64 parameter vector for the refiner.

"""
    SeedSource

Where a refinement starts. [`TableSeed`](@ref) reads a stored `Float64` table;
[`MultistartSeed`](@ref) needs only the orbit structure.
"""
abstract type SeedSource end

"""
    TableSeed()

Use the stored seed table shipped in `src/data`.
"""
struct TableSeed <: SeedSource end

"""
    MultistartSeed(; nstarts = 2000, rng_seed = 0x5eed)

Recover the rule from its orbit structure alone: Levenberg–Marquardt from `nstarts`
random interior starting points (a fixed-seed RNG, so the search is repeatable), keeping
the first solution with positive weights and interior nodes. This is the licence-free
route of PLAN §0.2.
"""
Base.@kwdef struct MultistartSeed <: SeedSource
    nstarts::Int = 2000
    rng_seed::UInt64 = 0x5eed
end

"""
    ExplicitSeed(θ; source = "user-supplied")

Start from a given `Float64` parameter vector (same layout as the stored table), e.g. one
transcribed from a published table whose licence permits it. The source is recorded in
the provenance.
"""
struct ExplicitSeed <: SeedSource
    θ::Vector{Float64}
    source::String
end
ExplicitSeed(θ::AbstractVector; source = "user-supplied") = ExplicitSeed(Float64.(θ), String(source))

describe(s::TableSeed) = "stored table"
describe(s::ExplicitSeed) = "explicit parameter vector (" * s.source * ")"
describe(s::MultistartSeed) = "orbit structure + multistart ($(s.nstarts) starts, rng seed $(repr(s.rng_seed)))"

"""
    multistart(structure, n; nstarts, rng_seed, basis, tol = 1e-13, first_only = true)

Search for fully symmetric simplex rules of degree `n` with the given orbit structure.
Start `i` draws from its own RNG seeded with `(rng_seed, i)`, and starts run in fixed
batches of 64 across threads, so the result does not depend on the number of threads.
Returns a vector of `(θ, min_weight, min_barycentric, start_index)` in canonical form,
deduplicated, ordered by start index. With `first_only`, stops after the first batch that
contains a valid rule.
"""
function multistart(structure::SymmetricStructure, n::Integer; nstarts::Integer = 2000,
                    rng_seed::Integer = 0x5eed, basis::InvariantBasis = invariant_basis(structure.N, n),
                    tol::Float64 = 1e-13, first_only::Bool = true, cancel = nothing)
    found = Vector{Tuple{Vector{Float64},Float64,Float64,Int}}()
    batch = 64
    for lo in 1:batch:nstarts
        checkcancel(cancel)
        idx = lo:min(lo + batch - 1, nstarts)
        results = Vector{Any}(nothing, length(idx))
        Threads.@threads :static for t in 1:Threads.nthreads()
            sys = SymmetricMomentSystem(structure, n, Float64, basis)
            for j in t:Threads.nthreads():length(idx)
                results[j] = multistart_attempt(sys, structure, start_rng(rng_seed, idx[j]), tol)
            end
        end
        for (j, res) in enumerate(results)
            res === nothing && continue
            θc, wmin, λmin = res
            any(f -> maximum(abs, f[1] - θc) < 1e-8, found) && continue
            push!(found, (θc, wmin, λmin, idx[j]))
        end
        first_only && !isempty(found) && break
    end
    return found
end

"One multistart attempt: LM from a random interior point, polished by Gauss–Newton."
function multistart_attempt(sys::SymmetricMomentSystem{Float64}, structure, rng, tol)
    θ0 = random_parameters(rng, structure, 1 / factorial(structure.N - 1))
    inside(θ) = rule_margins(structure, θ)[2] > -0.05
    θ, nr = levenberg_marquardt(sys, θ0; maxiter = 300, tol = tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    wmin, λmin = rule_margins(structure, res.θ)
    (wmin > 0 && λmin > 0) || return nothing
    _, θc = canonicalize(structure, res.θ)
    return θc, wmin, λmin
end

"""
    candidate_structures(npts, n, N = 3; max_unknowns = typemax(Int))

Orbit structures on the `(N-1)`-simplex with exactly `npts` points and at least as many
unknowns as `S_N`-invariant equations for degree `n` (and at most `max_unknowns`), ordered
by number of unknowns (fewest first), then by the orbit counts in canonical order. The
one-point orbit `[N]` (the centroid) is used at most once.
"""
function candidate_structures(npts::Integer, n::Integer, N::Integer = 3; max_unknowns::Integer = typemax(Int))
    m = n_invariants(N, n)
    types = sort(orbit_pattern_types(N); by = p -> (-length(p.mult), p.mult))   # canonical order
    sizes = orbit_size.(types)
    unk = nunknowns.(types)
    out = SymmetricStructure[]
    counts = zeros(Int, length(types))
    function rec(i, pts, u)
        if i > length(types)
            (pts == npts && m <= u <= max_unknowns) || return
            pats = reduce(vcat, [fill(types[j].mult, counts[j]) for j in eachindex(types)]; init = Vector{Int}[])
            push!(out, SymmetricStructure(pats, N))
            return
        end
        cmax = sizes[i] == 1 ? 1 : (npts - pts) ÷ sizes[i]
        for c in 0:cmax
            counts[i] = c
            rec(i + 1, pts + c * sizes[i], u + c * unk[i])
        end
        counts[i] = 0
    end
    rec(1, 0, 0)
    key(s) = (nunknowns(s), [count(o -> o.mult == t.mult, s.orbits) for t in types])
    sort!(out; by = key)
    return out
end

# One independent stream per start index. (Julia does not promise identical RNG streams
# across versions, which is one reason shipped seeds are stored rather than regenerated.)
start_rng(seed::Integer, i::Integer) = Random.Xoshiro(UInt64(seed) ⊻ (UInt64(i) * 0x9e3779b97f4a7c15))
