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

Search for fully symmetric triangle rules of degree `n` with the given orbit structure.
Start `i` draws from its own RNG seeded with `(rng_seed, i)`, and starts run in fixed
batches of 64 across threads, so the result does not depend on the number of threads.
Returns a vector of `(θ, min_weight, min_barycentric, start_index)` in canonical form,
deduplicated, ordered by start index. With `first_only`, stops after the first batch that
contains a valid rule.
"""
function multistart(structure::SymmetricStructure, n::Integer; nstarts::Integer = 2000,
                    rng_seed::Integer = 0x5eed, basis::InvariantBasis = invariant_basis(3, n),
                    tol::Float64 = 1e-13, first_only::Bool = true, cancel = nothing)
    found = Vector{Tuple{Vector{Float64},Float64,Float64,Int}}()
    batch = 64
    for lo in 1:batch:nstarts
        checkcancel(cancel)
        idx = lo:min(lo + batch - 1, nstarts)
        results = Vector{Any}(nothing, length(idx))
        Threads.@threads :static for t in 1:Threads.nthreads()
            sys = TriangleMomentSystem(structure, n, Float64, basis)
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
function multistart_attempt(sys::TriangleMomentSystem{Float64}, structure, rng, tol)
    θ0 = random_parameters(rng, structure, 0.5)
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
    candidate_structures(npts, n, N = 3)

Orbit structures on the triangle with exactly `npts` points and at least as many unknowns
as `S_3`-invariant equations for degree `n`, ordered by number of unknowns (fewest first),
then by the orbit counts.
"""
function candidate_structures(npts::Integer, n::Integer, N::Integer = 3)
    N == 3 || throw(NotYetImplemented("orbit structure enumeration for S_$N", "v0.2"))
    m = n_invariants(3, n)
    out = SymmetricStructure[]
    for c in 0:1, a in 0:(npts ÷ 3)
        rest = npts - c - 3a
        (rest >= 0 && rest % 6 == 0) || continue
        b = rest ÷ 6
        unknowns = c + 2a + 3b
        unknowns >= m || continue
        # canonical orbit order (see `canonicalize`): general, vertex-type, centroid
        pats = vcat(fill([1, 1, 1], b), fill([2, 1], a), fill([3], c))
        push!(out, SymmetricStructure(pats, 3))
    end
    sort!(out; by = s -> (nunknowns(s), [count(o -> o.mult == p, s.orbits) for p in ([1, 1, 1], [2, 1], [3])]))
    return out
end

# One independent stream per start index. (Julia does not promise identical RNG streams
# across versions, which is one reason shipped seeds are stored rather than regenerated.)
start_rng(seed::Integer, i::Integer) = Random.Xoshiro(UInt64(seed) ⊻ (UInt64(i) * 0x9e3779b97f4a7c15))
