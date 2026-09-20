# Node elimination (PLAN §6 Tier 4), pulled forward by the v0.2 gate: a from-scratch
# multistart search for tetrahedral rules grew 3–5× more expensive per degree.
#
# Two operations, both working on orbit structures so that symmetry is kept exactly:
#
#   grow       start from a rule of the previous degree, which already integrates every
#              moment but the top ones, add orbits with small weights, and refit: a solution
#              of an underdetermined system near a known point is easy to find.
#   eliminate  repeatedly try to remove points — drop an orbit, or merge two of an orbit's
#              barycentric values so that it becomes a smaller orbit type — and refit the
#              remaining parameters onto the moment variety. A move is accepted if the fit
#              converges with positive weights and interior nodes. Stop when no move works.
#
# Everything here runs in Float64; the result is a seed for `refine_symmetric`.

"A candidate elimination move: the new structure, its starting parameters, points removed."
struct EliminationMove
    structure::SymmetricStructure
    θ0::Vector{Float64}
    removed::Int
    score::Float64
    what::String
end

# Orbit blocks as (mult, [w, v₁, …]) pairs, and back.
function orbit_blocks(s::SymmetricStructure, θ)
    return [(copy(o.mult), θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
end
function from_blocks(N, blocks)
    s = SymmetricStructure([b[1] for b in blocks], N)
    return canonicalize(s, reduce(vcat, (b[2] for b in blocks); init = Float64[]))
end

"All drop and merge moves for the rule `(s, θ)`, most promising first."
function elimination_moves(s::SymmetricStructure, θ::Vector{Float64})
    N = s.N
    blocks = orbit_blocks(s, θ)
    moves = EliminationMove[]
    for (i, (mult, p)) in enumerate(blocks)
        o = s.orbits[i]
        w = p[1]
        mass = w * orbit_size(o)
        # drop the orbit, spreading its mass over the others
        rest = [(m, copy(q)) for (j, (m, q)) in enumerate(blocks) if j != i]
        if !isempty(rest)
            total = sum(q[1] * orbit_size(OrbitPattern(m, N)) for (m, q) in rest)
            for (_, q) in rest
                q[1] *= (total + mass) / total
            end
            st, θ0 = from_blocks(N, rest)
            push!(moves, EliminationMove(st, θ0, orbit_size(o), mass, "drop orbit $(mult)"))
        end
        # merge two distinct values of the orbit into one
        vals = pattern_values(o, p[2:end])
        r = length(mult)
        has_centroid = any(o2 -> length(o2.mult) == 1, s.orbits)
        for a in 1:(r - 1), b in (a + 1):r
            (r == 2 && has_centroid) && continue          # would duplicate the centroid
            newmult = copy(mult)
            newvals = copy(vals)
            v = (mult[a] * vals[a] + mult[b] * vals[b]) / (mult[a] + mult[b])
            newmult[a] += mult[b]
            newvals[a] = v
            deleteat!(newmult, b)
            deleteat!(newvals, b)
            order = sortperm(collect(zip(-newmult, newvals)))
            newmult, newvals = newmult[order], newvals[order]
            newpat = OrbitPattern(newmult, N)
            neww = mass / orbit_size(newpat)
            merged = (newmult, vcat(neww, newvals[1:(end - 1)]))
            st, θ0 = from_blocks(N, [j == i ? merged : blocks[j] for j in eachindex(blocks)])
            push!(moves, EliminationMove(st, θ0, orbit_size(o) - orbit_size(newpat), abs(vals[a] - vals[b]),
                                         "merge $(mult) → $(newmult)"))
        end
    end
    # most points removed first; among equals, the least disruptive (small mass or close values)
    sort!(moves; by = m -> (-m.removed, m.score))
    return moves
end

"Fit `(s, θ0)` onto the degree-`n` moment variety; the canonical parameters if valid, else `nothing`."
function fit_symmetric(s::SymmetricStructure, θ0::Vector{Float64}, n::Integer, basis::InvariantBasis;
                       tol::Float64 = 1e-13, maxiter::Int = 400)
    sys = SymmetricMomentSystem(s, n, Float64, basis)
    inside(θ) = rule_margins(s, θ)[2] > -0.05
    θ, nr = levenberg_marquardt(sys, θ0; maxiter, tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    wmin, λmin = rule_margins(s, res.θ)
    (wmin > 0 && λmin > 0) || return nothing
    return canonicalize(s, res.θ)[2]
end

"""
    eliminate(s, θ, n; basis, log = nothing, rng = nothing) -> (s, θ)

Greedy node elimination on a valid degree-`n` rule: try the moves of
[`elimination_moves`](@ref) in order, in parallel batches, accept the first (in that order)
whose refit is valid, and repeat until no move succeeds. With `rng`, moves removing the
same number of points are tried in random order, so that repeated chains explore different
local minima.
"""
function eliminate(s::SymmetricStructure, θ::Vector{Float64}, n::Integer;
                   basis::InvariantBasis = invariant_basis(s.N, n), log = nothing, cancel = nothing,
                   rng = nothing)
    m = size(basis.Q, 2)
    while true
        checkcancel(cancel)
        moves = filter(mv -> nunknowns(mv.structure) >= m - 2, elimination_moves(s, θ))
        rng === nothing || sort!(moves; by = mv -> (-mv.removed, rand(rng)))
        accepted = nothing
        batch = max(Threads.nthreads(), 1)
        for lo in 1:batch:length(moves)
            idx = lo:min(lo + batch - 1, length(moves))
            results = Vector{Any}(nothing, length(idx))
            Threads.@threads :static for j in eachindex(idx)
                move = moves[idx[j]]      # a name not used outside the threaded body
                results[j] = fit_symmetric(move.structure, move.θ0, n, basis)
            end
            k = findfirst(!isnothing, results)
            if k !== nothing
                accepted = (moves[idx[k]], results[k])
                break
            end
        end
        accepted === nothing && return s, θ
        mv, θn = accepted
        log === nothing || println(log, "    eliminate: ", mv.what, " → ", npoints(mv.structure), " points, ",
                                   nunknowns(mv.structure), " unknowns")
        s, θ = mv.structure, θn
    end
end

"""
    grow(s, θ, n; basis, ntries = 64, rng_seed) -> (s, θ) or nothing

Turn a rule of degree `< n` into a valid degree-`n` rule by adding orbits with small
weights until the system has comfortably more unknowns than equations, then refitting from
`ntries` random placements of the new orbits. Several recipes for the added orbits are
tried in turn (one of each smaller type plus general orbits; general orbits only; a larger
surplus), since which one works depends on the degree.
"""
function grow(s::SymmetricStructure, θ::Vector{Float64}, n::Integer;
              basis::InvariantBasis = invariant_basis(s.N, n), ntries::Int = 64, rng_seed::Integer = 0x9e0)
    N = s.N
    m = size(basis.Q, 2)
    types = sort(orbit_pattern_types(N); by = p -> (-length(p.mult), p.mult))
    general, small = types[1], types[2:(end - 1)]            # never a second centroid
    recipes = [(small, m + max(2, m ÷ 4)), (OrbitPattern[], m + max(2, m ÷ 4)),
               (small, m + max(4, m ÷ 2)), (vcat(small, small), m + max(4, m ÷ 2))]
    for (k, (extra, target)) in enumerate(recipes)
        added = copy(extra)
        unknowns = nunknowns(s) + sum(nunknowns, added; init = 0)
        while unknowns < target
            push!(added, general)
            unknowns += nunknowns(general)
        end
        got = _grow_try(s, θ, n, basis, added, ntries, rng_seed + k)
        got === nothing || return got
    end
    return nothing
end

function _grow_try(s, θ, n, basis, added, ntries, rng_seed)
    N = s.N
    blocks = orbit_blocks(s, θ)
    mass = 1 / factorial(N - 1)
    results = Vector{Any}(nothing, ntries)
    Threads.@threads :static for t in 1:ntries
        rng = start_rng(rng_seed, t)
        new = Tuple{Vector{Int},Vector{Float64}}[]
        for p in added
            u = -log.(rand(rng, length(p.mult)))
            u ./= sum(u)
            push!(new, (copy(p.mult), vcat(0.02 * mass / orbit_size(p), [u[i] / p.mult[i] for i in 1:(length(p.mult) - 1)])))
        end
        st, θ0 = from_blocks(N, vcat(blocks, new))
        θn = fit_symmetric(st, θ0, n, basis)
        results[t] = θn === nothing ? nothing : (st, θn)
    end
    k = findfirst(!isnothing, results)
    return k === nothing ? nothing : results[k]
end

"""
    grow_and_eliminate(s, θ, n; chains = 16, basis, rng_seed) -> (s, θ, npoints per chain) or nothing

Run `chains` independent grow → eliminate chains from the lower-degree rule `(s, θ)`,
each with its own random placement of the added orbits and its own move order, and keep
the rule with the fewest points (ties: the largest minimum barycentric coordinate).
Deterministic for a given `rng_seed`.
"""
function grow_and_eliminate(s::SymmetricStructure, θ::Vector{Float64}, n::Integer; chains::Int = 16,
                            basis::InvariantBasis = invariant_basis(s.N, n), rng_seed::Integer = 0xe11,
                            cancel = nothing)
    best = nothing
    counts = Int[]
    for c in 1:chains
        checkcancel(cancel)
        g = grow(s, θ, n; basis, ntries = 16, rng_seed = rng_seed + 1000c)
        g === nothing && continue
        se, θe = eliminate(g[1], g[2], n; basis, rng = start_rng(rng_seed, c), cancel)
        push!(counts, npoints(se))
        key = (npoints(se), -rule_margins(se, θe)[2])
        if best === nothing || key < best[3]
            best = (se, θe, key)
        end
    end
    return best === nothing ? nothing : (best[1], best[2], counts)
end
