# Seeds for fully symmetric rules on the square and the cube, and on the disk and the ball
# (PLAN §6 Tier 4), all in Float64: multistart over the orbit structures of a given point count
# at low degree, and grow → eliminate above it, as for the simplex (refine/elimination.jl), on
# the orbits of symmetry/box.jl. The orbits are the same on the box and the ball; what differs
# is the moment system and where an orbit's values may lie, which an `OrbitGeometry` supplies.
# The results are seeds for `refine_box` and `refine_round`.

abstract type OrbitGeometry end
"The box `[-1, 1]^D`, with the symmetrised Legendre moments of [`BoxMomentSystem`](@ref)."
struct BoxGeometry <: OrbitGeometry end
"The unit disk or ball, with the moments of [`RoundMomentSystem`](@ref)."
struct RoundGeometry <: OrbitGeometry end
orbit_system(::BoxGeometry, s, n, ::Type{S}) where {S} = BoxMomentSystem(s, n, S)
orbit_system(::RoundGeometry, s, n, ::Type{S}) where {S} = RoundMomentSystem(s, n, S)
orbit_equations(::BoxGeometry, D, n) = length(box_exponents(D, n))
orbit_equations(::RoundGeometry, D, n) = length(round_terms(D, n))
orbit_margins(::BoxGeometry, s, θ) = box_margins(s, θ)
orbit_margins(::RoundGeometry, s, θ) = round_margins(s, θ)
orbit_volume(::BoxGeometry, D) = 2.0^D
orbit_volume(::RoundGeometry, D) = D == 2 ? Float64(π) : 4π / 3
# random values for an orbit with these multiplicities, inside the domain
random_values(::BoxGeometry, mult, rng) = 0.03 .+ 0.94 .* rand(rng, length(mult))
random_values(::RoundGeometry, mult, rng) = (0.03 .+ 0.94 .* rand(rng, length(mult))) ./ sqrt(max(1, sum(mult; init = 0)))
# while fitting, values may wander through zero (the orbit is symmetric) but not far outside
loosely_inside(::BoxGeometry, o, v) = all(x -> abs(x) < 1.05, v)
loosely_inside(::RoundGeometry, o, v) = sum((o.mult[i] * v[i]^2 for i in eachindex(v)); init = 0.0) < 1.1

"Fit `(s, θ0)` onto the degree-`n` moment variety; the canonical `(s, θ)` if it is a positive interior rule, else `nothing`."
function box_fit(s::BoxStructure, θ0::Vector{Float64}, n::Integer; tol::Float64 = 1e-13, maxiter::Int = 400,
                 geom::OrbitGeometry = BoxGeometry())
    sys = orbit_system(geom, s, n, Float64)
    inside(θ) = all(o_off -> loosely_inside(geom, o_off[1], [θ[o_off[2] + 1 + i] for i in 1:nparams(o_off[1])]),
                    zip(s.orbits, param_offsets(s)))
    θ, nr = levenberg_marquardt(sys, θ0; maxiter, tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    sc, θc = canonicalize(s, res.θ)
    wmin, cmin, gmin = orbit_margins(geom, sc, θc)
    (wmin > 0 && cmin > 1e-6 && gmin > 1e-6) || return nothing
    return sc, θc
end

function random_box_start(s::BoxStructure, rng; geom::OrbitGeometry = BoxGeometry())
    V, N = orbit_volume(geom, s.D), npoints(s)
    θ = Float64[]
    for o in s.orbits
        push!(θ, V / N * (0.5 + rand(rng)))
        append!(θ, random_values(geom, o.mult, rng))
    end
    return θ
end

"""
    box_structures(npts, n, D; slack = 2)

Every fully symmetric orbit structure on `[-1, 1]^D` with exactly `npts` points and between
`m` and `m + slack` unknowns, `m` the number of degree-`n` moment equations.
"""
function box_structures(npts::Integer, n::Integer, D::Integer; slack::Integer = 2, geom::OrbitGeometry = BoxGeometry())
    types = box_orbit_types(D)
    m = orbit_equations(geom, D, n)
    out = BoxStructure[]
    function rec(i, left, counts)
        if i > length(types)
            left == 0 || return
            s = BoxStructure(D, reduce(vcat, [fill(types[k], counts[k]) for k in eachindex(types)]; init = BoxOrbit[]))
            m <= nunknowns(s) <= m + slack && push!(out, s)
            return
        end
        t = types[i]
        maxc = isempty(t.mult) ? min(1, left) : left ÷ orbit_size(t)
        for c in 0:maxc
            rec(i + 1, left - c * orbit_size(t), vcat(counts, c))
        end
    end
    rec(1, Int(npts), Int[])
    return out
end

"""
    box_multistart(n, D; npts, nstarts = 200, slack = 2, rng_seed) -> (s, θ) or nothing

The fully symmetric rule of degree `n` on `[-1, 1]^D` with the fewest points in `npts`: for
each count in turn, `nstarts` random starts on every orbit structure with that many points.
Among the rules found at the first count that has any, the one whose nodes keep farthest
from the boundary.
"""
function box_multistart(n::Integer, D::Integer; npts, nstarts::Integer = 200, slack::Integer = 2,
                        rng_seed::Integer = 0xb0c5, cancel = nothing, geom::OrbitGeometry = BoxGeometry())
    for N in npts
        found = Any[]
        for s in box_structures(N, n, D; slack, geom)
            results = Vector{Any}(nothing, nstarts)
            Threads.@threads :static for t in 1:nstarts
                checkcancel(cancel)
                rng = start_rng(rng_seed + N, t)
                results[t] = box_fit(s, random_box_start(s, rng; geom), n; geom)
            end
            append!(found, filter(!isnothing, results))
        end
        isempty(found) && continue
        return found[argmax([orbit_margins(geom, f...)[2] for f in found])]
    end
    return nothing
end

# --- grow and eliminate ------------------------------------------------------------------------

# an orbit and its values from group sizes and values in any order
function make_orbit(mult::Vector{Int}, vals::Vector{Float64}, D)
    order = sortperm(mult; rev = true)
    return BoxOrbit(mult[order], D), vals[order]
end
function from_box_blocks(D, blocks)
    s = BoxStructure(D, [b[1] for b in blocks])
    return canonicalize(s, reduce(vcat, (b[2] for b in blocks); init = Float64[]))
end

"All drop, merge and zero moves on the rule `(s, θ)`: (structure, θ0, points removed, score)."
function box_moves(s::BoxStructure, θ::Vector{Float64})
    D = s.D
    blocks = [(o, θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
    has_centre = any(b -> isempty(b[1].mult), blocks)
    moves = Any[]
    for (i, (o, p)) in enumerate(blocks)
        mass = p[1] * orbit_size(o)
        v = p[2:end]
        replaced(no, nv) = from_box_blocks(D, [j == i ? (no, vcat(mass / orbit_size(no), nv)) : blocks[j]
                                               for j in eachindex(blocks)])
        # drop the orbit, spreading its mass over the others
        rest = [(b[1], copy(b[2])) for (j, b) in enumerate(blocks) if j != i]
        if !isempty(rest)
            total = sum(b[2][1] * orbit_size(b[1]) for b in rest)
            foreach(b -> b[2][1] *= (total + mass) / total, rest)
            push!(moves, (from_box_blocks(D, rest)..., orbit_size(o), mass))
        end
        # merge two of the orbit's values into one
        for k in 1:(length(v) - 1), l in (k + 1):length(v)
            mult = copy(o.mult)
            nv = copy(v)
            nv[k] = (mult[k] * v[k] + mult[l] * v[l]) / (mult[k] + mult[l])
            mult[k] += mult[l]
            deleteat!(mult, l)
            deleteat!(nv, l)
            no, nvs = make_orbit(mult, nv, D)
            push!(moves, (replaced(no, nvs)..., orbit_size(o) - orbit_size(no), abs(v[k] - v[l])))
        end
        # send one of its values to zero
        for k in eachindex(v)
            mult = deleteat!(copy(o.mult), k)
            (isempty(mult) && has_centre) && continue           # never a second centre
            no, nvs = make_orbit(mult, deleteat!(copy(v), k), D)
            push!(moves, (replaced(no, nvs)..., orbit_size(o) - orbit_size(no), abs(v[k])))
        end
    end
    return moves
end

"""
    box_eliminate(s, θ, n; rng) -> (s, θ)

Greedy node elimination on a valid degree-`n` box rule: try the moves of [`box_moves`](@ref)
that remove the most points first, accept the first whose refit is valid, repeat until none
is. With `rng`, moves removing the same number of points are tried in random order.
"""
function box_eliminate(s::BoxStructure, θ::Vector{Float64}, n::Integer; rng = nothing, cancel = nothing,
                       geom::OrbitGeometry = BoxGeometry())
    m = orbit_equations(geom, s.D, n)
    while true
        checkcancel(cancel)
        moves = filter(mv -> nunknowns(mv[1]) >= m, box_moves(s, θ))
        sort!(moves; by = mv -> (-mv[3], rng === nothing ? mv[4] : rand(rng)))
        accepted = nothing
        batch = max(Threads.nthreads(), 1)
        for lo in 1:batch:length(moves)
            idx = lo:min(lo + batch - 1, length(moves))
            results = Vector{Any}(nothing, length(idx))
            Threads.@threads :static for j in eachindex(idx)
                mv = moves[idx[j]]
                results[j] = box_fit(mv[1], mv[2], n; geom)
            end
            k = findfirst(!isnothing, results)
            k === nothing || (accepted = results[k]; break)
        end
        accepted === nothing && return s, θ
        s, θ = accepted
    end
end

"""
    box_grow(s, θ, n; ntries = 32, rng_seed) -> (s, θ) or nothing

From a rule of lower degree, add general orbits with small weights until the unknowns
exceed the degree-`n` equations comfortably, and refit.
"""
function box_grow(s::BoxStructure, θ::Vector{Float64}, n::Integer; ntries::Int = 32, rng_seed::Integer = 0x96,
                  geom::OrbitGeometry = BoxGeometry())
    D = s.D
    m = orbit_equations(geom, D, n)
    general = last(box_orbit_types(D))
    blocks = [(o, θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
    for (k, target) in enumerate((m + max(2, m ÷ 4), m + max(4, m ÷ 2)))
        added = 0
        while nunknowns(s) + added * nunknowns(general) < target
            added += 1
        end
        results = Vector{Any}(nothing, ntries)
        Threads.@threads :static for t in 1:ntries
            rng = start_rng(rng_seed + k, t)
            new = [(general, vcat(0.02 * orbit_volume(geom, D) / orbit_size(general), random_values(geom, general.mult, rng)))
                   for _ in 1:added]
            st, θ0 = from_box_blocks(D, vcat(blocks, new))
            results[t] = box_fit(st, θ0, n; geom)
        end
        i = findfirst(!isnothing, results)
        i === nothing || return results[i]
    end
    return nothing
end

"""
    box_grow_and_eliminate(s, θ, n; chains = 8, rng_seed) -> (s, θ, points per chain) or nothing

Independent grow → eliminate chains from the lower-degree rule `(s, θ)`; the rule with the
fewest points, ties broken by the largest distance of the nodes from the boundary.
"""
function box_grow_and_eliminate(s::BoxStructure, θ::Vector{Float64}, n::Integer; chains::Int = 8,
                                rng_seed::Integer = 0xe11, cancel = nothing, geom::OrbitGeometry = BoxGeometry())
    best, counts = nothing, Int[]
    for c in 1:chains
        checkcancel(cancel)
        g = box_grow(s, θ, n; rng_seed = rng_seed + 1000c, geom)
        g === nothing && continue
        se, θe = box_eliminate(g[1], g[2], n; rng = start_rng(rng_seed, c), cancel, geom)
        push!(counts, npoints(se))
        key = (npoints(se), -orbit_margins(geom, se, θe)[2])
        (best === nothing || key < best[3]) && (best = (se, θe, key))
    end
    return best === nothing ? nothing : (best[1], best[2], counts)
end
