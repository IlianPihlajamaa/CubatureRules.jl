# Growing and eliminating O_h-symmetric sphere rules (PLAN §6 Tier 4), the Lebedev
# counterpart of `refine/elimination.jl`.
#
# Multistart alone reproduces Lebedev's counts to degree 17 and then stalls. Degree 19 needs
# twelve invariant equations, so the structures that can meet it have twelve-dimensional
# parameter spaces, and random starts land in the right basin too rarely to be worth waiting
# for. The simplex tables hit the same wall at degree 20, and the answer there was to stop
# starting cold:
#
#   grow       take the rule two degrees down — it already integrates every invariant but
#              the top ones — add light orbits, and refit. Solving an underdetermined system
#              near a known point is easy.
#   eliminate  repeatedly try to remove points: drop an orbit, or collapse it onto a smaller
#              orbit type, and refit the rest. Accept whatever still lands on a positive,
#              non-degenerate rule. Stop when nothing does.
#
# Everything here is Float64; the result is a seed for `refine_octahedral`.

"A candidate elimination move: the new structure, its starting parameters, points removed."
struct OctahedralMove
    structure::OctahedralStructure
    θ0::Vector{Float64}
    removed::Int
    what::String
end

octahedral_blocks(s::OctahedralStructure, θ) =
    [(o.kind, Float64.(θ[(off + 1):(off + nunknowns(o))])) for (o, off) in zip(s.orbits, param_offsets(s))]

function from_octahedral_blocks(blocks)
    s = OctahedralStructure([b[1] for b in blocks])
    return s, reduce(vcat, (b[2] for b in blocks); init = Float64[])
end

"The orbit types an orbit can collapse onto, each of them smaller."
octahedral_degradations(kind::Symbol) = kind === :d ? (:b, :c) :
                                        kind in (:b, :c) ? (:a3, :a1, :a2) : ()

"Parameters for `kind → to`, keeping the orbit's mass and a plausible coordinate."
function degrade_params(kind::Symbol, to::Symbol, p::Vector{Float64})
    w = p[1] * orbit_size(OctahedralOrbit(kind)) / orbit_size(OctahedralOrbit(to))
    to in (:a1, :a2, :a3) && return [w]
    return [w, clamp(p[2], 0.05, to === :b ? 0.68 : 0.95)]
end

"""
    octahedral_elimination_moves(s, θ)

Every way to make the rule one orbit smaller: drop the orbit outright, or collapse it onto a
smaller type. Largest saving first. The parameterless orbits are unique, so a collapse onto
one already present is not offered.
"""
function octahedral_elimination_moves(s::OctahedralStructure, θ::Vector{Float64})
    blocks = octahedral_blocks(s, θ)
    present = Set(o.kind for o in s.orbits)
    moves = OctahedralMove[]
    for (i, (kind, p)) in enumerate(blocks)
        o = s.orbits[i]
        mass = p[1] * orbit_size(o)
        rest = [(k, copy(q)) for (j, (k, q)) in enumerate(blocks) if j != i]
        if !isempty(rest)
            total = sum(q[1] * orbit_size(OctahedralOrbit(k)) for (k, q) in rest)
            for (_, q) in rest
                q[1] *= (total + mass) / total            # the dropped orbit's mass, spread
            end
            st, θ0 = from_octahedral_blocks(rest)
            push!(moves, OctahedralMove(st, θ0, orbit_size(o), "drop $(kind)"))
        end
        for to in octahedral_degradations(kind)
            # only a₁, a₂ and a₃ are unique; several b or c orbits differing in their
            # parameter are perfectly ordinary, so blocking those would discard real moves
            (to in (:a1, :a2, :a3) && to in present) && continue
            blocks2 = [(j == i ? (to, degrade_params(kind, to, p)) : (k, copy(q)))
                       for (j, (k, q)) in enumerate(blocks)]
            st, θ0 = from_octahedral_blocks(blocks2)
            push!(moves, OctahedralMove(st, θ0, orbit_size(o) - orbit_size(OctahedralOrbit(to)),
                                        "$(kind) to $(to)"))
        end
    end
    return sort!(moves; by = mv -> -mv.removed)
end

"""
    fit_octahedral(s, θ0, n; tol, maxiter)

Land `θ0` on the degree-`n` moment variety for structure `s`, or `nothing` if it does not
reach a rule that is positive, non-degenerate and has distinct nodes.
"""
function fit_octahedral(s::OctahedralStructure, θ0::Vector{Float64}, n::Integer;
                        tol::Float64 = 1e-13, maxiter::Int = 400)
    sys = OctahedralMomentSystem(s, n, Float64)
    inside(θ) = octahedral_margins(s, θ)[2] > -0.02
    θ, nr = levenberg_marquardt(sys, θ0; maxiter, tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    wmin, dmin = octahedral_margins(s, res.θ)
    (wmin > 0 && dmin > 1e-6) || return nothing
    node_separation(s, res.θ) > 1e-6 || return nothing
    return canonical_octahedral(s, res.θ)
end

"""
    eliminate_octahedral(s, θ, n; log, cancel, rng) -> (s, θ)

Apply the largest elimination move that still fits, and repeat until none does.
"""
function eliminate_octahedral(s::OctahedralStructure, θ::Vector{Float64}, n::Integer;
                              log = nothing, cancel = nothing, rng = nothing)
    m = length(invariant_exponents(n))
    while true
        checkcancel(cancel)
        moves = filter(mv -> nunknowns(mv.structure) >= m - 2, octahedral_elimination_moves(s, θ))
        rng === nothing || sort!(moves; by = mv -> (-mv.removed, rand(rng)))
        accepted = nothing
        batch = max(Threads.nthreads(), 1)
        for lo in 1:batch:length(moves)
            idx = lo:min(lo + batch - 1, length(moves))
            results = Vector{Any}(nothing, length(idx))
            Threads.@threads :static for j in eachindex(idx)
                move = moves[idx[j]]          # a name used nowhere else in this function
                results[j] = fit_octahedral(move.structure, move.θ0, n)
            end
            k = findfirst(!isnothing, results)
            if k !== nothing
                accepted = (moves[idx[k]], results[k])
                break
            end
        end
        accepted === nothing && return s, θ
        mv, θn = accepted
        log === nothing || println(log, "    eliminate: ", mv.what, " -> ", npoints(mv.structure),
                                   " points, ", nunknowns(mv.structure), " unknowns")
        s, θ = mv.structure, θn
    end
end

"""
    grow_octahedral(s, θ, n; ntries, rng_seed, cancel) -> (s, θ) or nothing

Add orbits to a lower-degree rule until the system has enough unknowns for degree `n`, then
fit.

The added orbits start at half the rule's mean weight, not at a token fraction of it. The
simplex version can afford a light touch because it grows one degree at a time; sphere
degrees step by two, so the rule below misses the top invariants by O(0.3), and an orbit
carrying 2% of a weight cannot absorb that — the fit stalls around 1e-3 and the new orbit
slides into an existing one instead. At half the mean weight it converges.
"""
function grow_octahedral(s::OctahedralStructure, θ::Vector{Float64}, n::Integer;
                         ntries::Int = 16, rng_seed::Integer = 0x9e0, cancel = nothing)
    for (k, added) in enumerate(grow_recipes(s, n))
        checkcancel(cancel)
        got = _grow_octahedral_try(s, θ, n, added, ntries, rng_seed + k)
        got === nothing || return got
    end
    return nothing
end

"""
    grow_recipes(s, n)

Sets of orbits to add to `s` so that the system has at least as many unknowns as degree `n`
has invariant equations, cheapest in added points first. A square system is tried before an
underdetermined one: it is where the minimal rules live, and it is what adding a single
24-point orbit to the rule two degrees down usually gives.
"""
function grow_recipes(s::OctahedralStructure, n::Integer)
    m = length(invariant_exponents(n))
    present = Set(o.kind for o in s.orbits)
    singles = [k for k in (:a1, :a2, :a3, :b, :c, :d) if !(k in (:a1, :a2, :a3) && k in present)]
    cands = [[k] for k in singles]
    for k1 in singles, k2 in singles
        (k1 in (:a1, :a2, :a3) && k1 == k2) && continue    # each parameterless orbit is unique
        push!(cands, sort([k1, k2]))
    end
    unique!(cands)
    enough(c) = nunknowns(s) + sum(nunknowns(OctahedralOrbit(k)) for k in c) >= m
    ok = filter(enough, cands)
    return sort!(ok; by = c -> (sum(orbit_size(OctahedralOrbit(k)) for k in c), length(c), c))
end

function _grow_octahedral_try(s, θ, n, added, ntries, rng_seed)
    blocks = octahedral_blocks(s, θ)
    wscale = 0.5 * 4 * pi / npoints(s)
    results = Vector{Any}(nothing, ntries)
    Threads.@threads :static for t in 1:ntries
        rng = start_rng(rng_seed, t)
        new = Tuple{Symbol,Vector{Float64}}[]
        for kind in added
            push!(new, (kind, random_octahedral_parameters(rng, OctahedralStructure([kind]), wscale)))
        end
        st, θ0 = from_octahedral_blocks(vcat(blocks, new))
        θn = fit_octahedral(st, θ0, n)
        results[t] = θn === nothing ? nothing : (st, θn)
    end
    k = findfirst(!isnothing, results)
    return k === nothing ? nothing : results[k]
end

"""
    grow_and_eliminate_octahedral(s, θ, n; chains, rng_seed, cancel) -> (s, θ, counts) or nothing

Run `chains` independent grow-then-eliminate chains from the lower-degree rule `(s, θ)` and
keep the smallest result, ties broken by the larger minimum weight. Deterministic for a
given `rng_seed`.
"""
function grow_and_eliminate_octahedral(s::OctahedralStructure, θ::Vector{Float64}, n::Integer;
                                       chains::Int = 16, rng_seed::Integer = 0xe11, cancel = nothing)
    best = nothing
    counts = Int[]
    for c in 1:chains
        checkcancel(cancel)
        g = grow_octahedral(s, θ, n; ntries = 16, rng_seed = rng_seed + 1000c, cancel)
        g === nothing && continue
        se, θe = eliminate_octahedral(g[1], g[2], n; rng = start_rng(rng_seed, c), cancel)
        push!(counts, npoints(se))
        key = (npoints(se), -octahedral_margins(se, θe)[1])
        if best === nothing || key < best[3]
            best = (se, θe, key)
        end
    end
    return best === nothing ? nothing : (best[1], best[2], counts)
end
