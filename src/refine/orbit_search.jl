# Seeds for fully symmetric rules on domains whose orbits are not those of the box: the
# pyramid and the wedge. The search is that of refine/box_search.jl (multistart over the orbit
# structures of a point count at low degree, grow → eliminate above), written once against a
# small interface that each structure type implements:
#
#     orbit_types(T), make_structure(T, orbits), general_orbit(T), structure_volume(T)
#     n_moment_equations(T, n), orbit_system(s, n, S), structure_margins(s, θ)
#     canonicalize(s, θ), random_params(o, rng), loosely_inside(o, v), structure_moves(s, θ)
#
# `structure_moves` returns (structure, starting parameters, points removed, score) tuples.

"Fit `(s, θ0)` onto the degree-`n` moment variety; the canonical `(s, θ)` if it is a positive interior rule, else `nothing`."
function orbit_fit(s, θ0::Vector{Float64}, n::Integer; tol::Float64 = 1e-13, maxiter::Int = 400)
    sys = orbit_system(s, n, Float64)
    inside(θ) = all(o_off -> loosely_inside(o_off[1], θ[(o_off[2] + 2):(o_off[2] + nunknowns(o_off[1]))]),
                    zip(s.orbits, param_offsets(s)))
    θ, nr = levenberg_marquardt(sys, θ0; maxiter, tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    sc, θc = canonicalize(s, res.θ)
    wmin, cmin, gmin = structure_margins(sc, θc)
    (wmin > 0 && cmin > 1e-6 && gmin > 1e-6) || return nothing
    return sc, θc
end

function random_orbit_start(s, rng)
    V, N = structure_volume(typeof(s)), npoints(s)
    θ = Float64[]
    for o in s.orbits
        push!(θ, V / N * (0.5 + rand(rng)))
        append!(θ, random_params(o, rng))
    end
    return θ
end

"""
    orbit_structures(T, npts, n; slack = 2)

Every orbit structure of type `T` with exactly `npts` points and between `m` and `m + slack`
unknowns, `m` the number of degree-`n` moment equations.
"""
function orbit_structures(::Type{T}, npts::Integer, n::Integer; slack::Integer = 2) where {T}
    types = orbit_types(T)
    m = n_moment_equations(T, n)
    out = T[]
    function rec(i, left, counts)
        if i > length(types)
            left == 0 || return
            s = make_structure(T, reduce(vcat, [fill(types[k], counts[k]) for k in eachindex(types)]; init = eltype(types)[]))
            m <= nunknowns(s) <= m + slack && push!(out, s)
            return
        end
        t = types[i]
        for c in 0:min(left ÷ orbit_size(t), nparams(t) == 0 ? 1 : typemax(Int))     # a fixed orbit only once
            # unknowns cannot fall: stop once a count alone exceeds the window
            c * nunknowns(t) > m + slack && break
            rec(i + 1, left - c * orbit_size(t), vcat(counts, c))
        end
    end
    rec(1, Int(npts), Int[])
    return out
end

"""
    orbit_multistart(T, n; npts, nstarts = 64, slack = 2, rng_seed) -> (s, θ) or nothing

The rule of type `T` and degree `n` with the fewest points in `npts`: for each count in turn,
`nstarts` random starts on every orbit structure with that many points. Among the rules found
at the first count that has any, the one whose nodes keep farthest from the boundary.
"""
function orbit_multistart(::Type{T}, n::Integer; npts, nstarts::Integer = 64, slack::Integer = 2,
                          rng_seed::Integer = 0x5e7, cancel = nothing) where {T}
    for N in npts
        found = Any[]
        for s in orbit_structures(T, N, n; slack)
            results = Vector{Any}(nothing, nstarts)
            Threads.@threads :static for t in 1:nstarts
                checkcancel(cancel)
                results[t] = orbit_fit(s, random_orbit_start(s, start_rng(rng_seed + N, t)), n)
            end
            append!(found, filter(!isnothing, results))
        end
        isempty(found) && continue
        return found[argmax([structure_margins(f...)[2] for f in found])]
    end
    return nothing
end

"""
    orbit_eliminate(s, θ, n; rng) -> (s, θ)

Greedy node elimination on a valid degree-`n` rule: try the moves of `structure_moves` that
remove the most points first, accept the first whose refit is valid, and repeat until none is.
With `rng`, moves removing the same number of points are tried in random order.
"""
function orbit_eliminate(s, θ::Vector{Float64}, n::Integer; rng = nothing, cancel = nothing)
    m = n_moment_equations(typeof(s), n)
    while true
        checkcancel(cancel)
        moves = filter(mv -> nunknowns(mv[1]) >= m, structure_moves(s, θ))
        sort!(moves; by = mv -> (-mv[3], rng === nothing ? mv[4] : rand(rng)))
        accepted = nothing
        batch = max(Threads.nthreads(), 1)
        for lo in 1:batch:length(moves)
            idx = lo:min(lo + batch - 1, length(moves))
            results = Vector{Any}(nothing, length(idx))
            Threads.@threads :static for j in eachindex(idx)
                mv = moves[idx[j]]
                results[j] = orbit_fit(mv[1], mv[2], n)
            end
            k = findfirst(!isnothing, results)
            k === nothing || (accepted = results[k]; break)
        end
        accepted === nothing && return s, θ
        s, θ = accepted
    end
end

"""
    orbit_grow(s, θ, n; ntries = 32, rng_seed) -> (s, θ) or nothing

From a rule of lower degree, add orbits with small weights until the unknowns exceed the
degree-`n` equations, and refit. Recipes for the added orbits are tried in turn, as in the
simplex search ([`grow`](@ref)): general orbits to a modest and then a larger surplus; one
of each smaller type besides; and one smaller orbit with a single unknown to spare, for a
tight rule below, whose fits otherwise land outside or on negative weights.
"""
function orbit_grow(s, θ::Vector{Float64}, n::Integer; ntries::Int = 32, rng_seed::Integer = 0x96)
    T = typeof(s)
    m = n_moment_equations(T, n)
    general = general_orbit(T)
    small = [t for t in orbit_types(T) if t != general && nparams(t) > 0]    # never a second fixed point
    recipes = vcat([(empty(small), m + max(2, m ÷ 4)), (empty(small), m + max(4, m ÷ 2)), (small, m + max(2, m ÷ 4))],
                   [([t], m + 1) for t in sort(small; by = orbit_size)])
    blocks = [(o, θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
    for (k, (extra, target)) in enumerate(recipes)
        added = copy(extra)
        unknowns = nunknowns(s) + sum(nunknowns, added; init = 0)
        while unknowns < target
            push!(added, general)
            unknowns += nunknowns(general)
        end
        results = Vector{Any}(nothing, ntries)
        Threads.@threads :static for t in 1:ntries
            rng = start_rng(rng_seed + k, t)
            new = [(o, vcat(0.02 * structure_volume(T) / orbit_size(o), random_params(o, rng))) for o in added]
            all_blocks = vcat(blocks, new)
            st, θ0 = canonicalize(make_structure(T, [b[1] for b in all_blocks]),
                                  reduce(vcat, (b[2] for b in all_blocks); init = Float64[]))
            results[t] = orbit_fit(st, θ0, n)
        end
        i = findfirst(!isnothing, results)
        i === nothing || return results[i]
    end
    return nothing
end

"""
    orbit_grow_and_eliminate(s, θ, n; chains = 8, rng_seed) -> (s, θ, points per chain) or nothing

Independent grow → eliminate chains from the lower-degree rule `(s, θ)`; the rule with the
fewest points, ties broken by the largest distance of the nodes from the boundary.
"""
function orbit_grow_and_eliminate(s, θ::Vector{Float64}, n::Integer; chains::Int = 8,
                                  rng_seed::Integer = 0xe11, cancel = nothing)
    best, counts = nothing, Int[]
    for c in 1:chains
        checkcancel(cancel)
        g = orbit_grow(s, θ, n; rng_seed = rng_seed + 1000c)
        g === nothing && continue
        se, θe = orbit_eliminate(g[1], g[2], n; rng = start_rng(rng_seed, c), cancel)
        push!(counts, npoints(se))
        key = (npoints(se), -structure_margins(se, θe)[2])
        (best === nothing || key < best[3]) && (best = (se, θe, key))
    end
    return best === nothing ? nothing : (best[1], best[2], counts)
end

# --- the pyramid ---------------------------------------------------------------------------------

orbit_types(::Type{PyramidStructure}) = pyramid_orbit_types()
make_structure(::Type{PyramidStructure}, orbits) = PyramidStructure(collect(PyramidOrbit, orbits))
general_orbit(::Type{PyramidStructure}) = PyramidOrbit([1, 1])
structure_volume(::Type{PyramidStructure}) = 4 / 3
n_moment_equations(::Type{PyramidStructure}, n) = length(pyramid_terms(n))
orbit_system(s::PyramidStructure, n, ::Type{S}) where {S} = PyramidMomentSystem(s, n, S)
structure_margins(s::PyramidStructure, θ) = pyramid_margins(s, θ)
random_params(o::PyramidOrbit, rng) = vcat(0.03 .+ 0.94 .* rand(rng, nparams(o.base)), 0.03 + 0.94 * rand(rng))
loosely_inside(o::PyramidOrbit, v) = all(x -> abs(x) < 1.05, v[1:(end - 1)]) && -0.05 < v[end] < 1.05

"Drop, merge and zero moves on a pyramid rule: the square orbit shrinks, the height stays."
function structure_moves(s::PyramidStructure, θ::Vector{Float64})
    blocks = [(o, θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
    rebuild(bs) = canonicalize(PyramidStructure([b[1] for b in bs]), reduce(vcat, (b[2] for b in bs); init = Float64[]))
    moves = Any[]
    for (i, (o, p)) in enumerate(blocks)
        mass = p[1] * orbit_size(o)
        v, z = p[2:(end - 1)], p[end]
        replaced(no, nv) = rebuild([j == i ? (no, vcat(mass / orbit_size(no), nv, z)) : blocks[j] for j in eachindex(blocks)])
        rest = [(b[1], copy(b[2])) for (j, b) in enumerate(blocks) if j != i]
        if !isempty(rest)
            total = sum(b[2][1] * orbit_size(b[1]) for b in rest)
            foreach(b -> b[2][1] *= (total + mass) / total, rest)
            push!(moves, (rebuild(rest)..., orbit_size(o), mass))
        end
        mult = o.base.mult
        for k in 1:(length(v) - 1), l in (k + 1):length(v)            # merge two values
            nm, nv = copy(mult), copy(v)
            nv[k] = (mult[k] * v[k] + mult[l] * v[l]) / (mult[k] + mult[l])
            nm[k] += nm[l]
            deleteat!(nm, l)
            deleteat!(nv, l)
            order = sortperm(nm; rev = true)
            no = PyramidOrbit(nm[order])
            push!(moves, (replaced(no, nv[order])..., orbit_size(o) - orbit_size(no), abs(v[k] - v[l])))
        end
        for k in eachindex(v)                                        # send a value to zero
            nm = deleteat!(copy(mult), k)
            no = PyramidOrbit(nm)
            push!(moves, (replaced(no, deleteat!(copy(v), k))..., orbit_size(o) - orbit_size(no), abs(v[k])))
        end
    end
    return moves
end

# --- the wedge -----------------------------------------------------------------------------------

orbit_types(::Type{WedgeStructure}) = wedge_orbit_types()
make_structure(::Type{WedgeStructure}, orbits) = WedgeStructure(collect(WedgeOrbit, orbits))
general_orbit(::Type{WedgeStructure}) = WedgeOrbit([1, 1, 1], true)
structure_volume(::Type{WedgeStructure}) = 1
n_moment_equations(::Type{WedgeStructure}, n) = n_wedge_equations(n)
orbit_system(s::WedgeStructure, n, ::Type{S}) where {S} = WedgeMomentSystem(s, n, S)
structure_margins(s::WedgeStructure, θ) = wedge_margins(s, θ)
function random_params(o::WedgeOrbit, rng)
    r = length(o.pattern.mult)
    u = -log.(rand(rng, r))                                   # Dirichlet(1, …, 1), as on the triangle
    u ./= sum(u)
    v = [u[i] / o.pattern.mult[i] for i in 1:(r - 1)]
    return o.mirrored ? vcat(v, 0.03 + 0.94 * rand(rng)) : v
end
function loosely_inside(o::WedgeOrbit, v)
    nv = ncoords(o.pattern)
    all(>(-0.05), pattern_values(o.pattern, v[1:nv])) || return false
    return !o.mirrored || -0.05 < v[end] < 1.05
end

"Drop, merge and flatten moves on a wedge rule: the triangle orbit shrinks, or the pair ±z joins at z = 0."
function structure_moves(s::WedgeStructure, θ::Vector{Float64})
    blocks = [(o, θ[(off + 1):(off + nunknowns(o))]) for (o, off) in zip(s.orbits, param_offsets(s))]
    rebuild(bs) = canonicalize(WedgeStructure([b[1] for b in bs]), reduce(vcat, (b[2] for b in bs); init = Float64[]))
    moves = Any[]
    for (i, (o, p)) in enumerate(blocks)
        mass = p[1] * orbit_size(o)
        nv = ncoords(o.pattern)
        vals = pattern_values(o.pattern, p[2:(1 + nv)])
        z = o.mirrored ? [p[end]] : Float64[]
        replaced(no, nv_, nz) = rebuild([j == i ? (no, vcat(mass / orbit_size(no), nv_, nz)) : blocks[j]
                                         for j in eachindex(blocks)])
        rest = [(b[1], copy(b[2])) for (j, b) in enumerate(blocks) if j != i]
        if !isempty(rest)                                            # drop the orbit
            total = sum(b[2][1] * orbit_size(b[1]) for b in rest)
            foreach(b -> b[2][1] *= (total + mass) / total, rest)
            push!(moves, (rebuild(rest)..., orbit_size(o), mass))
        end
        mult = o.pattern.mult
        if mult == [1, 1, 1]                                         # two values meet: (m, m, 1 − 2m)
            for k in 1:2, l in (k + 1):3
                no = WedgeOrbit([2, 1], o.mirrored)
                push!(moves, (replaced(no, [(vals[k] + vals[l]) / 2], z)..., orbit_size(o) - orbit_size(no),
                              abs(vals[k] - vals[l])))
            end
        elseif mult == [2, 1]                                        # to the centroid
            no = WedgeOrbit([3], o.mirrored)
            push!(moves, (replaced(no, Float64[], z)..., orbit_size(o) - orbit_size(no), abs(vals[1] - 1 / 3)))
        end
        if o.mirrored                                                # the pair ±z joins at z = 0
            no = WedgeOrbit(o.pattern, false)
            push!(moves, (replaced(no, vals[1:nv], Float64[])..., orbit_size(o) - orbit_size(no), abs(z[1])))
        end
    end
    return moves
end
