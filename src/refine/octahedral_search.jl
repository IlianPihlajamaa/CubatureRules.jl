# Searching for O_h-symmetric sphere rules (PLAN §6 Tier 3, v0.4): the Lebedev counterpart
# of `refine/seeds.jl`, which does the same for simplices.
#
# The shape of the search is identical — enumerate orbit structures with the right point
# count, start Levenberg–Marquardt from random admissible parameters, polish with
# Gauss–Newton, keep what converges to a positive interior rule — and only the notion of
# "admissible" differs. On a simplex it is the interior of the simplex; here it is the set
# where no orbit has collapsed onto a smaller one:
#
#     b(l)    degenerates at l = 0 (→ a₁), 2l² = 1 (→ a₂), 3l² = 1 (→ a₃)
#     c(p)    degenerates at p = 0 (→ a₁) and 2p² = 1 (→ a₂)
#     d(r,s)  degenerates whenever a coordinate vanishes (→ c) or two coincide (→ b)
#
# A degenerate structure is not wrong, merely misdescribed: it is a rule with fewer points
# than its structure claims, which would then be counted wrongly. The search rejects those
# and lets the smaller structure be found on its own terms.

"How far the orbit's parameters sit from a degeneracy, in coordinate units."
function orbit_margin(o::OctahedralOrbit, p)
    k = o.kind
    k in (:a1, :a2, :a3) && return 1.0
    if k === :b
        l = p[1]
        m2 = 1 - 2l^2
        m2 > 0 || return -1.0
        m = sqrt(m2)
        return min(l, m, abs(l - m))
    elseif k === :c
        q = p[1]
        r2 = 1 - q^2
        r2 > 0 || return -1.0
        r = sqrt(r2)
        return min(q, r, abs(q - r))
    else
        r, s = p[1], p[2]
        u2 = 1 - r^2 - s^2
        u2 > 0 || return -1.0
        u = sqrt(u2)
        return min(r, s, u, abs(r - s), abs(s - u), abs(r - u))
    end
end

"""
    octahedral_margins(structure, θ) -> (min weight, min coordinate margin)

The two quantities that decide whether a solution is a rule of the shape it claims: the
smallest weight, and the smallest distance to an orbit degeneracy.
"""
function octahedral_margins(s::OctahedralStructure, θ::AbstractVector)
    offs = param_offsets(s)
    wmin = Inf
    dmin = Inf
    for (i, o) in enumerate(s.orbits)
        wmin = min(wmin, Float64(θ[offs[i] + 1]))
        dmin = min(dmin, orbit_margin(o, Float64.(θ[(offs[i] + 2):(offs[i] + nunknowns(o))])))
    end
    return wmin, dmin
end

"""
    node_separation(structure, θ)

Smallest distance between two distinct nodes of the expanded rule. Two orbits of the same
kind at nearly equal parameters pass [`octahedral_margins`](@ref) and still describe one
orbit twice, which this catches.
"""
function node_separation(s::OctahedralStructure, θ::AbstractVector)
    xs, _ = expand(s, Float64.(θ))
    d = Inf
    for i in eachindex(xs), j in (i + 1):length(xs)
        d = min(d, sqrt(sum(abs2, xs[i] - xs[j])))
    end
    return d
end

"Random admissible parameters for a structure, with weights around `wscale` per point."
function random_octahedral_parameters(rng, s::OctahedralStructure, wscale::Float64)
    θ = Float64[]
    for o in s.orbits
        push!(θ, wscale * (0.5 + rand(rng)))
        k = o.kind
        if k === :b
            # l ∈ (0, 1/√2), away from the three degeneracies
            for _ in 1:100
                l = 0.02 + 0.68 * rand(rng)
                orbit_margin(o, [l]) > 0.03 && (push!(θ, l); break)
            end
        elseif k === :c
            for _ in 1:100
                p = 0.02 + 0.96 * rand(rng)
                orbit_margin(o, [p]) > 0.03 && (push!(θ, p); break)
            end
        elseif k === :d
            for _ in 1:200
                r, t = rand(rng), rand(rng)
                orbit_margin(o, [r, t]) > 0.03 && (append!(θ, (r, t)); break)
            end
        end
    end
    length(θ) == nunknowns(s) || throw(RefinementError("Lebedev", "could not draw admissible parameters"))
    return θ
end

"""
    multistart_octahedral(structure, n; nstarts, rng_seed, tol, first_only, cancel)

Random restarts for a degree-`n` rule of the given structure, returning
`(θ, min weight, min margin, start index)` for each distinct solution found.
"""
function multistart_octahedral(s::OctahedralStructure, n::Integer; nstarts::Integer = 2000,
                               rng_seed::Integer = 0x5eed, tol::Float64 = 1e-13,
                               first_only::Bool = true, cancel = nothing)
    found = Vector{Tuple{Vector{Float64},Float64,Float64,Int}}()
    wscale = 4π / npoints(s)
    batch = 64
    for lo in 1:batch:nstarts
        checkcancel(cancel)
        idx = lo:min(lo + batch - 1, nstarts)
        results = Vector{Any}(nothing, length(idx))
        Threads.@threads :static for t in 1:Threads.nthreads()
            sys = OctahedralMomentSystem(s, n, Float64)
            for j in t:Threads.nthreads():length(idx)
                results[j] = octahedral_attempt(sys, s, start_rng(rng_seed, idx[j]), wscale, tol)
            end
        end
        for (j, res) in enumerate(results)
            res === nothing && continue
            θc, wmin, dmin = res
            any(f -> maximum(abs, f[1] - θc) < 1e-8, found) && continue
            push!(found, (θc, wmin, dmin, idx[j]))
        end
        first_only && !isempty(found) && break
    end
    return found
end

"One attempt: Levenberg–Marquardt from a random admissible point, polished by Gauss–Newton."
function octahedral_attempt(sys::OctahedralMomentSystem{Float64}, s::OctahedralStructure, rng,
                            wscale::Float64, tol::Float64)
    θ0 = random_octahedral_parameters(rng, s, wscale)
    inside(θ) = octahedral_margins(s, θ)[2] > 0.0
    θ, nr = levenberg_marquardt(sys, θ0; maxiter = 300, tol = tol, accept = inside)
    nr <= 1e-10 || return nothing
    res = gauss_newton(sys, θ; step_tol = 1e-15, res_floor = tol, rank_rtol = 1e-13, maxiter = 10)
    res.residual <= 1e-12 || return nothing
    wmin, dmin = octahedral_margins(s, res.θ)
    (wmin > 0 && dmin > 1e-6) || return nothing
    node_separation(s, res.θ) > 1e-6 || return nothing
    return canonical_octahedral(s, res.θ), wmin, dmin
end

"""
    canonical_octahedral(structure, θ)

`θ` with each orbit's parameters put in a canonical form, so that solutions differing only
by a relabelling compare equal: `d`-orbits get their two parameters sorted, and repeated
orbits of one kind are ordered by their first parameter.
"""
function canonical_octahedral(s::OctahedralStructure, θ::AbstractVector)
    offs = param_offsets(s)
    blocks = [collect(Float64.(θ[(offs[i] + 1):(offs[i] + nunknowns(o))])) for (i, o) in enumerate(s.orbits)]
    for (i, o) in enumerate(s.orbits)
        o.kind === :d && (blocks[i][2:3] = sort(blocks[i][2:3]))
    end
    # order repeats of the same kind by their leading parameter, keeping kinds in place
    for k in (:b, :c, :d)
        pos = findall(o -> o.kind === k, s.orbits)
        length(pos) > 1 && (blocks[pos] = blocks[pos[sortperm([b[2] for b in blocks[pos]])]])
    end
    return reduce(vcat, blocks; init = Float64[])
end

"""
    octahedral_candidate_structures(npts, n; max_unknowns)

Every `O_h` structure with exactly `npts` points and at least as many unknowns as there are
invariant equations for degree `n`, fewest unknowns first. The parameterless orbits `a₁`,
`a₂` and `a₃` can each appear at most once, since each is a single orbit.
"""
function octahedral_candidate_structures(npts::Integer, n::Integer; max_unknowns::Integer = typemax(Int))
    m = length(invariant_exponents(n))
    out = OctahedralStructure[]
    for na1 in 0:1, na2 in 0:1, na3 in 0:1
        fixed = 6na1 + 12na2 + 8na3
        rest = npts - fixed
        rest >= 0 || continue
        for nd in 0:(rest ÷ 48)
            left = rest - 48nd
            left % 24 == 0 || continue
            nbc = left ÷ 24
            for nb in 0:nbc
                nc = nbc - nb
                u = na1 + na2 + na3 + 2nb + 2nc + 3nd
                (m <= u <= max_unknowns) || continue
                kinds = Symbol[]
                na1 == 1 && push!(kinds, :a1)
                na2 == 1 && push!(kinds, :a2)
                na3 == 1 && push!(kinds, :a3)
                append!(kinds, fill(:b, nb))
                append!(kinds, fill(:c, nc))
                append!(kinds, fill(:d, nd))
                push!(out, OctahedralStructure(kinds))
            end
        end
    end
    return sort!(out; by = s -> (nunknowns(s), [o.kind for o in s.orbits]))
end
