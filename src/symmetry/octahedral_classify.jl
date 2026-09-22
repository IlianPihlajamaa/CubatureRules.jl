# Recovering an O_h orbit structure from a bare point set (v0.4).
#
# The rest of the sphere machinery works in orbit-parameter space: a structure plus one
# weight and up to two coordinates per orbit. Anything arriving from outside — a table from
# another package, a file of numbers someone generated years ago — arrives instead as a list
# of points and weights. This turns the second into the first, which is what makes such a
# rule refinable, comparable and describable rather than just usable.
#
# Classification is by the shape of the sorted absolute coordinates, which is exactly what
# the orbit types are:
#
#     (1,0,0)      a₁        (t,t,0)  with 2t² = 1   a₂        (t,t,t)  a₃
#     (l,l,m)      b         (p,q,0)                 c         (r,s,u)  d
#
# Points are then grouped by their parameters, and the group sizes must match the orbit
# sizes exactly — 6, 12, 8, 24, 24, 48. A rule that does not decompose this way is not
# O_h-symmetric in the sense the rest of the code means, and saying so is more useful than
# guessing.

"""
    octahedral_orbit_of(x; tol) -> (kind, parameters)

The orbit type of a single point on the unit sphere, and the parameters that identify its
orbit: `Float64[]` for the three parameterless types, `[l]` for `b`, `[p]` for `c` (the
larger of the two non-zero coordinates), and `[r, s]` for `d` (sorted ascending).
"""
function octahedral_orbit_of(x, tol::Float64 = 1e-8)
    v = sort!(abs.(Float64.(collect(x))))          # ascending: v1 ≤ v2 ≤ v3
    iszero3 = v[1] <= tol
    iszero2 = v[2] <= tol
    eq12 = abs(v[1] - v[2]) <= tol
    eq23 = abs(v[2] - v[3]) <= tol
    if iszero2                                      # two zeros: (0,0,1)
        return :a1, Float64[]
    elseif iszero3                                  # one zero
        return eq23 ? (:a2, Float64[]) : (:c, [v[3]])
    elseif eq12 && eq23
        return :a3, Float64[]
    elseif eq12 || eq23
        # (l,l,m): the repeated value is l, whichever end it sits at
        return :b, [eq12 ? v[1] : v[2]]
    else
        return :d, [v[1], v[2]]
    end
end

"""
    classify_octahedral(nodes, weights; tol) -> (structure, θ) or nothing

Recover the `O_h` orbit structure and parameters `θ = [w₁, p₁…, w₂, p₂…, …]` of a rule given
only its nodes and weights. `nothing` when the points do not decompose into whole orbits
with a single weight each — that is, when the rule is not `O_h`-symmetric to `tol`.

The weights returned are per point, as everywhere else in this package.
"""
function classify_octahedral(nodes, weights; tol::Float64 = 1e-8)
    length(nodes) == length(weights) || throw(ArgumentError("nodes and weights differ in length"))
    # Group by proximity, not by rounding onto a grid: two points of one orbit can land on
    # opposite sides of a grid line and split it in two. There are at most a few dozen
    # orbits, so a scan over the open groups costs nothing.
    kinds = Symbol[]
    params = Vector{Vector{Float64}}()
    members = Vector{Vector{Int}}()
    for i in eachindex(nodes)
        length(nodes[i]) == 3 || return nothing
        kind, p = octahedral_orbit_of(nodes[i], tol)
        g = findfirst(j -> kinds[j] === kind && maximum(abs, params[j] - p; init = 0.0) <= tol,
                      eachindex(kinds))
        if g === nothing
            push!(kinds, kind); push!(params, p); push!(members, [i])
        else
            push!(members[g], i)
        end
    end
    order = sortperm(collect(zip(string.(kinds), [isempty(p) ? 0.0 : p[1] for p in params])))
    orbits = OctahedralOrbit[]
    blocks = Vector{Vector{Float64}}()
    for g in order
        kind, idx = kinds[g], members[g]
        o = OctahedralOrbit(kind)
        length(idx) == orbit_size(o) || return nothing        # not a whole orbit
        w = Float64.(weights[idx])
        (maximum(w) - minimum(w)) <= tol * max(1.0, maximum(abs, w)) || return nothing
        # average the parameters over the orbit rather than trusting one representative
        ps = [octahedral_orbit_of(nodes[i], tol)[2] for i in idx]
        p = isempty(first(ps)) ? Float64[] : [sum(q[j] for q in ps) / length(ps) for j in eachindex(first(ps))]
        push!(orbits, o)
        push!(blocks, vcat(sum(w) / length(w), p))
    end
    isempty(orbits) && return nothing
    s = OctahedralStructure(orbits)
    return s, reduce(vcat, blocks; init = Float64[])
end
