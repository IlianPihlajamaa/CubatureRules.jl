# Invariance of a rule under the symmetry group it claims, by orbits, for every group alike.
# Two nodes lie in the same orbit exactly when their canonical keys agree (the coordinates the
# group permutes, sorted, and those it reflects, made positive), so grouping the nodes by key
# finds the candidate orbits in O(N log N), and the rule is invariant when the nodes of each
# key carry every weight on whole orbits: each weight on `k` copies of the key's orbit (most
# rules have `k = 1` and one weight per orbit).

"""
    orbit_key(g, x, tol) -> (key, size) or nothing

The canonical form of the node `x` under the group `g`, and the size of its orbit, coordinates
within `tol` of each other counted as equal; `nothing` when `g` does not act on points of
`x`'s dimension.
"""
orbit_key(::SymmetryGroup, x, tol) = nothing
orbit_key(::Reflection, x, tol) = length(x) == 1 ? ([abs(x[1])], abs(x[1]) > tol ? 2 : 1) : nothing
function orbit_key(::SymmetricGroup{N}, x, tol) where {N}
    length(x) == N - 1 || return nothing
    λ = sort!(vcat(1 - sum(x), collect(x)); rev = true)
    return λ, permutation_count(λ, tol)
end
function orbit_key(g::Union{Oh,D4,B4}, x, tol)
    length(x) == signed_dimension(g) || return nothing
    k = sort!(abs.(collect(x)); rev = true)
    return k, permutation_count(k, tol) * 2^count(>(tol), k)
end
function orbit_key(::C4v, x, tol)                    # the square's group on (x, y), the height fixed
    length(x) == 3 || return nothing
    k = sort!(abs.([x[1], x[2]]); rev = true)
    return vcat(k, x[3]), permutation_count(k, tol) * 2^count(>(tol), k)
end
function orbit_key(::D3h, x, tol)                    # the triangle's group on (x, y), times z ↦ −z
    length(x) == 3 || return nothing
    λ = sort!([1 - x[1] - x[2], x[1], x[2]]; rev = true)
    return vcat(λ, abs(x[3])), permutation_count(λ, tol) * (abs(x[3]) > tol ? 2 : 1)
end
signed_dimension(::Oh) = 3
signed_dimension(::D4) = 2
signed_dimension(::B4) = 4

# a node as a vector of its coordinates: a number in one dimension (where `collect` would
# give a zero-dimensional array, which Julia 1.11 cannot sort)
node_coordinates(x::Number) = [x]
node_coordinates(x) = collect(x)

"`n!/∏ mⱼ!` for the runs of equal entries (within `tol`) of the sorted vector `k`: its distinct permutations."
function permutation_count(k, tol)
    total, run = factorial(length(k)), 1
    for i in 2:length(k)
        if abs(k[i] - k[i - 1]) <= tol
            run += 1
        else
            total ÷= factorial(run)
            run = 1
        end
    end
    return total ÷ factorial(run)
end

"""
    group_orbits(g, xs, ws, tol) -> (keys, totals) or nothing

The nodes `xs` grouped into the orbits of `g`: each orbit's key and total weight, or
`nothing` if the nodes of some key do not carry each weight on whole orbits — the rule is
then not invariant under `g` to tolerance `tol` (or its nodes are too close to tell), or `g`
does not act on its domain.

Keys within `tol` of each other have first coordinates within `tol`, so the nodes are sorted by
the first coordinate of their key and grouped within each window of such first coordinates:
sorting whole keys lexicographically would interleave two orbits whose first coordinates
agree to rounding, as `(½, ½, 0)` and `(½, ¼, ¼)` do on a triangle.
"""
function group_orbits(g::SymmetryGroup, xs, ws, tol)
    kss = [orbit_key(g, x, tol) for x in xs]
    any(isnothing, kss) && return nothing
    keys = [k[1] for k in kss]
    order = sortperm(keys; by = first)
    wscale = isempty(ws) ? zero(eltype(ws)) : maximum(abs, ws)
    reps, totals = similar(keys, 0), similar(ws, 0)
    taken = falses(length(xs))
    for a in eachindex(order)
        i = order[a]
        taken[i] && continue
        grp = [i]
        taken[i] = true
        b = a + 1
        while b <= length(order) && first(keys[order[b]]) - first(keys[i]) <= tol
            j = order[b]
            if !taken[j] && maximum(abs, keys[j] .- keys[i]) <= tol
                push!(grp, j)
                taken[j] = true
            end
            b += 1
        end
        complete_orbits(xs, ws, grp, kss[i][2], tol, wscale) || return nothing
        push!(reps, keys[i])
        push!(totals, sum(ws[q] for q in grp))
    end
    return reps, totals
end

# Whether the nodes `grp`, all with one key, carry every weight on whole orbits of size
# `size`: the nodes of each weight are `k` copies of the orbit, each point `k` times. A rule
# may place one orbit at several weights, or repeat a point (Grundmann–Möller does both),
# and is still invariant.
function complete_orbits(xs, ws, grp, size, tol, wscale)
    byweight = sort(grp; by = q -> ws[q])
    a = 1
    while a <= length(byweight)
        b = a
        while b < length(byweight) && abs(ws[byweight[b + 1]] - ws[byweight[a]]) <= tol * wscale
            b += 1
        end
        pts = sort([node_coordinates(xs[q]) for q in byweight[a:b]])
        runs = Int[1]
        for m in 2:length(pts)
            maximum(abs, pts[m] .- pts[m - 1]) <= tol ? (runs[end] += 1) : push!(runs, 1)
        end
        length(runs) == size && all(==(runs[1]), runs) || return false
        a = b + 1
    end
    return true
end

"""
    check_invariance(g, xs, ws, tol) -> Bool or nothing

Whether the nodes `xs` and weights `ws` are invariant under the group `g`, to `tol` in the
coordinates and relative to the largest weight in the weights (see [`group_orbits`](@ref));
`false` also when `g` does not act on the nodes' dimension, `nothing` for [`NoSymmetry`](@ref).
"""
check_invariance(::NoSymmetry, xs, ws, tol) = nothing
check_invariance(g::SymmetryGroup, xs, ws, tol) = group_orbits(g, xs, ws, tol) !== nothing
