# Fully symmetric rules on the wedge (PLAN §6 Tier 3): the reference wedge, the reference
# triangle times [-1, 1], is invariant under the 6 permutations of the triangle's barycentric
# coordinates (affine maps of the triangle onto itself) and the reflection z ↦ −z; 12 in all,
# the group D₃ₕ of the prism with an equilateral base.
#
# An orbit is a triangle orbit (symmetry/orbits.jl) either at z = 0 or at the pair ±z, so its
# parameters are the triangle orbit's free barycentric values and, for a pair, z in (0, 1).
#
# The group is a product acting on separate coordinates, so its invariant polynomials are
# the products of the triangle's invariants (the columns of `invariant_basis(3, n)` over the
# orthonormal Dubiner basis) and the even polynomials in z. With orthonormal even Legendre
# polynomials in z the products are orthonormal on the wedge, and only the constant has a
# nonzero integral.
#
# `Q` is computed in Float64, so its columns are invariant only to rounding; each orbit sums
# them over all of its triangle images, which makes the equations exactly those of invariant
# polynomials (their Reynolds averages) at any precision.

"""
    WedgeOrbit(pattern, mirrored)

An orbit type of the wedge's symmetry group: the triangle orbit `pattern` (an
[`OrbitPattern`](@ref) of `S₃`) at `z = 0`, or, if `mirrored`, at the pair `±z`.
"""
struct WedgeOrbit
    pattern::OrbitPattern
    mirrored::Bool
end
WedgeOrbit(mult::AbstractVector{<:Integer}, mirrored::Bool) = WedgeOrbit(OrbitPattern(mult, 3), mirrored)
Base.:(==)(a::WedgeOrbit, b::WedgeOrbit) = a.pattern == b.pattern && a.mirrored == b.mirrored
Base.hash(o::WedgeOrbit, h::UInt) = hash(o.pattern, hash(o.mirrored, hash(:WedgeOrbit, h)))
nparams(o::WedgeOrbit) = ncoords(o.pattern) + o.mirrored       # barycentric values, then z
nunknowns(o::WedgeOrbit) = 1 + nparams(o)
orbit_size(o::WedgeOrbit) = orbit_size(o.pattern) * (o.mirrored ? 2 : 1)
wedge_orbit_types() = [WedgeOrbit(m, z) for m in ([3], [2, 1], [1, 1, 1]) for z in (false, true)]

"""
    WedgeStructure(orbits)

A fully symmetric rule on the reference wedge as a list of [`WedgeOrbit`](@ref)s. Its
parameters hold, per orbit, the weight of each point, the free barycentric values of the
triangle orbit and, for a mirrored orbit, the height `z > 0`.
"""
struct WedgeStructure
    orbits::Vector{WedgeOrbit}
end
nunknowns(s::WedgeStructure) = sum(nunknowns, s.orbits; init = 0)
npoints(s::WedgeStructure) = sum(orbit_size, s.orbits; init = 0)
function param_offsets(s::WedgeStructure)
    offs, k = Int[], 0
    for o in s.orbits
        push!(offs, k)
        k += nunknowns(o)
    end
    return offs
end
"The orbits as stored in a table: the triangle pattern's multiplicities, then 0 (at `z = 0`) or 2 (the pair `±z`)."
wedge_structure_lists(s::WedgeStructure) = [vcat(o.pattern.mult, o.mirrored ? 2 : 0) for o in s.orbits]
WedgeStructure(lists::AbstractVector{<:AbstractVector{<:Integer}}) =
    WedgeStructure([WedgeOrbit(l[1:(end - 1)], l[end] == 2) for l in lists])
Base.show(io::IO, s::WedgeStructure) = print(io, "WedgeStructure(", wedge_structure_lists(s), ")")

"""
    expand(s::WedgeStructure, θ) -> (nodes, weights)

The nodes (Cartesian, on the reference wedge) and weights of the rule `(s, θ)`.
"""
function expand(s::WedgeStructure, θ::AbstractVector{S}) where {S}
    nodes, ws = Vector{Vector{S}}(), S[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        nv = ncoords(o.pattern)
        vals = pattern_values(o.pattern, θ[(off + 2):(off + 1 + nv)])
        zs = o.mirrored ? (θ[off + 2 + nv], -θ[off + 2 + nv]) : (zero(S),)
        for z in zs, lab in o.pattern.labels
            push!(nodes, [vals[lab[2]], vals[lab[3]], z])
            push!(ws, θ[off + 1])
        end
    end
    return nodes, ws
end

"""
    WedgeMomentSystem(structure, n, S)

The fully symmetric moment system for exactness to degree `n` on the reference wedge, in the
orthonormal invariant polynomials described at the top of this file: for each degree `k` of
the triangle's invariants and each even degree `2m ≤ n − k` in z. Callable: `sys(θ) -> (r, J)`,
or `sys(θ; jacobian = false)`.
"""
struct WedgeMomentSystem{S,B<:SimplexBasis{2,S}}
    structure::WedgeStructure
    n::Int
    basis::B
    blocks::Vector{Tuple{UnitRange{Int},UnitRange{Int},Matrix{S}}}   # Q by degree, as in SymmetricMomentSystem
    rows::Vector{Tuple{Int,Int,Int}}       # equation → (triangle invariant c, its degree k, m)
    mass::S                                # ∫ of the constant invariant, (Q[1, 1] / √2) · √2
end
function WedgeMomentSystem(s::WedgeStructure, n::Integer, ::Type{S}) where {S}
    Q = invariant_basis(3, n)
    b = SimplexBasis{2,S}(Int(n))
    blocks = invariant_blocks(Q, b, S)
    rows = Tuple{Int,Int,Int}[]
    for (k, rk) in enumerate(Q.ranks), c in (sum(Q.ranks[1:(k - 1)]; init = 0) + 1):(sum(Q.ranks[1:k]))
        for m in 0:((n - (k - 1)) ÷ 2)
            push!(rows, (c, k - 1, m))
        end
    end
    sort!(rows; by = t -> (t[2] + 2t[3], t[1], t[3]))
    mass = S(Q.Q[1, 1]) * dubiner_mass(S) * sqrt(S(2))
    return WedgeMomentSystem{S,typeof(b)}(s, Int(n), b, blocks, rows, mass)
end
n_wedge_equations(n::Integer) = sum(rk * ((n - k) ÷ 2 + 1) for (k, rk) in zip(0:n, invariant_basis(3, n).ranks))
n_equations(sys::WedgeMomentSystem) = length(sys.rows)
n_unknowns(sys::WedgeMomentSystem) = nunknowns(sys.structure)

function (sys::WedgeMomentSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    s = sys.structure
    L = basis_length(sys.basis)
    mT = isempty(sys.blocks) ? 0 : last(sys.blocks[end][2])
    M = sys.n ÷ 2 + 1
    r = zeros(S, length(sys.rows))
    J = jacobian ? zeros(S, length(sys.rows), nunknowns(s)) : zeros(S, 0, 0)
    x = Vector{S}(undef, 2)
    for (o, off) in zip(s.orbits, param_offsets(s))
        w = S(θ[off + 1])
        nv = ncoords(o.pattern)
        rr = length(o.pattern.mult)
        vals = pattern_values(o.pattern, S[θ[off + 1 + i] for i in 1:nv])
        # the Dubiner basis summed over the triangle images, and its derivatives in the values
        A = zeros(S, L)
        dA = zeros(S, L, nv)
        for lab in o.pattern.labels
            x[1], x[2] = vals[lab[2]], vals[lab[3]]
            φ, G = evaluate!(sys.basis, x; gradient = jacobian)
            A .+= φ
            jacobian || continue
            for i in 1:nv, j in 1:2
                l = lab[j + 1]
                d = l == i ? one(S) : l == rr ? -S(o.pattern.mult[i]) / o.pattern.mult[rr] : zero(S)
                iszero(d) || (dA[:, i] .+= d .* view(G, :, j))
            end
        end
        # the triangle invariants, a, and their derivatives
        a = zeros(S, mT)
        da = zeros(S, mT, nv)
        for (rows, cols, B) in sys.blocks
            a[cols] = transpose(B) * A[rows]
            jacobian && (da[cols, :] = transpose(B) * dA[rows, :])
        end
        # the even Legendre polynomials in z, summed over the orbit's heights
        if o.mirrored
            z = S(θ[off + 2 + nv])
            P, dP = legendre_orthonormal(2M, z)
            Z = [2P[2m + 1] for m in 0:(M - 1)]
            dZ = [2dP[2m + 1] for m in 0:(M - 1)]
        else
            P, _ = legendre_orthonormal(2M, zero(S))
            Z = [P[2m + 1] for m in 0:(M - 1)]
            dZ = zeros(S, M)
        end
        for (row, (c, _, m)) in enumerate(sys.rows)
            r[row] += w * a[c] * Z[m + 1]
            jacobian || continue
            J[row, off + 1] += a[c] * Z[m + 1]
            for i in 1:nv
                J[row, off + 1 + i] += w * da[c, i] * Z[m + 1]
            end
            o.mirrored && (J[row, off + 2 + nv] += w * a[c] * dZ[m + 1])
        end
    end
    r[findfirst(t -> t[1] == 1 && t[3] == 0, sys.rows)] -= sys.mass
    return r, J
end

"`(smallest weight, smallest barycentric value or distance of z from 0 and 1, smallest gap between an orbit's values)`"
function wedge_margins(s::WedgeStructure, θ::AbstractVector{S}) where {S}
    wmin, cmin, gmin = typemax(S), typemax(S), typemax(S)
    for (o, off) in zip(s.orbits, param_offsets(s))
        wmin = min(wmin, θ[off + 1])
        nv = ncoords(o.pattern)
        vals = pattern_values(o.pattern, θ[(off + 2):(off + 1 + nv)])
        cmin = min(cmin, minimum(vals))
        for i in eachindex(vals), j in (i + 1):length(vals)
            gmin = min(gmin, abs(vals[i] - vals[j]))
        end
        o.mirrored && (z = θ[off + 2 + nv]; cmin = min(cmin, z, 1 - z))
    end
    return wmin, cmin, gmin
end

"`(s, θ)` with each orbit's values in a canonical order (as for the triangle), `z ≥ 0`, and the orbits sorted."
function canonicalize(s::WedgeStructure, θ::AbstractVector)
    S = eltype(θ)
    blocks = Tuple{WedgeOrbit,Vector{S}}[]
    for (o, off) in zip(s.orbits, param_offsets(s))
        nv = ncoords(o.pattern)
        st, θt = canonicalize(SymmetricStructure([o.pattern], 3), θ[(off + 1):(off + 1 + nv)])
        no = WedgeOrbit(st.orbits[1], o.mirrored)
        push!(blocks, (no, o.mirrored ? vcat(θt, abs(θ[off + 2 + nv])) : collect(θt)))
    end
    sort!(blocks; by = b -> (orbit_size(b[1]), b[1].pattern.mult, b[1].mirrored, Float64.(b[2][2:end])))
    return WedgeStructure([b[1] for b in blocks]), reduce(vcat, (b[2] for b in blocks); init = S[])
end
