# Moment equations for fully symmetric rules: invariant basis × orbits → residual and
# Jacobian (PLAN §6 Tier 3). The same code serves the triangle (S₃) and the tetrahedron (S₄).
#
# Unknowns θ: per orbit [w, v₁, …] (weight per point, free barycentric values).
# Full residual  r_k = Σ_nodes w φ_k(x) − ∫ φ_k   over the orthonormal Dubiner basis to degree n.
# Reduced system Qᵀ r with Q the invariant basis — as many equations as the Molien count.

"""
    SymmetricMomentSystem(structure, n, S, basis = invariant_basis(structure.N, n))

The `S_N`-invariant moment system for exactness to degree `n` of a fully symmetric rule on
the reference `(N-1)`-simplex with the given orbit structure, evaluated in number type `S`.
Callable: `sys(θ) -> (r, J)`.
"""
struct SymmetricMomentSystem{S,D,B<:SimplexBasis{D,S}}
    structure::SymmetricStructure
    n::Int
    Q::Matrix{S}
    basis::B
    mass::S
end

function SymmetricMomentSystem(structure::SymmetricStructure, n::Integer, ::Type{S},
                               basis::InvariantBasis = invariant_basis(structure.N, n)) where {S}
    basis.n == n || throw(ArgumentError("invariant basis degree mismatch"))
    basis.N == structure.N || throw(ArgumentError("invariant basis is for S_$(basis.N), structure for S_$(structure.N)"))
    D = structure.N - 1
    b = SimplexBasis{D,S}(Int(n))
    mass = D == 2 ? dubiner_mass(S) : simplex_basis_mass(b)
    return SymmetricMomentSystem{S,D,typeof(b)}(structure, Int(n), S.(basis.Q), b, mass)
end

"The triangle case of [`SymmetricMomentSystem`](@ref)."
TriangleMomentSystem(structure::SymmetricStructure, n::Integer, ::Type{S}, args...) where {S} =
    SymmetricMomentSystem(structure, n, S, args...)

n_equations(sys::SymmetricMomentSystem) = size(sys.Q, 2)
n_unknowns(sys::SymmetricMomentSystem) = nunknowns(sys.structure)

"Full orthonormal-basis residual and Jacobian (before projection onto invariants)."
function full_residual(sys::SymmetricMomentSystem{S,D}, θ::AbstractVector; jacobian::Bool = true) where {S,D}
    L = basis_length(sys.basis)
    p = n_unknowns(sys)
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, p) : zeros(S, 0, 0)
    x = Vector{S}(undef, D)
    d = Vector{S}(undef, D)
    for (o, off) in zip(sys.structure.orbits, param_offsets(sys.structure))
        w = S(θ[off + 1])
        nv = ncoords(o)
        v = S[θ[off + 1 + i] for i in 1:nv]
        vals = pattern_values(o, v)
        rr = length(o.mult)
        for lab in o.labels
            for j in 1:D
                x[j] = vals[lab[j + 1]]
            end
            φ, G = evaluate!(sys.basis, x; gradient = jacobian)
            for k in 1:L
                r[k] += w * φ[k]
            end
            jacobian || continue
            for k in 1:L
                J[k, off + 1] += φ[k]
            end
            for i in 1:nv
                # ∂(value with label l)/∂vᵢ: 1 if l == i, −mᵢ/m_r if l is the last value
                for j in 1:D
                    l = lab[j + 1]
                    d[j] = l == i ? one(S) : l == rr ? -S(o.mult[i]) / o.mult[rr] : zero(S)
                end
                all(iszero, d) && continue
                for k in 1:L
                    g = G[k, 1] * d[1]
                    for j in 2:D
                        g += G[k, j] * d[j]
                    end
                    J[k, off + 1 + i] += w * g
                end
            end
        end
    end
    r[1] -= sys.mass
    return r, J
end

function (sys::SymmetricMomentSystem)(θ::AbstractVector)
    r, J = full_residual(sys, θ)
    Qt = transpose(sys.Q)
    return Qt * r, Qt * J
end

# ---------------------------------------------------------------------------------------
# Validity and canonical form of a parameter vector.

"Smallest weight and smallest barycentric coordinate over all nodes."
function rule_margins(structure::SymmetricStructure, θ::AbstractVector)
    λs, ws = expand(structure, θ)
    return minimum(ws), minimum(minimum, λs)
end

"""
    canonicalize(structure, θ) -> (structure, θ)

Canonical parameters: within each orbit, values of equal multiplicity sorted ascending;
orbits of the same pattern sorted by their parameters; orbits ordered by pattern type.
Two parameter vectors describing the same rule canonicalise identically.
"""
function canonicalize(structure::SymmetricStructure, θ::AbstractVector{S}) where {S}
    blocks = Tuple{Vector{Int},Vector{S}}[]
    for (o, off) in zip(structure.orbits, param_offsets(structure))
        w = θ[off + 1]
        vals = pattern_values(o, θ[(off + 2):(off + nunknowns(o))])
        # sort values within groups of equal multiplicity; the pattern's mult vector is
        # non-increasing, so groups are contiguous
        mult = o.mult
        order = sortperm(collect(zip(-mult, vals)))
        m2 = mult[order]
        v2 = vals[order]
        push!(blocks, (m2, vcat(w, v2[1:(end - 1)])))
    end
    sort!(blocks; by = b -> (-length(b[1]), b[1], b[2][2:end]...), lt = _lexless)
    patterns = [OrbitPattern(b[1], structure.N) for b in blocks]
    return SymmetricStructure(structure.N, patterns), reduce(vcat, (b[2] for b in blocks); init = S[])
end

_lexless(a::Tuple, b::Tuple) = isless(a, b)

"A uniformly random interior parameter vector for `structure`, weights near equal."
function random_parameters(rng::AbstractRNG, structure::SymmetricStructure, mass::Real)
    θ = Float64[]
    npts = npoints(structure)
    for o in structure.orbits
        push!(θ, mass / npts * (0.5 + rand(rng)))
        r = length(o.mult)
        u = -log.(rand(rng, r))            # Dirichlet(1, …, 1)
        u ./= sum(u)
        for i in 1:(r - 1)
            push!(θ, u[i] / o.mult[i])
        end
    end
    return θ
end
