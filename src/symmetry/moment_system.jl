# Moment equations for fully symmetric rules: invariant basis × orbits → residual and
# Jacobian (PLAN §6 Tier 3).
#
# Unknowns θ: per orbit [w, v₁, …] (weight per point, free barycentric values).
# Full residual  r_k = Σ_nodes w φ_k(x) − ∫ φ_k   over the orthonormal Dubiner basis to degree n.
# Reduced system Qᵀ r with Q the invariant basis — as many equations as the Molien count.

"""
    TriangleMomentSystem(structure, n, S)

The `S_3`-invariant moment system for exactness to degree `n` of a triangle rule with the
given orbit structure, evaluated in number type `S`. Callable: `sys(θ) -> (r, J)`.
"""
struct TriangleMomentSystem{S}
    structure::SymmetricStructure
    n::Int
    Q::Matrix{S}
    ws::DubinerWorkspace{S}
    φ::Vector{S}
    gx::Vector{S}
    gy::Vector{S}
    mass::S
end

function TriangleMomentSystem(structure::SymmetricStructure, n::Integer, ::Type{S},
                              basis::InvariantBasis = invariant_basis(3, n)) where {S}
    structure.N == 3 || throw(ArgumentError("TriangleMomentSystem needs a triangle structure"))
    basis.n == n || throw(ArgumentError("invariant basis degree mismatch"))
    L = dubiner_length(n)
    return TriangleMomentSystem{S}(structure, Int(n), S.(basis.Q), DubinerWorkspace{S}(n),
                                   zeros(S, L), zeros(S, L), zeros(S, L), dubiner_mass(S))
end

n_equations(sys::TriangleMomentSystem) = size(sys.Q, 2)
n_unknowns(sys::TriangleMomentSystem) = nunknowns(sys.structure)

"Full orthonormal-basis residual and Jacobian (before projection onto invariants)."
function full_residual(sys::TriangleMomentSystem{S}, θ::AbstractVector; jacobian::Bool = true) where {S}
    L = dubiner_length(sys.n)
    p = n_unknowns(sys)
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, p) : zeros(S, 0, 0)
    φ, gx, gy = sys.φ, sys.gx, sys.gy
    for (o, off) in zip(sys.structure.orbits, param_offsets(sys.structure))
        w = S(θ[off + 1])
        nv = ncoords(o)
        v = S[θ[off + 1 + i] for i in 1:nv]
        vals = pattern_values(o, v)
        rr = length(o.mult)
        for lab in o.labels
            x, y = vals[lab[2]], vals[lab[3]]
            if jacobian
                dubiner!(φ, gx, gy, sys.ws, x, y)
            else
                dubiner!(φ, sys.ws, x, y)
            end
            for k in 1:L
                r[k] += w * φ[k]
            end
            jacobian || continue
            for k in 1:L
                J[k, off + 1] += φ[k]
            end
            for i in 1:nv
                # ∂(value with label l)/∂vᵢ: 1 if l == i, −mᵢ/m_r if l is the last value
                dx = lab[2] == i ? one(S) : lab[2] == rr ? -S(o.mult[i]) / o.mult[rr] : zero(S)
                dy = lab[3] == i ? one(S) : lab[3] == rr ? -S(o.mult[i]) / o.mult[rr] : zero(S)
                (iszero(dx) && iszero(dy)) && continue
                for k in 1:L
                    J[k, off + 1 + i] += w * (gx[k] * dx + gy[k] * dy)
                end
            end
        end
    end
    r[1] -= sys.mass
    return r, J
end

function (sys::TriangleMomentSystem)(θ::AbstractVector)
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
