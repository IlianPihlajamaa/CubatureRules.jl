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
Callable: `sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)` for the residual alone.
"""
struct SymmetricMomentSystem{S,D,B<:SimplexBasis{D,S}}
    structure::SymmetricStructure
    n::Int
    # Q, one diagonal block per degree: (basis functions of degree k, their invariants, block).
    # The invariants of degree k combine only basis functions of degree k, so Q is block
    # diagonal and 96% zeros at degree 30; multiplying by it densely was 40% of a build.
    blocks::Vector{Tuple{UnitRange{Int},UnitRange{Int},Matrix{S}}}
    m::Int
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
    return SymmetricMomentSystem{S,D,typeof(b)}(structure, Int(n), invariant_blocks(basis, b, S),
                                                size(basis.Q, 2), b, mass)
end

"The diagonal blocks of `basis.Q`, by degree, in type `S`; checks that nothing lies outside them."
function invariant_blocks(basis::InvariantBasis, b::SimplexBasis, ::Type{S}) where {S}
    out = Tuple{UnitRange{Int},UnitRange{Int},Matrix{S}}[]
    c = 0
    for (rows, rk) in zip(degree_blocks(b), basis.ranks)
        cols = (c + 1):(c + rk)
        c += rk
        rk == 0 && continue
        push!(out, (rows, cols, S.(basis.Q[rows, cols])))
    end
    for (rows, cols, _) in out, j in cols, i in axes(basis.Q, 1)
        i in rows || iszero(basis.Q[i, j]) || error("the invariant basis is not block diagonal by degree")
    end
    return out
end

"`Qᵀ v` for a vector or matrix `v`, one degree block at a time."
function project(sys::SymmetricMomentSystem{S}, v::AbstractVector) where {S}
    out = zeros(S, sys.m)
    for (rows, cols, B) in sys.blocks
        out[cols] = transpose(B) * view(v, rows)
    end
    return out
end
function project(sys::SymmetricMomentSystem{S}, J::AbstractMatrix) where {S}
    out = zeros(S, sys.m, size(J, 2))
    for (rows, cols, B) in sys.blocks
        out[cols, :] = transpose(B) * view(J, rows, :)
    end
    return out
end

"The triangle case of [`SymmetricMomentSystem`](@ref)."
TriangleMomentSystem(structure::SymmetricStructure, n::Integer, ::Type{S}, args...) where {S} =
    SymmetricMomentSystem(structure, n, S, args...)

n_equations(sys::SymmetricMomentSystem) = sys.m
n_unknowns(sys::SymmetricMomentSystem) = nunknowns(sys.structure)

# The three accumulations of `full_residual`: r += w φ, J[:, col] += φ and
# J[:, col] += w (G d). For BigFloat they add in place (see core/mpfr.jl), which halves the
# time of a Jacobian at high degree; for other types they are the plain expressions.
function add_scaled!(r, w, φ, L)
    @inbounds for k in 1:L
        r[k] += w * φ[k]
    end
end
function add_scaled!(r::Vector{BigFloat}, w::BigFloat, φ::AbstractVector{BigFloat}, L)
    @inbounds for k in 1:L
        mp_fma!(r[k], w, φ[k], r[k])
    end
end
function add_column!(J, col, φ, L)
    @inbounds for k in 1:L
        J[k, col] += φ[k]
    end
end
function add_column!(J::Matrix{BigFloat}, col, φ::AbstractVector{BigFloat}, L)
    @inbounds for k in 1:L
        mp_add!(J[k, col], J[k, col], φ[k])
    end
end
function add_gradient!(J, col, G, d, w, L, D)
    @inbounds for k in 1:L
        g = G[k, 1] * d[1]
        for j in 2:D
            g += G[k, j] * d[j]
        end
        J[k, col] += w * g
    end
end
function add_gradient!(J::Matrix{BigFloat}, col, G::AbstractMatrix{BigFloat}, d::Vector{BigFloat}, w::BigFloat, L, D)
    wd = [iszero(d[j]) ? nothing : w * d[j] for j in 1:D]
    @inbounds for k in 1:L, j in 1:D
        wd[j] === nothing || mp_fma!(J[k, col], G[k, j], wd[j], J[k, col])
    end
end

"Full orthonormal-basis residual and Jacobian (before projection onto invariants)."
function full_residual(sys::SymmetricMomentSystem{S,D}, θ::AbstractVector; jacobian::Bool = true) where {S,D}
    L = basis_length(sys.basis)
    p = n_unknowns(sys)
    r = bigfloats(S, L)              # distinct entries: the BigFloat path adds in place
    J = jacobian ? bigfloats(S, L, p) : zeros(S, 0, 0)
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
            add_scaled!(r, w, φ, L)
            jacobian || continue
            add_column!(J, off + 1, φ, L)
            for i in 1:nv
                # ∂(value with label l)/∂vᵢ: 1 if l == i, −mᵢ/m_r if l is the last value
                for j in 1:D
                    l = lab[j + 1]
                    d[j] = l == i ? one(S) : l == rr ? -S(o.mult[i]) / o.mult[rr] : zero(S)
                end
                all(iszero, d) && continue
                add_gradient!(J, off + 1 + i, G, d, w, L, D)
            end
        end
    end
    r[1] -= sys.mass
    return r, J
end

function (sys::SymmetricMomentSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    r, J = full_residual(sys, θ; jacobian)
    return project(sys, r), jacobian ? project(sys, J) : zeros(S, 0, 0)
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
