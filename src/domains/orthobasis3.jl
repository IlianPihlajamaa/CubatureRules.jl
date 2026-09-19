# Proriol–Koornwinder–Dubiner basis on the reference tetrahedron 0, e₁, e₂, e₃, and a common
# interface over the simplex bases used by the moment systems, the invariant bases and
# verification.
#
#     φ_pqr = c · L_p(s₁, t₁) · Q_q^(2p+1)(s₂, t₂) · P_r^(2p+2q+2, 0)(2z − 1)
#
#     s₁ = 2x + y + z − 1,  t₁ = 1 − y − z       L_p(s, t) = t^p P_p(s/t)
#     s₂ = 2y + z − 1,      t₂ = 1 − z           Q_q^α(s, t) = t^q P_q^(α,0)(s/t)
#
# Both homogenised factors are evaluated by recurrences in (s, t), so there is no division
# anywhere and the basis is polynomial on the whole closed tetrahedron. With the collapsed
# coordinates a = s₁/t₁, b = s₂/t₂, c = 2z − 1 the Jacobian is (1/8)((1−b)/2)((1−c)/2)², which
# gives ∫ φ_pqr² = c² / (2 (2p+1)(p+q+1)(2p+2q+2r+3)); c is chosen to make that 1. At
# p = q = r = 0 the norm is 1/6, the volume.
#
# Ordering: by total degree k = p + q + r, then lexicographically in (p, q); degree blocks
# are contiguous, as the invariant-basis construction requires.

"Indices (p, q, r) of the tetrahedral basis up to degree n, in basis order."
function tet_indices(n::Int)
    out = NTuple{3,Int}[]
    for k in 0:n, p in 0:k, q in 0:(k - p)
        push!(out, (p, q, k - p - q))
    end
    return out
end
tet_length(n::Int) = ((n + 1) * (n + 2) * (n + 3)) ÷ 6
tet_block(k::Int) = (tet_length(k - 1) + 1):tet_length(k)

"""
    TetWorkspace{S}(n; normalize = true)

Scratch space for the degree-`n` tetrahedral Dubiner basis and its gradient in type `S`.
"""
struct TetWorkspace{S}
    n::Int
    idx::Vector{NTuple{3,Int}}
    c::Vector{S}
    L::Vector{S}; Ls::Vector{S}; Lt::Vector{S}
    Q::Matrix{S}; Qs::Matrix{S}; Qt::Matrix{S}      # Q[q+1, p+1]
    R::Matrix{S}; Rz::Matrix{S}                      # R[r+1, p+q+1]
end

function TetWorkspace{S}(n::Int; normalize::Bool = true) where {S}
    idx = tet_indices(n)
    c = [normalize ? sqrt(S(2 * (2p + 1) * (p + q + 1) * (2p + 2q + 2r + 3))) : one(S) for (p, q, r) in idx]
    v() = Vector{S}(undef, n + 1)
    m() = Matrix{S}(undef, n + 1, n + 1)
    return TetWorkspace{S}(n, idx, c, v(), v(), v(), m(), m(), m(), m(), m())
end

# Homogenised Jacobi Q_n^(α,0)(s, t) = t^n P_n^(α,0)(s/t), n = 0:m, with partials.
function _homog_jacobi!(Q, Qs, Qt, col, α, m, s, t, grad, ::Type{S}) where {S}
    Q[1, col] = one(S); Qs[1, col] = zero(S); Qt[1, col] = zero(S)
    m >= 1 || return
    Q[2, col] = (α + 1) * t + (α + 2) * (s - t) / 2
    Qs[2, col] = S(α + 2) / 2
    Qt[2, col] = S(α) / 2
    for n in 2:m
        a1 = 2n * (n + α) * (2n + α - 2)
        a2 = (2n + α - 1) * α^2
        a3 = (2n + α - 2) * (2n + α - 1) * (2n + α)
        a4 = 2 * (n + α - 1) * (n - 1) * (2n + α)
        f = a2 * t + a3 * s
        Q[n + 1, col] = (f * Q[n, col] - a4 * t^2 * Q[n - 1, col]) / a1
        if grad
            Qs[n + 1, col] = (a3 * Q[n, col] + f * Qs[n, col] - a4 * t^2 * Qs[n - 1, col]) / a1
            Qt[n + 1, col] = (a2 * Q[n, col] + f * Qt[n, col] - a4 * (2t * Q[n - 1, col] + t^2 * Qt[n - 1, col])) / a1
        end
    end
end

# Plain Jacobi P_n^(α,0)(x), n = 0:m, and derivative, into column `col`.
function _jacobi_col!(P, Pd, col, α, m, x, grad, ::Type{S}) where {S}
    P[1, col] = one(S); Pd[1, col] = zero(S)
    m >= 1 || return
    P[2, col] = (α + 1) + (α + 2) * (x - 1) / 2
    Pd[2, col] = S(α + 2) / 2
    for n in 2:m
        a1 = 2n * (n + α) * (2n + α - 2)
        a2 = (2n + α - 1) * α^2
        a3 = (2n + α - 2) * (2n + α - 1) * (2n + α)
        a4 = 2 * (n + α - 1) * (n - 1) * (2n + α)
        P[n + 1, col] = ((a2 + a3 * x) * P[n, col] - a4 * P[n - 1, col]) / a1
        grad && (Pd[n + 1, col] = ((a2 + a3 * x) * Pd[n, col] + a3 * P[n, col] - a4 * Pd[n - 1, col]) / a1)
    end
end

"""
    tet_dubiner!(φ, G, ws, x, y, z)

Evaluate the tetrahedral basis at `(x, y, z)` into `φ`; if `G` (an `L × 3` matrix) is not
`nothing`, also the gradient.
"""
function tet_dubiner!(φ, G, ws::TetWorkspace{S}, x, y, z) where {S}
    n = ws.n
    grad = G !== nothing
    s1 = 2x + y + z - 1; t1 = 1 - y - z
    s2 = 2y + z - 1;     t2 = 1 - z
    c3 = 2z - 1
    L, Ls, Lt = ws.L, ws.Ls, ws.Lt
    L[1] = one(S); Ls[1] = zero(S); Lt[1] = zero(S)
    if n >= 1
        L[2] = S(s1); Ls[2] = one(S); Lt[2] = zero(S)
    end
    for p in 1:(n - 1)
        L[p + 2] = ((2p + 1) * s1 * L[p + 1] - p * t1^2 * L[p]) / (p + 1)
        if grad
            Ls[p + 2] = ((2p + 1) * (L[p + 1] + s1 * Ls[p + 1]) - p * t1^2 * Ls[p]) / (p + 1)
            Lt[p + 2] = ((2p + 1) * s1 * Lt[p + 1] - p * (2t1 * L[p] + t1^2 * Lt[p])) / (p + 1)
        end
    end
    for p in 0:n
        _homog_jacobi!(ws.Q, ws.Qs, ws.Qt, p + 1, 2p + 1, n - p, s2, t2, grad, S)
    end
    for pq in 0:n
        _jacobi_col!(ws.R, ws.Rz, pq + 1, 2pq + 2, n - pq, c3, grad, S)
    end
    for (i, (p, q, r)) in enumerate(ws.idx)
        c = ws.c[i]
        l, q_, r_ = L[p + 1], ws.Q[q + 1, p + 1], ws.R[r + 1, p + q + 1]
        φ[i] = c * l * q_ * r_
        if grad
            ls, lt = Ls[p + 1], Lt[p + 1]
            qs, qt = ws.Qs[q + 1, p + 1], ws.Qt[q + 1, p + 1]
            rz = ws.Rz[r + 1, p + q + 1]
            G[i, 1] = c * 2ls * q_ * r_
            G[i, 2] = c * ((ls - lt) * q_ * r_ + l * 2qs * r_)
            G[i, 3] = c * ((ls - lt) * q_ * r_ + l * (qs - qt) * r_ + l * q_ * 2rz)
        end
    end
    return φ
end

# ---------------------------------------------------------------------------------------
# A common interface over the orthonormal simplex bases.

"""
    SimplexBasis{D,S}(n; normalize = true)

The orthonormal Dubiner basis of degree `≤ n` on the reference `D`-simplex (`D = 2, 3`), in
type `S`. `evaluate!(b, x)` fills and returns `(φ, G)`, the values and the `L × D` gradient.
"""
struct SimplexBasis{D,S,W}
    n::Int
    ws::W
    φ::Vector{S}
    G::Matrix{S}
    gx::Vector{S}     # scratch for the 2D evaluator
    gy::Vector{S}
end

function SimplexBasis{D,S}(n::Int; normalize::Bool = true) where {D,S}
    L = simplex_basis_length(D, n)
    ws = D == 2 ? DubinerWorkspace{S}(n; normalize) :
         D == 3 ? TetWorkspace{S}(n; normalize) :
         throw(NotYetImplemented("orthonormal bases on $D-simplices", "a later release"))
    return SimplexBasis{D,S,typeof(ws)}(n, ws, zeros(S, L), zeros(S, L, D), zeros(S, L), zeros(S, L))
end

simplex_basis_length(D::Int, n::Int) = D == 2 ? dubiner_length(n) : D == 3 ? tet_length(n) :
                                       binomial(n + D, D)
basis_length(b::SimplexBasis{D}) where {D} = simplex_basis_length(D, b.n)

"Index range of the degree-`k` block."
function degree_block(D::Int, k::Int)
    D == 2 && return dubiner_index(0, k):dubiner_index(k, 0)
    D == 3 && return tet_block(k)
    throw(NotYetImplemented("orthonormal bases on $D-simplices", "a later release"))
end
degree_blocks(b::SimplexBasis{D}) where {D} = [degree_block(D, k) for k in 0:(b.n)]

function evaluate!(b::SimplexBasis{2}, x; gradient::Bool = true)
    if gradient
        dubiner!(b.φ, b.gx, b.gy, b.ws, x[1], x[2])
        b.G[:, 1] .= b.gx
        b.G[:, 2] .= b.gy
    else
        dubiner!(b.φ, b.ws, x[1], x[2])
    end
    return b.φ, b.G
end

function evaluate!(b::SimplexBasis{3}, x; gradient::Bool = true)
    tet_dubiner!(b.φ, gradient ? b.G : nothing, b.ws, x[1], x[2], x[3])
    return b.φ, b.G
end

"`∫ φ₁` over the reference simplex: `c₀ / D!` (all other basis functions integrate to 0)."
simplex_basis_mass(b::SimplexBasis{D,S}) where {D,S} = b.ws.c[1] / S(factorial(D))
