# Orthonormal polynomial bases, used both to write moment equations and to verify rules.
#
# Monomials are deliberately not used: their Gram matrix is so ill-conditioned that
# verification in them produces false failures above degree ~15 (PLAN §8).

# ---------------------------------------------------------------------------------------
# Jacobi recurrence coefficients for the orthonormal polynomials of (1-x)^α (1+x)^β:
#
#     b_{k+1} p_{k+1}(x) = (x - a_k) p_k(x) - b_k p_{k-1}(x)

function jacobi_a(k::Integer, α, β)
    s = 2k + α + β
    if k == 0
        return (β - α) / (α + β + 2)
    end
    return (β^2 - α^2) / (s * (s + 2))
end

function jacobi_b(k::Integer, α, β)
    k >= 1 || throw(ArgumentError("b_k is defined for k ≥ 1"))
    if k == 1
        # the general formula has a removable 0/0 when α + β = -1
        return sqrt(4 * (1 + α) * (1 + β) / ((2 + α + β)^2 * (3 + α + β)))
    end
    s = 2k + α + β
    return sqrt(4k * (k + α) * (k + β) * (k + α + β) / (s^2 * (s + 1) * (s - 1)))
end

"""
    jacobi_mass(α, β)

`∫₋₁¹ (1-x)^α (1+x)^β dx = 2^(α+β+1) Γ(α+1) Γ(β+1) / Γ(α+β+2)`. Exact for integer
parameters.
"""
function jacobi_mass(α::Integer, β::Integer)
    return big(2)^(α + β + 1) * factorial(big(α)) * factorial(big(β)) // factorial(big(α + β + 1))
end
function jacobi_mass(α, β)
    # BigFloat at the ambient precision, which callers set explicitly (see `with_bits`)
    a, b = big(α), big(β)
    return 2^(a + b + 1) * SpecialFunctions.gamma(a + 1) * SpecialFunctions.gamma(b + 1) /
           SpecialFunctions.gamma(a + b + 2)
end

"""
    jacobi_orthonormal!(p, x, α, β)

Fill `p[k+1]` with the orthonormal Jacobi polynomial `p_k(x)`, `k = 0:length(p)-1`,
normalised against the weight `(1-x)^α (1+x)^β` on `[-1, 1]`.
"""
function jacobi_orthonormal!(p::AbstractVector{S}, x, α, β) where {S}
    n = length(p) - 1
    μ = jacobi_mass(α, β)
    p[1] = one(S) / sqrt(S(μ))
    n == 0 && return p
    Sα, Sβ = S(α), S(β)
    p[2] = (x - jacobi_a(0, Sα, Sβ)) * p[1] / jacobi_b(1, Sα, Sβ)
    for k in 1:(n - 1)
        p[k + 2] = ((x - jacobi_a(k, Sα, Sβ)) * p[k + 1] - jacobi_b(k, Sα, Sβ) * p[k]) /
                   jacobi_b(k + 1, Sα, Sβ)
    end
    return p
end

"""
    legendre_orthonormal!(p, x)

Orthonormal Legendre polynomials on `[-1, 1]`: `p[k+1] = √((2k+1)/2) P_k(x)`. In exact
arithmetic (`Rational`) the normalisation is omitted and `P_k` is returned.
"""
function legendre_orthonormal!(p::AbstractVector{S}, x; normalize::Bool = true) where {S}
    n = length(p) - 1
    p[1] = one(S)
    n >= 1 && (p[2] = S(x))
    for k in 1:(n - 1)
        p[k + 2] = ((2k + 1) * x * p[k + 1] - k * p[k]) / (k + 1)
    end
    if normalize
        for k in 0:n
            p[k + 1] *= sqrt(S(2k + 1) / 2)
        end
    end
    return p
end

# ---------------------------------------------------------------------------------------
# Proriol–Koornwinder–Dubiner basis on the reference triangle (0,0), (1,0), (0,1).
#
#     φ_pq(x, y) = c_pq · L_p(s, t) · P_q^(2p+1, 0)(2y - 1),
#     s = 2x + y - 1,  t = 1 - y,  L_p(s, t) = t^p P_p(s / t)
#
# L_p is evaluated by the homogenised Legendre recurrence, so there is no division by t and
# the basis is a polynomial everywhere including the collapsed vertex. The orthonormal
# constant is c_pq = √(2 (2p+1) (p+q+1)); the reference triangle has area 1/2.
#
# Ordering: by total degree k = p + q, then by p ascending. `dubiner_index(p, q)` gives the
# 1-based position.

dubiner_index(p::Int, q::Int) = ((p + q) * (p + q + 1)) ÷ 2 + p + 1
dubiner_length(n::Int) = ((n + 1) * (n + 2)) ÷ 2

"""
    DubinerWorkspace{S}(n)

Scratch space for evaluating the degree-`n` Dubiner basis (and gradient) in type `S`.
"""
struct DubinerWorkspace{S}
    n::Int
    L::Vector{S}
    Ls::Vector{S}
    Lt::Vector{S}
    J::Vector{S}
    Jb::Vector{S}
    c::Vector{S}
end

function DubinerWorkspace{S}(n::Int; normalize::Bool = true) where {S}
    c = Vector{S}(undef, dubiner_length(n))
    for k in 0:n, p in 0:k
        q = k - p
        c[dubiner_index(p, q)] = normalize ? sqrt(S(2 * (2p + 1) * (p + q + 1))) : one(S)
    end
    z() = Vector{S}(undef, n + 1)
    return DubinerWorkspace{S}(n, z(), z(), z(), z(), z(), c)
end

"""
    dubiner!(φ, ws, x, y)
    dubiner!(φ, gx, gy, ws, x, y)

Evaluate the Dubiner basis up to degree `ws.n` at `(x, y)` into `φ`; optionally the partial
derivatives into `gx`, `gy`.
"""
dubiner!(φ, ws::DubinerWorkspace, x, y) = _dubiner!(φ, nothing, nothing, ws, x, y)
dubiner!(φ, gx, gy, ws::DubinerWorkspace, x, y) = _dubiner!(φ, gx, gy, ws, x, y)

function _dubiner!(φ, gx, gy, ws::DubinerWorkspace{S}, x, y) where {S}
    n = ws.n
    grad = gx !== nothing
    L, Ls, Lt, J, Jb = ws.L, ws.Ls, ws.Lt, ws.J, ws.Jb
    s = 2x + y - 1
    t = 1 - y
    b = 2y - 1
    # homogenised Legendre L_p(s, t) and its partials
    L[1] = one(S); Ls[1] = zero(S); Lt[1] = zero(S)
    if n >= 1
        L[2] = S(s); Ls[2] = one(S); Lt[2] = zero(S)
    end
    for p in 1:(n - 1)
        L[p + 2] = ((2p + 1) * s * L[p + 1] - p * t^2 * L[p]) / (p + 1)
        if grad
            Ls[p + 2] = ((2p + 1) * (L[p + 1] + s * Ls[p + 1]) - p * t^2 * Ls[p]) / (p + 1)
            Lt[p + 2] = ((2p + 1) * s * Lt[p + 1] - p * (2t * L[p] + t^2 * Lt[p])) / (p + 1)
        end
    end
    for p in 0:n
        α = 2p + 1
        m = n - p                          # maximal q for this p
        J[1] = one(S); Jb[1] = zero(S)
        if m >= 1
            J[2] = (α + 1) + (α + 2) * (b - 1) / 2
            Jb[2] = S(α + 2) / 2
        end
        for q in 2:m
            a1 = 2q * (q + α) * (2q + α - 2)
            a2 = (2q + α - 1) * α^2
            a3 = (2q + α - 2) * (2q + α - 1) * (2q + α)
            a4 = 2 * (q + α - 1) * (q - 1) * (2q + α)
            J[q + 1] = ((a2 + a3 * b) * J[q] - a4 * J[q - 1]) / a1
            if grad
                Jb[q + 1] = ((a2 + a3 * b) * Jb[q] + a3 * J[q] - a4 * Jb[q - 1]) / a1
            end
        end
        for q in 0:m
            i = dubiner_index(p, q)
            c = ws.c[i]
            φ[i] = c * L[p + 1] * J[q + 1]
            if grad
                gx[i] = c * 2 * Ls[p + 1] * J[q + 1]
                gy[i] = c * ((Ls[p + 1] - Lt[p + 1]) * J[q + 1] + 2 * L[p + 1] * Jb[q + 1])
            end
        end
    end
    return φ
end

"""
    dubiner_mass(S)

`∫ φ_00` over the reference triangle for the orthonormal basis: `√2 · ½ = 1/√2`. All
other basis functions integrate to zero.
"""
dubiner_mass(::Type{S}) where {S} = one(S) / sqrt(S(2))
