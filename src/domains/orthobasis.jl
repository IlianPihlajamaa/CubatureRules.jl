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
    # The recurrence coefficients as numbers of type S, for the in-place BigFloat kernel,
    # which must not convert an integer on every step: L[p+2] = A_p s L[p+1] − B_p t² L[p],
    # J[q+1] = (C2 + C3 b) J[q] − C4 J[q-1] (column p + 1), J[2] = D0 + D1 b.
    A::Vector{S}
    B::Vector{S}
    C2::Matrix{S}
    C3::Matrix{S}
    C4::Matrix{S}
    D0::Vector{S}
    D1::Vector{S}
    tmp::Vector{S}
end

"Coefficients of the Jacobi `P^(α,0)` recurrence in the kernels: `(a2/a1, a3/a1, a4/a1)`."
function jacobi_step_coefficients(q::Int, α::Int, ::Type{S}) where {S}
    a1 = 2q * (q + α) * (2q + α - 2)
    a2 = (2q + α - 1) * α^2
    a3 = (2q + α - 2) * (2q + α - 1) * (2q + α)
    a4 = 2 * (q + α - 1) * (q - 1) * (2q + α)
    return S(a2) / S(a1), S(a3) / S(a1), S(a4) / S(a1)
end

function DubinerWorkspace{S}(n::Int; normalize::Bool = true) where {S}
    c = bigfloats(S, dubiner_length(n))
    for k in 0:n, p in 0:k
        q = k - p
        c[dubiner_index(p, q)] = normalize ? sqrt(S(2 * (2p + 1) * (p + q + 1))) : one(S)
    end
    z() = bigfloats(S, n + 1)
    A = [S(2p + 1) / S(p + 1) for p in 0:n]
    B = [S(p) / S(p + 1) for p in 0:n]
    C2, C3, C4 = bigfloats(S, n + 1, n + 1), bigfloats(S, n + 1, n + 1), bigfloats(S, n + 1, n + 1)
    for p in 0:n, q in 2:(n - p)
        C2[q + 1, p + 1], C3[q + 1, p + 1], C4[q + 1, p + 1] = jacobi_step_coefficients(q, 2p + 1, S)
    end
    D0 = [S(2p + 1) / 2 for p in 0:n]
    D1 = [S(2p + 3) / 2 for p in 0:n]
    return DubinerWorkspace{S}(n, z(), z(), z(), z(), z(), c, A, B, C2, C3, C4, D0, D1, bigfloats(S, 8))
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
    homog_legendre_mp!(L, Ls, Lt, A, B, s, t, tt, n, grad, u, v)

In place, for BigFloat: the homogenised Legendre `L_p(s, t)`, `p = 0:n`, and its partials,
with `tt = t²`, the coefficients `A_p = (2p+1)/(p+1)`, `B_p = p/(p+1)` and scratch `u`, `v`.
Shared by the triangle and tetrahedron kernels.
"""
function homog_legendre_mp!(L, Ls, Lt, A, B, s, t, tt, n, grad, u, v)
    mp_set_si!(L[1], 1); mp_set_si!(Ls[1], 0); mp_set_si!(Lt[1], 0)
    if n >= 1
        mp_set!(L[2], s); mp_set_si!(Ls[2], 1); mp_set_si!(Lt[2], 0)
    end
    @inbounds for p in 1:(n - 1)
        a, b = A[p + 1], B[p + 1]
        mp_mul!(u, s, L[p + 1]); mp_mul!(u, u, a)                  # A s L[p+1]
        mp_mul!(v, tt, L[p]); mp_mul!(v, v, b)                     # B t² L[p]
        mp_sub!(L[p + 2], u, v)
        grad || continue
        mp_fma!(u, s, Ls[p + 1], L[p + 1]); mp_mul!(u, u, a)       # A (L[p+1] + s Ls[p+1])
        mp_mul!(v, tt, Ls[p]); mp_mul!(v, v, b)
        mp_sub!(Ls[p + 2], u, v)
        mp_mul!(u, s, Lt[p + 1]); mp_mul!(u, u, a)                 # A s Lt[p+1]
        mp_mul!(v, t, L[p]); mp_twice!(v, v)                       # B (2t L[p] + t² Lt[p])
        mp_fma!(v, tt, Lt[p], v); mp_mul!(v, v, b)
        mp_sub!(Lt[p + 2], u, v)
    end
    return L
end

# The same evaluation for BigFloat, allocation-free (see core/mpfr.jl): 8× faster at degree
# 40, and the results agree with the generic method to the working precision (the
# coefficients are rounded ratios here, where the generic method divides at every step).
function _dubiner!(φ::AbstractVector{BigFloat}, gx, gy, ws::DubinerWorkspace{BigFloat}, x::BigFloat, y::BigFloat)
    n = ws.n
    grad = gx !== nothing
    prec = precision(ws.c[1])
    unshare!(φ, prec)
    grad && (unshare!(gx, prec); unshare!(gy, prec))
    L, Ls, Lt, J, Jb = ws.L, ws.Ls, ws.Lt, ws.J, ws.Jb
    s, t, b, tt, t1, t2, _, f = ws.tmp
    mp_twice!(s, x); mp_add!(s, s, y); mp_sub_si!(s, s, 1)       # s = 2x + y − 1
    mp_si_sub!(t, 1, y)                                            # t = 1 − y
    mp_twice!(b, y); mp_sub_si!(b, b, 1)                           # b = 2y − 1
    mp_mul!(tt, t, t)
    homog_legendre_mp!(L, Ls, Lt, ws.A, ws.B, s, t, tt, n, grad, t1, t2)
    @inbounds for p in 0:n
        m = n - p
        mp_set_si!(J[1], 1); mp_set_si!(Jb[1], 0)
        if m >= 1
            mp_fma!(J[2], ws.D1[p + 1], b, ws.D0[p + 1])
            mp_set!(Jb[2], ws.D1[p + 1])
        end
        for q in 2:m
            C2, C3, C4 = ws.C2[q + 1, p + 1], ws.C3[q + 1, p + 1], ws.C4[q + 1, p + 1]
            mp_fma!(f, C3, b, C2)                                   # f = C2 + C3 b
            mp_mul!(t2, C4, J[q - 1])
            mp_fms!(J[q + 1], f, J[q], t2)                          # f J[q] − C4 J[q−1]
            grad || continue
            mp_mul!(t1, C3, J[q]); mp_fma!(t1, f, Jb[q], t1)        # f Jb[q] + C3 J[q]
            mp_mul!(t2, C4, Jb[q - 1])
            mp_sub!(Jb[q + 1], t1, t2)
        end
        for q in 0:m
            i = dubiner_index(p, q)
            c = ws.c[i]
            mp_mul!(t1, c, L[p + 1]); mp_mul!(φ[i], t1, J[q + 1])     # c L J
            grad || continue
            mp_mul!(t1, c, Ls[p + 1]); mp_twice!(t1, t1)
            mp_mul!(gx[i], t1, J[q + 1])                              # 2c Ls J
            mp_sub!(t1, Ls[p + 1], Lt[p + 1]); mp_mul!(t1, t1, J[q + 1])
            mp_mul!(t2, L[p + 1], Jb[q + 1]); mp_twice!(t2, t2)
            mp_add!(t1, t1, t2); mp_mul!(gy[i], c, t1)               # c ((Ls − Lt) J + 2 L Jb)
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
