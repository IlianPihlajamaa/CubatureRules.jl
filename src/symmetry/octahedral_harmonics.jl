# The O_h moment system in an orthonormal invariant basis (PLAN §6 Tier 3).
#
# OctahedralMomentSystem writes the equations in the invariants p₄^a p₆^b, which integrate to
# exact rationals but behave like monomials: the system's condition number is 1e54 at degree
# 133 and 1e81 at 201. The precision that asks for is not what costs; the solver is. Above a
# condition number of about 1e10 a least-squares step cannot use a Float64 factorisation with
# BigFloat iterative refinement and falls back to a BigFloat QR: at degree 133, 10 s a step
# against a fraction of a second.
#
# This system is written in the O_h-invariant spherical harmonics instead, which are
# orthonormal. A combination of P̄_ℓ^m(cos θ) cos mφ with ℓ even and m ≡ 0 (mod 4) is invariant
# under D4h, the 16 symmetries that keep the z-axis, exactly and whatever its coefficients.
# Averaged over the rotations x → y → z → x it is O_h-invariant, again exactly: D4h, D4h·R
# and D4h·R² are the three cosets of D4h in O_h, so the average over {I, R, R²} of a
# D4h-invariant function is its average over all of O_h. The coefficients can therefore be
# computed once, in Float64 — the nullspace of f(Rx) − f(x) over random points — and the
# functions
#
#     F(x) = (f(x) + f(Rx) + f(R²x)) / 3
#
# are exactly invariant and orthonormal to Float64 rounding. The equations are exact whatever
# that rounding is; only the conditioning depends on it. Every F of degree ℓ > 0 integrates to
# zero, so there are no moments to compute. The unknowns are those of OctahedralMomentSystem,
# so the same seeds refine in either.

"The number of `O_h`-invariant harmonics of degree `ℓ`: the solutions of `4a + 6b = ℓ`."
invariant_harmonic_count(ℓ::Integer) = count(b -> (ℓ - 6b) % 4 == 0, 0:(ℓ ÷ 6))

# the rotation x → y → z → x and its square, on points and on tangent vectors alike
rotate3(v::SVector{3}, j::Integer) = j == 0 ? v : j == 1 ? SVector(v[2], v[3], v[1]) : SVector(v[3], v[1], v[2])

# s_m P̄_ℓ^m(cos θ) cos mφ for m = 0, 4, 8, …, ℓ, with s_0 = 1 and s_m = √2, in Float64
function d4h_harmonics(ℓ::Integer, x, P::Vector{Float64})
    z = x[3]
    ρ = sqrt(max(0.0, (1 - z) * (1 + z)))
    legendre_normalised!(P, ℓ, z, ρ)
    φ = atan(x[2], x[1])
    return [(m == 0 ? 1.0 : sqrt(2.0)) * P[legendre_index(ℓ, m)] * cos(m * φ) for m in 0:4:ℓ]
end

# Householder QR → the full orthogonal factor, by applying the reflectors to the identity
function householder_q(F::PivotedQR64)
    m = size(F.QR, 1)
    Q = Matrix{Float64}(I, m, m)
    for k in length(F.τ):-1:1, j in 1:m
        w = Q[k, j]
        for i in (k + 1):m
            w += F.QR[i, k] * Q[i, j]
        end
        w *= F.τ[k]
        Q[k, j] -= w
        for i in (k + 1):m
            Q[i, j] -= w * F.QR[i, k]
        end
    end
    return Q
end

const OH_HARMONIC_COEFFS = Dict{Int,Matrix{Float64}}()
const OH_HARMONIC_LOCK = ReentrantLock()

"""
    oh_harmonic_coefficients(ℓ) -> Matrix{Float64}

Orthonormal coefficient vectors, over `s_m P̄_ℓ^m(cos θ) cos mφ` with `m = 0, 4, …, ℓ`, of the
`O_h`-invariant harmonics of degree `ℓ`: the nullspace of `f(Rx) − f(x)` at random points,
from the package's own pivoted QR, so that it is the same on every platform. Cached.
"""
function oh_harmonic_coefficients(ℓ::Integer)
    cached = lock(() -> get(OH_HARMONIC_COEFFS, Int(ℓ), nothing), OH_HARMONIC_LOCK)
    cached === nothing || return cached
    ms = 0:4:ℓ
    M, d = length(ms), invariant_harmonic_count(ℓ)
    C = if d == 0 || isodd(ℓ)
        zeros(M, 0)
    elseif d == M                                             # ℓ = 0: the constant is invariant already
        Matrix{Float64}(I, M, M)
    else
        rng = Random.Xoshiro(0x0c7a + ℓ)
        K = 2M + 8
        A = zeros(2K, M)
        P = zeros(legendre_length(ℓ))
        for k in 1:K
            x = SVector{3,Float64}(normalize(randn(rng, 3)))
            e0 = d4h_harmonics(ℓ, x, P)
            A[2k - 1, :] = d4h_harmonics(ℓ, rotate3(x, 1), P) .- e0
            A[2k, :] = d4h_harmonics(ℓ, rotate3(x, 2), P) .- e0
        end
        F = pivoted_qr64(Matrix(A'))                       # Aᵀ = Q R: the last d columns of Q span null(A)
        r = M - d
        big, small = r == 0 ? Inf : abs(F.QR[r, r]), abs(F.QR[min(r + 1, M), min(r + 1, M)])
        (big > 1e-6 * abs(F.QR[1, 1]) && small < 1e-10 * abs(F.QR[1, 1])) ||
            error("O_h-invariant harmonics of degree $ℓ: no clear rank gap ($big, $small)")
        householder_q(F)[:, (r + 1):M]
    end
    lock(() -> (OH_HARMONIC_COEFFS[Int(ℓ)] = C), OH_HARMONIC_LOCK)
    return C
end

# ∂x/∂(parameter) at the orbit's representative, one tangent vector per free parameter
function representative_derivatives(o::OctahedralOrbit, p, ::Type{S}) where {S}
    k = o.kind
    k in (:a1, :a2, :a3) && return SVector{3,S}[]
    if k === :b                                                # (l, l, √(1−2l²))
        l = S(p[1])
        return [SVector{3,S}(one(S), one(S), -2l / sqrt(one(S) - 2l^2))]
    elseif k === :c                                            # (q, √(1−q²), 0)
        q = S(p[1])
        return [SVector{3,S}(one(S), -q / sqrt(one(S) - q^2), zero(S))]
    else                                                       # (r, s, √(1−r²−s²))
        r, s = S(p[1]), S(p[2])
        u = sqrt(one(S) - r^2 - s^2)
        return [SVector{3,S}(one(S), zero(S), -r / u), SVector{3,S}(zero(S), one(S), -s / u)]
    end
end

"""
    OctahedralHarmonicSystem(structure, n, S)

The `O_h`-invariant moment system for exactness to degree `n` on the sphere, in the
orthonormal invariant harmonics: the same unknowns and the same solutions as
[`OctahedralMomentSystem`](@ref), with a condition number that does not grow like the
monomials'. Callable: `sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)`.
"""
struct OctahedralHarmonicSystem{S}
    structure::OctahedralStructure
    n::Int
    ℓs::Vector{Int}                 # the degrees that carry invariant harmonics
    coefs::Vector{Matrix{S}}        # per degree: (m = 0, 4, …, ℓ) × invariant harmonics
    neq::Int
    rec_a::Matrix{S}                # the normalised Legendre recurrence, a[l+1, m+1] and b[l+1, m+1]
    rec_b::Matrix{S}
    rec_d::Vector{S}                # P̄_m^m = rec_d[m] sin θ P̄_{m−1}^{m−1}, P̄_0^0 = 1/√(4π) = rec_d[n+1]
    rec_s::Vector{S}                # P̄_{m+1}^m = rec_s[m+1] cos θ P̄_m^m
    dl::Vector{S}                   # per (ℓ ∈ ℓs, m = 0, 4, …, ℓ) in that order, √((2ℓ+1)(ℓ²−m²)/(2ℓ−1))
end

function OctahedralHarmonicSystem(structure::OctahedralStructure, n::Integer, ::Type{S}) where {S}
    ℓs = [ℓ for ℓ in 0:2:n if invariant_harmonic_count(ℓ) > 0]
    coefs = [S.(oh_harmonic_coefficients(ℓ)) for ℓ in ℓs]
    a, b = zeros(S, n + 2, n + 2), zeros(S, n + 2, n + 2)
    for m in 0:(n + 1), l in (m + 2):(n + 1)
        a[l + 1, m + 1] = sqrt(S(4l^2 - 1) / S(l^2 - m^2))
        b[l + 1, m + 1] = sqrt(S((l - 1)^2 - m^2) / S(4 * (l - 1)^2 - 1))
    end
    d = [[-sqrt(S(2m + 1) / S(2m)) for m in 1:n]; one(S) / sqrt(4 * S(π))]
    s = [sqrt(S(2m + 3)) for m in 0:n]
    dl = [ℓ == 0 ? zero(S) : sqrt(S((2ℓ + 1) * (ℓ^2 - m^2)) / S(2ℓ - 1)) for ℓ in ℓs for m in 0:4:ℓ]
    return OctahedralHarmonicSystem{S}(structure, Int(n), ℓs, coefs, sum(size.(coefs, 2)), a, b, d, s, dl)
end
n_equations(sys::OctahedralHarmonicSystem) = sys.neq
n_unknowns(sys::OctahedralHarmonicSystem) = nunknowns(sys.structure)

# P̄_l^m(z) for 0 ≤ m ≤ l ≤ N, into P[l+1, m+1], by the recurrence of legendre_normalised!
# with its coefficients precomputed. Only the columns the system reads are filled, m ≡ 0
# (mod 4); the diagonal is needed throughout to start them.
function legendre_table!(P::Matrix{S}, N::Integer, z::S, ρ::S, sys::OctahedralHarmonicSystem{S}) where {S}
    P[1, 1] = sys.rec_d[end]
    for m in 1:N
        P[m + 1, m + 1] = sys.rec_d[m] * ρ * P[m, m]
    end
    for m in 0:4:N
        m < N && (P[m + 2, m + 1] = sys.rec_s[m + 1] * z * P[m + 1, m + 1])
        for l in (m + 2):N
            P[l + 1, m + 1] = sys.rec_a[l + 1, m + 1] * (z * P[l, m + 1] - sys.rec_b[l + 1, m + 1] * P[l - 1, m + 1])
        end
    end
    return P
end

# The same in BigFloat, in place (core/mpfr.jl): the table is filled three times per orbit,
# and allocating every operation made it most of a system evaluation. `P` must hold
# distinct numbers at the working precision (`bigfloats`).
function legendre_table!(P::Matrix{BigFloat}, N::Integer, z::BigFloat, ρ::BigFloat,
                         sys::OctahedralHarmonicSystem{BigFloat})
    mp_set!(P[1, 1], sys.rec_d[end])
    for m in 1:N
        mp_mul!(P[m + 1, m + 1], sys.rec_d[m], P[m, m])
        mp_mul!(P[m + 1, m + 1], P[m + 1, m + 1], ρ)
    end
    t = BigFloat(0)
    for m in 0:4:N
        if m < N
            mp_mul!(P[m + 2, m + 1], sys.rec_s[m + 1], P[m + 1, m + 1])
            mp_mul!(P[m + 2, m + 1], P[m + 2, m + 1], z)
        end
        for l in (m + 2):N
            mp_mul!(t, sys.rec_b[l + 1, m + 1], P[l - 1, m + 1])     # b P̄_{l−2}
            mp_fms!(t, z, P[l, m + 1], t)                            # z P̄_{l−1} − b P̄_{l−2}
            mp_mul!(P[l + 1, m + 1], sys.rec_a[l + 1, m + 1], t)
        end
    end
    return P
end

# v[i] += x y, and d = x y or d = x y − a: in place for BigFloat (core/mpfr.jl), whose
# allocating arithmetic was most of an evaluation; by value otherwise. Use the returned d.
@inline fma_into!(v::Vector{BigFloat}, i, x::BigFloat, y::BigFloat) = (mp_fma!(v[i], x, y, v[i]); nothing)
@inline fma_into!(v::Vector, i, x, y) = (v[i] = muladd(x, y, v[i]); nothing)
@inline mul_to(d::BigFloat, x::BigFloat, y::BigFloat) = mp_mul!(d, x, y)
@inline mul_to(_, x, y) = x * y
@inline fms_to(d::BigFloat, x::BigFloat, y::BigFloat, a::BigFloat) = mp_fms!(d, x, y, a)
@inline fms_to(_, x, y, a) = x * y - a
zero_out!(v::Vector{BigFloat}) = foreach(x -> mp_set_si!(x, 0), v)
zero_out!(v::Vector) = fill!(v, 0)

# Σ_m C[m, k] v[off + m]
function project(acc::BigFloat, C::Matrix{BigFloat}, k, v::Vector{BigFloat}, off)
    mp_set_si!(acc, 0)
    for m in axes(C, 1)
        mp_fma!(acc, C[m, k], v[off + m], acc)
    end
    return acc
end
project(_, C, k, v, off) = sum(C[m, k] * v[off + m] for m in axes(C, 1))

# The equations are accumulated in the D4h harmonics s_m P̄_ℓ^m cos mφ, summed over the three
# rotations, and projected onto the invariant harmonics once per orbit: projecting each
# rotation and each derivative separately was most of the cost.
function (sys::OctahedralHarmonicSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    L, p, n = sys.neq, n_unknowns(sys), sys.n
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, p) : zeros(S, 0, 0)
    P = bigfloats(S, n + 2, n + 2)
    raw = [bigfloats(S, length(sys.dl)) for _ in 1:3]     # Σ_j f(R^j x), then its two derivatives
    dbuf, acc = bigfloats(S, 2)
    s2 = sqrt(S(2))
    nm = n ÷ 4 + 1
    ccm, ssm = zeros(S, nm), zeros(S, nm)                 # s_m cos mφ and −s_m m sin mφ / sin θ
    A, B = [zeros(S, nm) for _ in 1:2], [zeros(S, nm) for _ in 1:2]
    offs = param_offsets(sys.structure)
    for (i, o) in enumerate(sys.structure.orbits)
        w = S(θ[offs[i] + 1])
        par = θ[(offs[i] + 2):(offs[i] + nunknowns(o))]
        n_o = S(orbit_size(o))
        x = orbit_representative(o, par, S)
        dx = jacobian ? representative_derivatives(o, par, S) : SVector{3,S}[]
        foreach(zero_out!, raw)
        for j in 0:2
            y = rotate3(x, j)
            z = y[3]
            ρ = sqrt(max(zero(S), (one(S) - z) * (one(S) + z)))
            legendre_table!(P, n, z, ρ, sys)
            cφ, sφ = iszero(ρ) ? (one(S), zero(S)) : (y[1] / ρ, y[2] / ρ)
            eθ, eφ = SVector{3,S}(z * cφ, z * sφ, -ρ), SVector{3,S}(-sφ, cφ, zero(S))
            dθs = [dot(eθ, rotate3(d, j)) for d in dx]          # how the image moves in θ and φ
            dφs = [dot(eφ, rotate3(d, j)) for d in dx]
            c2, s2φ = cφ^2 - sφ^2, 2cφ * sφ
            c4, s4 = c2^2 - s2φ^2, 2c2 * s2φ                    # cos 4φ, sin 4φ; then mφ by rotation
            cm, sm = one(S), zero(S)
            for t in 1:nm
                m = 4(t - 1)
                t > 1 && ((cm, sm) = (cm * c4 - sm * s4, sm * c4 + cm * s4))
                k = m == 0 ? one(S) : s2
                ccm[t] = k * cm
                ssm[t] = iszero(ρ) || m == 0 ? zero(S) : -k * m * sm / ρ
                for q in eachindex(dx)                          # ∂/∂par = fθ dθ + (fφ / sin θ) dφ
                    A[q][t] = ssm[t] * dφs[q]
                    B[q][t] = iszero(ρ) ? zero(S) : -ccm[t] * dθs[q] / ρ
                end
            end
            idx = 0
            for ℓ in sys.ℓs
                zℓ = ℓ * z
                for (t, m) in enumerate(0:4:ℓ)
                    idx += 1
                    Pm = P[ℓ + 1, m + 1]
                    fma_into!(raw[1], idx, Pm, ccm[t])
                    isempty(dx) && continue
                    # d = −sin θ dP̄_ℓ^m/dθ = dl P̄_{ℓ−1}^m − ℓ cos θ P̄_ℓ^m, from the column already
                    # in the table (dl = 0 for ℓ = 0, and P̄_{ℓ−1}^ℓ = 0 for m = ℓ)
                    d = mul_to(dbuf, zℓ, Pm)
                    d = fms_to(dbuf, sys.dl[idx], P[max(ℓ, 1), m + 1], d)
                    for q in eachindex(dx)
                        fma_into!(raw[1 + q], idx, Pm, A[q][t])
                        fma_into!(raw[1 + q], idx, d, B[q][t])
                    end
                end
            end
        end
        c = n_o / 3                                             # F = (f(x) + f(Rx) + f(R²x)) / 3
        row, off = 0, 0
        for C in sys.coefs
            for k in axes(C, 2)
                row += 1
                v = project(acc, C, k, raw[1], off)
                r[row] += w * c * v
                jacobian || continue
                J[row, offs[i] + 1] = c * v
                for q in eachindex(dx)
                    J[row, offs[i] + 1 + q] = w * c * project(acc, C, k, raw[1 + q], off)
                end
            end
            off += size(C, 1)
        end
    end
    r[1] -= sqrt(4 * S(π))                                      # ∫ F₀ = ∫ 1/√(4π); the rest integrate to 0
    return r, J
end
