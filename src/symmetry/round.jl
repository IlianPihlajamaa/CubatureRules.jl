# Fully symmetric rules on the disk and the ball (PLAN §6 Tier 3): orbits of the signed
# permutations of the coordinates, the groups B₂ (the 8 symmetries of the square) and B₃ = O_h,
# which are those of symmetry/box.jl. An orbit has the same types as on the square and the cube,
# and only the constraint on its values changes: the representative inside the unit ball.
#
# The moment equations are written in orthonormal invariant polynomials of the disk and the
# ball, so that, as on the box, nothing has to be orthogonalised and only the constant has a
# nonzero integral.
#
# On the disk, the Zernike polynomials p_k(2r² − 1) r^m cos mθ with p_k orthonormal for the
# weight (1 + t)^m are orthogonal, and the invariant ones are those with m ≡ 0 (mod 4):
# r^m cos mθ = Re (x + iy)^m, and cos 4kθ is unchanged by θ → −θ, π − θ and π/2 − θ.
#
# On the ball, p_k(2r² − 1) r^ℓ F(x̂) with F an O_h-invariant spherical harmonic of degree ℓ
# (symmetry/octahedral_harmonics.jl) and p_k orthonormal for (1 + t)^(ℓ + 1/2): r^ℓ F(x̂) is a
# harmonic polynomial, so these are polynomials of degree 2k + ℓ, and orthonormal on the ball.

"Orthonormal Jacobi polynomials `p₀ … p_K` of `(1 − t)^α (1 + t)^β` on `[−1, 1]` at `t`, and their derivatives."
function jacobi_orthonormal_d(K::Integer, α, β, t::S) where {S}
    p, dp = zeros(S, K + 1), zeros(S, K + 1)
    p[1] = one(S) / sqrt(S(jacobi_mass(α, β)))
    Sα, Sβ = S(α), S(β)
    for k in 0:(K - 1)
        a, bnext = jacobi_a(k, Sα, Sβ), jacobi_b(k + 1, Sα, Sβ)
        bk = k == 0 ? zero(S) : jacobi_b(k, Sα, Sβ)
        prev, dprev = k == 0 ? (zero(S), zero(S)) : (p[k], dp[k])
        p[k + 2] = ((t - a) * p[k + 1] - bk * prev) / bnext
        dp[k + 2] = (p[k + 1] + (t - a) * dp[k + 1] - bk * dprev) / bnext
    end
    return p, dp
end

"""
    round_terms(D, n)

The orthonormal invariant polynomials of degree at most `n` on the unit disk (`D = 2`) or ball
(`D = 3`), by degree: `(k, m)` for `p_k(2r² − 1) Re (x + iy)^m`, or `(k, ℓ, j)` for
`p_k(2r² − 1) r^ℓ F_{ℓj}(x̂)`.
"""
function round_terms(D::Integer, n::Integer)
    if D == 2
        terms = [(k, m) for m in 0:4:n for k in 0:((n - m) ÷ 2)]
        return sort!(terms; by = t -> (2t[1] + t[2], t[2]))
    elseif D == 3
        terms = [(k, ℓ, j) for ℓ in 0:2:n for j in 1:invariant_harmonic_count(ℓ) for k in 0:((n - ℓ) ÷ 2)]
        return sort!(terms; by = t -> (2t[1] + t[2], t[2], t[3]))
    end
    throw(ArgumentError("fully symmetric rules on the ball are for D = 2 and 3, not $D"))
end

"""
    oh_invariant_harmonics(ℓs, coefs, x̂, tangents) -> (values, derivatives)

The `O_h`-invariant spherical harmonics of the degrees `ℓs` at the unit vector `x̂` (one
vector per degree), and their derivatives along each tangent vector, as in
[`OctahedralHarmonicSystem`](@ref): `F = Cᵀ (f(x̂) + f(Rx̂) + f(R²x̂)) / 3` over the D4h
harmonics `f`, so exactly invariant whatever the rounding of `C`.
"""
function oh_invariant_harmonics(ℓs, coefs, x̂::SVector{3,S}, tangents) where {S}
    nℓ = maximum(ℓs; init = 0)
    P = zeros(S, legendre_length(nℓ))
    vals = [zeros(S, size(C, 2)) for C in coefs]
    ders = [[zeros(S, size(C, 2)) for _ in tangents] for C in coefs]
    s2 = sqrt(S(2))
    for j in 0:2
        y = rotate3(x̂, j)
        z = y[3]
        ρ = sqrt(max(zero(S), (one(S) - z) * (one(S) + z)))
        legendre_normalised!(P, nℓ, z, ρ)
        φ = iszero(ρ) ? zero(S) : atan(y[2], y[1])
        cφ, sφ = iszero(ρ) ? (one(S), zero(S)) : (y[1] / ρ, y[2] / ρ)
        eθ, eφ = SVector{3,S}(z * cφ, z * sφ, -ρ), SVector{3,S}(-sφ, cφ, zero(S))
        dθs = [dot(eθ, rotate3(t, j)) for t in tangents]
        dφs = [dot(eφ, rotate3(t, j)) for t in tangents]
        for (i, (ℓ, C)) in enumerate(zip(ℓs, coefs))
            M = size(C, 1)
            f = zeros(S, M)
            fd = [zeros(S, M) for _ in tangents]
            for (t, m) in enumerate(0:4:ℓ)
                sm = m == 0 ? one(S) : s2
                Pm = P[legendre_index(ℓ, m)]
                cm, sn = cos(m * φ), sin(m * φ)
                f[t] = sm * Pm * cm
                iszero(ρ) && continue                            # every gradient vanishes at the poles
                # ρ dP̄_ℓ^m/dθ = ℓ cos θ P̄_ℓ^m − √((2ℓ+1)(ℓ²−m²)/(2ℓ−1)) P̄_{ℓ−1}^m
                e = ℓ == 0 ? zero(S) : sqrt(S((2ℓ + 1) * (ℓ^2 - m^2)) / S(2ℓ - 1))
                Pl = ℓ >= 1 && m <= ℓ - 1 ? P[legendre_index(ℓ - 1, m)] : zero(S)
                dPdθ = (ℓ * z * Pm - e * Pl) / ρ
                for q in eachindex(tangents)
                    fd[q][t] = sm * (dPdθ * cm * dθs[q] - Pm * m * sn * dφs[q] / ρ)
                end
            end
            vals[i] .+= (C' * f) ./ 3
            for q in eachindex(tangents)
                ders[i][q] .+= (C' * fd[q]) ./ 3
            end
        end
    end
    return vals, ders
end

"""
    RoundMomentSystem(structure, n, S)

The fully symmetric moment system for exactness to degree `n` on the unit disk or ball, for
an orbit structure of [`BoxOrbit`](@ref)s, in the orthonormal invariant polynomials of
[`round_terms`](@ref). Callable: `sys(θ) -> (r, J)`, or `sys(θ; jacobian = false)`.
"""
struct RoundMomentSystem{S}
    structure::BoxStructure
    n::Int
    terms::Vector{Any}
    ℓs::Vector{Int}                 # the ball: the degrees with invariant harmonics, and their coefficients
    coefs::Vector{Matrix{S}}
    mass0::S                        # ∫ of the constant term
end
function RoundMomentSystem(structure::BoxStructure, n::Integer, ::Type{S}) where {S}
    D = structure.D
    terms = round_terms(D, n)
    ℓs = D == 3 ? [ℓ for ℓ in 0:2:n if invariant_harmonic_count(ℓ) > 0] : Int[]
    coefs = [S.(oh_harmonic_coefficients(ℓ)) for ℓ in ℓs]
    # the constant term is c p₀ (times F₀ = 1/√(4π) on the ball) and the domain's volume
    mass0 = D == 2 ? sqrt(S(π)) : sqrt(4 * S(π) / 3)
    return RoundMomentSystem{S}(structure, Int(n), Any[t for t in terms], ℓs, coefs, mass0)
end
n_equations(sys::RoundMomentSystem) = length(sys.terms)
n_unknowns(sys::RoundMomentSystem) = nunknowns(sys.structure)

# the invariant polynomials at x and their derivatives along each direction in dxs
function round_values(sys::RoundMomentSystem{S}, x::SVector{D,S}, dxs) where {S,D}
    n = sys.n
    r2 = dot(x, x)
    t = 2r2 - 1
    vals = zeros(S, length(sys.terms))
    ders = [zeros(S, length(sys.terms)) for _ in dxs]
    if D == 2
        z = Complex{S}(x[1], x[2])
        for m in 0:4:n
            c = sqrt(ldexp(one(S), m + 2) / (m == 0 ? 2 * S(π) : S(π)))
            p, dp = jacobi_orthonormal_d((n - m) ÷ 2, 0, m, t)
            Cm = real(z^m)
            gC = m == 0 ? SVector{2,S}(zero(S), zero(S)) : SVector{2,S}(m * real(z^(m - 1)), -m * imag(z^(m - 1)))
            for (i, (k, mm)) in enumerate(sys.terms)
                mm == m || continue
                vals[i] = c * p[k + 1] * Cm
                for (q, dx) in enumerate(dxs)
                    ders[q][i] = c * (dp[k + 1] * 4 * dot(x, dx) * Cm + p[k + 1] * dot(gC, dx))
                end
            end
        end
    else
        r = sqrt(r2)
        if iszero(r)                                       # the centre: only ℓ = 0 survives
            for (i, (k, ℓ, _)) in enumerate(sys.terms)
                ℓ == 0 || continue
                p, _ = jacobi_orthonormal_d(k, 0, 1 // 2, t)
                vals[i] = sqrt(ldexp(one(S), 2) * sqrt(S(2))) * p[k + 1] / sqrt(4 * S(π))
            end
            return vals, ders
        end
        x̂ = x / r
        tangents = [(dx - x̂ * dot(x̂, dx)) / r for dx in dxs]
        F, dF = oh_invariant_harmonics(sys.ℓs, sys.coefs, x̂, tangents)
        for (iℓ, ℓ) in enumerate(sys.ℓs)
            c = sqrt(ldexp(one(S), ℓ + 2) * sqrt(S(2)))       # √(2^(ℓ + 5/2))
            p, dp = jacobi_orthonormal_d((n - ℓ) ÷ 2, 0, ℓ + 1 // 2, t)
            rl = r^ℓ
            for (i, (k, ll, j)) in enumerate(sys.terms)
                ll == ℓ || continue
                vals[i] = c * p[k + 1] * rl * F[iℓ][j]
                for (q, dx) in enumerate(dxs)
                    radial = dot(x̂, dx)
                    ders[q][i] = c * (dp[k + 1] * 4 * dot(x, dx) * rl * F[iℓ][j] +
                                      p[k + 1] * (ℓ == 0 ? zero(S) : ℓ * r^(ℓ - 1) * radial) * F[iℓ][j] +
                                      p[k + 1] * rl * dF[iℓ][q][j])
                end
            end
        end
    end
    return vals, ders
end

function (sys::RoundMomentSystem{S})(θ::AbstractVector; jacobian::Bool = true) where {S}
    s = sys.structure
    D, L = s.D, length(sys.terms)
    r = zeros(S, L)
    J = jacobian ? zeros(S, L, nunknowns(s)) : zeros(S, 0, 0)
    for (o, off) in zip(s.orbits, param_offsets(s))
        w = S(θ[off + 1])
        sz = orbit_size(o)
        x = SVector{D,S}(representative(o, S[θ[off + 1 + i] for i in 1:nparams(o)]))
        # ∂x/∂vᵢ: the coordinates that carry the orbit's i-th value
        dxs = jacobian ? [SVector{D,S}(ntuple(j -> o.slots[j] == i ? one(S) : zero(S), D)) for i in 1:nparams(o)] :
              SVector{D,S}[]
        vals, ders = round_values(sys, x, dxs)
        r .+= (w * sz) .* vals
        jacobian || continue
        J[:, off + 1] .+= sz .* vals
        for i in 1:nparams(o)
            J[:, off + 1 + i] .+= (w * sz) .* ders[i]
        end
    end
    r[1] -= sys.mass0                                       # only the constant integrates to nonzero
    return r, J
end

"""
    round_margins(s, θ) -> (smallest weight, smallest distance of a value from 0 and of a node from the boundary, smallest gap)

As [`box_margins`](@ref), for the unit disk or ball.
"""
function round_margins(s::BoxStructure, θ::AbstractVector{S}) where {S}
    wmin, cmin, gmin = typemax(S), typemax(S), typemax(S)
    for (o, off) in zip(s.orbits, param_offsets(s))
        wmin = min(wmin, θ[off + 1])
        v = θ[(off + 2):(off + nunknowns(o))]
        isempty(v) && continue
        cmin = min(cmin, minimum(v), 1 - sqrt(sum(o.mult[i] * v[i]^2 for i in eachindex(v))))
        for i in eachindex(v), j in (i + 1):length(v)
            gmin = min(gmin, abs(v[i] - v[j]))
        end
    end
    return wmin, cmin, gmin
end
