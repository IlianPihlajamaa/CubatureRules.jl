# Real spherical harmonics on S², and the Fourier basis on the circle (PLAN §8, v0.4).
#
# These are the orthonormal bases in which a rule on a sphere is verified. Harmonics rather
# than monomials, for two reasons: monomials restricted to a sphere are linearly dependent
# (Σxᵢ² = 1), so they are not a basis at all, and the harmonics of degree ≤ d span exactly
# the polynomials of degree ≤ d restricted to the sphere, which is the claim being checked.
#
# Everything is evaluated by the fully normalised recurrence, never by factorials: (l+m)!
# overflows Float64 by l = 85 and costs precision long before that, while the recurrence
# below is stable to the thousands.

"""
    legendre_normalised!(P, n, z, ρ)

The fully normalised associated Legendre functions `P̄_l^m(z)` for `0 ≤ m ≤ l ≤ n`, with
`ρ = √(1 - z²)`, written into `P` indexed by [`legendre_index`](@ref). Normalised so that
`P̄_l^0` is the l-th orthonormal harmonic on the sphere:

    P̄_0^0 = 1/√(4π),  P̄_m^m = -√((2m+1)/(2m)) ρ P̄_{m-1}^{m-1},  P̄_{m+1}^m = √(2m+3) z P̄_m^m
    P̄_l^m = √((4l²-1)/(l²-m²)) ( z P̄_{l-1}^m - √(((l-1)²-m²)/(4(l-1)²-1)) P̄_{l-2}^m )
"""
function legendre_normalised!(P::AbstractVector{S}, n::Integer, z::S, ρ::S) where {S}
    P[legendre_index(0, 0)] = one(S) / sqrt(4 * S(π))
    for m in 1:n
        # the diagonal, carrying the Condon–Shortley phase
        P[legendre_index(m, m)] = -sqrt(S(2m + 1) / S(2m)) * ρ * P[legendre_index(m - 1, m - 1)]
    end
    for m in 0:(n - 1)
        P[legendre_index(m + 1, m)] = sqrt(S(2m + 3)) * z * P[legendre_index(m, m)]
    end
    for m in 0:n, l in (m + 2):n
        a = sqrt(S(4l^2 - 1) / S(l^2 - m^2))
        b = sqrt(S((l - 1)^2 - m^2) / S(4 * (l - 1)^2 - 1))
        P[legendre_index(l, m)] = a * (z * P[legendre_index(l - 1, m)] - b * P[legendre_index(l - 2, m)])
    end
    return P
end

"Index of `P̄_l^m` in the packed triangle, `m ≤ l`."
legendre_index(l::Integer, m::Integer) = (l * (l + 1)) ÷ 2 + m + 1
legendre_length(n::Integer) = ((n + 1) * (n + 2)) ÷ 2

"Index of the real harmonic `Y_l^m`, `-l ≤ m ≤ l`, in the degree-graded ordering."
harmonic_index(l::Integer, m::Integer) = l^2 + l + m + 1
harmonic_length(n::Integer) = (n + 1)^2
"Index range of the degree-`l` block: the `2l+1` harmonics of exact degree `l`."
harmonic_block(l::Integer) = (l^2 + 1):((l + 1)^2)

"""
    real_harmonics!(Y, P, n, x)

The real orthonormal spherical harmonics up to degree `n` at the point `x` on the unit
sphere in `R³`, written into `Y` by [`harmonic_index`](@ref):

    Y_l^0 = P̄_l^0,  Y_l^m = √2 P̄_l^m cos(mφ),  Y_l^{-m} = √2 P̄_l^m sin(mφ)

`P` is scratch of length [`legendre_length`](@ref)`(n)`.
"""
function real_harmonics!(Y::AbstractVector{S}, P::AbstractVector{S}, n::Integer, x) where {S}
    z = S(x[3])
    ρ = sqrt(max(zero(S), (one(S) - z) * (one(S) + z)))       # sin θ, without cancelling near the poles
    legendre_normalised!(P, n, z, ρ)
    φ = atan(S(x[2]), S(x[1]))
    s2 = sqrt(S(2))
    for l in 0:n
        Y[harmonic_index(l, 0)] = P[legendre_index(l, 0)]
        for m in 1:l
            p = s2 * P[legendre_index(l, m)]
            Y[harmonic_index(l, m)] = p * cos(m * φ)
            Y[harmonic_index(l, -m)] = p * sin(m * φ)
        end
    end
    return Y
end

"""
    harmonic_gradient_bound(l, S)

A rigorous bound on `|∇_s Y|` for any orthonormal harmonic of degree `l`, used only to size
the verification tolerance: `√(l(l+1))` times the sup norm `√((2l+1)/4π)` of the harmonic
itself. It is a bound rather than the gradient at the point, so the tolerance it produces
is conservative — by a factor that grows like `l`, which costs about one digit of slack at
degree 20 and never lets a wrong rule pass, only a marginal one fail less often.
"""
harmonic_gradient_bound(l::Integer, ::Type{S}) where {S} =
    sqrt(S(l) * S(l + 1)) * sqrt(S(2l + 1) / (4 * S(π)))

"""
    fourier_circle!(F, n, x)

The orthonormal Fourier basis on the unit circle — `1/√(2π)`, then `cos(kφ)/√π` and
`sin(kφ)/√π` for `k = 1..n` — at the point `x = (cos φ, sin φ)`. Block `k` holds the two
functions of exact degree `k`, so the indexing matches the harmonic one with `2k+1`
replaced by `2`.
"""
function fourier_circle!(F::AbstractVector{S}, n::Integer, x) where {S}
    φ = atan(S(x[2]), S(x[1]))
    F[1] = one(S) / sqrt(2 * S(π))
    c = one(S) / sqrt(S(π))
    for k in 1:n
        F[2k] = c * cos(k * φ)
        F[2k + 1] = c * sin(k * φ)
    end
    return F
end

circle_length(n::Integer) = 2n + 1
circle_block(k::Integer) = k == 0 ? (1:1) : (2k:(2k + 1))
