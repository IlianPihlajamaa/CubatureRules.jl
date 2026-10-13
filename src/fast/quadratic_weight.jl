# The Gauss rule on [0, 1] of the weight s^j / ((s − s₀)² + κ²)^{p/2}, j ∈ {0, 1}, p ∈ {1, 3},
# in Float64 and in microseconds: the one-dimensional step of every fast kernel rule
# (notes/fast-kernel-rules.md). In units of a ray or an edge, it is the radial weight of
# DuffySinh (j = 1) and the angular weight of DuffyGauss and DuffyFinitePart (j = 0).
#
# With s = s₀ + κ sinh z the weight's singularity is gone: ds/√Q = dz and ds/Q^{3/2} =
# dz/(κ cosh z)². Gauss–Legendre in z is a discretisation, and the Stieltjes procedure on it
# gives the recurrence; QL gives the nodes, a Newton step on the recurrence refines them, and
# the weights are Christoffel numbers, a sum of positive terms.
#
# Measured against the certified rules, everything may be Float64 except the z-coordinates.
# s grows like e^z on an interval up to |z| ≈ 35, so a rounding error eps·|z| in z is a
# relative error eps·|z| in s, which moments of degree k amplify k-fold; and the ends must
# satisfy κ (sinh z₁ − sinh z₀) = 1 exactly, or every high moment is off by (1 − s₀) eps. So the
# ends, the Gauss–Legendre nodes and z = z_m + h x are double-double, and sinh and cosh take a
# first-order correction from the low part. The moment errors are then 4e-15 at worst, against
# 1e-15 for the certified rule rounded to Float64 (the full Float64x2 pipeline matches the
# latter, at five to ten times the cost).

# Gauss–Legendre on [−1, 1], nodes in double-double, for the grid sizes the fast rules use.
# This caches an intermediate, not a rule (see "No caching" in the design notes).
const FAST_GL_SIZES = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 20, 24, 28, 32, 40, 48, 56, 64, 80, 96, 112, 128,
                       160, 192, 224, 256)
const FAST_GL = Dict{Int,NTuple{3,Vector{Float64}}}()
const FAST_GL_LOCK = ReentrantLock()

const FAST_GL_READY = Threads.Atomic{Bool}(false)

"Compute every table in `FAST_GL_SIZES` (once; a builder calls this when it is made)."
function fast_gl_init()
    FAST_GL_READY[] && return nothing
    lock(FAST_GL_LOCK) do
        FAST_GL_READY[] && return nothing
        for M in FAST_GL_SIZES
            haskey(FAST_GL, M) && continue
            x, w = setprecision(() -> gauss_jacobi_work(M, 0, 0, 160)[1:2], BigFloat, 160)
            xh = Float64.(x)
            FAST_GL[M] = (xh, Float64.(x .- BigFloat.(xh)), Float64.(w))
        end
        FAST_GL_READY[] = true
    end
    return nothing
end

"The `M`-point Gauss–Legendre rule on `[−1, 1]` as `(x_hi, x_lo, w)`, `M` one of `FAST_GL_SIZES`."
function fast_gl(M::Int)
    FAST_GL_READY[] || fast_gl_init()
    return FAST_GL[M]                                 # read-only once ready: no lock
end
"The smallest size in `FAST_GL_SIZES` that is at least `M` (the largest if none is)."
fast_gl_size(M) = (i = findfirst(>=(M), FAST_GL_SIZES); i === nothing ? last(FAST_GL_SIZES) : FAST_GL_SIZES[i])

"""
    QuadraticWeightWork(mmax)

Buffers for [`quadratic_weight_rule!`](@ref): the grid, the recurrence, and the rule, in
`x[1:m]`, `w[1:m]`.
"""
struct QuadraticWeightWork
    s::Vector{Float64}; λ::Vector{Float64}; p0::Vector{Float64}; p1::Vector{Float64}
    α::Vector{Float64}; β::Vector{Float64}; d::Vector{Float64}; e::Vector{Float64}
    x::Vector{Float64}; w::Vector{Float64}
end
function QuadraticWeightWork(mmax::Integer)
    M = last(FAST_GL_SIZES)
    return QuadraticWeightWork(zeros(M), zeros(M), zeros(M), zeros(M), zeros(mmax), zeros(mmax),
                               zeros(mmax), zeros(mmax), zeros(mmax), zeros(mmax))
end

"Half the length of the interval in `z`, `s = s₀ + κ sinh z` (or `s − s₀ = e^z` for `κ = 0`)."
quadratic_weight_halflength(s0, κ) = κ > 0 ? (asinh((1 - s0) / κ) - asinh(-s0 / κ)) / 2 : log1p(-1 / s0) / 2

"""
    quadratic_weight_gridsize(m, L, j, p) -> M

The Gauss–Legendre grid for an `m`-point rule, `L` half the length of the interval in `z`:
calibrated on the weight's moments over `s₀ ∈ [−30, 4]`, `κ ∈ [1e-14, 30]` and `m ≤ 10`, where
`M − 8 − 2m` was at most `10.5 L` for three of the weights and `17.4 L` for `q^{-3/2}` in the
angular direction, and taken with a margin (scripts/validate_fast_kernel_rules.jl).
"""
function quadratic_weight_gridsize(m, L, j, p)
    c = (j, p) == (0, 3) ? 20.0 : 12.0
    return fast_gl_size(ceil(Int, 16 + 2m + c * min(L, 6) + (c / 3) * max(L - 6, 0)))
end

# sinh x and cosh x from one expm1, both to a few ulps for |x| ≤ 700 (no cancellation: for x < 0
# both terms of the sinh are negative)
@inline function sinh_cosh(x::Float64)
    em = expm1(x)
    ep = em + 1
    return (em + em / ep) / 2, (ep + 1 / ep) / 2
end

function quadratic_weight_grid!(W::QuadraticWeightWork, M, s0, κ, j, p)
    xh, xl, gw = fast_gl(M)
    s, λ = W.s, W.λ
    if κ > 0
        a0 = dd_div((-s0, 0.0), κ)                        # sinh z₀ and sinh z₁ = sinh z₀ + 1/κ
        a1 = dd_add(a0, dd_div((1.0, 0.0), κ))
        z0 = dd_asinh(a0)
        dz = dd_add(dd_asinh(a1), dd_neg(z0))
        hz = (dz[1] / 2, dz[2] / 2)
        @inbounds for i in 1:M
            o = two_sum(1.0, xh[i])
            v = dd_mul(hz, (o[1] / 2, (o[2] + xl[i]) / 2))  # (z − z₀)/2
            u = dd_add(z0, v)                               # (z + z₀)/2
            su, cu = sinh_cosh(u[1])
            sv, cv = sinh_cosh(v[1])
            si = 2κ * (cu + su * u[2]) * (sv + cv * v[2])
            f = hz[1] * gw[i] * (j == 1 ? si : 1.0)
            if p == 3
                z = dd_add(u, v)
                sz, cz = sinh_cosh(z[1])
                cz = κ * (cz + sz * z[2])
                f /= cz * cz
            end
            s[i], λ[i] = si, f
        end
    else                                                # κ = 0, s₀ < 0: s − s₀ = e^z
        z0 = dd_log((-s0, 0.0))
        dz = dd_add(dd_log(two_sum(1.0, -s0)), dd_neg(z0))
        hz = (dz[1] / 2, dz[2] / 2)
        @inbounds for i in 1:M
            o = two_sum(1.0, xh[i])
            v = dd_mul(hz, (o[1] / 2, (o[2] + xl[i]) / 2))
            em = expm1(2v[1]) + exp(2v[1]) * 2v[2]          # e^{z − z₀} − 1
            si = -s0 * em
            f = hz[1] * gw[i] * (j == 1 ? si : 1.0)
            if p == 3
                ez = -s0 * (1 + em)
                f /= ez * ez
            end
            s[i], λ[i] = si, f
        end
    end
    return W
end

function quadratic_weight_stieltjes!(W::QuadraticWeightWork, M, m)
    s, λ, p0, p1, α, β = W.s, W.λ, W.p0, W.p1, W.α, W.β
    t0, t1 = 0.0, 0.0
    @inbounds for i in 1:M
        p0[i] = 0.0
        p1[i] = 1.0
        t0 += λ[i]
        t1 += λ[i] * s[i]
    end
    β[1], α[1] = t0, t1 / t0
    nrm = t0
    @inbounds for k in 1:(m - 1)
        a, b = α[k], k == 1 ? 0.0 : β[k]
        n2, xn2 = 0.0, 0.0
        for i in 1:M
            pn = (s[i] - a) * p1[i] - b * p0[i]
            p0[i] = p1[i]
            p1[i] = pn
            l = λ[i] * pn * pn
            n2 += l
            xn2 += l * s[i]
        end
        α[k + 1], β[k + 1] = xn2 / n2, n2 / nrm
        nrm = n2
    end
    return W
end

# Eigenvalues of the symmetric tridiagonal matrix with diagonal d[1:n] and off-diagonal
# e[2:n], by the implicit QL algorithm (Numerical Recipes' tqli, without eigenvectors).
function tridiagonal_eigenvalues!(d, e, n)
    @inbounds begin
        for i in 2:n
            e[i - 1] = e[i]
        end
        e[n] = 0.0
        for l in 1:n
            iter = 0
            while true
                mm = l
                while mm < n
                    abs(e[mm]) <= eps() * (abs(d[mm]) + abs(d[mm + 1])) && break
                    mm += 1
                end
                mm == l && break
                (iter += 1) > 60 && throw(RefinementError("fast kernel rules", "the tridiagonal QL iteration did not converge"))
                g = (d[l + 1] - d[l]) / (2e[l])
                r = hypot(g, 1.0)
                g = d[mm] - d[l] + e[l] / (g + copysign(r, g))
                s, c, p = 1.0, 1.0, 0.0
                i = mm - 1
                while i >= l
                    f = s * e[i]
                    b = c * e[i]
                    r = hypot(f, g)
                    e[i + 1] = r
                    if r == 0
                        d[i + 1] -= p
                        e[mm] = 0.0
                        break
                    end
                    s = f / r
                    c = g / r
                    g = d[i + 1] - p
                    r = (d[i] - g) * s + 2c * b
                    p = s * r
                    d[i + 1] = g + p
                    g = c * r - b
                    i -= 1
                end
                r == 0 && i >= l && continue
                d[l] -= p
                e[l] = g
                e[mm] = 0.0
            end
        end
    end
    return d
end

function quadratic_weight_gauss!(W::QuadraticWeightWork, m)
    α, β, d, e = W.α, W.β, W.d, W.e
    @inbounds for k in 1:m
        d[k] = α[k]
        e[k] = k == 1 ? 0.0 : sqrt(β[k])
    end
    tridiagonal_eigenvalues!(d, e, m)
    @inbounds for i in 1:m
        x = d[i]
        # Newton on the monic πₘ, then the Christoffel number from the orthonormal p̂ₖ
        for _ in 1:2
            pm1, pp, dpm1, dp = 0.0, 1.0, 0.0, 0.0
            for k in 1:m
                b = k == 1 ? 0.0 : β[k]
                pn = (x - α[k]) * pp - b * pm1
                dpn = pp + (x - α[k]) * dp - b * dpm1
                pm1, pp, dpm1, dp = pp, pn, dp, dpn
            end
            x -= pp / dp
        end
        q0, q1, acc = 0.0, 1.0, 1.0
        for k in 1:(m - 1)
            b = k == 1 ? 0.0 : sqrt(β[k])
            qn = ((x - α[k]) * q1 - b * q0) / sqrt(β[k + 1])
            q0, q1 = q1, qn
            acc += qn * qn
        end
        W.x[i] = x
        W.w[i] = β[1] / acc
    end
    return W
end

"""
    quadratic_weight_rule!(W, s0, κ, m; j = 1, p = 1) -> W

The `m`-point Gauss rule on `[0, 1]` of `s^j / ((s − s₀)² + κ²)^{p/2}` into `W.x[1:m]` and
`W.w[1:m]` (unordered), in Float64; `κ ≥ 0`, and `s₀ < 0` when `κ = 0`. The grid size comes
from [`quadratic_weight_gridsize`](@ref). Allocation-free.
"""
function quadratic_weight_rule!(W::QuadraticWeightWork, s0::Float64, κ::Float64, m::Int; j::Int = 1, p::Int = 1)
    M = max(quadratic_weight_gridsize(m, quadratic_weight_halflength(s0, κ), j, p), fast_gl_size(m + 1))
    quadratic_weight_grid!(W, M, s0, κ, j, p)
    quadratic_weight_stieltjes!(W, M, m)
    quadratic_weight_gauss!(W, m)
    return W
end
