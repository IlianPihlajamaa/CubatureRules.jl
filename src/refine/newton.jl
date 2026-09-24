# Generic refinement drivers (PLAN §1, §2.7, §6 Tier 3).
#
# Minimal symmetric rules sit at or near the degeneracy locus of the moment map, so the
# linear solve is rank-revealing by default (truncated SVD in Float64, column-pivoted QR in
# extended precision) — Gauss–Newton with a pseudoinverse — not as a fallback. The condition number of the Jacobian is measured at
# every iteration and drives the guard digits.

"""
    RefinementError(family, message)

Thrown when refinement does not converge. The registry never falls back to returning an
unrefined seed.
"""
struct RefinementError <: Exception
    family::String
    msg::String
end
Base.showerror(io::IO, e::RefinementError) = print(io, "refinement failed for ", e.family, ": ", e.msg)

"""
    RefineResult

Outcome of a refinement: the parameters, the final residual norm, the number of
iterations, the condition-number estimate of the final Jacobian (and the largest seen),
whether it converged, and the per-iteration history `(residual, step, cond)`.
"""
struct RefineResult{S}
    θ::Vector{S}
    residual::S
    iterations::Int
    cond::Float64
    cond_max::Float64
    converged::Bool
    history::Vector{NTuple{3,Float64}}
end

"""
    estimate_cond(J)

Condition number `σ₁/σ_min` of `J`. Computed from a Float64 copy (cheap, and ample for a
number whose only uses are choosing guard digits and being reported); if that copy is
too ill-conditioned for Float64 to resolve, it is estimated in 256-bit arithmetic from a
pivoted QR, by power and inverse iteration on `R` ([`triangular_cond_mp`](@ref)).
"""
function estimate_cond(J::AbstractMatrix)
    σ = svdvals(Float64.(J))
    κ = iszero(σ[end]) ? Inf : σ[1] / σ[end]
    if !(κ < 1e12) && !(eltype(J) <: Float64)
        # R of a pivoted QR has the singular values of J; iterating on it costs O(n²) a step,
        # where a BigFloat SVD cost O(n³) and was most of a large Lebedev build
        κ = with_bits(256) do
            triangular_cond_mp(pivoted_qr_mp(J)[1])
        end
    end
    return κ
end

"""
    lsq_step(J, r; rank_rtol) -> (Δ, cond, rank)

Rank-revealing least-squares step for `J Δ = r`. In `Float64`: truncated SVD, minimum-norm
solution. In extended precision: Householder QR with column pivoting (pure Julia, so the
result is identical on every platform), truncated where `|R_kk| < rank_rtol · |R_11|`, with
the free variables set to zero.

`cond` differs between the two. In `Float64` the SVD is already in hand, so it is
`σ₁ / σ_min` over all singular values. In extended precision it comes from the pivoted
diagonal instead — free, since that is what the rank test already reads, but a lower bound
rather than the true ratio. Either way a rank-deficient Jacobian reports a large number,
which is the honest answer: the extended-precision estimate reaches one entry past the
truncation, the entry that failed the rank test. Callers needing the true ratio ask
[`estimate_cond`](@ref) directly; [`gauss_newton`](@ref) does, at the first and last iterate.
"""
function lsq_step(J::AbstractMatrix{Float64}, r::AbstractVector{Float64}; rank_rtol)
    F = svd(J)
    σ = F.S
    keep = σ .> rank_rtol * σ[1]
    c = F.U' * r
    y = [keep[i] ? c[i] / σ[i] : 0.0 for i in eachindex(σ)]
    κ = iszero(σ[end]) ? Inf : σ[1] / σ[end]
    return F.V * y, κ, count(keep)
end

function lsq_step(J::AbstractMatrix{S}, r::AbstractVector{S}; rank_rtol) where {S}
    m, n = size(J)
    if S === BigFloat && m == n
        fast = mixed_precision_step(J, r)
        fast === nothing || return fast
    end
    if S === BigFloat
        R, τ, perm = pivoted_qr_mp(J)                 # R in the upper triangle
        qt = b -> apply_qt_mp(R, τ, b)
    else
        F = qr(J, ColumnNorm())
        R, perm = F.R, F.p
        qt = b -> F.Q' * b
    end
    k = min(m, n)
    rank = 0
    for i in 1:k
        abs(R[i, i]) > rank_rtol * abs(R[1, 1]) || break
        rank = i
    end
    c = qt(r)[1:rank]
    z = UpperTriangular(R[1:rank, 1:rank]) \ c
    Δ = zeros(S, n)
    Δ[perm[1:rank]] .= z
    # The condition number here is the free one: column pivoting leaves |R_ii| non-increasing,
    # so the ratio of its ends measures the solve that was actually performed. It is a lower
    # bound on σ₁/σ_min, and deliberately so — `estimate_cond` costs a second factorisation in
    # extended precision (on a 352 × 352 system, 1.6× the QR above), which is too much to pay
    # every iteration for a diagnostic. `gauss_newton` pays it twice per run instead.
    # A rank-deficient Jacobian must still report a large number rather than the flattering
    # conditioning of its retained columns, so the estimate reaches one past the truncation:
    # that entry is what failed the rank test, and its ratio is the honest signal.
    κ = if rank == 0
        Inf
    elseif rank < k
        d = abs(R[rank + 1, rank + 1])
        iszero(d) ? Inf : Float64(abs(R[1, 1]) / d)
    else
        Float64(abs(R[1, 1]) / abs(R[rank, rank]))
    end
    return Δ, κ, rank
end

# --- the step in mixed precision ----------------------------------------------------------
#
# A column-pivoted QR in BigFloat costs O(n³) allocating operations: at triangle degree 40
# it was 58% of a build. A Newton step does not need it when the Jacobian is well
# conditioned. Factor a Float64 copy, solve, and refine the step in BigFloat: each
# correction costs one matrix-vector product and gains about 16 − log₁₀ κ digits, so the
# step reaches working precision in a few O(n²) passes. The result solves J Δ = r to working
# precision, as the QR does.
#
# The Float64 factorisation is written out in plain loops, with no BLAS and no SIMD
# reductions, so its bits are the same on every machine and the refinement stays
# reproducible, as the pure-Julia QR was chosen to be. Square, full-rank and κ < 1e10 only;
# anything else (Lebedev at degree 125 has κ ≈ 1e55) takes the BigFloat QR.

"""
    PivotedQR64

Householder QR with column pivoting of a `Float64` matrix, computed in plain loops (see
above): `R` in the upper triangle, the Householder vectors below it with an implicit unit
first entry, their coefficients `τ`, and the column permutation `p`.
"""
struct PivotedQR64
    QR::Matrix{Float64}
    τ::Vector{Float64}
    p::Vector{Int}
end

function pivoted_qr64(A::AbstractMatrix{Float64})
    m, n = size(A)
    QR = Matrix{Float64}(A)
    τ = zeros(min(m, n))
    p = collect(1:n)
    for k in 1:min(m, n)
        # pivot: the remaining column of largest norm
        best, jbest = -1.0, k
        for j in k:n
            s = 0.0
            for i in k:m
                s += QR[i, j] * QR[i, j]
            end
            s > best && ((best, jbest) = (s, j))
        end
        if jbest != k
            for i in 1:m
                QR[i, k], QR[i, jbest] = QR[i, jbest], QR[i, k]
            end
            p[k], p[jbest] = p[jbest], p[k]
        end
        normx = sqrt(best)
        normx == 0 && continue
        x1 = QR[k, k]
        β = x1 >= 0 ? -normx : normx
        τ[k] = (β - x1) / β
        scale = 1 / (x1 - β)
        for i in (k + 1):m
            QR[i, k] *= scale
        end
        QR[k, k] = β
        for j in (k + 1):n
            w = QR[k, j]
            for i in (k + 1):m
                w += QR[i, k] * QR[i, j]
            end
            w *= τ[k]
            QR[k, j] -= w
            for i in (k + 1):m
                QR[i, j] -= w * QR[i, k]
            end
        end
    end
    return PivotedQR64(QR, τ, p)
end

"The solution of `A x = b` for square, full-rank `A` from its [`PivotedQR64`](@ref)."
function solve64(F::PivotedQR64, b::Vector{Float64})
    QR, τ = F.QR, F.τ
    n = size(QR, 2)
    c = copy(b)
    for k in 1:n                                   # c = Qᵀ b
        w = c[k]
        for i in (k + 1):length(c)
            w += QR[i, k] * c[i]
        end
        w *= τ[k]
        c[k] -= w
        for i in (k + 1):length(c)
            c[i] -= w * QR[i, k]
        end
    end
    z = zeros(n)                                   # R z = c
    for k in n:-1:1
        s = c[k]
        for j in (k + 1):n
            s -= QR[k, j] * z[j]
        end
        z[k] = s / QR[k, k]
    end
    x = zeros(n)
    x[F.p] = z
    return x
end

"""
    mixed_precision_step(J, r) -> (Δ, cond, rank) or nothing

The step of [`lsq_step`](@ref) for a square BigFloat system, from a `Float64` factorisation
and iterative refinement; `nothing` when that does not apply (not finite in `Float64`,
`cond ≥ 1e10`, or no convergence), and the caller falls back to the BigFloat QR.
"""
function mixed_precision_step(J::AbstractMatrix{BigFloat}, r::AbstractVector{BigFloat})
    n = size(J, 2)
    J64 = Float64.(J)
    all(isfinite, J64) || return nothing
    F = pivoted_qr64(J64)
    d1, dn = abs(F.QR[1, 1]), abs(F.QR[n, n])
    (d1 > 0 && dn > 1e-10 * d1) || return nothing
    κ = d1 / dn
    prec = precision(BigFloat)
    # the corrections level off near κ u ‖Δ‖, which is what any backward-stable solve
    # (the BigFloat QR included) reaches; asking for less would never finish
    tol = 64 * BigFloat(κ) * ldexp(BigFloat(1), -prec)
    Δ = zeros(BigFloat, n)
    res = collect(r)
    for _ in 1:(2 + ceil(Int, prec * log10(2) / max(16 - log10(κ), 1)))
        s = maximum(abs, res)
        iszero(s) && return Δ, κ, n
        c = solve64(F, Float64.(res ./ s))         # scaled, so a tiny residual cannot underflow
        δ = BigFloat.(c) .* s
        Δ .+= δ
        maximum(abs, δ) <= tol * maximum(abs, Δ) && return Δ, κ, n
        res = r - J * Δ
    end
    return nothing
end

# --- QR in BigFloat, in place -------------------------------------------------------------
#
# For systems too ill-conditioned for the Float64 path (Lebedev reaches κ ≈ 1e55), the step
# and the condition estimate both need a factorisation in extended precision. The generic
# column-pivoted QR recomputes every column norm at every step and allocates on every
# operation, and the condition estimate was a full BigFloat SVD: together 73% of a
# degree-65 Lebedev build, and O(n³) with a large constant. This is Householder QR with
# column pivoting written with the in-place MPFR operations of core/mpfr.jl, in fixed loop
# order, so it is deterministic like the generic one.

"""
    pivoted_qr_mp(J) -> (QR, τ, p)

Householder QR with column pivoting of a BigFloat matrix, on a copy at the current
precision: `R` in the upper triangle, the Householder vectors below it (unit first entry
implicit), their coefficients `τ`, and the column permutation `p`.
"""
function pivoted_qr_mp(J::AbstractMatrix)
    m, n = size(J)
    prec = precision(BigFloat)
    # a genuine copy: `BigFloat(x; precision)` returns `x` itself when the precision already
    # matches, and writing into it in place would overwrite the caller's Jacobian
    A = [mp_set!(BigFloat(0; precision = prec), BigFloat(x)) for x in J]
    τ = bigfloats(BigFloat, min(m, n))
    p = collect(1:n)
    nrm = bigfloats(BigFloat, n)
    w, t, β, x1 = bigfloats(BigFloat, 4)
    for k in 1:min(m, n)
        best = k
        for j in k:n
            mp_set_si!(nrm[j], 0)
            for i in k:m
                mp_fma!(nrm[j], A[i, j], A[i, j], nrm[j])
            end
            nrm[j] > nrm[best] && (best = j)
        end
        if best != k
            for i in 1:m
                A[i, k], A[i, best] = A[i, best], A[i, k]
            end
            p[k], p[best] = p[best], p[k]
            nrm[k], nrm[best] = nrm[best], nrm[k]
        end
        iszero(nrm[k]) && continue
        mp_sqrt!(β, nrm[k])
        mp_set!(x1, A[k, k])
        signbit(x1) || mp_sub!(β, zero(BigFloat), β)                 # β = −sign(x₁) ‖x‖
        mp_sub!(t, β, x1); mp_div!(τ[k], t, β)                         # τ = (β − x₁)/β
        mp_sub!(t, x1, β)
        for i in (k + 1):m
            mp_div!(A[i, k], A[i, k], t)                               # v = x / (x₁ − β)
        end
        mp_set!(A[k, k], β)
        for j in (k + 1):n
            mp_set!(w, A[k, j])
            for i in (k + 1):m
                mp_fma!(w, A[i, k], A[i, j], w)
            end
            mp_mul!(w, w, τ[k])
            mp_sub!(A[k, j], A[k, j], w)
            for i in (k + 1):m
                mp_mul!(t, w, A[i, k]); mp_sub!(A[i, j], A[i, j], t)
            end
        end
    end
    return A, τ, p
end

"`Qᵀ b` for the Householder factors of [`pivoted_qr_mp`](@ref)."
function apply_qt_mp(A, τ, b::AbstractVector)
    m = size(A, 1)
    c = [BigFloat(x) for x in b]
    for k in eachindex(τ)
        w = c[k] + sum((A[i, k] * c[i] for i in (k + 1):m); init = zero(BigFloat))
        w *= τ[k]
        c[k] -= w
        for i in (k + 1):m
            c[i] -= w * A[i, k]
        end
    end
    return c
end

"""
    triangular_cond_mp(R; iterations = 12) -> Float64

`σ₁ / σ_min` of the upper-triangular `R`, by power iteration on `RᵀR` for `σ₁` and inverse
iteration, two triangular solves a step, for `σ_min`. O(n²) a step, where an SVD is O(n³).
The estimates converge from below, so after a fixed number of steps the result is a lower
bound, in practice within a small factor, which is what choosing guard digits needs.
"""
function triangular_cond_mp(R::AbstractMatrix; iterations::Int = 12)
    n = size(R, 2)
    any(i -> iszero(R[i, i]), 1:n) && return Inf
    U = UpperTriangular(R[1:n, 1:n])
    start() = [BigFloat(1) + BigFloat(i) / (2n) for i in 1:n]      # fixed, so deterministic
    x = start()
    σ1² = zero(BigFloat)
    for _ in 1:iterations
        y = U' * (U * x)
        σ1² = norm(y) / norm(x)
        x = y / norm(y)
    end
    x = start()
    σn⁻² = zero(BigFloat)
    for _ in 1:iterations
        y = U \ (U' \ x)
        σn⁻² = norm(y) / norm(x)
        x = y / norm(y)
    end
    return Float64(sqrt(σ1² * σn⁻²))
end

"""
    seed_cond(κ64, bigjac) -> Float64

The condition number to choose guard digits from at a seed. `κ64` is the estimate from the
`Float64` Jacobian, which is cheap but saturates near 1e16. When it is too large for `Float64`
to resolve, `bigjac()` supplies the Jacobian in `BigFloat` and the accurate estimate is taken
from that.

Choosing the guard from a saturated estimate is not safe, only slow: the refinement notices
and runs again with more guard. At Lebedev degree 125 the `Float64` estimate was 5.2e16
against a true 9.0e54, and the first attempt — 45% of the refinement — was thrown away.
One extended-precision factorisation at the seed is much cheaper than that.
"""
function seed_cond(κ64::Real, bigjac)
    κ64 < 1e12 && return Float64(κ64)
    return with_bits(256) do
        estimate_cond(bigjac())
    end
end

"""
    gauss_newton(F, θ0; step_tol, res_floor, rank_rtol, maxiter = 60, cancel = nothing,
                 verbose = 0, initial_cond = nothing)

Gauss–Newton with backtracking line search and a rank-revealing solve, in the arithmetic of
`θ0`. `F(θ)` returns `(r, J)`. Converged when the step falls below `step_tol` or the
residual below `res_floor`. Checks the cancellation token once per iteration.
"""
function gauss_newton(F, θ0::AbstractVector{S}; step_tol, res_floor, rank_rtol,
                      maxiter::Int = 60, cancel = nothing, verbose::Integer = 0,
                      initial_cond::Union{Nothing,Real} = nothing) where {S}
    θ = collect(θ0)
    r, J = F(θ)
    nr = norm(r)
    history = NTuple{3,Float64}[]
    # The accurate condition number is sampled at the two ends of the run and estimated
    # cheaply in between (see `lsq_step`). It chooses guard digits, so it has to be the real
    # σ₁/σ_min and not the QR lower bound; but the Jacobian of a Gauss–Newton run varies
    # smoothly, so the seed and the solution bracket it, and `refine_octahedral` re-runs with
    # more guard if the number that comes back asks for it. A caller that has already
    # measured it at the seed passes `initial_cond` and saves the factorisation.
    κ = initial_cond === nothing ? estimate_cond(J) : Float64(initial_cond)
    κmax = κ
    converged = nr <= res_floor
    it = 0
    if verbose >= 1
        @info @sprintf("Gauss–Newton: %d equations, %d unknowns, %d bits, cond %.2e, residual %.3e",
                       length(r), length(θ), precision(S), κ, Float64(nr))
    end
    while !converged && it < maxiter
        it += 1
        checkcancel(cancel)
        titer = time()
        Δ, κ, _ = lsq_step(J, r; rank_rtol)
        κmax = max(κmax, κ)
        t = one(S)
        accepted = false
        local θn, rn, Jn, nrn
        for halving in 0:29
            θn = θ - t * Δ
            rn, Jn = F(θn)
            nrn = norm(rn)
            if nrn <= (1 - t / 10_000) * nr || nrn <= res_floor
                accepted = true
                break
            end
            verbose >= 2 && @info @sprintf("    backtracking %2d: t = 2^-%d, residual %.3e > %.3e",
                                           halving + 1, halving, Float64(nrn), Float64(nr))
            t /= 2
        end
        step = norm(Δ) * t
        push!(history, (Float64(accepted ? nrn : nr), Float64(step), κ))
        if verbose >= 1
            @info @sprintf("  iter %2d: residual %.3e  step %.3e  cond %.2e  %6.1f s%s",
                           it, Float64(accepted ? nrn : nr), Float64(step), κ, time() - titer,
                           accepted ? "" : "  (no descent)")
        end
        if !accepted
            # no descent possible: converged only if we are already at the noise floor
            converged = norm(Δ) <= step_tol || nr <= 16res_floor
            verbose >= 1 && @info(converged ? "  stopped at the noise floor" :
                                              "  stopped without descent and above the noise floor")
            break
        end
        θ, r, J, nr = θn, rn, Jn, nrn
        converged = step <= step_tol || nr <= res_floor
    end
    # `estimate_cond`, not another `lsq_step`: the step it would compute is thrown away, and
    # on a large system that discarded factorisation costs as much as an iteration.
    κ = estimate_cond(J)
    if verbose >= 1
        @info @sprintf("Gauss–Newton %s after %d iterations: residual %.3e, cond %.2e",
                       converged ? "converged" : "gave up", it, Float64(nr), κ)
    end
    return RefineResult{S}(θ, nr, it, κ, max(κmax, κ), converged, history)
end

"""
    levenberg_marquardt(F, θ0; maxiter = 400, tol = 1e-14)

Levenberg–Marquardt in `Float64`, used for multistart seed search where starting points
are far from any solution. Finishes with Gauss–Newton steps.
"""
function levenberg_marquardt(F, θ0::AbstractVector{Float64}; maxiter::Int = 400, tol::Float64 = 1e-14,
                             accept = θ -> true)
    θ = copy(θ0)
    r, J = F(θ)
    nr = norm(r)
    μ = 1e-3 * maximum(sum(abs2, J; dims = 1))
    for _ in 1:maxiter
        nr <= tol && break
        A = J' * J
        g = J' * r
        Δ = (A + μ * Diagonal(diag(A) .+ 1e-12)) \ g
        θn = θ - Δ
        if !all(isfinite, θn) || !accept(θn)
            μ *= 8
            μ > 1e20 && break
            continue
        end
        rn, Jn = F(θn)
        nrn = norm(rn)
        if nrn < nr
            θ, r, J, nr = θn, rn, Jn, nrn
            μ = max(μ / 3, 1e-15)
        else
            μ *= 4
            μ > 1e20 && break
        end
    end
    return θ, nr
end

"""
    guard_bits_from_cond(κ)

Guard bits for a Newton solve whose Jacobian has condition number `κ`: the digits it
destroys, rounded up to a multiple of 8 bits, plus a fixed margin (PLAN §2.7 — measured,
not assumed). The rounding keeps the working precision identical across platforms even if
the Float64 estimate of κ differs in its last bits.
"""
guard_bits_from_cond(κ::Real) = 32 + 8cld(ceil(Int, log2(clamp(Float64(κ), 1.0, 2.0^1000))), 8)
