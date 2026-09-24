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
too ill-conditioned for Float64 to resolve, the SVD is redone in 256-bit arithmetic.
"""
function estimate_cond(J::AbstractMatrix)
    σ = svdvals(Float64.(J))
    κ = iszero(σ[end]) ? Inf : σ[1] / σ[end]
    if !(κ < 1e12) && !(eltype(J) <: Float64)
        σb = with_bits(256) do
            GenericLinearAlgebra.svdvals!(BigFloat.(J))
        end
        κ = iszero(σb[end]) ? Inf : Float64(σb[1] / σb[end])
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
    F = qr(J, ColumnNorm())
    R = F.R
    k = min(m, n)
    rank = 0
    for i in 1:k
        abs(R[i, i]) > rank_rtol * abs(R[1, 1]) || break
        rank = i
    end
    c = (F.Q' * r)[1:rank]
    z = UpperTriangular(R[1:rank, 1:rank]) \ c
    Δ = zeros(S, n)
    Δ[F.p[1:rank]] .= z
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
